# GAM movement models — design decisions & notes

Notes from building the GAM path of the deer movement pipeline (the SSF-as-GAM
analogue of the existing iSSF path). Captures the *why* behind the choices so
they don't have to be re-derived. Companion file: `kernel_exponent.html` (a
visual walk-through of the redistribution-kernel exponent).

Reference: Klappstein et al. (2024), *Step selection functions with non-linear
and random effects*, Methods Ecol Evol — fit SSFs as a stratified Cox PH model
in `mgcv` (`family = cox.ph`).

---

## 1. Fitting

An (integrated) SSF is fit as a stratified Cox PH model:
`gam(cbind(times, stratum) ~ <terms>, family = cox.ph, weights = obs)`, where
`times` is a constant event time, `stratum` groups each observed step with its
random points, and `obs` is the used/available (1/0) indicator. Fit via REML.
See `prepare_gam_data()` / `fit_gam_mod()` in `helper_functions.R`.

---

## 2. The redistribution kernel for GAMs

**`amt::redistribution_kernel()` does not work on GAMs** — it asserts a
`fit_clogit` object and reads fixed coefficients via `coef()` + `model.matrix()`.
Our models are `gam`/`cox.ph` objects whose effects are penalised smooths
evaluated through `predict.gam()`. So we wrote our own, mirroring amt's
machinery exactly and swapping only the one load-bearing step.

**The math (parametric / gamma–von Mises design).** For each candidate cell the
unnormalised weight is

```
w(cell) = h(sl_, ta_) · exp(η) / sl_
```

- `η = predict(gam, type = "link")` — the cox.ph linear predictor (fitted
  movement corrections + habitat). This replaces amt's `model.matrix %*% coefs`.
- `h` = the **tentative** gamma × von Mises the random steps were sampled from
  (the proposal cloud). It rides on the data as `attr(stp.var, "sl_")` /
  `attr(stp.var, "ta_")`, so **no change to fit was needed** to recover it.
- `/ sl_` = the polar → planar Jacobian.

This is exactly `amt::ssf_weights(..., compensate.movement = TRUE)`. The fitted
GAM coefficients are *tentative* (deviations from `h`); re-adding `h` and the
Jacobian reconstructs the true 2-D movement kernel × habitat — equivalent to
Klappstein Table 1's corrections, done numerically. Normalising over the disc
(`/ Σw`) gives a probability per cell (this is what the log-score reads).

Implemented as `redistribution_kernel_gam()` (+ `gam_kernel_weights()`,
`gam_movement_kernel()`, `gam_cov_fun()`) in `helper_functions.R`, with the same
signature shape and return object as the amt version.

**Validated:** for a movement-only GAM the kernel's implied mean step length
matched the tentative gamma's mean to within ~1 m (134.7 vs 135.6) — confirming
the `h` and Jacobian terms are correct. Cost ≈ 25 ms/kernel (≈925-cell disc);
`predict.gam` is ~5 ms of that, `extract_covariates` ~11 ms is the real cost, so
a GAM kernel is only ~25 % slower than the iSSF equivalent.

---

## 3. Movement-kernel parameterisation

Movement covariates enter **both** as parametric main effects
(`sl_ + log(sl_) + cos(ta_)`) **and** as cyclic-spline interactions with time of
day (the zebra model). The main effects carry the baseline correction
(gamma rate/shape, von Mises concentration) and are never penalised away; the
`s(tod_, bs='cc', by=…)` smooths add the (optional) time-of-day modulation. This
is the smooth-cyclic analogue of the iSSF baseline
`(sl_ + log(sl_) + cos(ta_)):tod_start_`.

Note: the `:`-only iSSF form is **not** under-specified — it is the *identical*
model (same column space, rank-verified) to `main + interaction`, just
relabelled (day/night values vs day + night−day contrast).

---

## 4. Basis dimension `k`

`k` is **not** model-selected like a covariate — in a penalised GAM it only caps
wiggliness; REML's smoothing parameter chooses the actual smoothness
(Klappstein §3.2.1; Wood 2017). So no AIC-over-`k` search. Instead:

- **`K_TOD = 10`** for the cyclic tod smooths, **capped per deer** at
  `(distinct tod − 1)`. **Raised from 8:** at 8 the tod smooth was k-bound
  (edf ~6 vs k'=7) on 84 deer; a `K_TOD` sweep showed edf **stabilises at ~6 by
  k=10** (the flag clears), while k≥12 over-parameterises and triggers
  convergence step-failures. `tod_` is the decimal *start hour*, and the per-deer
  cap uses the **raw** distinct-tod count (~24–50) — inflated because non-4h data
  gaps (every deer has them, up to ~600 h) shift the fix phase off the 4-h grid,
  and timestamps carry ~1-min jitter. The **meaningful** resolution is only ~12
  (≈6 base positions × the phase shifts), so the cap effectively never binds, k=10
  sits safely under it, and higher k does not (it exceeds the real tod support).
- **`K_NDVI = 5`** (per landcover class) for the `fs` smooths — wide headroom
  (max edf/k' ~0.34). The 1-D smooths `s(HR_edge_end)`, `s(HR_center_end)`,
  `s(ndvi_end)` keep the mgcv default `k=10`: none are k-bound (max ratios
  0.71–0.78), and a sweep showed their edf does **not** stabilise as k rises (it
  drifts 6→7→8 from k10→k20) — the regularisation regime, where k is a sensible
  cap and REML's penalty controls smoothness. **Rule of thumb:** raise k only when
  a smooth is *k-bound and its edf stabilises* once given room (tod); leave it when
  it is *not k-bound and edf keeps drifting* (HR / ndvi).
- A **post-fit k-audit** (`gam_smooth_diag()` → `k.check` edf vs k') reports any
  smooth pressed against its ceiling. The k-index / p-value columns are ignored:
  they are residual-based and unreliable for cox.ph SSFs (the same reason
  Klappstein advises against PH residuals). The trustworthy signal is edf vs k'.

---

## 5. Habitat smooths: Model GS (global + group)

NDVI-by-landcover and HR-center-by-landcover use a **global smooth + hierarchical
`fs` deviations** (Pedersen et al. 2019 "Model GS"): `s(ndvi_end) +
s(ndvi_end, wiscland_end, bs='fs')`. The global smooth mirrors the iSSF main
effect and lets sparse landcover classes pool toward the shared response instead
of toward zero. The `fs` already carries per-class intercepts (covering the
`wiscland_end` main effect), so it is not listed separately. (See §10 for why
this plain GS was chosen over a veg-restricted hybrid, a per-class random slope,
and Model GI.)

---

## 6. Home-range centre: `_end` only, not `_start`

We model attraction to the home-range centre with `s(HR_center_end)` (a
selection term) and **dropped** `s(HR_center_start, by=…)` (a movement term).
They are *not* redundant twins — and `_end` is the one that matters:

- `s(HR_center_end)` is **selection**: cells nearer the centre get higher weight.
  Because candidate cells toward the centre are systematically *closer* to it
  than cells away, this produces a **drift** toward the centre (mean reversion) —
  the thing that bounds a home range. Quantified: ~+54 m/step toward centre,
  **at every heading**.
- `s(HR_center_start, by=…)` only modulates step length/persistence by current
  position; its angles are relative to the deer's **heading, not the centre**, so
  it produces **no net homing** (heading-averaged ≈ 0; it just amplifies the
  current direction — diffusion, not drift).

So `_end` → bounded home range, `_start` → diffusive walk. For UD/log-score
comparison they are not interchangeable; `_end` is kept.

---

## 7. NDVI is unusable in winter → season-conditional models

**Finding:** 107 of 372 deer have NDVI NAs — and they are *exactly* the 107
non-breeding (`nb`, winter) deer (0 in `fa`/`pf`; every `nb` deer affected,
35 of them 100 % NA). Cause: MODIS NDVI is masked under snow / dormancy. So the
NDVI models can't be fit for `nb` deer (this caused the old errors).

**Resolution:** `make_formulas()` builds **6 models per deer, named by model
number `"1","2","3","4","5","6"`**:

| slot | non-`nb` (fa/pf) | `nb` (winter) |
|------|------------------|---------------|
| 1,2,3 | move / +HR-edge / +HR-center | same |
| **4** | NDVI (with HR-center) | landcover (with HR-center) |
| **5** | NDVI (no HR-center) | landcover (no HR-center) |
| 6 | HR-center × landcover | same |

Slots 4 & 5 always mean "the resource-selection model" (NDVI where it exists, a
landcover-only substitute where it doesn't). The numbers are contiguous `"1"…"6"`,
so the results object is a **named list of 6** in which position = model number.

`fit_GAM.R` uses this `1–6` scheme (the HR-center × landcover model is **6**);
treat it as canonical. `fit_issf.R` now uses the **same** contiguous
`"1","2","3","4","5","6"` numbering (its HR-center × landcover model is likewise
**6**), so GAM and iSSF model numbers correspond one-to-one.

---

## 8. Convergence & stability

A diagnostic sweep (4 model variants × 3 deer) isolated **two** sources of the
"iteration limit" / "step failure" warnings:

1. **The three cyclic move smooths** (`by=sl_`, `by=log(sl_)`, `by=cos(ta_)`):
   `sl_` and `log(sl_)` are monotonically related, so their tod-smooths are
   near-collinear → REML can't separate the smoothing parameters → iteration
   limit, ~12× slower. The audit also showed `select` was shrinking most of these
   away anyway. **Recommended fix (deferred): reduce to the single `by=sl_`
   smooth (zebra form).**
2. **The unpenalised 10-level `wiscland_end` factor:** for a deer using only a
   few classes, sparse levels are near-separated in cox.ph → step failure. Fixed
   by `s(wiscland_end, bs='re')` — a **random effect, which is a penalised smooth**
   (the Klappstein equivalence): per-class intercepts shrunk toward the mean.
   Validated ~17× faster, no warnings. Applied to the `nb` landcover substitutes.
   (Not available in the iSSF — `clogit` has no penalised RE without `glmmTMB`/
   `coxme` — so the iSSF `nb` models keep the fixed factor.)

### `select = TRUE` → dropped for production

`select` (Marra & Wood double penalty) can shrink a smooth past linear to zero.
Useful for exploration, but dropped for the production run because:
- With a fixed, deliberate model set, **between-model filtering is the selection
  step**; we don't also want per-deer *within*-model term removal, which makes
  "model 4" a different effective structure on different deer and muddies the
  comparison.
- It's the wrong level of selection (within- vs between-model) and redundant with
  AIC/log-score filtering.
- Bonus: lighter optimisation → faster and more reliable convergence.
- Regularisation is retained — REML smoothing penalties still shrink unsupported
  smooths to a line/constant; they just can't be zeroed entirely. (`fs`/`re`
  penalisation is independent of `select`.)

---

## 9. Parallelisation & performance

`fit_GAM.R` parallelises **across deer** (`furrr::future_map`, `detectCores()−1`
workers), not across models. fit is low-memory tabular work (no rasters), 372
independent deer ≫ cores, so deer-level parallelism gives the best
load-balancing (a free worker grabs the next deer rather than waiting on a deer's
slowest model) and avoids re-shipping each deer's table to model-level workers.
Each worker fits one deer and returns its status + k-audit; the main process
reduces them.

---

## 10. Choosing the NDVI × landcover structure (GS vs the alternatives)

*Records the reasoning behind the §5 form, after evaluating the hierarchical
alternatives. Reference: Pedersen et al. (2019), "Hierarchical generalized
additive models in ecology," PeerJ — model taxonomy G / GS / GI / S / I.*

**Requirements.** Models 4/5 (resource selection) need: per-class landcover
selection; an NDVI response that can be nonlinear and differ by class; pooling so
sparse classes borrow strength; and **prediction for a landcover class a deer
never visited** (for simulation / scoring).

**Alternatives evaluated and rejected:**

- **Veg-restricted hybrid** — `s(wiscland_end, bs='re') + s(ndvi_veg, by=is_veg)
  + s(ndvi_veg, veg_class, bs='fs', by=is_veg)`, with a 0/1 `is_veg` indicator
  zeroing the NDVI smooths on classes where NDVI is biologically meaningless
  (forest / wetlands / developed) and a separate all-class `re` for intercepts.
  **Problem:** this bolts a GI-style explicit `re` onto a GS smoother, but the
  `fs` already carries per-class intercepts (Pedersen p17) → veg-class intercepts
  are double-counted and REML drains the `re`. Diagnosed: the landcover `re` was
  "removed" (edf≈0) in **296/523 breeding** fits vs only **14/214 winter** (winter
  has no NDVI block to compete) — i.e. the collapse was a structural artifact, not
  biology.
- **Per-class random slope** (`s(ndvi_slope, veg_class, bs='re')` replacing the
  `fs`) — motivated by edf showing the per-class `fs` is mostly flat (median edf
  **0.85** summed over 6 classes; global smooth median **2.0**, ≈ near-linear).
  **Tested head-to-head on 18 breeding deer:** AIC tie (median ΔAIC 0.5), the
  landcover `re` rescued in only **3/10** collapsed deer, and it **discards
  genuine per-class curvature** where present (`8125_fa_2021`: `fs` edf 2.93
  capturing non-monotonic corn / grassland NDVI responses → the linear slope
  shrank to 0 and lost them). Rejected.
- **Canonical GI** (`s(ndvi_end) + s(ndvi_end, by=wiscland_end, m=1)
  + s(wiscland_end, bs='re')`) — cleanest intercept/shape separation and
  individual per-class wiggliness. **Rejected:** Pedersen p21 — **GI cannot
  predict unobserved group levels**, which is our hard requirement; also fragile
  on sparse classes (one smoothing parameter each, no shape-pooling) and the
  individual wiggliness buys little given how weak the per-class deviations are.

**Chosen — plain Model GS** (§5): `s(ndvi_end) + s(ndvi_end, wiscland_end,
bs='fs')`, all 10 classes, no `is_veg` filter, no separate `re`. It

1. **resolves the intercept problem by construction** — the `fs` is the single
   home for per-class landcover selection, so nothing to double-count or collapse;
2. **predicts unobserved classes** (GS property, p21) via
   `drop.unused.levels = FALSE` (an unvisited class keeps a ~0-shrunk coefficient
   = the population-average fallback);
3. **pools** sparse classes toward the global response;
4. is the **simplest**, canonical form (and what the original code used).

**Trade-off accepted.** GS swaps the hard "NDVI ≡ 0 on non-veg classes" prior for
the `fs`'s soft shrinkage (an uninformative class's deviation shrinks toward the
global). Residual risk: the *shared* global smooth is influenced by the
data-heavy non-veg classes (forest dominates the rows), so if NDVI there carried
*systematic* (not random) structure it could leak into the veg response. Judged
acceptable under a benign-noise assumption; the hard veg filter is the fallback
only if that assumption is ever in doubt.

**Why AIC didn't decide it.** Decomposing the `fs`-vs-`re` AIC gap
(`ΔAIC = Δdeviance + 2·Δedf`) showed `fs` spends its extra edf right at AIC's
break-even rate — each extra edf buys ≈2 log-likelihood, exactly what AIC charges
— so AIC is ~indifferent between the two. The choice therefore rests on structure
and the prediction property, not on fit. (Cf. §8: model-set selection is
between-model, hence `select = FALSE`.)

---

## Key files / functions

| file | what |
|------|------|
| `helper_functions.R` | `redistribution_kernel_gam()`, `gam_kernel_weights()`, `gam_movement_kernel()`, `gam_cov_fun()`; `fit_gam_mod()` (+ `select`, `gam_smooth_diag()` k-audit); `prepare_gam_data()` |
| `fit_GAM.R` | per-deer `make_formulas(k_tod, k_ndvi, season)`; parallel loop; config `K_TOD`/`K_NDVI`/`SELECT`; end-of-run k-audit → `results/k_audit_gam.rds` |
| `fit_issf.R` | season-aware `make_formulas(season)` mirroring the GAM |
| `docs/kernel_exponent.html` | visual walk-through of the kernel exponent (h · exp(η) / sl) |
