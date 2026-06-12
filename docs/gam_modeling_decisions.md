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

- **`K_TOD = 8`** for the cyclic tod smooths, **capped per deer** at
  `(distinct tod − 1)`. With 4-hour fixes the day is sampled at ~12 two-hour
  positions (min ~10), so the cap is normally a no-op at 8.
- **`K_NDVI = 5`** for the `fs` smooths — NDVI is continuous with rich support,
  not tod-constrained; the audit shows wide headroom (edf ratio ~0.15).
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
`wiscland_end` main effect), so it is not listed separately.

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

**Resolution:** `make_formulas()` (in both `fit_GAM.R` and `fit_issf.R`) builds
**6 models per deer, named by model number `"1","2","3","4","5","8"`**:

| slot | non-`nb` (fa/pf) | `nb` (winter) |
|------|------------------|---------------|
| 1,2,3 | move / +HR-edge / +HR-center | same |
| **4** | NDVI model 4 | landcover model 6 |
| **5** | NDVI model 5 | landcover model 7 |
| 8 | HR-center × landcover | same |

Slots 4 & 5 always mean "the resource-selection model" (NDVI where it exists,
landcover-only substitute where it doesn't). Models 6 & 7 never appear under
their own numbers; model 8 keeps its number (hence the 1,2,3,4,5,8 gap). The
results object is therefore a **named list of 6** — downstream code that indexes
by position must read `names()` (model 8 is at position 6).

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

## Key files / functions

| file | what |
|------|------|
| `helper_functions.R` | `redistribution_kernel_gam()`, `gam_kernel_weights()`, `gam_movement_kernel()`, `gam_cov_fun()`; `fit_gam_mod()` (+ `select`, `gam_smooth_diag()` k-audit); `prepare_gam_data()` |
| `fit_GAM.R` | per-deer `make_formulas(k_tod, k_ndvi, season)`; parallel loop; config `K_TOD`/`K_NDVI`/`SELECT`; end-of-run k-audit → `results/k_audit_gam.rds` |
| `fit_issf.R` | season-aware `make_formulas(season)` mirroring the GAM |
| `docs/kernel_exponent.html` | visual walk-through of the kernel exponent (h · exp(η) / sl) |
