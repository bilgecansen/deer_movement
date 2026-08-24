# GAM path — decision inventory and behavioural spec

Two review aids for checking whether the code matches what you actually intended.

Read this **instead of** the code comments, not alongside them. The comments in the
source were written by me, so if I misunderstood something they state the
misunderstanding confidently — they are the least independent part of the
artifact. Everything below is phrased as a claim about the *world* or about
*outputs*, so you can judge it true or false without reading any R.

**Who decided** column: `you` = you asked for it explicitly; `me` = I chose it
without asking; `inherited` = it was already in the codebase before this work.
The `me` rows are where miscommunication concentrates — those are the ones to
read first.

---

## Part 1 — Decision inventory

### 1.1 Model set and structure

| # | Decision | Setting | Alternative not taken | Who |
|---|---|---|---|---|
| 1 | Null model | movement + `s(HR_center_end)`, fit and stored separately per deer | Keep it as a numbered slot | you |
| 2 | HR-edge model | Retired from the candidate set; `HR_edge_end` still computed in the wrangle | Keep as a candidate | you |
| 3 | Numbered models | 4: (1) movement only, (2) HR-center + resource, (3) resource only, (4) HR-center × landcover | — | you |
| 4 | Winter resource models | `nb` substitutes `s(wiscland_end, bs='re')` for NDVI in slots 2 & 3 | Skip winter deer, or fit NDVI anyway | inherited |
| 5 | Landcover interactions | Hierarchical "GS": global smooth + shrunk per-class `bs="fs"` deviations | Independent `by=` smooths per class | inherited |
| 6 | Movement kernel | Parametric: `sl_ + log(sl_) + cos(ta_)`, plus `s(tod_, bs='cc', by=sl_)` | Non-parametric `s(sl_)`/`s(ta_)` on the uniform-disc design | inherited |
| 7 | Landcover reference | `forest`; `open_water` excluded entirely as a step mask | Any other reference class | inherited |

### 1.2 Fitting

| # | Decision | Setting | Alternative not taken | Who |
|---|---|---|---|---|
| 8 | Smoothing selection | `method = "REML"` | GCV/UBRE, ML | inherited |
| 9 | Extra shrinkage | `select = FALSE` | `TRUE` (double penalty, per-deer term removal) | inherited |
| 10 | tod basis | `K_TOD = 10`, capped per deer at `distinct(tod) − 1` | 8 (was k-bound on 84 deer); ≥12 over-parameterises | inherited |
| 11 | fs basis | `K_NDVI = 5` per landcover class | Higher; the audit shows wide headroom | inherited |
| 12 | Unvisited landcover classes | `drop.unused.levels = FALSE` — a class the deer never used gets ~0 deviation, i.e. the population-average response, so prediction on it works | Drop the level (prediction would fail) | inherited |
| 13 | Design used for fitting | `stp.var` (gamma/von Mises, 25 random points/step) | `stp.var.nonp` (uniform disc, 50 points) | inherited |
| 14 | Per-model failure | Caught and stored as a status object; the deer continues | Abort the deer | inherited |
| 15 | Steps actually fit | Only steps with a turning angle. First-in-burst steps and whole 2-step bursts are dropped by `random_steps` | — (amt behaviour) | inherited |

**#15 is worth a hard look.** On the test deer, models are fit on **42** steps but
log-scored on **77**. That is not a bug, but "in-sample log score" includes 35
steps the model never saw, and `n_steps` in the logscore files is not the fitting
sample size.

### 1.3 Simulation

| # | Decision | Setting | Alternative not taken | Who |
|---|---|---|---|---|
| 16 | Paths per model | `n_sim = 10` | More (tighter UD/ES estimates, linearly more time) | inherited |
| 17 | What gets simulated | Numbered models only; the null is never simulated | Simulate the null too | you |
| 18 | Path length | Per burst, matches the observed number of steps in that burst | Fixed length | inherited |
| 19 | Start point | First **observed** location of each burst; month-chunks continue from the previous chunk's end | Random start, or observed start per step | inherited |
| 20 | Step interval | `dt = 4 hours`, hard-coded | Read from the data | inherited |
| 21 | Movement compensation | `compensate.movement = TRUE` — kernel includes `phi − log(sl_)` | FALSE (no tentative-kernel correction) | inherited |
| 22 | Candidate disc radius | `max.dist` = 0.99 quantile of the fitted gamma | Wider (slower), narrower (truncates long steps) | inherited |
| 23 | Off-map tolerance | If >1% of candidate endpoints fall outside the map, the kernel returns NULL and the path stops | Clamp to the edge; larger tolerance | inherited |
| 24 | Raster crop | `CROP_BUFFER_M = 3000` around the track | 5000 (previous value) | you |
| 25 | NDVI resampling | `method = 'near'` when resampling NDVI onto the landcover grid | bilinear | inherited |

### 1.4 Scoring

| # | Decision | Setting | Alternative not taken | Who |
|---|---|---|---|---|
| 26 | Log-score steps | Every observed step in `stp`, including those excluded from fitting (see #15) | Restrict to fitted steps | inherited |
| 27 | Unscoreable steps | `sum(logp, na.rm = TRUE)`; `n_steps` counts only non-NA. A step whose kernel fails is **silently dropped**, not counted as a failure | Propagate NA; count and report drops | inherited |
| 28 | Null scoring | Scored in the **same run** as the numbered models, same rasters, same steps | Separate invocation | me |
| 29 | delta_logp storage | Lives on the numbered rows; the null file holds raw scores only | Store delta in both; compute downstream | me |
| 30 | UD overlap statistic | `ctmm::overlap(...)$CI[1,2,2]` — the point estimate, CI discarded | Keep the interval | inherited |
| 31 | Simulated-path ctmm fits | `ctmm.select` on the observed track; each simulated path fit with that structure as a guess; **failed fits are dropped**, and the remainder averaged | Refit each freely; fail the deer | inherited |
| 32 | SVF agreement | `1 − ∫\|γ_a−γ_b\| / max(∫γ_a, ∫γ_b)`, integrated in log-Δt | Linear-Δt; data-only normalisation | inherited |
| 33 | Energy score matching | Observed and simulated times rounded to the nearest hour to pair them | Exact timestamp matching | inherited |

**#27 is the one I would most want you to look at.** A model whose kernel fails on
some steps gets a *better-looking* total log score than one that succeeds
everywhere, because failures are dropped rather than penalised. The
"n_steps identical across models" check now guards this, and it currently passes
on all 372 GAM deer — but the guard is new, not the behaviour.

### 1.5 Filtering

| # | Decision | Setting | Alternative not taken | Who |
|---|---|---|---|---|
| 34 | Gate order | Sequential: `bat_uds` → `svf_agree` → `delta_logp` → `p_excd`; only survivors flow on | Apply all four independently | inherited |
| 35 | `ud_min`, `svf_min` | 0.8, 0.8 | — | inherited |
| 36 | `delta_logp_min` | **+3** — a model must *beat* the null by 3 | −3 ("don't lose by more than 3"), which is what the HR-edge null used | me |
| 37 | `p_excd_min` | 0.5, i.e. energy score below the deer's median step length | — | inherited |
| 38 | Observed `sl_` for `p_excd` | From the track (**all** steps), not the fitted model frame | Fitted-model steps only (what iSSF does) | me |
| 39 | Null and the gates | The null never enters the gates; there is no null-as-fallback branch | Keep the fallback | you |
| 40 | "Environmental" models | 2, 3, 4; model 1 (movement-only) is the non-env comparator | — | me |

**#36 is a live decision, not a default.** It changed the meaning of gate 3 from
"don't be much worse than the null" to "be meaningfully better". I inferred it
from your earlier null-2/null-3 comparison. If that inference is wrong, this
single number changes which deer pass.

### 1.6 Data preparation

| # | Decision | Setting | Alternative not taken | Who |
|---|---|---|---|---|
| 41 | Random points | 25 per step (`stp.var`) | 50 (`stp.var.nonp`) | inherited |
| 42 | Water | `open_water` cells excluded from random-point sampling | Allow, with an indicator | inherited |
| 43 | HR-centre transform | `log1p(HR_center)` available as a separate layer | Raw metres only | inherited |
| 44 | NDVI layer times | **Mid-month** (changed this session from 1st-of-month) | 1st of month — caused a fit/predict mismatch on 46.7% of days | you |
| 45 | NDVI lookup window | `max_time = 31 days`, `when = "any"` (nearest in time) | `when = "before"`; a shorter window | inherited |

---

## Part 2 — Behavioural spec

Statements about what the pipeline **does**, derived from the code. Each should be
readable as true or false against your intent.

### Fitting

1. For each deer-season-year, five models are fit: one null and four numbered.
2. The null is `movement + s(HR_center_end)` and goes into its own file.
3. A deer's fitted models see only steps that have a turning angle — the first step
   of every burst is excluded, and bursts of two steps are excluded entirely.
4. Every fit uses 25 random points per observed step, drawn from a gamma/von Mises
   proposal, with open-water cells excluded.
5. Each model is a stratified Cox proportional-hazards fit: one stratum per
   observed step, one "event" per stratum, constant event time, used/available
   passed as the censoring weight.
6. In winter (`nb`), models 2 and 3 use a penalised landcover random effect instead
   of NDVI. In every other season they use NDVI plus a per-landcover-class NDVI
   deviation.
7. A model that fails to fit is recorded as a failure and the other models continue.
   Five such failures exist across the current 372 deer.

### Simulation

8. Only the four numbered models are simulated. The null never is.
9. Ten paths are simulated per model per deer.
10. Each simulated path starts at the deer's **observed** first location for that
    burst, and runs for exactly as many steps as the burst has.
11. Each step samples one endpoint from a redistribution kernel evaluated at every
    map cell within the 99th-percentile step length of the current position.
12. The kernel weight at a cell is `exp(habitat effect + log tentative movement
    density − log(step length))`. The final term converts the polar movement
    density to planar coordinates. *(Verified analytically — check K2.)*
13. If more than 1% of candidate endpoints fall off the map, the path stops early
    and that simulation returns nothing.
14. NDVI is swapped to the layer for the step's calendar month whenever the month
    changes within a burst.

### Scoring

15. The log score is evaluated at **every observed step**, including the ones no
    model was fit on.
16. For each step, the kernel is normalised over the disc and the probability at
    the deer's actual endpoint is read off; the score is the sum of those logs.
17. Steps where that fails contribute nothing and are not counted — they lower
    `n_steps` rather than lowering the score.
18. `delta_logp` is each numbered model's total minus the null's total, both
    computed in the same run over the same steps.
19. The null's own scores are written to a separate file and are never read by the
    filtering step.

### Metrics

20. UD overlap is the Bhattacharyya coefficient between the observed track's AKDE
    and the average AKDE of the ten simulated paths.
21. Any simulated path whose ctmm fit fails is dropped, and the metric is computed
    from those that remain — with no record of how many were dropped.
22. SVF agreement compares the model-implied semi-variance curves of observed and
    simulated tracks, integrated over log time-lag; 1 means identical.
23. The energy score pairs observed and simulated positions by rounding both to the
    nearest hour.

### Filtering

24. Only numbered models are filtered. The null is never a candidate.
25. The four gates run in order and each sees only the previous gate's survivors.
26. A model passes gate 3 only if it beats the null by at least 3 log units.
27. A deer with no model passing all four gates contributes nothing — there is no
    fallback.
28. A deer "passed" if at least one of models 2, 3, 4 survives; surviving only with
    model 1 counts as failed.

---

## Things this document deliberately does not claim

- That the modelling choices are *correct*. It records what they are.
- That the checks in `scripts/checks/` are exhaustive. They cover the kernel
  mathematics, the data contract, and cross-output sanity — not the science.
- That passing checks means the code is right. A check that has never been seen to
  fail is weak evidence; only K2 currently ships with a deliberate-break control.
