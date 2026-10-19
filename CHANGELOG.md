# Changelog

## 0.8.5 (YYYY-MM-DD)

Found by a statistical review of the package for the manuscript. No number in the six reproduced analyses, the side-by-side comparisons or the manuscript changes. The fixes do change results for the inputs they concern (zero prior weights, near-collinear designs, separated matched sets, fine grids of large values), which is the point of them.

**Fixed**

- `CrossPred` looked up a prediction value or lag with `np.isclose` at its default relative tolerance of 1e-5, so on a fine grid of large values (exposures near 1e5 at unit spacing, such as air pressure in pascals) `slice_var`, the slice plots, `overlay_slices` and `DLNM.rr_at` silently used the neighbouring row. Values now match to 1e-10 relative, nearest row first, as R's lookup by name does exactly.
- The link detection reported statsmodels' `LogC` link, log(1 - p), as a log link, so `crosspred` exponentiated it and labelled a ratio of (1 - p) a relative risk. It is now reported as `"logc"` and not exponentiated.
- `attr_table` (and so `DLNM.attributable`) accepted a model without a log or logit link and returned meaningless fractions, for example -23,826 for a Gaussian model, where `attrdl` refuses. It now refuses too.

- `fit_glm` counted rows with a zero prior weight (`var_weights` or `freq_weights`) in the residual degrees of freedom, so the quasi-likelihood dispersion and every standard error were too small (by 13% with a quarter of the rows weighted zero). R's `glm.fit` counts only rows with a non-zero weight (`n.ok`), and so does `fit_glm` now. `res.df_resid` is set to the same `n - p`.
- `fit_glm` tested for aliased columns with `lm.fit`'s tolerance, 1e-7. `glm.fit` uses `min(1e-7, epsilon / 1000)`, which is 1e-11 at its default `epsilon`, on the weighted design, so a column that R keeps (relative residual norm between 1e-11 and 1e-7 after orthogonalisation) was dropped, changing the coefficients, the residual degrees of freedom and the dispersion. `fit_glm` now uses `glm.fit`'s tolerance and prior weights. Exact duplicates are still dropped.
- `fit_pgam`/`fit_pglm` made the same zero-weight error in the scale: the starting REML scale and the Fletcher scale divided by a count that included rows with zero prior weight, and the Fletcher correction was averaged over them, so the reported scale and every standard error were too small (by about 14% with a quarter of the rows weighted zero). `nobs` now counts rows with a non-zero weight.
- `fit_clogit`'s step-halving path could stop at a point that is not a maximum, without saying so, when the likelihood had no finite maximum (separated data). It now raises a `ConvergenceWarning` that a coefficient may be infinite, as `coxph` does, and every result has `res.converged`.

**Changed**

- `model_grid` warns when the candidate models were fitted to different numbers of observations. Models with different maximum lags drop different numbers of initial days, so their information criteria are not comparable without refitting on the common rows, as R's `AIC()` also warns.

- `fit_pgam`/`fit_pglm` raise a `ConvergenceWarning` when the smoothing-parameter search does not converge, as mgcv warns when its outer iteration fails. Previously only `converged=False` said so.
- An offline test tier, `tests/test_reference_offline.py`, runs with the ordinary suite. It reruns the second stage of Gasparrini et al. (2015) and of Gasparrini and Armstrong (2013) from R's own first-stage reductions, replays R's simulated coefficients of Vicedo-Cabrera et al. (2019), and checks every `attr_table` component and its intervals. Until now the stored R results of the published analyses were read only by the example scripts.
- `tests/test_reproductions_slow.py` runs the six published reproductions and the side-by-side analyses under pytest (`DLNMPY_SLOW=1 pytest -m slow`, about an hour, needs network once).
- `tools/side_by_side.py` no longer calls its check against the simulated truth "coverage": it reports the share of one replicate's prediction grid at which the band contains the true curve.
- Targeted tests for three earlier defects that had none (D16, D32 and D38 in `docs/defects.md`), in `tests/test_penalized_guards.py`.

**Documentation**

- `docs/theory.md` states what the simulated intervals condition on (coefficient uncertainty only, the MMT as reference, Psi known for BLUPs, Vp for penalised fits, the MMT grid) and how to report their Monte Carlo error. `qaic` documents that it ignores prior weights, as `fqaic` does.
- The REML criterion of `fit_pgam` takes the log pseudo-determinant of the penalty directly from its eigenvalues, where mgcv uses Wood's (2011) stable transform. The docstring and `docs/penalized.md` no longer say the criterion is computed "exactly as" in mgcv, and state that this loses accuracy when smoothing parameters differ by a factor of about e^20 or more.

## 0.8.4 (2026-10-13)

Six published analyses are now rerun from their authors' own code, with the R side committed. Four of them (Gasparrini 2014 *Stat Med*, Gasparrini & Leone 2014 *BMC MRM*, Vicedo-Cabrera, Sera & Gasparrini 2019 *Epidemiology*, Gasparrini et al. 2017 *Biometrics*) were reproduced for 0.8.1 and 0.8.2 without the scripts being committed. They are now `examples/statmed2014.py`, `leone2014.py`, `vicedo2019.py` and `biomet2017.py`, each with a `tools/reproduce_*.R` that pins the authors' repository and writes `tests/reference/*.json`. The two earlier reproductions gained stage-by-stage comparisons with R, `examples/lancet2015_vs_r.py` and `examples/bmcmrm2013_vs_r.py`. There are 161,356 comparisons with R in all, with no failures. `docs/validation.md` has the results and `docs/defects.md` lists every defect found since 0.4.1, with how it was found and the test that guards it. `paper/` holds the manuscript, an executable Quarto document that reruns all of this.

**Fixed**

- `simulate_coef(normals=Z)` was documented to reproduce R's draws from R's normal deviates. It does not in general: an eigenvector is defined up to its sign and numpy's and R's eigenvalue routines choose differently (8 of 20 on the London model of Gasparrini & Leone 2014, moving the attributable-number intervals by up to 1.2%). `simulate_coef` now takes `eigvecs=` (R's `eigen(vcov)$vectors`), which reproduces R's draws exactly. Without it the sign of each eigenvector is now fixed (largest entry positive), so seeded draws are the same on every platform. **Seeded results change** relative to 0.8.3 for `attrdl(sim=True)`, `findmin`/`mmt`, `attr_table`, `bootstrap` and the `DLNM` methods that use them. The new draws are statistically equivalent to the old ones, but the individual numbers differ.
- `fit_clogit` lost the parameter names on its step-halving path (taken when statsmodels' Newton diverges), so `crosspred(cb, res, name=...)` raised "model parameters are not named".
- `int=`, R's partial match of `intercept` used in the authors' scripts, was rejected. It is accepted wherever `Bound=` is.

**Corrections to earlier entries**

- 0.5.0: the QAIC defect was fixed in 0.4.1. Gasparrini & Armstrong (2013) compare their three lag models graphically and do not use QAIC. The authors' script computes it, and the pre-0.4.1 criterion would have ranked the three in the wrong order.
- 0.6.0 and 0.7.0 said every defect was in the Python-only layers. Several were in ported functions (`qtest` in 0.6.0, and `attrdl`, `crossreduce`, `cr`/`cbpen`, `poly`, `crossbasis(group=)` and `MixMeta` in 0.7.0). See `docs/defects.md`.
- 0.6.1: the `__version__` fix first shipped in 0.7.0. The 0.6.1 wheel still reported 0.6.0.
- 0.8.0 was released (tag `v0.8.0`) and was wrongly headed "unreleased".
- 0.8.1: the entry said four papers but covers three (the fourth, *Biometrics* 2017, is 0.8.2). The "197 quantities to 1.8e-10" figure for the *Stat Med* paper cannot be traced. See `docs/validation.md` for what the committed script shows. `Bound=` being rejected was a defect, which the entry called a rough edge.

## 0.8.3 (2026-09-21)

This release changes citation metadata only. No numerical result changes.

- `CITATION.cff` carries the author's ORCID and affiliation, and `.zenodo.json` describes the archive, so each release is deposited in Zenodo with a DOI and can be cited properly.

## 0.8.2 (2026-09-21)

This release comes from reproducing Gasparrini, Scheipl, Armstrong & Kenward (2017) *Biometrics*, the penalised DLNM paper, against mgcv. The penalty matrices, their ranks and the REML/ML criterion reproduce mgcv to machine precision, and at matched smoothing parameters the coefficients agree to 1e-10. The outer smoothing-parameter search is weaker than mgcv's Newton method, however, and it failed in three ways.

**Fixed: a silent wrong answer**

- `fit_pgam`/`fit_pglm` could return all-NaN coefficients, covariance and criterion while reporting `converged=True`, after running the smoothing parameters off to ~1e19. Nothing in the result flagged it, and the failure surfaced later as `crosspred`'s "coef/vcov not consistent with basis matrix", which points at the basis when the fault is in the fit. A non-finite fit now reports `converged=False` and raises a `ConvergenceWarning` naming the smoothing parameters it reached.

**Fixed: fits that could not start**

- The starting smoothing parameters were derived from an *unpenalised* fit. When the penalty is the only thing making the model estimable (a rare-outcome Poisson against a wide cross-basis, which is the paper's second example), that fit diverges, every starting value was NaN, and the search aborted with a bare `LinAlgError: Eigenvalues did not converge`. Starting values now fall back to a lightly penalised fit and then to unweighted cross-products, so a finite start always exists.
- A trial value of the smoothing parameters that makes the penalised information singular no longer aborts the search. It is treated as an infinitely bad point, so the optimiser steps back out (as the meta-analysis fitter already did). If the final fit is singular the error now says so and suggests reducing `df` or passing `sp`.

**Changed**

- The README no longer claims smoothing parameters and coefficients agree with mgcv to 1e-4 end to end. That holds at matched smoothing parameters. The selection itself can differ on flat or near-singular problems. The corrected statement, and the advice to pass `sp=` when mgcv's exact answer is needed, are now in the README and `docs/penalized.md`.

## 0.8.1 (2026-09-21)

These changes were found by reproducing four more landmark papers against R. The DLNM core came through unchanged: the exposure-lag-response framework paper (Gasparrini 2014 *Stat Med*) agrees on 197 quantities to 1.8e-10, the attributable-risk paper (Gasparrini & Leone 2014) to 1e-13, and the projections tutorial (Vicedo-Cabrera, Sera & Gasparrini 2019) to 4.3e-12. The reproductions did find two defects in `fit_clogit`, which serves case-crossover and matched case-control designs, and two rough edges in porting the authors' own scripts.

**Fixed: a silent wrong answer**

- `fit_clogit` returned **all-NaN coefficients without raising** on designs whose columns are large, such as a cross-basis of cumulative exposure reaching 1e4. statsmodels' Newton takes the full step with no safeguard, so the first step overflows `exp(X @ beta)` and never recovers. The failure then surfaced downstream as the unrelated-sounding "coef/vcov not consistent with basis matrix". A step-halving Newton, as `survival::coxph` uses, now takes over when that happens. On the model that failed (Gasparrini 2014 *Stat Med* Table II, model 2, a piecewise-constant lag on cumulative radon exposure) coefficients now agree with `coxph` to 7.1e-10 and standard errors to 2.7e-12.
- `fit_clogit`'s covariance differenced the score at a fixed absolute step of 1e-4. When the design columns are large the coefficients are of order 1e-5, so that step measured the wrong curvature and **every standard error and confidence interval was up to 5e-3 out**, although the docstring claimed 1e-10. The step now follows each coefficient's own standard error, calibrated against R in both regimes: the large-column case improves from 5.6e-03 to 1.8e-10, and the nested case-control fit in `tools/side_by_side.py`, whose standard errors reach 6, improves from 2.1e-08 to 3.3e-09.

**Changed**

These were found while reproducing two more landmark papers against R (Gasparrini & Leone 2014 *BMC MRM*, and the Vicedo-Cabrera, Sera & Gasparrini 2019 *Epidemiology* projections tutorial). Both reproduce exactly (the attributable fractions of the 2014 paper to 1e-13, the 2019 projected fractions to 4.3e-12), but two things made porting the authors' own scripts harder than it should be.

- Basis arguments now accept any prefix of `Boundary.knots`, as R's partial matching does. The published code of the dlnm authors writes `argvar = list(fun = "ns", Bound = range(x))`, which previously raised `TypeError` on the first line of an otherwise faithful translation. `Bound=`, `Bo=`, `Boundary=` and the exact R spellings all work, in `onebasis`, `crossbasis`, `ns` and `bs`. Unknown arguments are still rejected, and giving the same argument under two spellings is an error.
- `simulate_coef(..., normals=Z)` takes a matrix of standard normal deviates instead of drawing its own, so a projection or attribution loop can be replayed exactly: feed in R's draws and the two languages agree without Monte Carlo noise. This was previously only possible through a private function. The docstring now also warns that the result is `(k, nsim)`, the transpose of `MASS::mvrnorm`.

## 0.8.0 (2026-09-21)

The multivariate meta-analysis is now a separate package, [mixmetapy](https://github.com/Logic06183/mixmetapy) 0.1.0, as `mixmeta` is separate from `dlnm` in R. The package had grown to cover mixed-effects meta-analysis, which is a methodology of its own with users who never fit a DLNM. Nothing numerical changed. The code moved unchanged with its R fixtures and history, and the three published analyses (`examples/lancet_2015.py`, `examples/bmcmrm_2013.py`, the Chicago attributable-risk example) give the same numbers to the last printed digit before and after.

- `mixmeta`, `MixMeta`, `vech` and `xpnd` are in mixmetapy. Install with `pip install "dlnmpy[twostage]"`.
- `stack_reduced` and `predict_reduced` stay in dlnmpy, now in `dlnmpy.twostage` (still `dl.stack_reduced`, `dl.predict_reduced`). They take and return plain arrays, so `predict_reduced` works with pooled coefficients from any source.
- dlnmpy no longer needs mixmetapy to import or to run a single-location analysis. The CI `minimal` job checks this without it installed.
- Plotting now needs matplotlib >= 3.6: `summary_figure` crashed on 3.5 (`pcolormesh` rejects non-finite coordinates there).
- Two tests no longer depend on library versions (`links.Probit` is lowercase before statsmodels 0.14, and `axvspan` returns a `Polygon` before matplotlib 3.9).
- A `.gitlab-ci.yml` now sits alongside the GitHub workflows, with the same checks and a PyPI trusted-publishing release job.

### Moving from 0.7

| 0.7 | 0.8 |
|---|---|
| `dl.mixmeta(y, S, ...)`, `dl.meta.mixmeta(...)` | `from mixmetapy import mixmeta` |
| `dl.MixMeta`, `dlnmpy.meta.vech`, `dlnmpy.meta.xpnd` | `mixmetapy.MixMeta`, `mixmetapy.vech`, `mixmetapy.xpnd` |
| `dlnmpy.meta.stack_reduced`, `dlnmpy.meta.predict_reduced` | `dl.stack_reduced`, `dl.predict_reduced` (unchanged) |

The old names keep working through 0.8 and 0.9: with mixmetapy installed they forward to it with a `DeprecationWarning`. Without it they raise an `ImportError` that says `pip install dlnmpy[twostage]`. They are removed in 1.0.

## 0.7.0 (2026-09-14)

This release follows a third audit. It consisted of three independent line-by-line reviews of the Python source against the R source (`dlnm` 2.4.10, `mixmeta`, mgcv, Gasparrini's `attrdl.R` and `findmin.R`), with every suspected defect confirmed by running it, after which the changed paths were compared numerically against R 4.3.3 with `dlnm` 2.4.7: `fit_glm` with an offset and an aliased column (coefficients to 1e-13, dispersion to 5e-10, intervals to 1e-11), `fit_clogit` with missing rows against `survival::clogit` (1e-13), `attrdl(group=)` against `attrdl.R` run per group (daily contributions to 1e-13), and the `dlnm()` one-call workflow, single and grouped, against the equivalent `glm` (all coefficients to 4e-13, attributable fractions to 1e-15). The numerical core again needed no change. Everything below is in the Python layers.

**Fixed: silent wrong answers**

- `fit_glm(..., exposure=)` on a design with an aliased column logged the exposure twice. The refit that drops the aliased column passed statsmodels' stored `exposure` (already `log(exposure)`) back as `exposure`, so every coefficient shifted (intercept 2.17 where R gives -9.98, and RR 1.1515 where R gives 1.1493). `offset=` was unaffected. The exposure is now carried over as an offset.
- `bootstrap_ci` on an uncentred prediction (`cen=False`) centred the simulated draws automatically while leaving the point estimate uncentred, so the fit fell outside its own interval. The draws now follow the prediction.
- `DLNM.rr_at()` and `DLNM.summary()` exponentiated the effect for any link, so a Gaussian model reported a "relative risk" of 8e9 for an additive effect of 26 deaths. The column is now `rr` for log and logit links and `effect`, on the linear-predictor scale, otherwise.
- `attrdl` and `attr_table` on a cross-basis built with `group=` formed the lagged exposures across group boundaries, so the last days of one series contributed to the first days of the next. Both take a `group` argument, and `DLNM.attributable()` passes it. Without it on a grouped basis a warning is raised. The reference `attrdl.R` has no such argument. The grouped total is one pooled fraction, as documented.
- `attr_table` with a reference below the 2.5th (or above the 97.5th) percentile produced a "moderate cold" range that ran backwards (silently 0) and an "extreme cold" range that contained the heat side. The cut-offs are now clamped to `cen`.
- `qaic()` returned a Poisson-likelihood number for a Gaussian or quasi-binomial fit. It now refuses anything but a (quasi-)Poisson family.
- `crossreduce(type="lag"|"var", value=[...])` took the first element of a vector, where R stops. It now raises.
- `cr(knots=)` with fewer than three knots returned a one-column basis with a zero penalty, and `cbpen` then divided by a zero eigenvalue (all-NaN penalty, no error). mgcv refuses this input, and so does dlnmpy now. `_rescale` refuses a penalty with no positive eigenvalue.
- `poly(degree=0)` returned a basis with no columns.
- `crossbasis(group=)` with a missing group label left those rows all-NaN and reported the wrong group count. Such a group is now refused, as is a `group` of the wrong length.

**Fixed: crashes and misleading state**

- `fit_clogit` failed with statsmodels' `MissingDataError` on the NaN rows every time-series cross-basis has. Rows with a missing outcome, design value or stratum are dropped (as `clogit`'s `na.action` does) and counted in `res.n_dropped`.
- `mixmeta` with a missing outcome died inside the moment estimator with `LinAlgError`. It now says what is missing.
- `fit_pglm(sp=1.0)` with a single penalty failed on a 0-d array. Scalars are accepted and the length is checked against the number of penalties.
- `converged` was never `False` on `MixMeta` or `PenalizedGLMResults` (it tested whether the objective was finite). It now reports the optimiser's own status. For penalised fits a negligible REML/ML gradient at the solution also counts, because Nelder-Mead can exhaust its iterations on a flat optimum (seen with scipy 1.18).
- `DLNM.figure()` on a fresh fit drew the curve without the MMT interval band, because it read the cached MMT before `predict()` had computed it.
- `DLNM.mmt()` claimed to cache and did not. Repeated calls with the same arguments now return the same object, and a call with other arguments replaces the default reference.

**Changed**

- `crosspred`/`crossreduce` warn (`CenteringWarning`) when the centring value is chosen automatically, the counterpart of R's `"centering value unspecified"` message. A value stored in `argvar` is the user's own and does not warn. Silence with `warnings.simplefilter("ignore", dlnmpy.predict.CenteringWarning)`.
- `dlnm(..., group=)` now fits a group intercept and a separate seasonal spline per group (block-diagonal `ns(time)`), which is the pooled model an analyst would write by hand. The previous formula fitted one intercept and one spline across the stacked series. The docstring says when to prefer the two-stage design instead.
- `fit_glm`'s quasi-likelihood dispersion with `freq_weights` follows R (`n - p` with `n` the number of rows) instead of statsmodels (`sum(weights) - p`). This was already the behaviour and is now documented.
- CI now covers Python 3.13, has a minimal-install job (numpy, scipy, pandas only), and has a build job that runs `python -m build` and `twine check` on every push. A `release.yml` publishes to PyPI by trusted publishing when a `v*` tag is pushed and the tag matches `pyproject.toml`.
- Two tests that needed statsmodels without saying so now skip when it is absent.
- The README gains a PyPI badge and status, the test count and a network note for the two paper reproductions. `CITATION.cff` and the conda recipe are at 0.7.0, and `docs/r-to-python.md` no longer recommends a private method.

## 0.6.1 (2026-09-02)

**Fixed**

- `dlnm(..., time=)` given a date column killed the interpreter. `datetime64` cast to float is nanoseconds since the epoch, so a fourteen-year daily series was read as 1.2e12 years, the seasonal spline was asked for 8.5e12 degrees of freedom, and `np.linspace` tried to allocate about 68,000 GB. The process was killed by the OS with no Python traceback. Date columns (and `datetime.date` objects) are now converted to day offsets from the first observation, and `time="date"` gives numerically identical results to the equivalent day index. Found by running the published 0.6.0 wheel through a full Chicago analysis in a notebook.
- The seasonal spline now raises `ValueError` when it would need at least as many degrees of freedom as there are observations, which catches any unit error of this kind before it exhausts memory.
- `__version__` was hardcoded in `__init__.py` and had drifted from `pyproject.toml`. It is now read from the installed package metadata, with a test to keep them tied.
- `basis.py` had a conditional whose branches were identical (`thr_value if side == "d" else thr_value`). It is simplified, with no change in behaviour.

**Changed**

- `ruff` is pinned in the `dev` extra and the lint rule set is selected explicitly. Ruff's default selection widens between releases, so a linter upgrade alone could turn CI red, with no change to this code.

## 0.6.0 (2026-09-02)

This release follows a second audit against R (`dlnm` 2.4.7 in a fresh R 4.3 install, `mixmeta` 1.2.2), this time on 64 edge cases the fixtures did not cover: negative and sub-period lag ranges, exposure-history matrices with `cumul`, non-integer `bylag` and `crossreduce(type="lag", value=2.5)`, `group`, every basis function with explicit knots or thresholds, the default and `by` prediction grids, `mkcen` corner cases, `logknots`/`equalknots`, `exphist` with `times` and `fill`, and `attrdl` in the forward, reduced-coefficient and daily forms. All 64 agree to 1e-6 or better (most to 1e-12), and `mixmeta` REML, BLUPs, predictions, Q/I² and univariate Q agree to 1e-8. The numerical core needed no change. The defects below were all in the Python-only layers.

**Fixed**

- `PenalizedGLMResults.edf_by("cb")` summed the effective df of every column whose name *started with* `cb`, so with a second basis called `cb_o3` it silently reported the two together (38.1 where the temperature basis had 9.4). `fit_pgam` used the same match to place the penalties and failed with a shape error in that case. Both now match a basis name the way `crosspred` does (name plus a `v1_l1`/`b1` label).
- `qaic()` raised `AttributeError` on a penalised fit. It now accepts them, using the total effective degrees of freedom as the parameter count (Gasparrini et al. 2017). `PenalizedGLMResults` gained a `fittedvalues` alias.
- `fit_pgam(..., offset=)` with a full-length offset mis-aligned it after rows with missing values were dropped. It now follows the same rows.
- `overlay_slices(labels=...)` with too few labels raised `IndexError` from inside matplotlib. It now says what is wrong.
- `qtest()` p-values were computed as `1 - cdf` and rounded to exactly 0 for large Q. The survival function is now used.
- `bootstrap_ci` is exported.
- `mkat` passes `min.n = nobs/2` to `pretty()` as R does (no case where it changes the grid was found in 20,000 random ranges, so this is a matter of fidelity to R and fixes no observed error).

**One R defect not reproduced**

- `ps(df=4)` (degree 3) builds in R and then fails inside `crosspred()` when the basis is rebuilt from its knots: a lag basis you can never predict from. `dlnmpy.ps` refuses it at construction with the minimum df in the message.

**Beyond R: the workflow**

- `dl.dlnm(data, outcome, exposure, lag, argvar, arglag, time, df_per_year, dow, controls, family, offset, group)` fits the whole model in one call and returns a `DLNM` bound to its data: `.mmt()`, `.rr_at(percentiles)`, `.attributable()`, `.predict(percentiles=...)`, `.reduce()`, `.qaic()`, `.summary()`, `.figure()`, `.plot()`, `.summary_figure()`. Defaults are the usual specification (natural spline with knots at the 10th, 75th and 90th percentiles, a lag spline with log-spaced knots, and 7 df per year of time). `ps`/`cr` bases route to the penalised fitter. It is tested to agree with the explicit pipeline to 1e-12.
- `CrossPred.to_dataframe()` returns the long ("tidy") table of lag-specific (or cumulative) effects.
- Two helpers are added: `percentile_knots(x, [10, 75, 90])` and `percentile_of(x, values)`.

**Beyond R: the figures**

- `plot_overall_risk` draws the exposure-response figure of a temperature-mortality paper (cold and heat sides in two hues, minimum-risk value with its interval, percentile axis along the top, exposure histogram underneath). `plot_attributable` draws attributable fractions or numbers by component with intervals.
- Styling is applied per figure through `plot.style()` (an `rc_context`), so plots look the same everywhere and the caller's rcParams are untouched. The 3D surface is coloured by value with a fine translucent mesh instead of the black wireframe, and no longer clips its z label. The contour's null-effect line is dashed and the frame closed.
- `examples/gallery.py` renders every figure in both themes.


## 0.5.1 (2026-08-30)

This release silences spurious floating-point `RuntimeWarning`s. numpy 2.0.x raises `divide by zero`, `overflow`, `invalid value` and `underflow` warnings from `matmul` on arrays that legitimately carry `NaN` (a cross-basis has `NaN` in its first `lag` rows by construction) and from the log-determinant the penalised REML optimiser evaluates while probing extreme smoothing parameters. The values were never affected, but a routine analysis printed a wall of warnings for anyone on Python 3.9 or 3.10, where numpy resolves to 2.0.x. A full Chicago analysis went from about 1400 warnings to 3, and the test suite from 52 to 1. The three that remain come from statsmodels' own IRLS calling `numpy.linalg.pinv` and are outside dlnmpy. Results are bit-identical: RR at 33 °C is 0.922319 and at -20 °C 1.329074 on numpy 2.0.2 and 2.5.2 alike.


## 0.5.0 (2026-08-30)

This release also reproduces Gasparrini & Armstrong (2013) *BMC Med Res Methodol* 13:1, the methodological paper behind the two-stage design (`examples/bmcmrm_2013.py`): minimum-mortality temperature 17.1 °C, pooled RR 1.101 (1.079-1.124) at 22 °C and 1.309 (1.245-1.376) at 0 °C, I² of 63.8/16.4/63.5%, and all three QAIC totals. Of 68 intermediates, 67 agree with R to 1e-5..1e-15. That paper selects its lag specification by QAIC, and the pre-0.5.0 implementation would have reversed its published conclusion. *(Corrected in 0.8.4: the paper compares its lag models graphically and does not use QAIC. The authors' script computes QAIC, and the pre-0.4.1 criterion would have ranked the three models in the wrong order.)*

This is a bug-fix release from a systematic audit against R (`dlnm` 2.4.10, `splines`, `mgcv`, `mixmeta`, `survival`), covering the basis functions, the prediction machinery, model fitting, meta-analysis, attribution and plotting. Everything below was verified against R and has a regression test that fails without the fix.

**Silent wrong answers**

- `fit_glm` ignored `offset`, `exposure`, `freq_weights` and `var_weights`. They are arguments of the model, but were forwarded to `GLM.fit()`, whose own `**kwargs` discards unknown names, so the caller silently got the un-offset fit. A Poisson model with a population offset returned coefficients `[0.850, 0.426]` where R's `offset()` gives `[0.497, 0.357]`. They are now explicit parameters, and survive both the NaN drop and the aliased-column refit.
- `crossbasis(argvar={..., "cen": v})` discarded the centring value. R keeps it as `attr(cb, "argvar")$cen` and uses it as the default reference for every later `crosspred`/`crossreduce`, so predictions were silently centred on the automatic value instead: a relative risk of 1.343 where R gives 1.166, 15% out.
- `strata` returned a rank-deficient basis when quantile breaks tied (a spike of zeros, say): an empty stratum, so an identically-zero column and a singular design, with nothing to signal it. R's `cut()` raises `'breaks' are not unique`, and so does this now. Values outside the outer edges are `NaN`, as `cut()` gives `NA`. They were previously a row of zeros.
- `plot(..., "slices", lag=L, cumul=True)` read the cumulative outcomes off `pred.predlag` instead of the integer lags they are defined on, so with `bylag != 1` it plotted a different lag: at `bylag=0.2`, `lag=2` drew the curve for lag 10, and `lag >= 4` raised `IndexError`. R sets `bylag` to 1 under `cumul` and indexes by name.
- `get_link` reported a complementary log-log link as `"log"`, because `"log"` is a substring of `"cloglog"` and was tested first, so `crosspred` exponentiated the fit and reported relative risks for a model that has none. The same applied to `loglog`.
- `cen=np.True_` / `np.False_` were coerced to 1.0 / 0.0: `np.bool_` is not a subclass of `bool`, so a numpy boolean fell through both branches of `mkcen`.

**Divergences from R**

- `bs` refused data whose interior knots all coincide with a boundary knot, which `splines::bs` accepts with a warning (only `ns` stops). This is reached by any exposure with a floor, such as zero-inflated rainfall or a detection limit. The port now warns and returns the basis, as R does, and emits R's "shoving 'interior' knots" warning when only some knots tie.
- `thr` computed its default threshold with NaN dropped. `dlnm` calls `median(x)` without `na.rm`, so a missing value makes the threshold, and the whole basis, `NA`. Dropping NaN silently gave a different threshold, and different results, from the same R code. **Behaviour change:** pass `thr.value` explicitly for data with missing values.

**Other**

- `summary_figure(pred, var=[<one value>])` raised `AttributeError`. `overlay_slices` now draws no legend for a single value.
- `overlay_slices(..., var=..., xlab=...)` titled the legend with the x-axis label, though the legend indexes the curves.
- `simulate_coef` is exported at top level, as the other attribution helpers are.


## 0.4.2 (2026-08-30)

This release reproduces Gasparrini et al. (2015) *Lancet* end to end on the author's published England & Wales data (`examples/lancet_2015.py`, [docs/validation.md](docs/validation.md)): identical minimum-mortality percentiles and attributable fractions, reduced coefficients to 1e-13.

Two fixes were found while doing it:

- `to_dataframe()` now carries the index of the pandas object the basis was built from, so `data.join(cb.to_dataframe("cb"))` aligns when `data` has been filtered or grouped and no longer has a default `RangeIndex`. Previously the result was silently all-NaN: the join succeeded and the row count was unchanged, so nothing signalled the problem. This is the central pattern of a two-stage multi-location analysis. Pass `index=` to override.
- Basis names are matched with an anchored separator, so a basis called `cb` no longer also matches the columns of a second basis called `cb2`. `extract_coef_vcov`/`crosspred` raised "found 4 coefficients matching basis 'cb'" for that perfectly reasonable naming.


## 0.4.1 (2026-08-30)

This release fixes `uncertainty.qaic`, which read the log-likelihood from `results.llf`. For a quasi family statsmodels divides the log-likelihood by the estimated dispersion, so the fit term was scaled by `1/phi`. Because phi differs between candidate models the error did not cancel in a comparison: it favoured the more over-dispersed, typically under-fitted, specification and could reverse the ranking returned by `model_grid`. The Poisson log-likelihood is now evaluated at the fitted values, as in the reference R implementations (R's `logLik()` is NA for quasipoisson for the same reason). New fixtures from R (`tools/make_fixtures_qaic.R`) pin both the QAIC values and the ranking of a six-model grid. The previous test asserted the implementation against itself and so could not detect this.

## 0.4.0 (2026-08-30)

The plotting is redesigned, with journal and colour themes (`plot.set_theme`), `overlay_slices`, `summary_figure` and R-consistent labels. `dlnmpy.uncertainty` adds parametric bootstrap helpers, QAIC and `model_grid`. The theory document is extended to cr, attribution, MMT, meta-analysis and penalised fitting. A cubic regression spline basis `cr` is added (a port of mgcv, matching R to 1e-14). Penalised fitting arrives as `fit_pgam`/`fit_pglm` with REML or ML smoothing-parameter selection reproducing `mgcv::gam(paraPen=)` (validated on Poisson, quasi-Poisson, ML, cr/cr and additional lag penalties). Plot defaults are restrained and journal-style, with `plot.set_style()`. `datasets.simulate_cities` is included (0.3.2).

## 0.3.2 (2026-08-29)

This release adds `datasets.simulate_cities()`, a pure-Python multi-location simulator, and the two-stage example no longer reads a file from the test fixtures. It also adds `environment.yml` and a conda-forge recipe skeleton in `conda/`.

## 0.3.1 (2026-08-29)

Plot labels now follow R: for `contour` and `3d`, `ylab` labels the lag axis and the new `zlab` labels the outcome (colour bar / z axis). `ns` and `bs` warn, as R does, when `df` is too small and the fallback is used.

## 0.3.0 (2026-08-29)

This release adds `dlnmpy.meta`: multivariate meta-analysis and meta-regression (`mixmeta`: REML/ML/fixed, unstr/diag/id), BLUPs, Q and I², predictions, plus `stack_reduced` and `predict_reduced` for two-stage designs. Validated against the R package mixmeta 1.2.2 (`tools/make_fixtures_meta.R`). The example is `examples/two_stage.py`.

## 0.2.0 (2026-08-29)

This release adds `dlnmpy.attribution`: `attrdl` (port of Gasparrini & Leone 2014 `attrdl.R`), `findmin` (port of Tobías et al. 2017 `findmin.R`), `mmt` and `attr_table`, validated against the reference R functions (`tools/make_fixtures_attr.R`). The example is `examples/attributable_risk.py`.

## 0.1.0 (2026-08-29)

This is the first release, a port of dlnm 2.4.10 (R) covering onebasis, crossbasis (series and exposure-history input, groups), the basis functions lin, poly, strata, thr, integer, ns, bs and ps, logknots/equalknots, exphist, crosspred (lag-specific, overall and cumulative effects, centring, sub-lag prediction, matrix `at`), crossreduce (overall, lag, var), cbpen, matplotlib plots, statsmodels fitting helper reproducing R's glm behaviour on rank-deficient designs, and the chicagoNMMAPS/drug/nested datasets. It is validated against R-generated fixtures and an end-to-end side-by-side comparison (tools/side_by_side.*) covering OLS, logistic, quasi-Poisson and conditional logistic models. It also adds fit_clogit.
