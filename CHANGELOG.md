# Changelog

## 0.6.1 (2026-09-02)

**Fixed**

- `dlnm(..., time=)` given a date column killed the interpreter. `datetime64` cast to float is nanoseconds since the epoch, so a fourteen-year daily series was read as 1.2e12 years, the seasonal spline was asked for 8.5e12 degrees of freedom, and `np.linspace` tried to allocate about 68,000 GB; the process was killed by the OS with no Python traceback. Date columns (and `datetime.date` objects) are now converted to day offsets from the first observation, and `time="date"` gives numerically identical results to the equivalent day index. Found by running the published 0.6.0 wheel through a full Chicago analysis in a notebook.
- The seasonal spline now raises `ValueError` when it would need at least as many degrees of freedom as there are observations, which catches any unit error of this kind rather than exhausting memory.
- `__version__` was hardcoded in `__init__.py` and had drifted from `pyproject.toml`; it is now read from the installed package metadata, with a test to keep them tied.
- `basis.py` had a conditional whose branches were identical (`thr_value if side == "d" else thr_value`); simplified, with no change in behaviour.

**Changed**

- `ruff` is pinned in the `dev` extra and the lint rule set is selected explicitly. Ruff's default selection widens between releases, so CI could turn red on a linter upgrade rather than on a change to this code.

## 0.6.0 (2026-09-02)

A second audit against R (`dlnm` 2.4.7 in a fresh R 4.3 install, `mixmeta` 1.2.2), this time on 64 edge cases the fixtures did not cover: negative and sub-period lag ranges, exposure-history matrices with `cumul`, non-integer `bylag` and `crossreduce(type="lag", value=2.5)`, `group`, every basis function with explicit knots or thresholds, the default and `by` prediction grids, `mkcen` corner cases, `logknots`/`equalknots`, `exphist` with `times` and `fill`, and `attrdl` in the forward, reduced-coefficient and daily forms. All 64 agree to 1e-6 or better (most to 1e-12); `mixmeta` REML, BLUPs, predictions, Q/I² and univariate Q agree to 1e-8. The numerical core needed no change. The defects below were all in the Python-only layers.

**Fixed**

- `PenalizedGLMResults.edf_by("cb")` summed the effective df of every column whose name *started with* `cb`, so with a second basis called `cb_o3` it silently reported the two together (38.1 where the temperature basis had 9.4). `fit_pgam` used the same match to place the penalties and failed with a shape error in that case. Both now match a basis name the way `crosspred` does (name plus a `v1_l1`/`b1` label).
- `qaic()` raised `AttributeError` on a penalised fit. It now accepts them, using the total effective degrees of freedom as the parameter count (Gasparrini et al. 2017); `PenalizedGLMResults` gained a `fittedvalues` alias.
- `fit_pgam(..., offset=)` with a full-length offset mis-aligned it after rows with missing values were dropped; it now follows the same rows.
- `overlay_slices(labels=...)` with too few labels raised `IndexError` from inside matplotlib; it now says what is wrong.
- `qtest()` p-values were computed as `1 - cdf` and rounded to exactly 0 for large Q; the survival function is used.
- `bootstrap_ci` is exported.
- `mkat` passes `min.n = nobs/2` to `pretty()` as R does (no case where it changes the grid was found in 20,000 random ranges, so this is fidelity rather than a fix).

**One R defect not reproduced**

- `ps(df=4)` (degree 3) builds in R and then fails inside `crosspred()` when the basis is rebuilt from its knots: a lag basis you can never predict from. `dlnmpy.ps` refuses it at construction with the minimum df in the message.

**Beyond R: the workflow**

- `dl.dlnm(data, outcome, exposure, lag, argvar, arglag, time, df_per_year, dow, controls, family, offset, group)` fits the whole model in one call and returns a `DLNM` bound to its data: `.mmt()`, `.rr_at(percentiles)`, `.attributable()`, `.predict(percentiles=...)`, `.reduce()`, `.qaic()`, `.summary()`, `.figure()`, `.plot()`, `.summary_figure()`. Defaults are the usual specification (natural spline with knots at the 10th, 75th and 90th percentiles; lag spline with log-spaced knots; 7 df per year of time). `ps`/`cr` bases route to the penalised fitter. Tested to agree with the explicit pipeline to 1e-12.
- `CrossPred.to_dataframe()` returns the long ("tidy") table of lag-specific (or cumulative) effects.
- `percentile_knots(x, [10, 75, 90])` and `percentile_of(x, values)`.

**Beyond R: the figures**

- `plot_overall_risk`: the exposure-response figure of a temperature-mortality paper (cold and heat sides in two hues, minimum-risk value with its interval, percentile axis along the top, exposure histogram underneath). `plot_attributable`: attributable fractions or numbers by component with intervals.
- Styling is applied per figure through `plot.style()` (an `rc_context`), so plots look the same everywhere and the caller's rcParams are untouched. The 3D surface is coloured by value with a fine translucent mesh instead of the black wireframe, and no longer clips its z label; the contour's null-effect line is dashed and the frame closed.
- `examples/gallery.py` renders every figure in both themes.


## 0.5.1 (2026-08-30)

Silences spurious floating-point `RuntimeWarning`s. numpy 2.0.x raises `divide by zero`, `overflow`, `invalid value` and `underflow` warnings from `matmul` on arrays that legitimately carry `NaN` -- a cross-basis has `NaN` in its first `lag` rows by construction -- and from the log-determinant the penalised REML optimiser evaluates while probing extreme smoothing parameters. The values were never affected, but a routine analysis printed a wall of warnings for anyone on Python 3.9 or 3.10, where numpy resolves to 2.0.x. A full Chicago analysis went from about 1400 warnings to 3, and the test suite from 52 to 1; the three that remain come from statsmodels' own IRLS calling `numpy.linalg.pinv`, not from dlnmpy. Results are bit-identical: RR at 33 °C is 0.922319 and at -20 °C 1.329074 on numpy 2.0.2 and 2.5.2 alike.


## 0.5.0 (2026-08-30)

Also reproduces Gasparrini & Armstrong (2013) *BMC Med Res Methodol* 13:1, the methodological paper behind the two-stage design (`examples/bmcmrm_2013.py`): minimum-mortality temperature 17.1 °C, pooled RR 1.101 (1.079-1.124) at 22 °C and 1.309 (1.245-1.376) at 0 °C, I² of 63.8/16.4/63.5%, and all three QAIC totals; 67 of 68 intermediates agree with R to 1e-5..1e-15. That paper selects its lag specification by QAIC, and the pre-0.5.0 implementation would have reversed its published conclusion.

Bug-fix release from a systematic audit against R (`dlnm` 2.4.10, `splines`, `mgcv`, `mixmeta`, `survival`), covering the basis functions, the prediction machinery, model fitting, meta-analysis, attribution and plotting. Everything below was verified against R and has a regression test that fails without the fix.

**Silent wrong answers**

- `fit_glm` ignored `offset`, `exposure`, `freq_weights` and `var_weights`. They are arguments of the model, but were forwarded to `GLM.fit()`, whose own `**kwargs` discards unknown names, so the caller silently got the un-offset fit. A Poisson model with a population offset returned coefficients `[0.850, 0.426]` where R's `offset()` gives `[0.497, 0.357]`. They are now explicit parameters, and survive both the NaN drop and the aliased-column refit.
- `crossbasis(argvar={..., "cen": v})` discarded the centring value. R keeps it as `attr(cb, "argvar")$cen` and uses it as the default reference for every later `crosspred`/`crossreduce`, so predictions were silently centred on the automatic value instead: a relative risk of 1.343 where R gives 1.166, 15% out.
- `strata` returned a rank-deficient basis when quantile breaks tied (a spike of zeros, say): an empty stratum, so an identically-zero column and a singular design, with nothing to signal it. R's `cut()` raises `'breaks' are not unique`, and so does this now. Values outside the outer edges are `NaN`, as `cut()` gives `NA`, rather than a row of zeros.
- `plot(..., "slices", lag=L, cumul=True)` read the cumulative outcomes off `pred.predlag` instead of the integer lags they are defined on, so with `bylag != 1` it plotted a different lag: at `bylag=0.2`, `lag=2` drew the curve for lag 10, and `lag >= 4` raised `IndexError`. R sets `bylag` to 1 under `cumul` and indexes by name.
- `get_link` reported a complementary log-log link as `"log"`, because `"log"` is a substring of `"cloglog"` and was tested first, so `crosspred` exponentiated the fit and reported relative risks for a model that has none. Same for `loglog`.
- `cen=np.True_` / `np.False_` were coerced to 1.0 / 0.0: `np.bool_` is not a subclass of `bool`, so a numpy boolean fell through both branches of `mkcen`.

**Divergences from R**

- `bs` refused data whose interior knots all coincide with a boundary knot, which `splines::bs` accepts with a warning (only `ns` stops). This is reached by any exposure with a floor -- zero-inflated rainfall, a detection limit. The port now warns and returns the basis, as R does, and emits R's "shoving 'interior' knots" warning when only some knots tie.
- `thr` computed its default threshold with NaN dropped; `dlnm` calls `median(x)` without `na.rm`, so a missing value makes the threshold, and the whole basis, `NA`. Dropping NaN silently gave a different threshold, and different results, from the same R code. **Behaviour change:** pass `thr.value` explicitly for data with missing values.

**Other**

- `summary_figure(pred, var=[<one value>])` raised `AttributeError`; `overlay_slices` draws no legend for a single value.
- `overlay_slices(..., var=..., xlab=...)` titled the legend with the x-axis label, though the legend indexes the curves.
- `simulate_coef` is exported at top level, as the other attribution helpers are.


## 0.4.2 (2026-08-30)

Reproduces Gasparrini et al. (2015) *Lancet* end to end on the author's published England & Wales data (`examples/lancet_2015.py`, [docs/validation.md](docs/validation.md)): identical minimum-mortality percentiles and attributable fractions, reduced coefficients to 1e-13.

Two fixes found while doing it:

- `to_dataframe()` now carries the index of the pandas object the basis was built from, so `data.join(cb.to_dataframe("cb"))` aligns when `data` has been filtered or grouped and no longer has a default `RangeIndex`. Previously the result was silently all-NaN — the join succeeded and the row count was unchanged, so nothing signalled the problem. This is the central pattern of a two-stage multi-location analysis. Pass `index=` to override.
- Basis names are matched with an anchored separator, so a basis called `cb` no longer also matches the columns of a second basis called `cb2`. `extract_coef_vcov`/`crosspred` raised "found 4 coefficients matching basis 'cb'" for that perfectly reasonable naming.


## 0.4.1 (2026-08-30)

Fixes `uncertainty.qaic`, which read the log-likelihood from `results.llf`. For a quasi family statsmodels divides the log-likelihood by the estimated dispersion, so the fit term was scaled by `1/phi`. Because phi differs between candidate models the error did not cancel in a comparison: it favoured the more over-dispersed, typically under-fitted, specification and could reverse the ranking returned by `model_grid`. The Poisson log-likelihood is now evaluated at the fitted values, as in the reference R implementations (R's `logLik()` is NA for quasipoisson for the same reason). New fixtures from R (`tools/make_fixtures_qaic.R`) pin both the QAIC values and the ranking of a six-model grid; the previous test asserted the implementation against itself and so could not detect this.

## 0.4.0 (2026-08-30)

Plotting redesign: journal and colour themes (`plot.set_theme`), `overlay_slices`, `summary_figure`, R-consistent labels. `dlnmpy.uncertainty`: parametric bootstrap helpers, QAIC and `model_grid`. Theory document extended to cr, attribution, MMT, meta-analysis and penalised fitting. Cubic regression spline basis `cr` (port of mgcv, matches R to 1e-14). Penalised fitting: `fit_pgam`/`fit_pglm` with REML or ML smoothing-parameter selection reproducing `mgcv::gam(paraPen=)` (validated on Poisson, quasi-Poisson, ML, cr/cr and additional lag penalties). Restrained journal-style plot defaults and `plot.set_style()`. `datasets.simulate_cities` (0.3.2).

## 0.3.2 (2026-08-29)

`datasets.simulate_cities()` (pure-Python multi-location simulator); the two-stage example no longer reads a file from the test fixtures. `environment.yml` and a conda-forge recipe skeleton in `conda/`.

## 0.3.1 (2026-08-29)

Plot labels now follow R: for `contour` and `3d`, `ylab` labels the lag axis and the new `zlab` labels the outcome (colour bar / z axis). `ns` and `bs` warn, as R does, when `df` is too small and the fallback is used.

## 0.3.0 (2026-08-29)

Adds `dlnmpy.meta`: multivariate meta-analysis and meta-regression (`mixmeta`: REML/ML/fixed, unstr/diag/id), BLUPs, Q and I², predictions, plus `stack_reduced` and `predict_reduced` for two-stage designs. Validated against the R package mixmeta 1.2.2 (`tools/make_fixtures_meta.R`). Example `examples/two_stage.py`.

## 0.2.0 (2026-08-29)

Adds `dlnmpy.attribution`: `attrdl` (port of Gasparrini & Leone 2014 `attrdl.R`), `findmin` (port of Tobías et al. 2017 `findmin.R`), `mmt` and `attr_table`, validated against the reference R functions (`tools/make_fixtures_attr.R`). Example `examples/attributable_risk.py`.

## 0.1.0 (2026-08-29)

First release. Port of dlnm 2.4.10 (R) covering onebasis, crossbasis (series and exposure-history input, groups), the basis functions lin, poly, strata, thr, integer, ns, bs and ps, logknots/equalknots, exphist, crosspred (lag-specific, overall and cumulative effects, centring, sub-lag prediction, matrix `at`), crossreduce (overall, lag, var), cbpen, matplotlib plots, statsmodels fitting helper reproducing R's glm behaviour on rank-deficient designs, and the chicagoNMMAPS/drug/nested datasets. Validated against R-generated fixtures and an end-to-end side-by-side comparison (tools/side_by_side.*) covering OLS, logistic, quasi-Poisson and conditional logistic models. Adds fit_clogit.
