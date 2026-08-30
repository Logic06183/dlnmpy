# Changelog

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
