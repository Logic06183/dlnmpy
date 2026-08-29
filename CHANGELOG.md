# Changelog

## 0.3.1 (2026-08-29)

Plot labels now follow R: for `contour` and `3d`, `ylab` labels the lag axis and the new `zlab` labels the outcome (colour bar / z axis). `ns` and `bs` warn, as R does, when `df` is too small and the fallback is used.

## 0.3.0 (2026-08-29)

Adds `dlnmpy.meta`: multivariate meta-analysis and meta-regression (`mixmeta`: REML/ML/fixed, unstr/diag/id), BLUPs, Q and I², predictions, plus `stack_reduced` and `predict_reduced` for two-stage designs. Validated against the R package mixmeta 1.2.2 (`tools/make_fixtures_meta.R`). Example `examples/two_stage.py`.

## 0.2.0 (2026-08-29)

Adds `dlnmpy.attribution`: `attrdl` (port of Gasparrini & Leone 2014 `attrdl.R`), `findmin` (port of Tobías et al. 2017 `findmin.R`), `mmt` and `attr_table`, validated against the reference R functions (`tools/make_fixtures_attr.R`). Example `examples/attributable_risk.py`.

## 0.1.0 (2026-08-29)

First release. Port of dlnm 2.4.10 (R) covering onebasis, crossbasis (series and exposure-history input, groups), the basis functions lin, poly, strata, thr, integer, ns, bs and ps, logknots/equalknots, exphist, crosspred (lag-specific, overall and cumulative effects, centring, sub-lag prediction, matrix `at`), crossreduce (overall, lag, var), cbpen, matplotlib plots, statsmodels fitting helper reproducing R's glm behaviour on rank-deficient designs, and the chicagoNMMAPS/drug/nested datasets. Validated against R-generated fixtures and an end-to-end side-by-side comparison (tools/side_by_side.*) covering OLS, logistic, quasi-Poisson and conditional logistic models. Adds fit_clogit.
