# Changelog

## 0.2.0 (2026-08-29)

Adds `dlnmpy.attribution`: `attrdl` (port of Gasparrini & Leone 2014 `attrdl.R`), `findmin` (port of Tobías et al. 2017 `findmin.R`), `mmt` and `attr_table`, validated against the reference R functions (`tools/make_fixtures_attr.R`). Example `examples/attributable_risk.py`.

## 0.1.0 (2026-08-29)

First release. Port of dlnm 2.4.10 (R) covering onebasis, crossbasis (series and exposure-history input, groups), the basis functions lin, poly, strata, thr, integer, ns, bs and ps, logknots/equalknots, exphist, crosspred (lag-specific, overall and cumulative effects, centring, sub-lag prediction, matrix `at`), crossreduce (overall, lag, var), cbpen, matplotlib plots, statsmodels fitting helper reproducing R's glm behaviour on rank-deficient designs, and the chicagoNMMAPS/drug/nested datasets. Validated against R-generated fixtures and an end-to-end side-by-side comparison (tools/side_by_side.*) covering OLS, logistic, quasi-Poisson and conditional logistic models. Adds fit_clogit.
