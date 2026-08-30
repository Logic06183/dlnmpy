# Penalised DLNMs

The R package supports penalised DLNMs (Gasparrini, Scheipl, Armstrong & Kenward 2017) through `ps`/`cr` bases and `cbPen`, with the fit delegated to `mgcv::gam(..., paraPen=list(cb=cbPen(cb)), method="REML")`. dlnmpy provides the same route end to end, without mgcv:

- `ps` (P-spline) and `cr` (cubic regression spline, a port of mgcv's `smooth.construct.cr.smooth.spec` and `crspl`) bases with their penalty matrices, matching R to 1e-14;
- `cbpen(cb, add_slag=...)`, the expanded and rescaled penalties for the cross-basis (matches `cbPen`);
- `fit_pgam(formula, data, penalties={"cb": cbpen(cb)}, family=..., method="reml"|"ml")`, a penalised IRLS fitter whose smoothing parameters (and, for quasi families, the scale) maximise the Laplace-approximate REML or ML criterion of Wood (2011), implemented exactly as in `mgcv::gam.fit3` (`fit_pglm` is the array-level version).

```python
cb = dl.crossbasis(temp, lag=21, argvar={"fun": "cr", "df": 8}, arglag={"fun": "ps", "df": 6})
fit = dl.fit_pgam("death ~ " + " + ".join(X.columns) + " + C(dow)", data,
                  penalties={"cb": dl.cbpen(cb)}, family="quasipoisson", method="reml")
fit.sp, fit.edf_by("cb_"), fit.scale
pred = dl.crosspred(cb, fit, cen=21, by=1, name="cb")   # works like any other model
```

## Agreement with mgcv

`tools/make_fixtures_pen.R` fits the same models with `gam` (mgcv 1.9). At fixed smoothing parameters the REML score, coefficients, covariance, scale, effective degrees of freedom and deviance agree to 1e-8 or better, which pins the criterion itself. With smoothing parameters selected (Poisson REML, quasi-Poisson REML with the scale estimated, Poisson ML, `cr`/`cr` quasi-Poisson, and a P-spline cross-basis with an additional ridge penalty on the lag coefficients), the optimised scores agree to 1e-5, the smoothing parameters to about 1e-4 relative and the cross-basis coefficients to 1e-4 or better; the residual differences come from the two optimisers stopping at slightly different points on a flat criterion. Two details worth knowing:

- mgcv reports, and scales the covariance by, the Fletcher (2012) corrected Pearson estimate of the scale rather than the REML scale; `fit_pgam` does the same (`fit.scale`; the REML value is `fit.reml_scale`).
- mgcv's ML criterion uses the log determinant of the penalised Hessian projected onto the range space of the penalty (the unpenalised coefficients are profiled, not integrated); `method="ml"` reproduces that.

The optimiser is a generic quasi-Newton on the log smoothing parameters with numerical derivatives, which is adequate for the two or three parameters of a cross-basis (a few seconds on the Chicago data) but slower than mgcv's Newton method with analytic derivatives; that is the natural place for a later optimisation.

## A practical note

A `ps` basis on a predictor with sparse tails (temperature) leaves the outer coefficients weakly penalised and can produce implausible estimates at the extremes: on Chicago, `ps(df=9)` for temperature gives RR 5.4 at 33 C, identically in R and Python. A `cr` basis for the predictor (knots at quantiles) or restricting predictions to the 1st to 99th percentiles avoids this; the 2017 paper discusses additional penalties (`add_slag`) for the lag dimension.
