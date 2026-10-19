# Penalised DLNMs

The R package supports penalised DLNMs (Gasparrini, Scheipl, Armstrong & Kenward 2017) through `ps`/`cr` bases and `cbPen`, with the fit delegated to `mgcv::gam(..., paraPen=list(cb=cbPen(cb)), method="REML")`. dlnmpy provides the same route end to end, without mgcv:

- `ps` (P-spline) and `cr` (cubic regression spline, a port of mgcv's `smooth.construct.cr.smooth.spec` and `crspl`) bases with their penalty matrices, matching R to 1e-14,
- `cbpen(cb, add_slag=...)`, the expanded and rescaled penalties for the cross-basis (matches `cbPen`), and
- `fit_pgam(formula, data, penalties={"cb": cbpen(cb)}, family=..., method="reml"|"ml")`, a penalised IRLS fitter whose smoothing parameters (and, for quasi families, the scale) maximise the Laplace-approximate REML or ML criterion of Wood (2011), with the terms of `mgcv::gam.fit3` (`fit_pglm` is the array-level version). One difference in how a term is evaluated: the log pseudo-determinant of the total penalty is taken directly from its eigenvalues, where mgcv uses Wood's (2011, Appendix B) stable similarity transform. The two agree to rounding error in ordinary use. The direct route loses accuracy when two smoothing parameters differ by a factor of about e^20 or more and the penalty with the smaller one penalises directions that the larger ones do not (an error of 3e-3 in the criterion at a log-ratio of 25 for the `ps` x `ps` cross-basis of the Chicago data, growing quickly beyond). When the penalties with large smoothing parameters are of full rank, as in the third London model of Gasparrini et al. (2017), where the ratio is about e^36, the computation is accurate to 1e-12. The outer search also evaluates the criterion from a warm-started P-IRLS, so repeated evaluations at the same smoothing parameters differ by about 5e-9, and `converged` can rest on a central-difference gradient of that slightly noisy function.

```python
cb = dl.crossbasis(temp, lag=21, argvar={"fun": "cr", "df": 8}, arglag={"fun": "ps", "df": 6})
fit = dl.fit_pgam("death ~ " + " + ".join(X.columns) + " + C(dow)", data,
                  penalties={"cb": dl.cbpen(cb)}, family="quasipoisson", method="reml")
fit.sp, fit.edf_by("cb_"), fit.scale
pred = dl.crosspred(cb, fit, cen=21, by=1, name="cb")   # works like any other model
```

## Agreement with mgcv

`tools/make_fixtures_pen.R` fits the same models with `gam` (mgcv 1.9). At fixed smoothing parameters the REML score, coefficients, covariance, scale, effective degrees of freedom and deviance agree to 1e-8 or better, which pins the criterion itself. With smoothing parameters selected (Poisson REML, quasi-Poisson REML with the scale estimated, Poisson ML, `cr`/`cr` quasi-Poisson, and a P-spline cross-basis with an additional ridge penalty on the lag coefficients), the tests in `tests/test_penalized.py` require the optimised scores to agree to 1e-4, the smoothing parameters to 1e-3 relative (1e-4 for the added ridge penalty) and the cross-basis coefficients to between 1e-6 and 5e-3 absolute, depending on the model. The residual differences come from the two optimisers stopping at slightly different points on a flat criterion. These are the fixture models. On harder problems the selected smoothing parameters can differ by much more (see the README and `docs/validation.md`), and passing mgcv's `sp=` reproduces its fit. Two details are worth knowing:

- mgcv reports, and scales the covariance by, the Fletcher (2012) corrected Pearson estimate of the scale in place of the REML scale. `fit_pgam` does the same (`fit.scale`, with the REML value in `fit.reml_scale`).
- mgcv's ML criterion uses the log determinant of the penalised Hessian projected onto the range space of the penalty (the unpenalised coefficients are profiled instead of integrated). `method="ml"` reproduces that.

The optimiser is a generic quasi-Newton on the log smoothing parameters with numerical derivatives, which is adequate for the two or three parameters of a cross-basis (a few seconds on the Chicago data) but slower than mgcv's Newton method with analytic derivatives. That is the natural place for a later optimisation.

## A practical note

A `ps` basis on a predictor with sparse tails (temperature) leaves the outer coefficients weakly penalised and can produce implausible estimates at the extremes: on Chicago, `ps(df=9)` for temperature gives RR 5.4 at 33 C, identically in R and Python. A `cr` basis for the predictor (knots at quantiles) or restricting predictions to the 1st to 99th percentiles avoids this. The 2017 paper discusses additional penalties (`add_slag`) for the lag dimension.


## Agreement with mgcv, and where it stops

The penalty matrices from `cbpen`, their ranks, and the REML/ML criterion
reproduce `mgcv` to machine precision. Given the same smoothing parameters,
the fitted coefficients agree to about 1e-10 and the criterion to 1e-9.

What differs is the *selection* of the smoothing parameters. `mgcv` uses a
Newton method on analytic derivatives of the criterion, while `fit_pgam` runs BFGS
then Nelder-Mead on numerical ones. On well-conditioned models both reach the
same optimum. On flat or near-singular ones they need not, and because the
criterion is flat there, a large difference in the smoothing parameters can sit
behind a tiny difference in the score: reproducing Gasparrini et al. (2017)
found the selected smoothing parameters differing on two of three penalised
models while the scores agreed to 1e-5 relative.

The practical advice is as follows.

- Check `.converged`. A fit that cannot produce finite estimates reports
  `converged=False` and warns. It no longer returns NaN silently.
- If you need `mgcv`'s exact answer, pass `sp=` explicitly. dlnmpy then
  reproduces its coefficients and surface to ~1e-10.
- If the fit is singular, the design is usually too rich for the data. Reduce
  `df` before reaching for a different optimiser.
