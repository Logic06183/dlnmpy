# Penalised DLNMs: status

The R package supports penalised DLNMs (Gasparrini et al. 2017) in two ways: `ps`/`cr` bases with `cbPen` to build penalty matrices for use in `mgcv::gam(..., paraPen=)`, and a purpose-built smooth constructor (`smooth.construct.cb.smooth.spec`) for `s(x, lag, bs="cb")`. Both rely on mgcv for the penalised fit (smoothing parameter selection by REML or GCV).

What `dlnmpy` provides now:

- the `ps` basis with its difference penalty matrix (matches R);
- `cbpen(cb, add_slag=...)`, giving the expanded and rescaled penalty matrices, their ranks and the smoothing parameter placeholders (matches R).

What is missing:

- the `cr` basis (mgcv's cubic regression spline construction, `smooth.construct.cr.smooth.spec`), which needs the natural cubic spline parameterisation with knots as parameters; not difficult, but it must be ported carefully for equivalence;
- a fitter. There is no mgcv in Python. The options are a penalised IRLS with the penalties from `cbpen` and smoothing parameters chosen by REML (Wood 2011), a mixed-model reparameterisation, or a Bayesian fit. Once the coefficients and their covariance are available, `crosspred` and `crossreduce` apply unchanged: the prediction algebra does not depend on how the coefficients were estimated.

An interim route is to fit in R (`gam` with `paraPen`) and bring `coef`/`vcov` into Python for prediction, or to fix smoothing parameters and fit a penalised GLM by hand: with penalty `S = sum_k sp_k S_k` the penalised IRLS step solves `(X'WX + S) beta = X'Wz`, and the Bayesian covariance is `(X'WX + S)^{-1} phi`.
