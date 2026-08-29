# dlnmpy

Distributed lag non-linear models (DLNMs) in Python.

`dlnmpy` is a from-first-principles port of the R package [`dlnm`](https://github.com/gasparrini/dlnm) by Antonio Gasparrini and Ben Armstrong. It reproduces the R package's cross-basis construction, prediction and reduction algorithms to machine precision, and it is built so that the same reference tests can be used to validate ports to other languages (Rust, Julia, JavaScript).

The package is validated against R on the `chicagoNMMAPS` examples from the `dlnm` vignettes: basis matrices, cross-basis matrices, predictions, confidence intervals and reduced coefficients agree with R to around 1e-12, and the statsmodels fitting path reproduces R's `glm()` coefficients to 1e-14 (see [Validation](#validation)).

## Why

DLNMs are the standard method for exposure-lag-response associations in environmental epidemiology, above all for temperature and mortality. The methodology lives in one R package. That is a problem for teams working in Python (Databricks, scikit-learn pipelines, deep learning workflows) and for anyone who wants to run the models in a compiled language at scale. This repository sets out to make the method itself portable: the algorithms are written down as a language-neutral specification, the numerical behaviour is pinned by fixtures generated from R, and the Python implementation is the first port.

## Installation

```bash
pip install git+https://github.com/Logic06183/dlnmpy.git
# or, for development
git clone https://github.com/Logic06183/dlnmpy.git && cd dlnmpy
pip install -e ".[dev]"
```

Requires Python 3.9+, numpy, scipy and pandas. `statsmodels` is needed to fit models through `fit_glm` and `matplotlib` for plotting; both are optional.

## Quick start

The workflow is the same as in R: build a cross-basis, include it in a regression model, predict.

```python
import numpy as np
import dlnmpy as dl

chicago = dl.datasets.chicago_nmmaps()

# cross-basis for temperature: natural splines in both dimensions, lag 0-21
cb = dl.crossbasis(chicago.temp, lag=21,
                   argvar={"fun": "ns", "df": 4},
                   arglag={"fun": "ns", "df": 4})

# design matrix: cross-basis + seasonality/trend spline + day of week
ns_time = dl.onebasis(chicago.time, "ns", df=7 * 14)
X = dl.design_matrix(chicago, ("cb", cb), ("ns_time", ns_time), intercept=False)
data = chicago.join(X)
model = dl.fit_glm("death ~ " + " + ".join(X.columns) + " + C(dow)", data,
                   family="quasipoisson")

# predictions centred at 21C, for every integer temperature
pred = dl.crosspred(cb, model, cen=21, by=1, name="cb")
pred.plot("overall", xlab="Temperature", ylab="RR")
pred.plot("slices", var=[-20, 33], lag=[0, 5])
pred.overall()          # DataFrame with var, fit, low, high (RR scale)
pred.slice_var(33)      # lag-response curve at 33C

# reduce to one-dimensional summaries (for meta-analysis, say)
red = dl.crossreduce(cb, model, cen=21, name="cb")
red.coef, red.vcov
```

Any fitting library can be used. `crosspred` and `crossreduce` only need the coefficients and covariance matrix of the cross-basis columns, so you can also pass them directly:

```python
pred = dl.crosspred(cb, coef=beta, vcov=V, model_link="log", cen=21, by=1)
```

The `examples/chicago_timeseries.py` script reproduces the four examples of the R vignette *Distributed lag linear and non-linear models for time series data* and writes the figures to `examples/figures/`.

## What is implemented

| R (`dlnm`) | Python (`dlnmpy`) | Notes |
|---|---|---|
| `onebasis(x, fun, ...)` | `onebasis(x, fun, **args)` | returns `OneBasis` |
| `crossbasis(x, lag, argvar, arglag, group)` | `crossbasis(x, lag, argvar, arglag, group)` | returns `CrossBasis`; `x` may be a series or a matrix of exposure histories |
| `lin`, `poly`, `strata`, `thr`, `integer` | same names in `dlnmpy.basis` | |
| `splines::ns`, `splines::bs` | `ns`, `bs` | literal ports, including extrapolation beyond the boundary knots |
| `ps` | `ps` | P-spline basis with difference penalty |
| `cr` | not yet | needs mgcv's cubic regression spline construction |
| `logknots`, `equalknots` | `logknots`, `equalknots` | |
| `exphist` | `exphist` | |
| `crosspred(...)` | `crosspred(...)` | returns `CrossPred`; `from`/`to` are `from_`/`to` |
| `crossreduce(...)` | `crossreduce(...)` | returns `CrossReduce` |
| `cbPen` | `cbpen` | penalty matrices; no penalised fitter yet (see `docs/penalized.md`) |
| `plot.crosspred` (`overall`, `slices`, `contour`, `3d`) | `CrossPred.plot(...)` | matplotlib |
| `plot.crossreduce` | `CrossReduce.plot(...)` | |
| `summary` methods | `.summary()` | |
| `chicagoNMMAPS`, `drug`, `nested` | `dlnmpy.datasets` | |

R argument spellings are accepted where they differ (`thr.value`/`thr`, `Boundary.knots`), so R code translates almost line by line. See `docs/r-to-python.md` for the full mapping and the differences that remain.

The cross-basis columns are named `v{i}.l{j}` as in R (or `{name}_v{i}_l{j}` in the DataFrame returned by `to_dataframe`, which is what `crosspred(..., name=...)` uses to find the right coefficients in a fitted model).

## Validation

`tools/make_fixtures.R` runs the R package on the vignette examples and writes every intermediate object to `tests/fixtures/` as CSV and JSON: 32 basis-function specifications (with and without values outside the fitting range), 9 cross-bases, 12 prediction objects, 7 reductions, the penalty matrices, and R's `pretty()`, `quantile()` and knot helpers. The Python test-suite compares against these numbers with absolute tolerances of 1e-10 to 1e-12.

```
pytest            # 63 tests
```

Two things worth knowing about equivalence with R:

- R's `glm` handles rank-deficient designs by dropping aliased columns (reported as `NA`); statsmodels uses a pseudo-inverse and a numerical rank. `fit_glm` reproduces R's behaviour (drops aliased columns using the same criterion as LINPACK's `dqrdc2`, solves the IRLS step by QR, and computes the quasi-likelihood dispersion with R's residual degrees of freedom). This is why the example-1 coefficients agree to 1e-14 rather than 1e-6.
- The reference level of categorical terms in patsy formulas is the first level in alphabetical order, as in R for character vectors; R factors with explicit levels may differ. This does not affect the cross-basis coefficients.

## Repository layout

```
src/dlnmpy/
  _rcompat.py    R's pretty(), quantile(type=7), median(), seq()   (ported from R source)
  _splines.py    splineDesign(), ns(), bs()                        (ported from R's splines)
  basis.py       lin, poly, strata, thr, integer, ns, bs, ps
  lag.py         mklag, seqlag, lag_matrix (tsModel::Lag), exphist
  knots.py       logknots, equalknots
  core.py        OneBasis, CrossBasis, onebasis(), crossbasis()
  predict.py     CrossPred, CrossReduce, crosspred(), crossreduce(), mkat, mkcen
  model.py       coefficient/vcov extraction, fit_glm (statsmodels), design_matrix
  penalty.py     cbpen
  plot.py        matplotlib plots
  datasets.py    chicagoNMMAPS, drug, nested
docs/
  theory.md          the mathematics, written as a language-neutral specification
  r-to-python.md     API mapping and behavioural differences
  porting-guide.md   how to port to another language and validate with the fixtures
  penalized.md       status of penalised DLNMs
tests/               pytest suite + R-generated fixtures
tools/               make_fixtures.R
examples/            vignette reproduction
```

## Roadmap

1. Cubic regression splines (`cr`) and a penalised fitter for `ps`/`cr` cross-bases (the R package delegates to `mgcv::gam`; the Python route is a REML or GCV fit with the penalty matrices from `cbpen`).
2. Attributable risk (`attrdl`, from Gasparrini & Leone 2014), not part of `dlnm` itself but used in nearly every application.
3. Multivariate meta-analysis of reduced coefficients (the `mixmeta` step of two-stage designs).
4. A Rust core with Python bindings, validated with the same fixtures.

## Prior work

Two earlier Python efforts exist: [aedessler/pydlnm](https://github.com/aedessler/pydlnm) (with meta-analysis and BLUPs, validated on multi-city data) and [`crossbasis`](https://pypi.org/project/crossbasis/) on PyPI. This repository differs in aiming for a literal, attribute-for-attribute port of the R algorithms (including R's own `pretty`, quantile and spline routines so defaults match), in shipping the R-generated fixtures as a language-neutral test oracle, and in documenting the method as a specification for further ports.

## References

- Gasparrini A, Armstrong B, Kenward MG. Distributed lag non-linear models. *Statistics in Medicine* 2010; 29(21):2224-2234.
- Gasparrini A. Distributed lag linear and non-linear models in R: the package dlnm. *Journal of Statistical Software* 2011; 43(8):1-20.
- Gasparrini A, Armstrong B. Reducing and meta-analysing estimates from distributed lag non-linear models. *BMC Medical Research Methodology* 2013; 13:1.
- Gasparrini A. Modeling exposure-lag-response associations with distributed lag non-linear models. *Statistics in Medicine* 2014; 33(5):881-899.
- Gasparrini A, Scheipl F, Armstrong B, Kenward MG. A penalized framework for distributed lag non-linear models. *Biometrics* 2017; 73(3):938-948.

## Licence

GPL-2.0-or-later. This is a derived work of the R package `dlnm` (GPL >= 2) and of R's `splines` and `base` packages (GPL >= 2), whose algorithms it reproduces. The `chicagoNMMAPS` data are redistributed from the `dlnm` package.
