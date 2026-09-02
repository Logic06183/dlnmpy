# dlnmpy

[![tests](https://github.com/Logic06183/dlnmpy/actions/workflows/ci.yml/badge.svg)](https://github.com/Logic06183/dlnmpy/actions/workflows/ci.yml)
![python](https://img.shields.io/badge/python-3.9%20%7C%203.11%20%7C%203.12-blue)
![licence](https://img.shields.io/badge/licence-GPL--2.0--or--later-green)

Distributed lag non-linear models (DLNMs) in Python. A port of the R package [`dlnm`](https://github.com/gasparrini/dlnm) by Antonio Gasparrini and Ben Armstrong, checked number for number against R, plus the parts of a temperature-mortality analysis that R leaves to loose scripts: the minimum mortality temperature, attributable fractions, two-stage meta-analysis and the figures.

Not affiliated with the authors of `dlnm`. The method is theirs; this repository makes it usable from Python.

## Does it give the same answer as R?

Yes, to the precision below. Every number here comes from the test-suite or from a script in `tools/` or `examples/` that you can run yourself.

| What | Agreement with R | Where |
|---|---|---|
| Basis and cross-basis matrices (32 basis specs, 9 cross-bases) | 1e-12 | `tests/`, fixtures from `tools/make_fixtures.R` |
| `crosspred` / `crossreduce` fits, standard errors, reduced coefficients | 1e-12 | `tests/` |
| `glm()` coefficients through statsmodels (`fit_glm`) | 1e-14 | `tests/test_model.py` |
| Five complete analyses end to end, 102 quantities (GLM, OLS, logistic, conditional logistic) | better than 1e-8, most 1e-12 | `tools/side_by_side.py` |
| 64 edge cases (negative lags, sub-periods, exposure histories, `bylag`, `group`, explicit knots, `attrdl` variants) | 1e-6 or better, most 1e-12 | audit for 0.6.0, see `CHANGELOG.md` |
| `attrdl.R`, `findmin.R` (attributable risk, MMT) | 1e-8 | `tests/test_attribution.py` |
| `mixmeta` (REML, BLUPs, predictions, Q, I²) | 1e-5 or better | `tests/test_meta.py` |
| `mgcv::gam` penalised DLNMs (scores, smoothing parameters, coefficients) | 1e-5, 1e-4, 1e-4 | `tests/test_penalized.py` |
| Gasparrini et al. 2015 *Lancet*, England and Wales, 10 regions | identical MMT percentiles; AF to 4e-5 points | `examples/lancet_2015.py` |
| Gasparrini and Armstrong 2013 *BMC MRM* | 67 of 68 intermediates to 1e-5..1e-15 | `examples/bmcmrm_2013.py` |

If you find a number that differs from R, open an issue with the R code and the Python code side by side. That is the most useful contribution this project can get.

## Who it is for

Time-series studies of temperature, air pollution or any lagged exposure against daily counts (quasi-Poisson), case-crossover and matched case-control designs (conditional logistic), cohort studies with exposure histories (`exphist`), and multi-location studies pooled with multivariate meta-analysis. Anything you would reach for `dlnm` in R to do, and the surrounding workflow: MMT with an interval, attributable numbers by cold, heat, extreme and moderate, the standard figures.

## Why

DLNMs are the standard method for exposure-lag-response associations in environmental epidemiology, above all for temperature and mortality. The methodology lives in one R package. That is a problem for teams working in Python (Databricks, scikit-learn pipelines, deep learning workflows) and for anyone who wants to run the models in a compiled language at scale. This repository sets out to make the method itself portable: the algorithms are written down as a language-neutral specification (`docs/theory.md`), the numerical behaviour is pinned by fixtures generated from R, and the Python implementation is the first port.

## Installation

```bash
pip install git+https://github.com/Logic06183/dlnmpy.git
# or, for development
git clone https://github.com/Logic06183/dlnmpy.git && cd dlnmpy
pip install -e ".[dev]"
```

With conda, create the environment from the file in the repo (conda-forge packages plus dlnmpy from GitHub):

```bash
conda env create -f environment.yml && conda activate dlnmpy
```

or add the pip line to an existing environment. A conda-forge recipe is in `conda/`; it can be submitted once the package is on PyPI.

Requires Python 3.9+, numpy, scipy and pandas. `statsmodels` is needed to fit models through `fit_glm` and `matplotlib` for plotting; both are optional. The example datasets (`chicagoNMMAPS`, `drug`, `nested`) ship inside the package, so nothing from R is needed at any point; `dlnmpy.datasets.simulate_cities()` generates multi-location data for the two-stage examples.

## Quick start

The whole temperature-mortality workflow is one call. `dlnm()` builds the cross-basis, adds the seasonal spline and day of week, fits the quasi-Poisson model, and the result knows its own data, so the minimum mortality temperature, the relative risks at chosen percentiles, the attributable fractions and the figure are each a method:

```python
import dlnmpy as dl

chicago = dl.datasets.chicago_nmmaps()

fit = dl.dlnm(chicago, outcome="death", exposure="temp", lag=21,
              argvar={"fun": "bs", "degree": 2, "knots": dl.percentile_knots(chicago.temp, [10, 75, 90])},
              arglag={"fun": "ns", "knots": dl.logknots(21, 3)},
              time="time", df_per_year=7, dow="dow")

fit.mmt()                    # MMTResult(mmt=25.50, 95% CI (25.20, 25.80), percentile=94.2)
fit.rr_at([1, 99])           # RR (95% CI) at the 1st and 99th percentiles, relative to the MMT
fit.attributable()           # total, cold, heat, extreme and moderate: AN and AF with intervals
fit.figure(xlab="Temperature (°C)")   # the curve, the MMT, a percentile axis, the exposure distribution
fit.predict(percentiles=[1, 99]).to_dataframe()   # tidy long table: var, lag, fit, se, low, high
```

![overall risk](examples/figures/overall_risk_colour.png)

Underneath, the workflow is the same as in R and every step is available on its own: build a cross-basis, include it in a regression model, predict.

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
| `ps`, `cr` | `ps`, `cr` | penalised bases with their penalty matrices (`cr` ports mgcv's construction) |
| `logknots`, `equalknots` | `logknots`, `equalknots` | |
| `exphist` | `exphist` | |
| `crosspred(...)` | `crosspred(...)` | returns `CrossPred`; `from`/`to` are `from_`/`to` |
| `crossreduce(...)` | `crossreduce(...)` | returns `CrossReduce` |
| `cbPen` + `mgcv::gam(paraPen=)` | `cbpen` + `fit_pgam` | penalised DLNMs with REML/ML smoothing, validated against mgcv (`docs/penalized.md`) |
| `plot.crosspred` (`overall`, `slices`, `contour`, `3d`) | `CrossPred.plot(...)` | matplotlib |
| `plot.crossreduce` | `CrossReduce.plot(...)` | |
| `summary` methods | `.summary()` | |
| `chicagoNMMAPS`, `drug`, `nested` | `dlnmpy.datasets` | |

R argument spellings are accepted where they differ (`thr.value`/`thr`, `Boundary.knots`), so R code translates almost line by line. See `docs/r-to-python.md` for the full mapping and the differences that remain.

The cross-basis columns are named `v{i}.l{j}` as in R (or `{name}_v{i}_l{j}` in the DataFrame returned by `to_dataframe`, which is what `crosspred(..., name=...)` uses to find the right coefficients in a fitted model).

## Figures

Plots are designed for print: neutral sans-serif, thin marks, hairline gridlines, no top or right spines, a light-grey interval band and a dashed null-effect line. The style is applied per figure (`plot.style()` is a context manager every plotting function uses), so a dlnmpy plot looks the same in a notebook, a script or inside someone else's figure, and your global matplotlib settings are left alone. `plot.set_theme("journal")` (greyscale, the default) or `plot.set_theme("colour")` (a colour-blind-safe palette validated for adjacent-pair separation: blue estimate, blue-to-red diverging surface with a neutral grey midpoint, blue for cold and red for heat, fixed-order hues for overlays).

Besides the four R plot types (`overall`, `slices`, `contour`, `3d`), there are the figures papers actually print: `plot_overall_risk` (the overall cumulative curve with cold and heat sides, the minimum-risk value and its interval, exposure percentiles along the top and the exposure distribution underneath, as in Gasparrini et al. 2015), `plot_attributable` (attributable fractions or numbers by component with intervals), `overlay_slices` (several curves on one axis with a legend, the R `plot()` + `lines()` idiom) and `summary_figure` (the four-panel figure below). All functions return matplotlib axes, so anything can be adjusted. `examples/gallery.py` renders every figure in both themes.

![summary figure](examples/figures/summary_colour.png)

![attributable](examples/figures/attributable_colour.png)

## Uncertainty beyond the delta method, and model selection

`dlnmpy.uncertainty` gives parametric-bootstrap intervals for any function of the coefficients (`bootstrap(fn, coef, vcov)`, `simulate_pred`, `empirical_ci`), the route used for the MMT and attributable numbers, and `qaic` with `model_grid` to compare cross-basis specifications by quasi-AIC (Gasparrini et al. 2010).

## Beyond the R package: attributable risk and the MMT

`dlnm` stops at prediction; the quantities papers actually report need Gasparrini's separate R scripts. Those are ported and validated here too (`dlnmpy.attribution`):

```python
res = dl.mmt(cb, model, x=chicago.temp, name="cb")       # minimum mortality temperature
res.mmt, res.low, res.high, res.percentile               # with empirical 95% CI and percentile

tab = dl.attr_table(chicago.temp, cb, chicago.death, model, cen=res.mmt, name="cb")
# component | range | an an_low an_high | af af_low af_high
# total, cold, heat, extreme cold, moderate cold, moderate heat, extreme heat

dl.attrdl(chicago.temp, cb, chicago.death, model, type="an", dir="forw", cen=res.mmt, name="cb")
```

`attrdl` is a port of `attrdl.R` (Gasparrini & Leone 2014): backward and forward perspectives, attributable numbers and fractions, daily or total, exposure ranges, reduced coefficients, and empirical intervals by simulation. `findmin` ports `findmin.R` (Tobías, Armstrong & Gasparrini 2017). Both match the R functions to 1e-8 or better on every deterministic quantity (`tests/test_attribution.py`); the simulation path is checked with the same coefficient draws in both languages. `mmt` and `attr_table` wrap them into the summaries used in practice, with one shared set of simulated coefficients so cold, heat and their extreme and moderate parts are mutually consistent. See `examples/attributable_risk.py`.

## Penalised DLNMs without mgcv

`fit_pgam` fits the penalised cross-basis models of Gasparrini et al. (2017) by penalised IRLS with smoothing parameters chosen by REML or ML (Wood 2011), reproducing `gam(..., paraPen=list(cb=cbPen(cb)), method="REML")`: scores to 1e-5, smoothing parameters to 1e-4, coefficients to 1e-4 or better against mgcv (`tests/test_penalized.py`). See `docs/penalized.md` and `examples/penalized_dlnm.py`.

## Two-stage designs: multivariate meta-analysis

Multi-location studies pool location-specific reduced coefficients with a multivariate random-effects meta-analysis; in R that is the `mixmeta` package. `dlnmpy.meta` implements the same model (single random level: REML, ML or fixed effects; unstructured, diagonal or identity between-location covariance; meta-regression; BLUPs with prediction intervals; Cochran's Q and I²), validated against `mixmeta` on a 12-location simulation.

```python
y, S = dl.stack_reduced([dl.crossreduce(cb_i, fit_i, cen=18, name="cb") for ...])
mm = dl.mixmeta(y, S, method="reml")                       # pooled coefficients, Psi, Q, I2
mr = dl.mixmeta(y, S, X=np.column_stack([np.ones(m), mean_temp]))   # meta-regression
pooled = dl.predict_reduced(cb_1, mm.coef_vec, mm.vcov, at=grid, cen=18)   # pooled curve (CrossPred)
blups = mm.blup(se=True)                                    # location-specific BLUPs and covariances
```

Where the two disagree it is because `mixmeta` stopped at its iteration limit: on the identity-structure fit R reports non-convergence and the Python optimum has a higher restricted log-likelihood by 0.76. For all converged fits coefficients, covariances, `Psi`, log-likelihood, AIC/BIC, Q, I², BLUPs and predictions agree to 1e-5 or better (1e-13 for fixed effects). See `examples/two_stage.py`.

## Validation

The table at the top summarises this section; the detail is here for anyone who wants to check the claims.

### A published analysis, reproduced whole

`examples/lancet_2015.py` reproduces the England & Wales analysis of

> Gasparrini A, Guo Y, Hashizume M, et al. Mortality risk attributable to high and low ambient temperature: a multicountry observational study. *The Lancet* 2015;**386**(9991):369-375.

from the author's [public code and data](https://github.com/gasparrini/2015_gasparrini_Lancet_Rcodedata) - 10 regions, 1993-2006, 7,573,716 deaths - running all five stages: a DLNM per region, reduction to the overall cumulative curve, multivariate meta-regression and BLUPs, the minimum-mortality temperature, and attributable deaths with empirical intervals.

|                              | dlnmpy                | R                 |
|------------------------------|-----------------------|-------------------|
| Minimum-mortality percentile | **89.5**              | 89.5              |
| Attributable fraction, total | **8.94%** (8.42-9.43) | 8.94% (8.42-9.43) |
| Attributable fraction, cold  | **8.63%** (8.13-9.12) | 8.63% (8.13-9.12) |
| Attributable fraction, heat  | **0.31%** (0.27-0.34) | 0.31% (0.27-0.34) |

These intervals come from replaying R's own simulated draws through `attrdl(..., coefsim=)`, so the comparison is deterministic; the example script draws its own and so differs in the third digit by Monte Carlo error. `examples/bmcmrm_2013.py` likewise reproduces Gasparrini & Armstrong (2013), the methodological paper behind this two-stage design: minimum-mortality temperature 17.1 °C, pooled RR 1.101 (1.079-1.124) at 22 °C and 1.309 (1.245-1.376) at 0 °C, I² of 63.8/16.4/63.5%, and all three QAIC model-comparison totals.

The minimum-mortality percentile matches as an exact integer in all 10 regions; the reduced coefficients agree to 1e-13 and the attributable fractions to 4e-05 percentage points. See [docs/validation.md](docs/validation.md) for the stage-by-stage agreement and for the two places where it is looser than machine precision (R's default `glm` convergence tolerance, and a deliberately overparameterised meta-regression).

### Unit fixtures

`tools/make_fixtures.R` runs the R package on the vignette examples and writes every intermediate object to `tests/fixtures/` as CSV and JSON: 32 basis-function specifications (with and without values outside the fitting range), 9 cross-bases, 12 prediction objects, 7 reductions, the penalty matrices, and R's `pretty()`, `quantile()` and knot helpers. The Python test-suite compares against these numbers with absolute tolerances of 1e-10 to 1e-12.

```
pytest            # 117 tests
```

A second, independent check lives in `tools/side_by_side.R` / `tools/side_by_side.py`: five complete analyses are run end to end in both languages, including the model fit, on data and specifications not used for the unit-test fixtures (the `drug` trial with OLS on exposure histories, the `nested` case-control study with conditional logistic regression, a logistic and a quasi-Poisson Chicago model with threshold, polynomial, strata and integer bases, 80% and 90% intervals, exposure-history matrices passed to `at`, and a four-city simulation with known truth). All 102 quantities compared agree to better than 1e-8; most to 1e-12. The report is printed by `python tools/side_by_side.py` and the test-suite runs it too.

Two things worth knowing about equivalence with R:

- R's `glm` handles rank-deficient designs by dropping aliased columns (reported as `NA`); statsmodels uses a pseudo-inverse and a numerical rank. `fit_glm` reproduces R's behaviour (drops aliased columns using the same criterion as LINPACK's `dqrdc2`, solves the IRLS step by QR, and computes the quasi-likelihood dispersion with R's residual degrees of freedom). This is why the example-1 coefficients agree to 1e-14 rather than 1e-6.
- R's `glm` stops iterating at a looser tolerance (`epsilon = 1e-8`) than `fit_glm` does; with R's default settings, standard errors differ from Python's at about 1e-7. Tightening R's `glm.control(epsilon = 1e-14)` brings the two to 1e-15, so the difference is R's convergence criterion, not the algorithms.
- Conditional logistic regression: statsmodels' default optimiser (BFGS) can stop on a flat likelihood with coefficients 0.02 away from `survival::clogit`, and its covariance is a numerical approximation. `fit_clogit` uses Newton-Raphson and an extrapolated finite-difference Hessian of the analytic score, which matches R to 1e-9.
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
  penalized.py   fit_pgam / fit_pglm: penalised IRLS with REML/ML smoothing (mgcv's criterion)
  plot.py        matplotlib plots (journal and colour themes; plot.style())
  workflow.py    dlnm(): one call from a data frame to MMT, RR table, attributable fractions, figure
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

1. Analytic derivatives for the penalised fitter (mgcv-speed smoothing parameter selection).
2. Parametric-bootstrap intervals for any derived quantity, and a Bayesian route (same design matrix in PyMC/numpyro).
3. Multilevel meta-analysis (nested random levels, as in `mixmeta`'s extended framework) and longitudinal/repeated-measures structures.
4. A Rust core with Python bindings, validated with the same fixtures.

## Status

Alpha. The numerical core has been stable since 0.4 and is pinned by the fixtures; the API of `dlnm()` and the plotting functions is new in 0.6.0 and may still change. Not on PyPI yet (install from GitHub above). Known gaps are in the roadmap below; the penalised fitter uses numerical derivatives and is slower than mgcv.

## Contributing

See `CONTRIBUTING.md`. The rule is numerical equivalence with R: changes to the core keep `pytest` green against the fixtures, and features that exist in R come with fixtures from `tools/make_fixtures.R`. Ports to Rust, Julia or JavaScript validated against the same fixtures are welcome.

## How to cite

There is no paper for `dlnmpy`. Cite the methods papers below for the models, and the software as:

> Parker C. dlnmpy: distributed lag non-linear models in Python (version 0.6.0). 2026. https://github.com/Logic06183/dlnmpy

A `CITATION.cff` is in the repository, so GitHub's "Cite this repository" button gives the same thing in BibTeX or APA. No DOI yet.

## Prior work

Two earlier Python efforts exist: [aedessler/pydlnm](https://github.com/aedessler/pydlnm) (with meta-analysis and BLUPs, validated on multi-city data) and [`crossbasis`](https://pypi.org/project/crossbasis/) on PyPI. This repository differs in aiming for a literal, attribute-for-attribute port of the R algorithms (including R's own `pretty`, quantile and spline routines so defaults match), in shipping the R-generated fixtures as a language-neutral test oracle, and in documenting the method as a specification for further ports.

## References

- Gasparrini A, Armstrong B, Kenward MG. Distributed lag non-linear models. *Statistics in Medicine* 2010; 29(21):2224-2234.
- Gasparrini A. Distributed lag linear and non-linear models in R: the package dlnm. *Journal of Statistical Software* 2011; 43(8):1-20.
- Gasparrini A, Armstrong B. Reducing and meta-analysing estimates from distributed lag non-linear models. *BMC Medical Research Methodology* 2013; 13:1.
- Gasparrini A. Modeling exposure-lag-response associations with distributed lag non-linear models. *Statistics in Medicine* 2014; 33(5):881-899.
- Gasparrini A, Armstrong B, Kenward MG. Multivariate meta-analysis for non-linear and other multi-parameter associations. *Statistics in Medicine* 2012; 31(29):3821-3839.
- Sera F, Armstrong B, Blangiardo M, Gasparrini A. An extended mixed-effects framework for meta-analysis. *Statistics in Medicine* 2019; 38(29):5429-5444.
- Gasparrini A, Leone M. Attributable risk from distributed lag models. *BMC Medical Research Methodology* 2014; 14:55.
- Tobías A, Armstrong B, Gasparrini A. Investigating uncertainty in the minimum mortality temperature. *Epidemiology* 2017; 28(1):72-76.
- Wood SN. Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. *JRSS B* 2011; 73(1):3-36.
- Gasparrini A, Scheipl F, Armstrong B, Kenward MG. A penalized framework for distributed lag non-linear models. *Biometrics* 2017; 73(3):938-948.

## Licence

GPL-2.0-or-later. This is a derived work of the R package `dlnm` (GPL >= 2) and of R's `splines` and `base` packages (GPL >= 2), whose algorithms it reproduces. The `chicagoNMMAPS` data are redistributed from the `dlnm` package.
