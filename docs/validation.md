# Validation

Two kinds of evidence that `dlnmpy` reproduces the R packages it ports.

1. **Unit fixtures** — `tools/make_fixtures*.R` run R (`dlnm` 2.4.10, `mgcv`,
   `mixmeta`, and Gasparrini's published `attrdl.R`/`findmin.R`) and write every
   intermediate to `tests/fixtures/`. The test suite compares against those
   numbers, and `tools/side_by_side.*` runs five complete analyses end to end.
2. **A published analysis, reproduced whole** — the section below.

## Gasparrini et al. (2015) *Lancet*, England & Wales

> Gasparrini A, Guo Y, Hashizume M, et al. Mortality risk attributable to high
> and low ambient temperature: a multicountry observational study.
> *The Lancet* 2015;**386**(9991):369-375.

The author's code and data are public at
[gasparrini/2015_gasparrini_Lancet_Rcodedata](https://github.com/gasparrini/2015_gasparrini_Lancet_Rcodedata).
`examples/lancet_2015.py` is a translation of its five stages, run on the 10
regions of England and Wales, 1993-2006, 5,113 days each, 7,573,716 deaths.
This exercises `crossbasis`, `fit_glm`, `crossreduce`, `mixmeta`, BLUPs,
`findmin` and `attrdl` in one pipeline.

### Result

|                                | dlnmpy                | R                     |
|--------------------------------|-----------------------|-----------------------|
| Minimum-mortality percentile   | **89.5**              | 89.5                  |
| Attributable fraction, total   | **8.94%** (8.42-9.43) | 8.94% (8.42-9.43)     |
| Attributable fraction, cold    | **8.63%** (8.13-9.12) | 8.63% (8.13-9.12)     |
| Attributable fraction, heat    | **0.31%** (0.27-0.34) | 0.31% (0.27-0.34)     |

The minimum-mortality percentile agrees as an exact integer in all 10 regions.
The intervals above were computed by replaying R's own simulated coefficient
draws through `attrdl(..., coefsim=)`, so they are a deterministic comparison;
`examples/lancet_2015.py` draws its own and so differs from R in the third
digit by Monte Carlo error.

### Agreement of each stage

| Stage | Quantity | max abs difference |
|---|---|---|
| 1. First stage | reduced coefficients | 1.2e-13 |
| 1. First stage | reduced covariance | 5.2e-12 |
| 2. Meta-regression | log-likelihood | 1.4e-08 |
| 2. Meta-regression | coefficients | 1.3e-05 |
| 2. Meta-regression | between-study covariance Psi | 3.1e-06 |
| 3. MMT | minimum-mortality percentile | 0 (exact, all regions) |
| 3. MMT | minimum-mortality temperature | 1e-10 |
| 4. Attributable | attributable deaths per region | 0.2 deaths (of thousands) |
| 4. Attributable | attributable fraction | 3.9e-05 percentage points |

Two notes on the two places agreement is looser than machine precision:

- **The first stage** matches to 1e-13 only when R's `glm` is given
  `glm.control(epsilon = 1e-14)`. With R's default `epsilon = 1e-8`, as in the
  published script, the reduced coefficients differ at 1e-7 — that is R's
  convergence tolerance, not a difference in the algorithms.
- **The meta-regression** is fitted by `mvmeta` in the published code and by a
  port of `mixmeta` here. The log-likelihood agrees to 1.4e-08, so both reach
  the same optimum, but the parameters are only pinned to about 1e-5 because
  the likelihood is flat: this model estimates 30 parameters from 10 studies,
  which the original script flags in a comment. For scale, R's own `mvmeta` and
  `mixmeta` differ from each other by 5.0e-06 on the same data, so `dlnmpy`
  sits inside the spread of two R implementations of the same model.

### Conventions that differ from `mvmeta`

Both are cosmetic, and both will bite a line-by-line translation:

- `mvmeta` and `mixmeta` report REML log-likelihoods that differ by a constant
  (15.85 on this fit). `dlnmpy` follows `mixmeta`.
- The flat coefficient vector and the rows and columns of `vcov()` are
  *predictor-major* in `mixmeta` and in `dlnmpy` (all intercepts, then all
  slopes), and *outcome-major* in `mvmeta`. The `(p, k)` coefficient matrix is
  the same in all three: `mvmeta`'s `$coefficients` matches `MixMeta.coef`.

## Missing data

Scattered missing values in the exposure and the outcome propagate exactly as
in R: the NaN pattern of the cross-basis is identical, the same rows are
dropped from the fit (3,208 of 5,114 in a test with 200 missing exposures and
100 missing outcomes), and the coefficients agree to 6.9e-15.
