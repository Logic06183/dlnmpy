# Validation

Two kinds of evidence that `dlnmpy` reproduces the R packages it ports.

1. **Unit fixtures** — `tools/make_fixtures*.R` run R (`dlnm` 2.4.10, `mgcv`,
   `mixmeta`, and Gasparrini's published `attrdl.R`/`findmin.R`) and write every
   intermediate to `tests/fixtures/`. The test suite compares against those
   numbers, and `tools/side_by_side.*` runs five complete analyses end to end.
2. **A published analysis, reproduced whole** — the section below.

The meta-analysis moved to the separate package
[mixmetapy](https://github.com/Logic06183/mixmetapy) in 0.8.0, and its own
evidence moved with it: the 12-location comparison with R's `mixmeta`, and
the conventions that differ from `mvmeta`, are in its `docs/validation.md`.
The two published analyses below stay here because they run both stages;
their meta-analysis rows are repeated there. The numbers below were produced
before the split and reproduced unchanged after it (18 September 2026).

## Gasparrini et al. (2015) *Lancet*, England & Wales

> Gasparrini A, Guo Y, Hashizume M, et al. Mortality risk attributable to high
> and low ambient temperature: a multicountry observational study.
> *The Lancet* 2015;**386**(9991):369-375.

The author's code and data are public at
[gasparrini/2015_gasparrini_Lancet_Rcodedata](https://github.com/gasparrini/2015_gasparrini_Lancet_Rcodedata).
`examples/lancet_2015.py` is a translation of its five stages, run on the 10
regions of England and Wales, 1993-2006, 5,113 days each, 7,573,716 deaths.
This exercises `crossbasis`, `fit_glm`, `crossreduce`, `mixmeta` (from mixmetapy), BLUPs,
`findmin` and `attrdl` in one pipeline.

### Result

Rerun with `python examples/lancet2015_vs_r.py` against `tests/reference/lancet2015.json`, which `tools/reproduce_lancet2015.R` writes from the authors' code (pinned commit 7c3e03c, R 4.3.3, dlnm 2.4.7, mvmeta 1.0.3, mixmeta 1.2.2). The R draws behind the intervals are recorded with their seed (20150725) and replayed in Python with R's normals and R's eigenvectors.

|                                | dlnmpy                     | R (authors' code, recorded seed) |
|--------------------------------|----------------------------|----------------------------------|
| Minimum-mortality percentile   | **89.5** (median)          | 89.5                             |
| Attributable fraction, total   | **8.94%** (8.40-9.48)      | 8.938% (8.404-9.477)             |
| Attributable fraction, cold    | **8.63%** (8.08-9.19)      | 8.630% (8.082-9.193)             |
| Attributable fraction, heat    | **0.31%** (0.27-0.35)      | 0.308% (0.266-0.348)             |

Earlier versions of this page gave intervals of 8.42-9.43 from an R run whose seed was not recorded; they cannot be reproduced and have been replaced. `examples/lancet_2015.py` draws its own samples and so differs from R in the second decimal of the intervals by Monte Carlo error.

### Agreement of each stage

Against the authors' code with the first-stage `glm` run to convergence (`glm.control(epsilon = 1e-14)`):

| Stage | Quantity | max abs difference |
|---|---|---|
| 1. First stage | reduced coefficients | 1.0e-13 |
| 1. First stage | reduced covariance | 5.2e-12 |
| 2. Meta-regression | log-likelihood | 1.3e-08 (Python higher) |
| 2. Meta-regression | coefficients | 1.5e-05 (3.0e-06 against R run to its optimum) |
| 2. Meta-regression | between-study covariance Psi | 3.5e-06 (3.5e-07 against R run to its optimum) |
| 3. MMT | minimum-mortality percentile | 0 (exact, all regions) |
| 3. MMT | minimum-mortality temperature | 4.3e-14 |
| 4. Attributable | attributable deaths, per region and total | 0.22 deaths (1.5e-05 relative) |
| 4. Attributable | attributable fraction | 3.4e-05 percentage points |

- **The first stage, published settings.** With R's default `epsilon = 1e-8`, as in the published script, the reduced coefficients differ by 4.8e-10 and the covariance by 1.9e-6 relative. (An earlier version of this page said 1e-7.) R's default tolerance acts mainly through the dispersion.
- **The meta-regression** estimates 30 parameters from 10 studies on a flat likelihood, so R's default stopping rule leaves it about 1e-5 from the optimum. `mvmeta` and `mixmeta` differ from each other by 5.0e-06 at their defaults and agree to 6e-8 when both are run with `reltol = 1e-14`. Against R at that optimum, Python's coefficients agree to 3.0e-6. (An earlier version of this page said dlnmpy "sits inside the spread of two R implementations"; the comparison at the optimum is the right one.)

### Conventions that differ from `mvmeta`

Both are cosmetic, and both will bite a line-by-line translation:

- `mvmeta` and `mixmeta` report REML log-likelihoods that differ by a constant
  (15.85 on this fit). `dlnmpy` follows `mixmeta`.
- The flat coefficient vector and the rows and columns of `vcov()` are
  *predictor-major* in `mixmeta` and in `dlnmpy` (all intercepts, then all
  slopes), and *outcome-major* in `mvmeta`. The `(p, k)` coefficient matrix is
  the same in all three: `mvmeta`'s `$coefficients` matches `MixMeta.coef`.

## Gasparrini & Armstrong (2013) *BMC Medical Research Methodology*

> Gasparrini A, Armstrong B. Reducing and meta-analysing estimates from
> distributed lag non-linear models. *BMC Medical Research Methodology*
> 2013;**13**:1.

The methodological paper behind the two-stage design, with code and data at
[gasparrini/2013_gasparrini_BMCmrm_Rcodedata](https://github.com/gasparrini/2013_gasparrini_BMCmrm_Rcodedata).
`examples/bmcmrm_2013.py` translates it. Same 10 regions as the Lancet example
but a different question, so it covers what that one does not: `crossreduce`
in its overall *and* predictor-specific (`type="var"`) forms, three competing
lag specifications compared by QAIC, random- and fixed-effects pooling,
Cochran's Q and I², meta-regression on latitude, and prediction from pooled
coefficients through a `onebasis`.

### Result

|                                   | dlnmpy                | R, authors' script  |
|-----------------------------------|-----------------------|-----------------------|
| Minimum-mortality temperature     | **17.1 °C**           | 17.1 °C               |
| Pooled RR at 22 °C                | **1.101** (1.079-1.124) | 1.101 (1.079-1.124) |
| Pooled RR at 0 °C                 | **1.309** (1.245-1.376) | 1.309 (1.245-1.376) |
| I², overall / at 22 °C / at 0 °C  | **63.8 / 16.4 / 63.5 %** | 63.8 / 16.4 / 63.5 % |
| QAIC, B-spline lag 0-21           | **402827.3**          | 402827.3              |
| QAIC, constant lag 0-3            | **406212.0**          | 406212.0              |
| QAIC, constant lag 0-21           | **405966.3**          | 405966.3              |

The R column is the authors' public script (updated December 2017), rerun on 24 Sep 2026 with R 4.3.3, dlnm 2.4.7 and mvmeta 1.0.3. The paper itself prints slightly different values from the original 2013 code: RR 1.101 (1.078-1.124) at 22 °C, 1.308 (1.245-1.375) at 0 °C, I² 63.7 / 16.0 / 63.4 %. The paper does not report QAIC; the script computes it.

`python examples/bmcmrm2013_vs_r.py` compares every intermediate against `tests/reference/bmcmrm2013.json`, which `tools/reproduce_bmcmrm2013.R` writes from the authors' code (pinned commit d542c93): the five reductions in all ten regions, the pooled and fixed-effects fits, `Psi`, log-likelihoods, Q, degrees of freedom and p-values, the meta-regression on latitude and its Wald tests, all twelve pooled curves and the QAIC per region. 22,771 comparisons, no failures. Against the authors' code run to convergence the reductions agree to 1e-13, their covariances to 2e-11, QAIC totals to 1e-6 and pooled relative risks to 2e-6. Against R's meta-analyses run to their optimum (`reltol = 1e-14`), pooled coefficients agree to 6e-10 and meta-regression coefficients to 3e-7. At its defaults R stops short of the optimum of the 22 °C pooled fit: log-likelihood 168.244063, where R at the tight tolerance gives 168.2440677686 and Python 168.2440677685. (The earlier "67 of 68 quantities" summary is superseded by this script.)

As in the Lancet example, comparing against the published script as written
shows larger differences that are not the port's: R's `glm` defaults to
`epsilon = 1e-8`, and the script uses `mvmeta`, which reports a REML
log-likelihood differing from `mixmeta`'s by a constant that depends on the
dimensions (4.6052 at k=4, 5.7565 at k=5 with an intercept-only model). With
`glm.control(epsilon = 1e-14)` and `mixmeta`, the 16 apparent discrepancies
reduce to the single one above.

### Why the QAIC fix mattered

This paper selects its lag specification by QAIC, so it is a direct test of
`uncertainty.qaic`, which before 0.4.1 read the log-likelihood from
statsmodels' `results.llf` — divided by the dispersion for a quasi family.
Rerunning this comparison with the old implementation:

| model | QAIC (0.4.1+) | QAIC (before) |
|---|---|---|
| B-spline of lag 0-21 (the paper's choice) | **402827.3** | 368280.5 |
| constant of lag 0-3 | 406212.0 | **350059.0** |
| constant of lag 0-21 | 405966.3 | 351008.5 |

The old criterion selected the constant lag 0-3 model, reversing the paper's
conclusion. Dividing the fit term by phi rewards the more over-dispersed, and
so the more under-fitted, specification.

## Missing data

Scattered missing values in the exposure and the outcome propagate exactly as
in R: the NaN pattern of the cross-basis is identical, the same rows are
dropped from the fit (3,208 of 5,114 in a test with 200 missing exposures and
100 missing outcomes), and the coefficients agree to 6.9e-15.

## The workflow layer against R (0.7.0)

The layers above `dlnm`'s own functions (`dlnm()`, `fit_glm`, `fit_clogit`, `attrdl(group=)`) were compared against hand-written R on 14 September 2026 (R 4.3.3, `dlnm` 2.4.7, `survival`, `mixmeta` 1.2.2). Scripts are not in the repository because they are one-off; the numbers are:

| Comparison | Max abs | Max rel |
|---|---|---|
| `fit_glm(exposure=pop)` with an aliased column vs `glm(offset=log(pop))`: 39 coefficients | 1.1e-13 | 6.3e-13 |
| same, standard errors and dispersion | 9.0e-12, 5.0e-10 | 1.3e-10, 3.4e-10 |
| same, `crosspred` RR and interval bounds on a 51-point grid | 2.1e-14, 1.1e-11 | 1.7e-14, 1.0e-11 |
| `fit_clogit` with NaN rows vs `clogit(na.action)`: coefficient, SE, log-likelihood | 1.0e-13, 2.1e-14, 2.8e-14 | 2.3e-13 |
| `attrdl(group=)` vs `attrdl.R` per group, 3 series: daily AN back and forward, daily AF | 2.7e-14 | 2.8e-13 |
| `dlnm(lag=21, time, dow)` vs `glm(death ~ cb + ns(time) + dow)`: 84 coefficients | 3.9e-13 | 4.8e-10 |
| same, `rr_at([1, 99])` RR and bounds, attributable fraction total/cold/heat | 5.1e-12, 1.1e-15 | 4.5e-12 |
| `dlnm(group=)` on two stacked series vs `glm` with `crossbasis(group=)`, group intercept and block-diagonal `ns(time)`: 86 coefficients | 9.4e-14 | 9.3e-13 |

One quirk of R worth knowing: `glm.control(epsilon = 1e-14)` also tightens the aliasing tolerance to `epsilon / 1000`, so an exact duplicate column is no longer detected and the fit does not converge. With `epsilon = 1e-10` R flags the alias and the two implementations agree as above.

One definitional point: the grouped total attributable number is one pooled fraction, `sum(an) / sum(cases)` over rows with a complete lag history, applied to all cases. It is not the sum of per-group totals, which weight each group by its own ratio of all cases to cases with a complete history; the two differ by a fraction of a percent when the first `lag` days of a series are unusual. Per-group totals are obtained by calling `attrdl` per group.

## Four further published analyses (0.8.4)

Each has an R script in `tools/` that pins the authors' repository, runs their code and writes `tests/reference/<name>.json` with the R and package versions, and a Python translation in `examples/` that compares every quantity and exits non-zero on any failure. Tolerances and every departure from the authors' code are stated at the top of each script.

| Paper | Authors' repository (pinned commit) | Python script | Comparisons | Result |
|---|---|---|---|---|
| Gasparrini (2014) *Stat Med*, uranium miners cohort, Cox models | `2014_gasparrini_StatMed_Rcodedata` (3ae639c) | `statmed2014.py` | 102,224 | 0 failures; coefficients and vcov to 1e-9, predictions within 4.4e-7 of R's SE; about 8 min |
| Gasparrini & Leone (2014) *BMC MRM*, attributable risk | `2014_gasparrini_BMCmrm_Rcodedata` (68a8d44) | `leone2014.py` | 22,690 | 0 failures; fractions to 1e-13 against R run to convergence; intervals from R's draws to 4e-6 relative |
| Vicedo-Cabrera, Sera & Gasparrini (2019) *Epidemiology*, projections | `2019_vicedo-cabrera_Epidem_Rcodedata` (313606c) | `vicedo2019.py` | 11,060 | 0 failures; fractions to 1e-12, intervals to 1e-10 |
| Gasparrini, Scheipl, Armstrong & Kenward (2017) *Biometrics*, penalised DLNMs | `2017_gasparrini_Biomet_Rcodedata` (531c5c8) | `biomet2017.py` | 66 | 0 failures at mgcv's smoothing parameters; dlnmpy's own selection never scored worse than mgcv |

These replace the figures quoted for the same papers in the 0.8.1 and 0.8.2 changelog entries, whose scripts were not committed. The 0.8.1 figure of "197 quantities to 1.8e-10" for the Stat Med paper cannot be traced: against the authors' default `coxph` fits, model 8 differs by 1.6e-6 standard errors, which is R's stopping rule. The R scripts therefore also store each fit refitted to convergence, and the tightest comparisons are against those.
