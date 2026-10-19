# Porting DLNMs to another language

The repository is organised so that a port to Rust, Julia, JavaScript or anything else can be built and validated without R installed. The three ingredients are the specification (`theory.md`), the fixtures (`tests/fixtures/`) and the Python implementation as a worked example.

## Order of work

The dependencies run in one direction, so this order keeps every step testable.

1. Port the R-compatible helpers in `_rcompat.py`: `quantile` type 7, `median`, `pretty` and `seq`. Test them against `helpers.json` (`pretty_cases`) and the quantile examples in `tests/test_basis.py`. Getting `pretty` right matters because it determines the default prediction grid and the automatic centring value.
2. Port the B-spline evaluation (`_splines.py::spline_design`) by translating `spline_basis` from R's `splines.c` as written. The edge cases (repeated knots, `x` equal to the last knot, `outer.ok`) are what make `ns`, `bs` and `ps` agree with R.
3. Port the basis functions (`basis.py`) and test each against `onebasis.json`. For every specification it holds the input `x`, the basis matrix, the attributes R stored, a vector `xpred` (including values outside the fitting range and a value beyond the boundary knots) and the basis on `xpred` rebuilt from the attributes. Both matrices must match to ~1e-12.
4. Port the lag functions (`lag.py`), `lag_matrix` with groups and `exphist`, and test them with `exphist_Q.csv` and the small cases in `cases.json` (`ex6.Q2`, `ex6.Q3`).
5. Build the cross-basis (`core.py`) and test it with `cb1_pm.csv` ... `cb8.csv` and the metadata in `cases.json` (`df`, `range`, `lag`, stored `argvar`/`arglag`).
6. Port prediction (`predict.py`), that is `mkat`, `mkcen`, `mkXpred` and `crosspred`, and test it with `mkat.json` and the `pred*` entries of `cases.json`, which carry R's own `coef` and `vcov` so the test isolates the DLNM algebra from the regression fitter.
7. Port the reduction, `crossreduce`, and test it with the `red*` entries.
8. Port the penalties, `cbpen`, and test them with `ex6.pen7`.
9. Model fitting is outside the DLNM core. Whatever fitter is used, the contract is: coefficient vector and covariance matrix for the cross-basis columns, plus the link. `cases.json` includes R's `glm` output (`coef_all`, `coef_names`, `dispersion`) for the vignette models so a fitter can be checked too.

## Fixture formats

- The `*.csv` files hold matrices as written by R's `write.csv` (header row, no row names, NA as empty). Row order is observation order.
- `onebasis.json` has the structure `{spec_name: {spec, x, basis, attributes, xpred, basis_pred}}`. Matrices are arrays of rows, and JSON `null` is NA.
- `cases.json` holds `ex1` ... `ex7`, each with cross-basis metadata (`cbmeta`), models (`modelmeta`: `coef_all`, `coef_names`, `coef`, `vcov`, `dispersion`, `link`), predictions (`predmeta`: `predvar`, `cen`, `lag`, `bylag`, `coef`, `vcov`, `matfit`, `matse`, `allfit`, `allse`, optional `cumfit`/`cumse`, and either `allRRfit`/`allRRlow`/`allRRhigh`/`matRRfit` or `alllow`/`allhigh`) and reductions (`redmeta`: `type`, `value`, `coef`, `vcov`, `basis`, `predvar`, `cen`, `lag`, `bylag`, `fit`, `se`, RR or plain limits).
- `helpers.json` and `mkat.json` hold small vectors.

`jsonlite` unboxes length-one vectors to scalars. Treat scalars and length-one arrays as equivalent.

## Regenerating the fixtures

```
Rscript tools/make_fixtures.R
```

requires R with `dlnm`, `splines` and `jsonlite`. The script is deterministic (no randomness is involved except the seed set for parity with the vignette). If a future `dlnm` release changes a default, re-run it and note the version in the commit.

## Conventions to keep across languages

- The cross-basis columns are ordered with the lag index fastest (`v1.l1, v1.l2, ...`).
- The stored attributes of a basis are exactly the arguments the basis function accepts (R's `formals()`), so prediction re-applies the same transformation. Do not store `df` for `ns`/`bs` (the knots are stored instead). Do store the full knot vector for `ps`.
- Missing values propagate as NaN rows and are never silently dropped by the basis code.
- Values outside the fitting range at prediction time are allowed and follow R's extrapolation rules (linear for `ns`, Taylor for `bs`, zero for `ps`, natural for the others).

## Suggested design for a compiled core

A Rust (or C++) core with bindings only needs to expose: `spline_design`, the basis functions, `crossbasis`, `mkXpred`, and the reduction matrix `M`. Everything else (the grid, centring, intervals, plotting) is cheap and can stay in the host language. The heavy operation in practice is the cross-basis for long series with many columns, and the prediction tensor for fine grids. Both are dense matrix products and vectorise well.
