# From R `dlnm` to `dlnmpy`

The public API mirrors the R package so that existing R analyses translate almost line by line. This page lists the mapping and the places where the two differ.

## Functions and objects

| R | Python |
|---|---|
| `onebasis(x, fun="ns", df=5)` | `dl.onebasis(x, "ns", df=5)` |
| `crossbasis(x, lag=21, argvar=list(fun="ns", df=4), arglag=list(fun="ns", df=4))` | `dl.crossbasis(x, lag=21, argvar={"fun": "ns", "df": 4}, arglag={"fun": "ns", "df": 4})` |
| `crossbasis(x, lag=c(2, 8), ...)` | `dl.crossbasis(x, lag=(2, 8), ...)` |
| `crossbasis(..., group=data$year)` | `dl.crossbasis(..., group=data.year)` |
| `argvar=list(fun="thr", thr=25)` | `argvar={"fun": "thr", "thr": 25}` (also `thr_value`, `thr.value`) |
| `argvar=list(fun="ns", knots=k, Boundary.knots=b)` | `argvar={"fun": "ns", "knots": k, "boundary_knots": b}` (also `Boundary.knots`) |
| `logknots(30, 3)`, `equalknots(x, fun="bs", df=5, degree=2)` | same |
| `exphist(exp, times, lag, fill)` | same |
| `glm(death ~ cb + ns(time, 98) + dow, quasipoisson(), data)` | `dl.fit_glm("death ~ " + " + ".join(cols) + " + C(dow)", data, family="quasipoisson")` where `cols` are the columns of `cb.to_dataframe("cb")` and `dl.onebasis(time, "ns", df=98).to_dataframe("ns_time")` |
| `clogit(case ~ cb + strata(riskset), nested)` | `dl.fit_clogit(nested.case, cb.to_dataframe("cb"), groups=nested.riskset)` |
| `crosspred(cb, model, at=0:20, bylag=0.2, cumul=TRUE)` | `dl.crosspred(cb, model, at=range(21), bylag=0.2, cumul=True, name="cb")` |
| `crosspred(cb, model, from=-10, to=30, by=1, cen=21)` | `dl.crosspred(cb, model, from_=-10, to=30, by=1, cen=21, name="cb")` |
| `crosspred(cb, coef=b, vcov=V, model.link="log")` | `dl.crosspred(cb, coef=b, vcov=V, model_link="log")` |
| `crossreduce(cb, model, type="var", value=33)` | `dl.crossreduce(cb, model, type="var", value=33, name="cb")` |
| `pred$allRRfit["10"]` | `pred.allRRfit[pred._var_index(10)]` or `pred.overall().query("var == 10")` |
| `pred$matRRfit`, `pred$matRRlow`, `pred$cumRRfit` | `pred.matRRfit`, `pred.matRRlow`, `pred.cumRRfit` |
| `pred$matfit`, `pred$matse`, `pred$allfit`, `pred$allse` | same names |
| `red$coefficients`, `vcov(red)` | `red.coef`, `red.vcov` |
| `plot(pred, "overall")` | `pred.plot("overall")` |
| `plot(pred, "slices", var=c(-20, 33), lag=c(0, 5))` | `pred.plot("slices", var=[-20, 33], lag=[0, 5])` |
| `plot(pred, "contour", xlab=, ylab=, key.title=title("RR"))` | `pred.plot("contour", xlab=, ylab=, zlab="RR")` (`ylab` is the lag axis, `zlab` the outcome, as in R) |
| `plot(pred, xlab=, ylab=, zlab=)` (3d) | `pred.plot("3d", xlab=, ylab=, zlab=)` |
| `plot(red)` | `red.plot()` |
| `summary(cb)`, `summary(pred)` | `print(cb.summary())`, `print(pred.summary())` |
| `cbPen(cb)` | `dl.cbpen(cb)` |
| `attrdl(x, cb, cases, model, cen=mmt)` (Gasparrini script) | `dl.attrdl(x, cb, cases, model, cen=mmt, name="cb")` |
| `findmin(cb, model, from, to, by)` (Gasparrini script) | `dl.findmin(cb, model, from_=, to=, by=, name="cb")`, or `dl.mmt(...)` with CI |
| `mixmeta(y ~ 1, S, method="reml")` (mixmeta package) | `dl.mixmeta(y, S, method="reml")` |
| `mixmeta(y ~ x, S)` | `dl.mixmeta(y, S, X=np.column_stack([np.ones(m), x]))` |
| `blup(mm, se=TRUE, pi=TRUE, vcov=TRUE)` | `mm.blup(se=True)` |
| `predict(mm, newdata, se=TRUE)` | `mm.predict(Xnew, se=True)` |
| `qtest(mm)`, `summary(mm)$i2stat` | `mm.qtest()` |
| `chicagoNMMAPS` | `dl.datasets.chicago_nmmaps()` |

## Linking a basis to the coefficients of a fitted model

R finds the coefficients of a cross-basis in a fitted model by matching the object name in `coef(model)`. Python has no equivalent of a formula picking up a matrix by name, so the link is explicit:

1. `cb.to_dataframe("cb")` gives a DataFrame with columns `cb_v1_l1, cb_v1_l2, ...`;
2. these columns go into the design (with `dl.design_matrix` or a plain `join`);
3. `crosspred(cb, model, name="cb")` finds the coefficients whose names match `cb` followed by `v{i}_l{j}` (or R's `v{i}.l{j}`) and orders them by `(i, j)`.

If the model was fitted elsewhere (R, Stan, a Rust fitter), pass `coef` and `vcov` directly.

## Behavioural differences

- **Rank deficiency.** R drops aliased design columns and reports `NA` coefficients. `fit_glm` does the same (the dropped names are in `results.aliased`); with `drop_aliased=False`, statsmodels' pseudo-inverse is used and standard errors can differ slightly on ill-conditioned designs.
- **Categorical reference level.** patsy's `C(dow)` uses the alphabetically first level as reference; R uses the first factor level. Coefficients of the cross-basis are unaffected.
- **Warnings and messages.** R prints a message when the centring value is chosen automatically; Python is silent (`pred.cen` holds the value).
- **Row names.** For a matrix `at`, R labels predictions with the row names; Python uses the row index `0..n-1` in `pred.predvar`.
- **`crosspred` for `mgcv::gam` smooths** (`type="gam"` in R) is not implemented, as there is no mgcv in Python. Penalised fitting is discussed in `penalized.md`.
- **Missing values in `x`** produce NaN rows in the basis and cross-basis, as in R. `fit_glm` drops them (`missing="drop"`).

## Objects

`OneBasis` and `CrossBasis` behave like numpy arrays (`np.asarray(cb)`, `cb.shape`, `cb[:, 0]`) and carry the attributes R stores on the matrix: `cb.df`, `cb.range`, `cb.lag`, `cb.argvar`, `cb.arglag`. `CrossPred` and `CrossReduce` are dataclasses whose fields follow the R list components; confidence limits and exponentiated versions are properties computed from `fit`, `se` and `ci_level`.
