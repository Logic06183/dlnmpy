# Distributed lag non-linear models: a language-neutral specification

This document states what the software computes, precisely enough that an implementation in any language can be written from it and checked against the fixtures in `tests/fixtures/`. Notation follows Gasparrini, Armstrong & Kenward (2010) and Gasparrini (2014), with the column ordering used by the R package since version 2.2.4.

## 1. The model

For observations `t = 1..n` with outcome `y_t`, exposure series `x_t` and lags `l = l0..L`, a DLNM is a generalised linear (or similar) model

    g(E[y_t]) = alpha + s(x_t, ...; eta) + sum_k gamma_k z_tk

where `s` is the exposure-lag-response function

    s(x, t; eta) = sum_{l=l0}^{L} f.w(x_{t-l}, l)

`f.w(x, l)` is a bi-dimensional function: `f` describes the exposure-response shape at a given lag and `w` the lag-response shape at a given exposure. Both are represented with basis functions. With `R(x)` an `n x vx` basis matrix for the predictor space and `C` an `(L-l0+1) x vl` basis matrix for the lag space (rows correspond to lags `l0..L`), the *cross-basis* is

    W = [ Q_1 C | Q_2 C | ... | Q_vx C ]           (n x vx*vl)

where `Q_v` is the `n x (L-l0+1)` matrix whose column `l` holds the values of the `v`-th predictor basis at time `t-l`, that is `Q_v[t, l] = R(x_{t-l})[v]`. The model is then linear in `W` with coefficient vector `eta` of length `vx*vl`, ordered with the lag index varying fastest: `v1.l1, v1.l2, ..., v1.lvl, v2.l1, ...`.

If instead of a time series the data are a matrix of exposure histories (rows = observations, columns = exposures at lags `l0..L`), `Q_v` is obtained by applying the `v`-th predictor basis to every cell of that matrix.

When a grouping variable is supplied, lags are computed within groups: the first `L` rows of each group are missing (NaN), never borrowing from the previous group.

## 2. Basis functions

Every basis function maps a vector `x` to a matrix and returns the attributes needed to apply the identical transformation to new values. Missing `x` give rows of NaN.

**lin**: `[x]`, or `[1, x]` with intercept.

**poly** (degree `d`, scale `s`, default `s = max|x|`): columns `(x/s)^k` for `k = 1..d`, or `k = 0..d` with intercept.

**strata** (breaks `b_1 < ... < b_m`, reference `ref`): indicators of the left-closed intervals `[min-0.0001, b_1), [b_1, b_2), ..., [b_m, max+0.0001)`; the `ref`-th indicator (1-based) is dropped unless `ref = 0`; with intercept a column of ones is prepended. Without explicit breaks, `df - intercept` breaks are placed at the quantiles `k/(df-intercept+1)` (R type 7 quantiles). With `df = 1` and intercept, the basis is a single column of ones.

**thr** (threshold values, side): `h` (high) gives `max(x - t, 0)`; `l` (low) gives `-min(x - t, 0)`; `d` (double) gives the two columns `-min(x - t_1, 0)` and `max(x - t_2, 0)`. Default threshold is the median of `x`; default side is `d` if two thresholds are given, else `h`.

**integer** (values): one indicator column per value; the first column is dropped unless intercept; `x` outside the values gives NaN.

**ns** (natural cubic spline, `splines::ns`): with boundary knots `B = (b_l, b_r)` (default range of `x`) and interior knots `k` (default: quantiles at `1/(m+1), ..., m/(m+1)` with `m = df - 1 - intercept`), form the cubic B-spline basis on the augmented knot sequence `sort(c(rep(B, 4), k))`, drop the first column unless intercept, then impose linearity beyond the boundary knots by projecting out the second-derivative constraints: with `D` the `2 x p` matrix of second derivatives of the B-spline basis at the boundary knots, take the full QR decomposition `D' = Q R` and use `basis %*% Q[, 3:p]`. Values outside the boundary knots are extrapolated linearly from the boundary (value and first derivative at the boundary knot). The result has `df` columns.

**bs** (B-spline, `splines::bs`, degree `d`, order `o = d+1`): cubic by default; interior knots as for ns with `m = df - o + (1 - intercept)`; augmented knots `sort(c(rep(B, o), k))`; drop the first column unless intercept. Values outside the boundary knots are extrapolated by a Taylor expansion of degree `d` around a pivot at `(3/4) b + (1/4) k_o` (left) or `(3/4) b + (1/4) k_{n-o}` (right), where the `k`'s are elements of the augmented sequence.

**ps** (P-spline): B-spline basis of degree `d` on `nik + 2d` equally spaced knots, `nik = df - d + 2 - intercept`, spanning `[min(x) - 0.001 range, max(x) + 0.001 range]` extended by `d` knot intervals on each side; values outside the outer knots give zero; the first column is dropped unless intercept. The penalty matrix is `S = D' D` with `D` the `diff`-th order difference matrix of size `(p+1) x (p+1)` (`p+1` being the number of columns before dropping the first), symmetrised, with the first row and column removed when the intercept is dropped.

**cr** (cubic regression spline, mgcv): with knots `x_1 < ... < x_k` (default: quantiles of the unique values of `x` at `0, 1/(k-1), ..., 1`), the spline is parameterised by its values at the knots. Let `h_j = x_{j+1} - x_j`, `D` the `(k-2) x k` second-difference matrix with rows `(1/h_j, -1/h_j - 1/h_{j+1}, 1/h_{j+1})` and `B` the `(k-2) x (k-2)` tridiagonal matrix with diagonal `(h_j + h_{j+1})/3` and off-diagonals `h_{j+1}/6`. Then `F = [0; B^-1 D; 0]'` maps knot values to second derivatives at the knots and the basis on `x in [x_j, x_{j+1}]` is `c_j^- F[:, j] + c_j^+ F[:, j+1] + a^- e_j + a^+ e_{j+1}` with `a^- = (x_{j+1} - x)/h_j`, `a^+ = (x - x_j)/h_j`, `c^- = ((x_{j+1}-x)^3/h_j - h_j (x_{j+1}-x))/6`, `c^+ = ((x-x_j)^3/h_j - h_j (x-x_j))/6` (Wood 2006, s4.1.2); beyond the outer knots the spline is linear. The wiggliness penalty is `S = D' B^-1 D`. The first column is dropped unless intercept, as is the first row/column of `S`.

B-spline evaluation follows de Boor's recursion exactly as implemented in R's `splines.c` (`spline_basis`), including the treatment of `x` equal to the last knot (assigned to the last interval).

## 3. Cross-basis construction

Inputs: `x` (vector, or matrix with `L-l0+1` columns), `lag = (l0, L)` (a scalar `L` means `(0, L)`; a negative scalar `-L` means `(-L, 0)`), `argvar`, `arglag`, optional `group`.

1. Predictor basis: `R = fun_var(as.vector(x), argvar)`. For a matrix `x`, the vector is the column-major flattening.
2. Lag basis: if `arglag` is empty or `l0 = L`, use `strata(df=1, intercept=TRUE)` (a column of ones, i.e. a moving sum). Otherwise `C = fun_lag(seq(l0, L), arglag)`; `intercept = TRUE` is added when the function accepts it and it was not specified; centring is never applied to the lag basis.
3. For each predictor basis column `v`, form `Q_v` (lagged copies for a series; reshape for a matrix) and set `W[:, vl*(v-1) + l] = Q_v C[:, l]`.
4. Store `df = (vx, vl)`, `range(x)`, `lag`, and the attributes of both bases (only the arguments that the basis function accepts, as R does through `formals()`), plus the centring value `cen` if one was supplied for the predictor.

## 4. Prediction (`crosspred`)

Inputs: the cross-basis (with its attributes), the coefficients `eta` and covariance `V` of its columns, the link, a grid of predictor values `at` (or `from`/`to`/`by`), a lag range (default the fitted one) and lag step `bylag`, a centring value `cen`, and the confidence level.

**Grid.** If `at` is not given: `lo = from or range[1]`, `hi = to or range[2]`; `pretty(c(lo, hi), n)` with `n = 50` if `by` is missing else `max(1, diff(range)/by)`; keep values within `[lo, hi]`; if `by` is given, return `seq(min(pretty), hi, by)`. R's `pretty` is reproduced exactly (`src/appl/pretty.c`). If `at` is a vector it is sorted and de-duplicated. If `at` is a matrix of exposure histories it must have `L-l0+1` columns and `bylag = 1`.

**Centring.** For `lin`, `strata`, `thr`, `integer` the value is used only if numeric (logical values are ignored). For other functions, `cen = NULL` or `TRUE` gives the automatic value `median(pretty(range))`; `FALSE` disables centring. Centring is always disabled when the predictor basis contains an intercept. Centring is applied marginally to the predictor basis: `R_c(x) = R(x) - R(cen)`.

**Lag-specific effects.** Let `predvar` be the grid (length `nv`) and `predlag = seq(l0, L, bylag)` (length `nl`). Form `varvec = rep(predvar, nl)` (or the column-major flattening of the matrix `at`) and `lagvec = rep(predlag, each = nv)`. Then

    Xpred = rowwise_tensor( R_c(varvec), C(lagvec) )        (nv*nl x vx*vl)

with the row-wise tensor product ordering columns as in the cross-basis (lag index fastest). Then

    matfit = reshape(Xpred eta, nv, nl)
    matse  = reshape( sqrt( rowSums( (Xpred V) * Xpred ) ), nv, nl )

**Overall cumulative effects.** Recompute `Xpred` for the integer lags `seq(l0, L)`; sum the `nv x p` blocks over lags to get `Xall`, then `allfit = Xall eta`, `allse = sqrt(rowSums((Xall V) * Xall))`. If cumulative effects are requested, the partial sums after each lag give `cumfit`/`cumse` (integer lags only).

**Intervals.** `z = Phi^{-1}(1 - (1 - level)/2)`; `low = fit - z se`, `high = fit + z se`. For log and logit links the exponentiated quantities are reported as relative risks (odds ratios).

For a one-dimensional basis (`onebasis`) the same applies with `lag = (0, 0)` and `Xpred = R_c(varvec)`.

## 5. Reduction (`crossreduce`)

Given `eta` (length `vx*vl`) and `V`, define a transformation `M` such that `theta = M eta` are the coefficients of a one-dimensional basis:

- overall cumulative exposure-response: `M = I_vx (x) (1' C)` where `1' C` is the column sum of the lag basis over the integer lags (a `1 x vl` row); the new basis is `R_c(at)` (`vx` columns).
- lag-specific exposure-response at lag `l*`: `M = I_vx (x) C(l*)`; new basis `R_c(at)`.
- predictor-specific lag-response at `x*`: `M = R_c(x*) (x) I_vl`; new basis `C(seq(l0, L, bylag))` (`vl` columns).

`(x)` is the Kronecker product. Then `theta = M eta`, `V_theta = M V M'`, `fit = basis theta`, `se = sqrt(rowSums((basis V_theta) * basis))`. The reduced fit for the overall cumulative summary equals `allfit` from `crosspred` exactly.

## 6. Exposure histories (`exphist`)

Given an exposure profile `e_1..e_m` at integer times and requested times `t`, the row for `t` holds `e_{t-l0}, e_{t-l0-1}, ..., e_{t-L}`, padding with `fill` where the index falls outside `1..m`.

## 7. Penalties (`cbpen`)

For `ps` (or `cr`) bases with penalty matrices `S_var` (`vx x vx`) and `S_lag` (`vl x vl`), the cross-basis penalties are `S_var (x) I_vl` and `I_vx (x) S_lag`, each divided by its largest eigenvalue. Additional lag penalties are `I_vx (x) S` with `S` likewise rescaled. Ranks are the numbers of eigenvalues above `10 * eps * max eigenvalue`.

## 8. Conventions that matter for numerical equivalence

- Quantiles are R type 7 (linear interpolation between order statistics), NaN removed.
- `median` is the mean of the two central order statistics.
- `pretty` follows `R_pretty` with `n = 5`, `min.n = n %/% 3`, `shrink.sml = 0.75`, `high.u.bias = 1.5`, `u5.bias = 0.5 + 1.5 * high.u.bias`; non-integer `n` is truncated.
- `seq(from, to, by)` uses `floor((to - from)/by + 1e-10)` steps.
- Column naming: `v{i}.l{j}` with 1-based indices.
- Lag basis values for non-integer lags are obtained by evaluating the lag basis function at those lags (splines and polynomials interpolate; `strata` steps; `integer` is not allowed).

## 9. Attributable risk (Gasparrini & Leone 2014)

With the fitted cross-basis coefficients `eta` and a counterfactual reference `x0` (typically the MMT), the risk on day `t` relative to `x0` is `RR_t = exp(sum_l [f.w(x_{t-l}, l) - f.w(x0, l)])` (backward perspective: the risk at `t` attributable to exposures over the past `L` lags). The attributable fraction is `AF_t = 1 - 1/RR_t = 1 - exp(-eta' Xall_t)` where `Xall_t` is the centred cross-basis row summed over lags, and the attributable number `AN_t = AF_t n_t` with `n_t` the observed cases. Totals: `AF = sum AN_t / sum n_t` over days with complete lagged exposure, `AN = AF * total cases`. The forward perspective attributes the exposure on day `t` to cases over the following `L` days, using the mean of future cases in place of `n_t`. Restricting `x` to a range (cold, heat, extreme) sets exposures outside it to `x0`, so they contribute no risk. Because `AF = 1 - exp(-sum)`, components do not add exactly. Empirical intervals come from `nsim` draws `eta* ~ N(eta, V)` pushed through the same formula; the draws use the eigen square root of `V`.

## 10. Minimum mortality exposure (Tobías, Armstrong & Gasparrini 2017)

The MMT is the `x` in a grid (by default the 1st to 99th percentiles at step 0.1) minimising the overall cumulative curve `Xall(x) eta`. Its uncertainty is the empirical distribution of the minimiser over draws `eta* ~ N(eta, V)`.

## 11. Multivariate meta-analysis (Gasparrini, Armstrong & Kenward 2012)

Location `i` supplies reduced coefficients `y_i` (k-vector) with covariance `S_i`. The model is `y_i ~ N(X_i beta, S_i + Psi)`, `X_i = x_i' (x) I_k`, with `Psi` the between-location covariance. Writing `Sigma_i = S_i + Psi` with Cholesky factor `U_i` (`Sigma_i = U_i' U_i`), the GLS estimate is `beta = (sum X_i' Sigma_i^-1 X_i)^-1 sum X_i' Sigma_i^-1 y_i`, computed by stacking `U_i^-T X_i` and `U_i^-T y_i` and solving least squares. The REML log-likelihood is

    l_R(Psi) = c - sum_i log|U_i| - 1/2 log|sum_i X_i' Sigma_i^-1 X_i| - 1/2 sum_i r_i' Sigma_i^-1 r_i

with `r_i = y_i - X_i beta(Psi)` and `c = -(n - p k)/2 log(2 pi) + log|chol(sum X_i'X_i)|`; ML drops the middle term and uses `c = -n/2 log(2 pi)`. `Psi` is parameterised through its Cholesky factor (unstructured), log variances (diagonal) or one log variance (identity), and maximised numerically. BLUPs: `y_i^ = X_i beta + Psi Sigma_i^-1 (y_i - X_i beta)` with covariance `X_i V_beta X_i' + Psi - Psi Sigma_i^-1 Psi`. Cochran's Q is the GLS residual sum of squares with `Psi = 0`, on `n - pk` degrees of freedom; `I^2 = max((Q - df)/Q, 0)`.

## 12. Penalised DLNMs (Gasparrini et al. 2017; Wood 2011)

For `ps`/`cr` bases the penalty on the cross-basis coefficients is `sum_k lambda_k eta' S_k eta` with `S_var (x) I` and `I (x) S_lag` (each rescaled by its largest eigenvalue), plus optional extra lag penalties. For fixed `lambda`, penalised IRLS solves `(X'WX + S_lambda) eta = X'W z` to convergence of the penalised deviance. The smoothing parameters minimise the Laplace-approximate negative REML

    V(rho) = D_p/(2 phi) - l_s(phi) + 1/2 log|X'WX + S_rho| - 1/2 log|S_rho|_+ - M_p/2 log(2 pi phi)

with `D_p` the penalised deviance, `l_s` the saturated log-likelihood, `|.|_+` the generalised determinant over the fixed rank of the total penalty, and `M_p` the null-space dimension; for quasi families `log phi` is optimised jointly and `l_s` takes the extended quasi-likelihood form. ML replaces `log|X'WX + S_rho|` by the log determinant of the Hessian projected onto the range space of the penalty and drops the `M_p` term. The posterior covariance is `(X'WX + S)^-1 phi`, effective degrees of freedom are the diagonal of `(X'WX + S)^-1 X'WX`, and the reported scale for quasi families is the Fletcher-corrected Pearson estimate, all as in `mgcv`.
