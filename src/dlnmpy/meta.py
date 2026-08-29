"""Multivariate meta-analysis and two-stage designs.

The second stage of a multi-location DLNM analysis pools the reduced
coefficients from each location (see :func:`dlnmpy.crossreduce`) with a
multivariate random-effects meta-analysis (Gasparrini, Armstrong & Kenward
2012, *Statistics in Medicine*; Sera, Armstrong, Blangiardo & Gasparrini
2019 for the extended framework). This module implements the model fitted
by the R package ``mixmeta`` (Gasparrini) for a single random level:

    y_i ~ N(X_i beta, S_i + Psi),   X_i = x_i' (x) I_k

with ``y_i`` the k outcomes of location i, ``S_i`` their within-location
covariance, ``x_i`` optional meta-predictors, and ``Psi`` the between-location
covariance. Estimation is by REML (default), ML or fixed effects, with an
unstructured, diagonal or identity ``Psi``. Best linear unbiased predictions
(BLUPs) of the location-specific coefficients, Cochran's Q, I² and
predictions at new covariate values follow ``mixmeta``'s formulas, and the
results are validated against it in ``tests/test_meta.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize
from scipy.stats import chi2, norm

from .predict import CrossPred, crosspred
from .core import CrossBasis, OneBasis

__all__ = ["mixmeta", "MixMeta", "vech", "xpnd", "predict_reduced", "stack_reduced"]


# ----------------------------------------------------------------------------
def vech(M: np.ndarray) -> np.ndarray:
    """Lower-triangular (column-major) half of a symmetric matrix, as R's
    ``vechMat``."""
    M = np.asarray(M)
    k = M.shape[0]
    return np.concatenate([M[i:, i] for i in range(k)])


def xpnd(v: np.ndarray) -> np.ndarray:
    """Inverse of :func:`vech`."""
    v = np.asarray(v, dtype=float).ravel()
    k = int(round((np.sqrt(1 + 8 * v.size) - 1) / 2))
    M = np.zeros((k, k))
    pos = 0
    for i in range(k):
        n = k - i
        M[i:, i] = v[pos:pos + n]
        pos += n
    return M + np.tril(M, -1).T


def _as_S_list(S, m: int, k: int) -> list[np.ndarray]:
    S = np.asarray(S, dtype=float)
    if S.ndim == 3:
        return [S[i] for i in range(m)]
    if S.ndim == 2 and S.shape == (m, k * (k + 1) // 2):
        return [xpnd(S[i]) for i in range(m)]
    if S.ndim == 2 and S.shape == (m, k):  # variances only
        return [np.diag(S[i]) for i in range(m)]
    if S.ndim == 1 and k == 1 and S.size == m:
        return [np.array([[s]]) for s in S]
    raise ValueError("'S' must be (m, k, k), (m, k*(k+1)/2) vech rows, (m, k) variances or a length-m vector")


# ----------------------------------------------------------------------------
class _Model:
    """Working quantities for the likelihoods (single random level)."""

    def __init__(self, y, S, X, bscov):
        self.y = np.asarray(y, dtype=float)
        if self.y.ndim == 1:
            self.y = self.y[:, None]
        self.m, self.k = self.y.shape
        self.Slist = _as_S_list(S, self.m, self.k)
        X = np.ones((self.m, 1)) if X is None else np.asarray(X, dtype=float)
        if X.ndim == 1:
            X = X[:, None]
        self.X = X
        self.p = X.shape[1]
        self.Xlist = [np.kron(X[i:i + 1, :], np.eye(self.k)) for i in range(self.m)]
        self.ylist = [self.y[i] for i in range(self.m)]
        self.nall = self.m * self.k
        self.bscov = bscov
        self.npar = {"unstr": self.k * (self.k + 1) // 2, "diag": self.k, "id": 1}[bscov]
        XtX = sum(Xi.T @ Xi for Xi in self.Xlist)
        self.const_reml = (-0.5 * (self.nall - self.p * self.k) * np.log(2 * np.pi)
                           + np.sum(np.log(np.diag(np.linalg.cholesky(XtX)))))
        self.const_ml = -0.5 * self.nall * np.log(2 * np.pi)

    # Psi parametrisation (as mixmeta: Cholesky for unstr, log-variances otherwise)
    def par2psi(self, par):
        k = self.k
        if self.bscov == "unstr":
            L = np.zeros((k, k))
            L[np.tril_indices(k)] = _vech_to_tril_order(par, k)
            return L @ L.T
        if self.bscov == "diag":
            return np.diag(np.exp(par))
        return np.eye(k) * np.exp(par[0])

    def psi2par(self, Psi):
        if self.bscov == "unstr":
            L = np.linalg.cholesky(Psi)
            return vech(L)
        if self.bscov == "diag":
            return np.log(np.diag(Psi))
        return np.array([np.log(np.mean(np.diag(Psi)))])

    def gls(self, Psi):
        Ulist, iUXl, iUyl = [], [], []
        for Xi, yi, Si in zip(self.Xlist, self.ylist, self.Slist):
            U = np.linalg.cholesky(Si + Psi).T  # upper, as R chol()
            Ulist.append(U)
            invU = np.linalg.solve(U, np.eye(U.shape[0]))
            iUXl.append(invU.T @ Xi)
            iUyl.append(invU.T @ yi)
        iUX = np.vstack(iUXl)
        iUy = np.concatenate(iUyl)
        coef, *_ = np.linalg.lstsq(iUX, iUy, rcond=None)
        return coef, Ulist, iUXl, iUX, iUy

    def loglik(self, par, method):
        Psi = self.par2psi(par)
        try:
            coef, Ulist, iUXl, iUX, iUy = self.gls(Psi)
        except np.linalg.LinAlgError:
            return -np.inf
        r = iUy - iUX @ coef
        res = -0.5 * float(r @ r)
        det1 = -sum(np.sum(np.log(np.diag(U))) for U in Ulist)
        if method == "ml":
            return self.const_ml + det1 + res
        tXWX = sum(A.T @ A for A in iUXl)
        det2 = -np.sum(np.log(np.diag(np.linalg.cholesky(tXWX))))
        return self.const_reml + det1 + det2 + res


def _vech_to_tril_order(par, k):
    """Map a column-major vech vector onto numpy's row-major tril_indices order."""
    L = np.zeros((k, k))
    pos = 0
    for j in range(k):
        n = k - j
        L[j:, j] = par[pos:pos + n]
        pos += n
    return L[np.tril_indices(k)]


# ----------------------------------------------------------------------------
@dataclass
class MixMeta:
    """Fitted multivariate meta-analysis / meta-regression.

    ``coef`` is a (p, k) matrix (rows: meta-predictors, columns: outcomes);
    ``coef_vec`` the flattened vector in ``mixmeta`` order (predictor-major);
    ``vcov`` its (pk, pk) covariance; ``Psi`` the (k, k) between-study
    covariance (zeros for fixed effects).
    """

    coef: np.ndarray
    vcov: np.ndarray
    Psi: np.ndarray
    method: str
    bscov: str
    loglik: float
    converged: bool
    y: np.ndarray
    S: list
    X: np.ndarray
    df_residual: int
    niter: int = 0

    @property
    def k(self) -> int:
        return self.coef.shape[1]

    @property
    def p(self) -> int:
        return self.coef.shape[0]

    @property
    def m(self) -> int:
        return self.y.shape[0]

    @property
    def coef_vec(self) -> np.ndarray:
        return self.coef.ravel()

    @property
    def npar(self) -> int:
        return self.p * self.k + (0 if self.method == "fixed" else {"unstr": self.k * (self.k + 1) // 2, "diag": self.k, "id": 1}[self.bscov])

    @property
    def aic(self) -> float:
        return -2 * self.loglik + 2 * self.npar

    @property
    def bic(self) -> float:
        n = self.m * self.k - (self.p * self.k if self.method == "reml" else 0)
        return -2 * self.loglik + np.log(n) * self.npar

    # --- predictions ----------------------------------------------------------
    def predict(self, X=None, se: bool = False, ci_level: float = 0.95):
        """Predicted outcomes at meta-predictor values ``X`` (default: the
        fitted rows). Returns ``fit`` (n, k) or, with ``se``, a dict with
        ``fit``, ``se``, ``low``, ``high`` and per-row covariance ``vcov``."""
        X = self.X if X is None else np.atleast_2d(np.asarray(X, dtype=float))
        if X.shape[1] != self.p:
            raise ValueError(f"'X' must have {self.p} columns")
        fit = np.zeros((X.shape[0], self.k))
        V = []
        for i in range(X.shape[0]):
            Xi = np.kron(X[i:i + 1], np.eye(self.k))
            fit[i] = Xi @ self.coef_vec
            V.append(Xi @ self.vcov @ Xi.T)
        if not se:
            return fit
        sd = np.sqrt(np.array([np.diag(v) for v in V]))
        z = norm.ppf(1 - (1 - ci_level) / 2)
        return {"fit": fit, "se": sd, "low": fit - z * sd, "high": fit + z * sd, "vcov": V}

    def blup(self, se: bool = False, pi_level: float = 0.95):
        """Best linear unbiased predictions of the location-specific outcomes.
        Returns ``blup`` (m, k) or, with ``se``, a dict with ``blup``, ``se``,
        ``pi_low``, ``pi_high`` and per-location covariance ``vcov``."""
        out = np.zeros((self.m, self.k))
        V = []
        for i in range(self.m):
            Xi = np.kron(self.X[i:i + 1], np.eye(self.k))
            pred = Xi @ self.coef_vec
            Vi = Xi @ self.vcov @ Xi.T
            if self.method != "fixed":
                Sigma = self.Psi + self.S[i]
                PinvS = self.Psi @ np.linalg.inv(Sigma)
                pred = pred + PinvS @ (self.y[i] - Xi @ self.coef_vec)
                Vi = Vi + self.Psi - PinvS @ self.Psi
            out[i] = pred
            V.append(Vi)
        if not se:
            return out
        sd = np.sqrt(np.array([np.diag(v) for v in V]))
        z = norm.ppf(1 - (1 - pi_level) / 2)
        return {"blup": out, "se": sd, "pi_low": out - z * sd, "pi_high": out + z * sd, "vcov": V}

    def qtest(self) -> dict:
        """Cochran's Q test of heterogeneity (multivariate and per outcome)
        and the corresponding I² statistics."""
        mod = _Model(self.y, self.S, self.X, self.bscov)
        coef, Ulist, iUXl, iUX, iUy = mod.gls(np.zeros((self.k, self.k)))
        r = iUy - iUX @ coef
        Q = [float(r @ r)]
        df = [self.m * self.k - self.p * self.k]
        if self.k > 1:
            cm = coef.reshape(self.p, self.k)
            for j in range(self.k):
                res = self.y[:, j] - self.X @ cm[:, j]
                Q.append(float(np.sum(res ** 2 / np.array([S[j, j] for S in self.S]))))
                df.append(self.m - self.p)
        Q, df = np.array(Q), np.array(df)
        pval = 1 - chi2.cdf(Q, df)
        i2 = np.maximum((Q - df) / Q * 100, 0)
        return {"Q": Q, "df": df, "pvalue": pval, "I2": i2}

    def summary(self) -> str:
        q = self.qtest()
        lines = [f"Multivariate {'meta-regression' if self.p > 1 else 'meta-analysis'}",
                 f"Dimension: {self.k}", f"Estimation method: {self.method.upper()}",
                 f"Studies: {self.m}, converged: {self.converged}",
                 f"logLik: {self.loglik:.4f}   AIC: {self.aic:.4f}   BIC: {self.bic:.4f}", "",
                 "Fixed-effects coefficients (rows: predictors, columns: outcomes):",
                 np.array2string(self.coef, precision=4), "",
                 f"Between-study covariance Psi ({self.bscov}):", np.array2string(self.Psi, precision=6), "",
                 f"Cochran Q: {q['Q'][0]:.3f} on {q['df'][0]} df (p = {q['pvalue'][0]:.4f}), I2 = {q['I2'][0]:.1f}%"]
        return "\n".join(lines)

    def __repr__(self):
        return f"MixMeta(method={self.method!r}, k={self.k}, p={self.p}, m={self.m}, converged={self.converged})"


def mixmeta(y, S, X=None, method: str = "reml", bscov: str = "unstr", init_psi=None,
            maxiter: int = 500, tol: float = 1e-10) -> MixMeta:
    """Fit a multivariate meta-analysis or meta-regression.

    Parameters
    ----------
    y : (m, k) array
        Outcomes (e.g. reduced coefficients) for m studies.
    S : (m, k, k) array, (m, k(k+1)/2) vech rows, (m, k) variances or list
        Within-study covariance matrices.
    X : (m, p) array, optional
        Meta-predictors, including a column of ones if an intercept is
        wanted. Default: intercept only.
    method : {"reml", "ml", "fixed"}
    bscov : {"unstr", "diag", "id"}
        Structure of the between-study covariance.
    """
    method = method.lower()
    if method not in ("reml", "ml", "fixed"):
        raise ValueError("'method' must be 'reml', 'ml' or 'fixed'")
    if bscov not in ("unstr", "diag", "id"):
        raise ValueError("'bscov' must be 'unstr', 'diag' or 'id'")
    mod = _Model(y, S, X, bscov)
    k, p = mod.k, mod.p

    if method == "fixed":
        Psi = np.zeros((k, k))
        coef, Ulist, iUXl, iUX, iUy = mod.gls(Psi)
        r = iUy - iUX @ coef
        ll = mod.const_ml - sum(np.sum(np.log(np.diag(U))) for U in Ulist) - 0.5 * float(r @ r)
        converged, niter = True, 0
    else:
        Psi0 = np.eye(k) * 0.001 if init_psi is None else np.asarray(init_psi, dtype=float)
        par0 = mod.psi2par(Psi0)

        def nll(par):
            v = mod.loglik(par, method)
            return 1e100 if not np.isfinite(v) else -v

        best = None
        for start in (par0, mod.psi2par(np.eye(k) * 0.1), mod.psi2par(_moment_psi(mod))):
            try:
                opt = optimize.minimize(nll, start, method="BFGS", options={"gtol": 1e-8, "maxiter": maxiter})
                opt = optimize.minimize(nll, opt.x, method="Nelder-Mead",
                                        options={"xatol": tol, "fatol": tol, "maxiter": 20000})
                opt = optimize.minimize(nll, opt.x, method="BFGS", options={"gtol": 1e-10, "maxiter": maxiter})
            except (np.linalg.LinAlgError, ValueError):
                continue
            if best is None or opt.fun < best.fun:
                best = opt
        if best is None:
            raise RuntimeError("meta-analysis did not converge")
        Psi = mod.par2psi(best.x)
        coef, Ulist, iUXl, iUX, iUy = mod.gls(Psi)
        ll = -float(best.fun)
        converged = bool(np.isfinite(ll))
        niter = int(best.nit)

    R = np.linalg.qr(iUX, mode="r")
    Rinv = np.linalg.solve(R, np.eye(R.shape[0]))
    vcov = Rinv @ Rinv.T
    df_res = mod.nall - p * k - (0 if method == "fixed" else mod.npar)
    return MixMeta(coef=coef.reshape(p, k), vcov=vcov, Psi=Psi, method=method, bscov=bscov,
                   loglik=ll, converged=converged, y=mod.y, S=mod.Slist, X=mod.X,
                   df_residual=df_res, niter=niter)


def _moment_psi(mod: _Model) -> np.ndarray:
    """Crude method-of-moments starting value for Psi."""
    coef, *_ = mod.gls(np.zeros((mod.k, mod.k)))
    res = np.array([mod.ylist[i] - mod.Xlist[i] @ coef for i in range(mod.m)])
    Sbar = sum(mod.Slist) / mod.m
    Psi = np.cov(res.T, ddof=1).reshape(mod.k, mod.k) - Sbar
    w, U = np.linalg.eigh((Psi + Psi.T) / 2)
    w = np.maximum(w, 1e-4)
    return U @ np.diag(w) @ U.T


# ----------------------------------------------------------------------------
# two-stage helpers
# ----------------------------------------------------------------------------
def stack_reduced(reduced: list) -> tuple[np.ndarray, list]:
    """Stack the coefficients and covariances of a list of ``CrossReduce``
    objects (one per location) into the ``y`` and ``S`` inputs of
    :func:`mixmeta`."""
    y = np.vstack([r.coef for r in reduced])
    S = [r.vcov for r in reduced]
    return y, S


def predict_reduced(basis, coef, vcov, at=None, from_=None, to=None, by=None, cen=None,
                    model_link: str = "log", ci_level: float = 0.95) -> CrossPred:
    """Exposure-response curve from reduced (one-dimensional) coefficients,
    e.g. pooled or BLUP coefficients from :func:`mixmeta`.

    ``basis`` is a ``CrossBasis`` whose predictor-space specification
    (``argvar``) defines the one-dimensional basis, or a ``OneBasis``.
    """
    if isinstance(basis, CrossBasis):
        av = dict(basis.argvar)
        fun = av.pop("fun")
        av.pop("cen", None)
        ob = OneBasis(matrix=np.zeros((1, basis.df[0])), fun=fun, attrs=av, range=basis.range)
    elif isinstance(basis, OneBasis):
        ob = basis
    else:
        raise TypeError("'basis' must be a CrossBasis or OneBasis")
    return crosspred(ob, coef=coef, vcov=vcov, model_link=model_link, at=at, from_=from_,
                     to=to, by=by, cen=cen, ci_level=ci_level)
