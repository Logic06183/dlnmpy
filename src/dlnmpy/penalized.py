"""Penalised GLMs with smoothing parameters chosen by REML or ML.

R users fit penalised DLNMs (Gasparrini et al. 2017) with ``mgcv::gam`` and
``paraPen`` penalties from ``cbPen``. Python has no mgcv, so this module
provides the piece that was missing: a penalised iteratively reweighted
least squares fitter whose smoothing parameters (and, for quasi families,
the scale) are chosen by maximising the Laplace-approximate restricted
marginal likelihood of Wood (2011), the criterion ``gam(method="REML")``
uses. The criterion is implemented exactly as in ``mgcv::gam.fit3``:

    V(rho) = D_p/(2 phi) - l_s(phi) + log|X'WX + S_rho|/2 - log|S_rho|_+/2
             - M_p log(2 pi phi)/2

with ``D_p`` the penalised deviance, ``l_s`` the saturated log-likelihood,
``S_rho = sum_k exp(rho_k) S_k``, ``|.|_+`` the generalised determinant and
``M_p`` the dimension of the penalty null space. Unlike mgcv it uses a
generic quasi-Newton optimiser on ``rho`` (numerical derivatives), which is
adequate for the two or three smoothing parameters of a cross-basis; the
result is validated against ``gam`` in ``tests/test_penalized.py``.

The fitted object exposes ``params`` (named), ``cov_params()`` (Bayesian
posterior covariance ``(X'WX + S)^-1 phi``), ``edf``, ``sp``, ``scale`` and
``link`` so that :func:`dlnmpy.crosspred` and friends accept it like a
statsmodels result.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import optimize
from scipy.special import gammaln

__all__ = ["fit_pglm", "fit_pgam", "PenalizedGLMResults"]


# ----------------------------------------------------------------------------
class _Family:
    def __init__(self, name: str):
        name = name.lower()
        self.name = name
        if name in ("poisson", "quasipoisson"):
            self.link, self.scale_known = "log", name == "poisson"
        elif name in ("gaussian", "normal"):
            self.link, self.scale_known = "identity", False
        elif name in ("binomial", "quasibinomial"):
            self.link, self.scale_known = "logit", name == "binomial"
        else:
            raise ValueError(f"unknown family '{name}'")

    def linkinv(self, eta):
        if self.link == "log":
            return np.exp(eta)
        if self.link == "identity":
            return eta
        return 1 / (1 + np.exp(-eta))

    def mu_eta(self, eta):  # d mu / d eta
        if self.link == "log":
            return np.exp(eta)
        if self.link == "identity":
            return np.ones_like(eta)
        p = 1 / (1 + np.exp(-eta))
        return p * (1 - p)

    def variance(self, mu):
        if self.link == "log":
            return mu
        if self.link == "identity":
            return np.ones_like(mu)
        return mu * (1 - mu)

    def dev_resids(self, y, mu, w):
        if self.link == "log":
            r = np.where(y > 0, y * np.log(np.where(y > 0, y / mu, 1.0)), 0.0) - (y - mu)
            return 2 * w * r
        if self.link == "identity":
            return w * (y - mu) ** 2
        y1 = np.where(y > 0, y * np.log(np.where(y > 0, y / mu, 1.0)), 0.0)
        y0 = np.where(y < 1, (1 - y) * np.log(np.where(y < 1, (1 - y) / (1 - mu), 1.0)), 0.0)
        return 2 * w * (y1 + y0)

    def initial_mu(self, y):
        if self.link == "log":
            return y + 0.1
        if self.link == "identity":
            return y.astype(float)
        return (y + 0.5) / 2

    def ls(self, y, w, scale):
        """Saturated log-likelihood (mgcv's ``family$ls``)."""
        nobs = int(np.sum(w > 0))
        if self.name == "poisson":
            return float(np.sum(w * (np.where(y > 0, y * np.log(np.where(y > 0, y, 1)) - y, 0.0) - gammaln(y + 1))))
        if self.name == "gaussian":
            return -nobs * np.log(2 * np.pi * scale) / 2 + np.sum(np.log(w[w > 0])) / 2
        if self.name == "binomial":
            y1 = np.where(y > 0, y * np.log(np.where(y > 0, y, 1)), 0.0)
            y0 = np.where(y < 1, (1 - y) * np.log(np.where(y < 1, 1 - y, 1)), 0.0)
            return float(np.sum(w * (y1 + y0)))
        # quasi families: extended quasi-likelihood form
        return -nobs * np.log(scale) / 2 + np.sum(np.log(w[w > 0])) / 2


# ----------------------------------------------------------------------------
@dataclass
class PenalizedGLMResults:
    """Result of :func:`fit_pglm` / :func:`fit_pgam`."""

    params: pd.Series
    vcov: np.ndarray
    sp: np.ndarray
    scale: float
    edf: np.ndarray
    deviance: float
    reml: float
    method: str
    family: str
    link: str
    converged: bool
    fitted_values: np.ndarray
    linear_predictor: np.ndarray
    nobs: int
    reml_scale: float = float("nan")
    penalty_names: list = field(default_factory=list)
    exog_names: list = field(default_factory=list)
    endog: np.ndarray | None = field(default=None, repr=False)
    exog: np.ndarray | None = field(default=None, repr=False)

    def cov_params(self) -> pd.DataFrame:
        return pd.DataFrame(self.vcov, index=self.params.index, columns=self.params.index)

    @property
    def bse(self) -> pd.Series:
        return pd.Series(np.sqrt(np.diag(self.vcov)), index=self.params.index)

    @property
    def edf_total(self) -> float:
        return float(np.sum(self.edf))

    @property
    def fittedvalues(self) -> np.ndarray:
        """Alias of ``fitted_values`` (statsmodels spelling), so helpers such as
        :func:`dlnmpy.qaic` accept either kind of result."""
        return self.fitted_values

    def edf_by(self, prefix: str) -> float:
        """Effective degrees of freedom of the basis named ``prefix`` (the
        ``name`` given to ``to_dataframe``)."""
        idx = _columns_of(list(self.params.index), prefix)
        return float(np.sum(self.edf[idx]))

    def summary(self) -> str:
        lines = [f"Penalised GLM ({self.family}, {self.link} link), smoothing by {self.method.upper()}",
                 f"observations: {self.nobs}   coefficients: {self.params.size}   total edf: {self.edf_total:.2f}",
                 "smoothing parameters: " + ", ".join(f"{n}={s:.4g}" for n, s in zip(self.penalty_names, self.sp)),
                 f"scale: {self.scale:.4f}   deviance: {self.deviance:.2f}   {self.method.upper()} score: {self.reml:.4f}",
                 f"converged: {self.converged}"]
        return "\n".join(lines)

    def __repr__(self):
        return (f"PenalizedGLMResults(family={self.family!r}, method={self.method!r}, "
                f"sp={np.round(self.sp, 4).tolist()}, edf={self.edf_total:.2f}, scale={self.scale:.4f})")


def _columns_of(names: list, prefix: str) -> list:
    """Indices of the design columns that belong to the basis called
    ``prefix``: the name followed by a column label of the form ``v1_l1`` /
    ``v1.l1`` (cross-basis) or ``b1`` (one-basis), as ``to_dataframe`` writes
    them. A bare ``startswith`` also picks up the columns of a second basis
    whose name extends the first (``cb`` and ``cb_o3``), and the penalty or the
    edf would silently be applied to both. Same rule as ``model._match``."""
    import re
    pat = re.compile(rf"^{re.escape(prefix)}[._]?(v\d+[._]l\d+|b\d+)$")
    return [i for i, n in enumerate(names) if n == prefix or pat.match(n)]


# ----------------------------------------------------------------------------
def _pirls(y, X, S, fam, w, offset, beta0=None, maxit=100, tol=1e-10):
    """Penalised IRLS for a fixed total penalty ``S``. Returns beta, W, dev."""
    n, p = X.shape
    if beta0 is None:
        mu = fam.initial_mu(y)
        eta = np.log(mu) if fam.link == "log" else (mu if fam.link == "identity" else np.log(mu / (1 - mu)))
    else:
        eta = X @ beta0 + offset
        mu = fam.linkinv(eta)
    dev_old = np.inf
    beta = beta0
    for it in range(maxit):
        mu_eta = fam.mu_eta(eta)
        var = fam.variance(mu)
        z = (eta - offset) + (y - mu) / mu_eta
        W = w * mu_eta ** 2 / var
        XtW = X.T * W
        A = XtW @ X + S
        beta_new = np.linalg.solve(A, XtW @ z)
        eta = X @ beta_new + offset
        mu = fam.linkinv(eta)
        dev = float(np.sum(fam.dev_resids(y, mu, w)))
        pen = float(beta_new @ S @ beta_new)
        beta = beta_new
        if abs(dev + pen - dev_old) < tol * (abs(dev + pen) + 0.1):
            break
        dev_old = dev + pen
    mu_eta = fam.mu_eta(eta)
    W = w * mu_eta ** 2 / fam.variance(mu)
    return beta, W, dev, mu, eta


def _logdet_plus_fixed_rank(S, rank):
    """Log generalised determinant of ``S`` over its (known) rank: the sum of
    the logs of the ``rank`` largest eigenvalues. Using a fixed rank keeps the
    criterion continuous as a smoothing parameter tends to zero."""
    ev = np.sort(np.linalg.eigvalsh((S + S.T) / 2))[::-1][:rank]
    return float(np.sum(np.log(np.maximum(ev, 1e-300))))


def _logdet_plus(S, tol=1e-10):
    ev = np.linalg.eigvalsh((S + S.T) / 2)
    ev = ev[ev > tol * max(ev.max(), 1e-300)]
    return float(np.sum(np.log(ev))), ev.size


def fit_pglm(y, X, penalties, family: str = "poisson", method: str = "reml", sp=None,
             offset=None, weights=None, exog_names=None, penalty_names=None,
             maxiter: int = 200) -> PenalizedGLMResults:
    """Fit a penalised GLM ``g(E y) = X beta`` with penalty ``sum_k sp_k beta' S_k beta``.

    Parameters
    ----------
    y : (n,) array
    X : (n, p) array or DataFrame (column names are used for ``params``)
    penalties : list of (p, p) arrays
        Penalty matrices, zero-padded to the full coefficient vector.
    family : {"poisson", "quasipoisson", "gaussian", "binomial", "quasibinomial"}
    method : {"reml", "ml"}
        Criterion for the smoothing parameters (ignored if ``sp`` is given).
    sp : array_like, optional
        Fixed smoothing parameters; if given no outer optimisation is done.
    """
    if isinstance(X, pd.DataFrame):
        exog_names = list(X.columns) if exog_names is None else exog_names
        X = X.to_numpy(dtype=float)
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).ravel()
    n, p = X.shape
    w = np.ones(n) if weights is None else np.asarray(weights, dtype=float)
    off = np.zeros(n) if offset is None else np.asarray(offset, dtype=float)
    fam = _Family(family)
    method = method.lower()
    if method not in ("reml", "ml"):
        raise ValueError("'method' must be 'reml' or 'ml'")
    reml_ind = 1 if method == "reml" else 0
    Slist = [np.asarray(S, dtype=float) for S in penalties]
    for S in Slist:
        if S.shape != (p, p):
            raise ValueError("each penalty must be a (p, p) matrix zero-padded to the full design")
    exog_names = exog_names or [f"x{i}" for i in range(p)]
    penalty_names = penalty_names or [f"S{i + 1}" for i in range(len(Slist))]
    nsp = len(Slist)
    # null-space dimension of the total penalty, and a basis of the null space
    # (unpenalised directions) used by the ML criterion
    if nsp:
        Stot = sum(Slist)
        _, rank_total = _logdet_plus(Stot)
        ev, U = np.linalg.eigh((Stot + Stot.T) / 2)
        Ur = U[:, ev > 1e-10 * max(ev.max(), 1e-300)]
    else:
        rank_total, Ur = 0, np.zeros((p, 0))
    Mp = p - rank_total

    state = {"beta": None}

    def fit_at(rho, log_phi=None):
        S = sum(np.exp(r) * Sk for r, Sk in zip(rho, Slist)) if nsp else np.zeros((p, p))
        beta, W, dev, mu, eta = _pirls(y, X, S, fam, w, off, beta0=state["beta"])
        state["beta"] = beta
        pen = float(beta @ S @ beta)
        XtWX = (X.T * W) @ X
        A = XtWX + S
        with np.errstate(divide="ignore", over="ignore", invalid="ignore", under="ignore"):
            _, logdetA = np.linalg.slogdet(A)
        if method == "ml" and Ur.shape[1]:
            # ML integrates out only the penalised coefficients: the log
            # determinant is that of the Hessian projected onto the range space
            # of the penalty (mgcv's MLpenalty1)
            with np.errstate(divide="ignore", over="ignore", invalid="ignore", under="ignore"):
                logdetA = np.linalg.slogdet(Ur.T @ A @ Ur)[1]
        logdetS = _logdet_plus_fixed_rank(S, rank_total) if nsp else 0.0
        if fam.scale_known:
            phi = 1.0
        elif log_phi is not None:
            phi = float(np.exp(log_phi))
        else:
            phi = float(np.sum(w * (y - mu) ** 2 / fam.variance(mu)) / (n - Mp))
        Dp = dev + pen
        crit = (Dp / (2 * phi) - fam.ls(y, w, phi)) + logdetA / 2 - logdetS / 2 - reml_ind * Mp / 2 * np.log(2 * np.pi * phi)
        return crit, dict(beta=beta, W=W, dev=dev, mu=mu, eta=eta, S=S, A=A, phi=phi, pen=pen)

    if sp is not None:
        rho = np.log(np.asarray(sp, dtype=float))
        if fam.scale_known:
            crit, out = fit_at(rho)
        else:
            phi0 = fit_at(rho)[1]["phi"]
            res = optimize.minimize_scalar(lambda lp: fit_at(rho, lp)[0], bracket=(np.log(phi0) - 1, np.log(phi0) + 1))
            crit, out = fit_at(rho, res.x)
        converged = True
    elif nsp == 0:
        crit, out = fit_at(np.array([]))
        rho, converged = np.array([]), True
    else:
        # starting values: scale each penalty so that sp*S has trace comparable to X'WX
        _, W0, *_ = _pirls(y, X, np.zeros((p, p)), fam, w, off)
        XtWX = (X.T * W0) @ X
        rho0 = np.array([np.log(np.trace(XtWX) / max(np.trace(Sk), 1e-12)) for Sk in Slist])
        if fam.scale_known:
            def obj(r):
                v, _ = fit_at(r)
                return v if np.isfinite(v) else 1e100
            x0 = rho0
        else:
            phi0 = fit_at(rho0)[1]["phi"]

            def obj(r):
                v, _ = fit_at(r[:-1], r[-1])
                return v if np.isfinite(v) else 1e100
            x0 = np.r_[rho0, np.log(phi0)]
        res = optimize.minimize(obj, x0, method="BFGS", options={"gtol": 1e-6, "maxiter": maxiter, "eps": 1e-5})
        res = optimize.minimize(obj, res.x, method="Nelder-Mead",
                                options={"xatol": 1e-5, "fatol": 1e-9, "maxiter": 2000, "initial_simplex": None})
        converged = bool(res.success or np.isfinite(res.fun))
        if fam.scale_known:
            rho = res.x
            crit, out = fit_at(rho)
        else:
            rho = res.x[:-1]
            crit, out = fit_at(rho, res.x[-1])

    A, W, beta, phi = out["A"], out["W"], out["beta"], out["phi"]
    Ainv = np.linalg.inv(A)
    F = Ainv @ ((X.T * W) @ X)
    edf = np.diag(F)
    if not fam.scale_known:
        # mgcv reports (and scales Vp by) the Fletcher (2012) corrected Pearson
        # estimate rather than the REML scale parameter
        mu = out["mu"]
        pearson = float(np.sum(w * (y - mu) ** 2 / fam.variance(mu)))
        scale_est = pearson / (n - float(np.sum(edf)))
        dvar = np.ones_like(mu) if fam.link == "log" else (np.zeros_like(mu) if fam.link == "identity" else 1 - 2 * mu)
        s_bar = max(-0.9, float(np.mean(dvar * (y - mu) / fam.variance(mu))))
        if np.isfinite(s_bar):
            scale_est = scale_est / (1 + s_bar)
        reml_scale, phi = phi, scale_est
    else:
        reml_scale = phi
    return PenalizedGLMResults(
        params=pd.Series(beta, index=exog_names), vcov=Ainv * phi, sp=np.exp(rho), scale=phi, reml_scale=reml_scale,
        edf=edf, deviance=out["dev"], reml=float(crit), method=method, family=fam.name, link=fam.link,
        converged=converged, fitted_values=out["mu"], linear_predictor=out["eta"], nobs=n,
        penalty_names=list(penalty_names), exog_names=list(exog_names), endog=y, exog=X)


# ----------------------------------------------------------------------------
def fit_pgam(formula: str, data, penalties: dict, family: str = "quasipoisson",
             method: str = "reml", sp=None, **kwargs) -> PenalizedGLMResults:
    """Penalised GLM from a formula, the Python counterpart of
    ``mgcv::gam(formula, family, paraPen=list(cb=cbPen(cb)), method="REML")``.

    ``penalties`` maps a design-column prefix (the ``name`` given to
    ``CrossBasis.to_dataframe``) to the output of :func:`dlnmpy.cbpen` or to a
    list of penalty matrices for those columns. Rows with missing values are
    dropped, as R does.
    """
    import statsmodels.formula.api as smf
    import statsmodels.api as sm

    fam = _Family(family)
    smfam = {"log": sm.families.Poisson(), "identity": sm.families.Gaussian(), "logit": sm.families.Binomial()}[fam.link]
    model = smf.glm(formula, data=data, family=smfam, missing="drop")
    X = pd.DataFrame(model.exog, columns=model.exog_names)
    y = model.endog
    # an offset or weights given for every row of ``data`` must follow the
    # rows that the NaN drop removed
    dropped = np.asarray(getattr(model.data, "missing_row_idx", []), dtype=int)
    if dropped.size:
        keep = np.setdiff1d(np.arange(len(data)), dropped)
        for key in ("offset", "weights"):
            v = kwargs.get(key)
            if v is not None and np.size(v) == len(data):
                kwargs[key] = np.asarray(v, dtype=float)[keep]
    p = X.shape[1]
    Slist, names = [], []
    for prefix, pen in penalties.items():
        cols = _columns_of(list(X.columns), prefix)
        mats = [v for k, v in pen.items() if isinstance(v, np.ndarray) and k not in ("rank", "sp")] if isinstance(pen, dict) else list(pen)
        labels = [k for k, v in pen.items() if isinstance(v, np.ndarray) and k not in ("rank", "sp")] if isinstance(pen, dict) else [f"{prefix}{i + 1}" for i in range(len(mats))]
        for lab, M in zip(labels, mats):
            M = np.asarray(M, dtype=float)
            if M.shape != (len(cols), len(cols)):
                raise ValueError(f"penalty '{lab}' is {M.shape} but '{prefix}' has {len(cols)} columns")
            full = np.zeros((p, p))
            full[np.ix_(cols, cols)] = M
            Slist.append(full)
            names.append(f"{prefix}:{lab}")
    return fit_pglm(y, X, Slist, family=family, method=method, sp=sp, penalty_names=names, **kwargs)
