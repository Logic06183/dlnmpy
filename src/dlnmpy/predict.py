"""Prediction from a fitted DLNM: ``crosspred`` and ``crossreduce``.

Both functions take a basis object and either a fitted model or the pair
``coef``/``vcov`` for the basis columns, and return the estimated
exposure-lag-response association on a grid of predictor values and lags.

``crosspred`` follows ``crosspred.R``: it rebuilds the (centred) basis
matrices on the prediction grid, forms the tensor product ``Xpred`` and
computes ``Xpred @ coef`` with standard errors from ``vcov``. Lag-specific
effects are returned as matrices (predictor x lag); the overall cumulative
effect sums ``Xpred`` across integer lags first.

``crossreduce`` follows ``crossreduce.R`` and Gasparrini & Armstrong (2013):
the bi-dimensional parameters are reduced to a one-dimensional basis
through a linear transformation ``M`` (overall cumulative or lag-specific
exposure-response, or predictor-specific lag-response).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.stats import norm

from ._rcompat import median, pretty
from .core import CrossBasis, OneBasis
from .lag import mklag, seqlag
from .model import extract_coef_vcov, get_link

__all__ = ["CrossPred", "CrossReduce", "crosspred", "crossreduce", "mkat", "mkcen"]

_NO_CEN_FUNS = ("thr", "strata", "integer", "lin")


# ----------------------------------------------------------------------------
# helpers shared by crosspred and crossreduce
# ----------------------------------------------------------------------------
def mkat(at, from_, to, by, rng, lag, bylag):
    """Prediction grid for the predictor (port of ``mkat.R``)."""
    if at is None:
        lo = rng[0] if from_ is None else from_
        hi = rng[1] if to is None else to
        nobs = 50 if by is None else max(1, (rng[1] - rng[0]) / by)
        pr = pretty([lo, hi], n=nobs)
        pr = pr[(pr >= lo) & (pr <= hi)]
        if by is None:
            return pr
        # seq(from=min(pretty), to=to, by=by)
        start = float(pr.min())
        n = int(np.floor((hi - start) / by + 1e-10))
        return start + np.arange(n + 1) * by
    at = np.asarray(at, dtype=float)
    if at.ndim == 2:
        if at.shape[1] != lag[1] - lag[0] + 1:
            raise ValueError("matrix in 'at' must have ncol=diff(lag)+1")
        if bylag != 1:
            raise ValueError("'bylag!=1 not allowed with 'at' in matrix form")
        return at
    return np.unique(at.ravel())


def mkcen(cen, basis, rng):
    """Centring value (port of ``mkcen.R``). Returns ``None`` for uncentred."""
    if isinstance(basis, CrossBasis):
        fun = basis.argvar.get("fun")
        stored = basis.argvar.get("cen")
        intercept = basis.argvar.get("intercept", False)
    else:
        fun = basis.fun
        stored = basis.cen
        intercept = basis.attrs.get("intercept", False)
    fname = fun if isinstance(fun, str) else getattr(fun, "__name__", "")
    if cen is None:
        cen = stored
    # np.bool_ is not a subclass of bool, so a numpy boolean would fall through
    # both branches and be coerced by float(cen) into 1.0 or 0.0
    isbool = isinstance(cen, (bool, np.bool_))
    if fname in _NO_CEN_FUNS:
        if isbool:
            cen = None
    else:
        if cen is None or (isbool and cen):
            cen = median(pretty(rng))
        elif isbool:
            cen = None
    if isinstance(intercept, (bool, np.bool_)) and intercept:
        cen = None
    return None if cen is None else float(cen)


def _basis_type(basis) -> str:
    if isinstance(basis, CrossBasis):
        return "cb"
    if isinstance(basis, OneBasis):
        return "one"
    raise TypeError("'basis' must be a CrossBasis or OneBasis object")


def _tensor(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Row-wise tensor product with the columns of ``b`` varying fastest
    (mgcv's ``tensor.prod.model.matrix`` ordering)."""
    return (a[:, :, None] * b[:, None, :]).reshape(a.shape[0], -1)


def mkXpred(basis, at, predvar, predlag, cen):
    """Basis matrix on the prediction grid (port of ``mkXpred.R``)."""
    at = np.asarray(at, dtype=float)
    if at.ndim == 2:
        varvec = at.ravel(order="F")
    else:
        varvec = np.tile(at, len(predlag))
    lagvec = np.repeat(predlag, len(predvar))
    if isinstance(basis, CrossBasis):
        basisvar = basis.basis_var(varvec)
        basislag = basis.basis_lag(lagvec)
        if cen is not None:
            basisvar = basisvar - basis.basis_var([cen])
        return _tensor(basisvar, basislag)
    basisvar = basis.transform(varvec)
    if cen is not None:
        basisvar = basisvar - basis.transform([cen])
    return basisvar


def _se(X: np.ndarray, vcov: np.ndarray) -> np.ndarray:
    return np.sqrt(np.maximum(0, np.sum((X @ vcov) * X, axis=1)))


def _resolve_coef(basis, model, coef, vcov, model_link, name, kind):
    ncol = basis.ncol
    if model is not None:
        coef, vcov = extract_coef_vcov(model, name, kind, ncol)
        model_link = get_link(model, model_link)
        model_class = type(model).__name__
    else:
        if coef is None or vcov is None:
            raise ValueError("At least 'model' or 'coef'-'vcov' must be provided")
        model_class = None
    coef = np.asarray(coef, dtype=float).ravel()
    vcov = np.asarray(vcov, dtype=float)
    if vcov.ndim == 0:
        vcov = vcov.reshape(1, 1)
    elif vcov.ndim == 1 and vcov.size == coef.size:
        vcov = np.diag(vcov)
    if coef.size != ncol or vcov.shape != (ncol, ncol) or np.any(np.isnan(coef)) or np.any(np.isnan(vcov)):
        raise ValueError("coef/vcov not consistent with basis matrix")
    return coef, vcov, model_link, model_class


# ----------------------------------------------------------------------------
@dataclass
class CrossPred:
    """Predictions from a DLNM (see :func:`crosspred`).

    The estimates are on the scale of the linear predictor. ``matfit`` and
    ``matse`` are (n predictor values x n lags) matrices of lag-specific
    effects; ``allfit``/``allse`` are the overall cumulative effects;
    ``cumfit``/``cumse`` (if ``cumul=True``) the incremental cumulative
    effects over integer lags. Exponentiated versions (``allRRfit``,
    ``matRRlow`` ...) are always available; ``is_exp`` says whether the
    model link (log or logit) makes them the natural scale.
    """

    predvar: np.ndarray
    lag: np.ndarray
    bylag: float
    coef: np.ndarray
    vcov: np.ndarray
    matfit: np.ndarray
    matse: np.ndarray
    allfit: np.ndarray
    allse: np.ndarray
    cen: float | None = None
    cumfit: np.ndarray | None = None
    cumse: np.ndarray | None = None
    ci_level: float = 0.95
    model_class: str | None = None
    model_link: str | None = None
    at_matrix: np.ndarray | None = field(default=None, repr=False)

    # --- derived quantities -------------------------------------------------
    @property
    def z(self) -> float:
        return float(norm.ppf(1 - (1 - self.ci_level) / 2))

    @property
    def predlag(self) -> np.ndarray:
        return seqlag(self.lag, self.bylag)

    @property
    def is_exp(self) -> bool:
        return self.model_link in ("log", "logit")

    @property
    def matlow(self): return self.matfit - self.z * self.matse
    @property
    def mathigh(self): return self.matfit + self.z * self.matse
    @property
    def alllow(self): return self.allfit - self.z * self.allse
    @property
    def allhigh(self): return self.allfit + self.z * self.allse
    @property
    def cumlow(self): return None if self.cumfit is None else self.cumfit - self.z * self.cumse
    @property
    def cumhigh(self): return None if self.cumfit is None else self.cumfit + self.z * self.cumse

    @property
    def matRRfit(self): return np.exp(self.matfit)
    @property
    def matRRlow(self): return np.exp(self.matlow)
    @property
    def matRRhigh(self): return np.exp(self.mathigh)
    @property
    def allRRfit(self): return np.exp(self.allfit)
    @property
    def allRRlow(self): return np.exp(self.alllow)
    @property
    def allRRhigh(self): return np.exp(self.allhigh)
    @property
    def cumRRfit(self): return None if self.cumfit is None else np.exp(self.cumfit)
    @property
    def cumRRlow(self): return None if self.cumfit is None else np.exp(self.cumlow)
    @property
    def cumRRhigh(self): return None if self.cumfit is None else np.exp(self.cumhigh)

    # --- convenience --------------------------------------------------------
    def _var_index(self, var) -> int:
        idx = np.nonzero(np.isclose(self.predvar, var))[0]
        if idx.size == 0:
            raise ValueError(f"'var'={var} must match values used for prediction")
        return int(idx[0])

    def _lag_index(self, lag) -> int:
        idx = np.nonzero(np.isclose(self.predlag, lag))[0]
        if idx.size == 0:
            raise ValueError(f"'lag'={lag} must match values used for prediction")
        return int(idx[0])

    def overall(self, exp: bool | None = None):
        """DataFrame of the overall cumulative association."""
        import pandas as pd
        e = self.is_exp if exp is None else exp
        f = np.exp if e else (lambda a: a)
        return pd.DataFrame({"var": self.predvar, "fit": f(self.allfit),
                             "low": f(self.alllow), "high": f(self.allhigh)})

    def slice_var(self, var, exp: bool | None = None, cumul: bool = False):
        """Lag-response curve at predictor value ``var`` (DataFrame)."""
        import pandas as pd
        i = self._var_index(var)
        e = self.is_exp if exp is None else exp
        f = np.exp if e else (lambda a: a)
        if cumul:
            if self.cumfit is None:
                raise ValueError("set cumul=True in crosspred() to get cumulative effects")
            fit, se, lags = self.cumfit[i], self.cumse[i], seqlag(self.lag)
        else:
            fit, se, lags = self.matfit[i], self.matse[i], self.predlag
        return pd.DataFrame({"lag": lags, "fit": f(fit), "low": f(fit - self.z * se),
                             "high": f(fit + self.z * se)})

    def slice_lag(self, lag, exp: bool | None = None):
        """Exposure-response curve at lag ``lag`` (DataFrame)."""
        import pandas as pd
        j = self._lag_index(lag)
        e = self.is_exp if exp is None else exp
        f = np.exp if e else (lambda a: a)
        fit, se = self.matfit[:, j], self.matse[:, j]
        return pd.DataFrame({"var": self.predvar, "fit": f(fit), "low": f(fit - self.z * se),
                             "high": f(fit + self.z * se)})

    def summary(self) -> str:
        lines = ["PREDICTIONS:", f"values: {len(self.predvar)}"]
        if self.cen is not None:
            lines.append(f"centered at: {self.cen:g}")
        lines += [f"range: {np.min(self.predvar):g} , {np.max(self.predvar):g}",
                  f"lag: {self.lag[0]} {self.lag[1]}",
                  f"exponentiated: {'yes' if self.is_exp else 'no'}",
                  f"cumulative: {'yes' if self.cumfit is not None else 'no'}",
                  "", "MODEL:", f"parameters: {self.coef.size}",
                  f"class: {self.model_class}", f"link: {self.model_link}"]
        return "\n".join(lines)

    def plot(self, *args, **kwargs):
        from .plot import plot_crosspred
        return plot_crosspred(self, *args, **kwargs)

    def __repr__(self):
        return (f"CrossPred(values={len(self.predvar)}, lag={tuple(int(v) for v in self.lag)}, bylag={self.bylag}, "
                f"cen={self.cen}, link={self.model_link!r}, cumul={self.cumfit is not None})")


def crosspred(basis, model=None, coef=None, vcov=None, model_link=None, at=None,
              from_=None, to=None, by=None, lag=None, bylag: float = 1,
              cen=None, ci_level: float = 0.95, cumul: bool = False,
              name: str = "cb") -> CrossPred:
    """Predict exposure-lag-response associations from a fitted DLNM.

    Parameters
    ----------
    basis : CrossBasis or OneBasis
    model : fitted model, optional
        A statsmodels results object (or any object with ``.params`` named
        as produced by ``basis.to_dataframe(name)`` and ``.cov_params()``).
    coef, vcov : array_like, optional
        Coefficients and covariance matrix for the basis columns, used
        instead of ``model``.
    model_link : str, optional
        Link function (``"log"``, ``"logit"``, ``"identity"`` ...). Extracted
        from ``model`` when possible.
    at : array_like, optional
        Predictor values (vector), or a matrix of exposure histories.
    from_, to, by : float, optional
        Alternative way to define the grid of predictor values.
    lag : int or (int, int), optional
        Lag range for prediction; defaults to the one used in ``basis``.
    bylag : float
        Step along lags (finer grids give smoother lag-response curves).
    cen : float or bool, optional
        Centring value (reference). ``None``/``True`` picks an automatic
        value for continuous functions; ``False`` disables centring.
    ci_level : float
    cumul : bool
        Also compute incremental cumulative effects over integer lags.
    name : str
        Name prefix of the basis columns in the model design (the ``name``
        given to ``to_dataframe``).
    """
    kind = _basis_type(basis)
    origlag = basis.lag if kind == "cb" else np.array([0, 0])
    lag = origlag if lag is None else mklag(lag)
    if not np.array_equal(lag, origlag) and cumul:
        raise ValueError("cumulative prediction not allowed for lag sub-period")
    lagfun = basis.arglag.get("fun") if kind == "cb" else None
    if bylag != 1 and lagfun == "integer":
        raise ValueError("prediction for non-integer lags not allowed for type 'integer'")
    if not (0 < ci_level < 1):
        raise ValueError("'ci_level' must be numeric and between 0 and 1")

    coef, vcov, model_link, model_class = _resolve_coef(basis, model, coef, vcov, model_link, name, kind)

    rng = basis.range
    at = mkat(at, from_, to, by, rng, lag, bylag)
    at_matrix = at if at.ndim == 2 else None
    predvar = np.arange(at.shape[0], dtype=float) if at.ndim == 2 else at
    predlag = seqlag(lag, bylag)
    cen = mkcen(cen, basis, rng)

    # lag-specific effects
    Xpred = mkXpred(basis, at, predvar, predlag, cen)
    matfit = (Xpred @ coef).reshape((len(predlag), len(predvar))).T
    matse = _se(Xpred, vcov).reshape((len(predlag), len(predvar))).T

    # overall cumulative (and incremental cumulative) effects over integer lags
    predlag_int = seqlag(lag)
    Xpred = mkXpred(basis, at, predvar, predlag_int, cen)
    nv = len(predvar)
    Xall = np.zeros((nv, Xpred.shape[1]))
    cumfit = cumse = None
    if cumul:
        cumfit = np.zeros((nv, len(predlag_int)))
        cumse = np.zeros((nv, len(predlag_int)))
    for i in range(len(predlag_int)):
        Xall = Xall + Xpred[i * nv:(i + 1) * nv]
        if cumul:
            cumfit[:, i] = Xall @ coef
            cumse[:, i] = _se(Xall, vcov)
    allfit = Xall @ coef
    allse = _se(Xall, vcov)

    return CrossPred(predvar=predvar, lag=lag, bylag=bylag, coef=coef, vcov=vcov,
                     matfit=matfit, matse=matse, allfit=allfit, allse=allse, cen=cen,
                     cumfit=cumfit, cumse=cumse, ci_level=ci_level,
                     model_class=model_class, model_link=model_link, at_matrix=at_matrix)


# ----------------------------------------------------------------------------
@dataclass
class CrossReduce:
    """Reduced one-dimensional summary of a DLNM (see :func:`crossreduce`)."""

    type: str
    value: float | None
    coef: np.ndarray
    vcov: np.ndarray
    basis: np.ndarray
    lag: np.ndarray
    bylag: float
    fit: np.ndarray
    se: np.ndarray
    predvar: np.ndarray | None = None
    cen: float | None = None
    ci_level: float = 0.95
    model_class: str | None = None
    model_link: str | None = None

    @property
    def z(self) -> float:
        return float(norm.ppf(1 - (1 - self.ci_level) / 2))

    @property
    def is_exp(self) -> bool:
        return self.model_link in ("log", "logit")

    @property
    def x(self) -> np.ndarray:
        """Horizontal axis: lags for ``type='var'``, predictor otherwise."""
        return seqlag(self.lag, self.bylag) if self.type == "var" else self.predvar

    @property
    def low(self): return self.fit - self.z * self.se
    @property
    def high(self): return self.fit + self.z * self.se
    @property
    def RRfit(self): return np.exp(self.fit)
    @property
    def RRlow(self): return np.exp(self.low)
    @property
    def RRhigh(self): return np.exp(self.high)

    def to_dataframe(self, exp: bool | None = None):
        import pandas as pd
        e = self.is_exp if exp is None else exp
        f = np.exp if e else (lambda a: a)
        col = "lag" if self.type == "var" else "var"
        return pd.DataFrame({col: self.x, "fit": f(self.fit), "low": f(self.low), "high": f(self.high)})

    def plot(self, *args, **kwargs):
        from .plot import plot_crossreduce
        return plot_crossreduce(self, *args, **kwargs)

    def __repr__(self):
        return (f"CrossReduce(type={self.type!r}, value={self.value}, params={self.coef.size}, "
                f"lag={tuple(int(v) for v in self.lag)}, cen={self.cen}, link={self.model_link!r})")


def crossreduce(basis: CrossBasis, model=None, type: str = "overall", value=None,
                coef=None, vcov=None, model_link=None, at=None, from_=None, to=None,
                by=None, lag=None, bylag: float = 1, cen=None, ci_level: float = 0.95,
                name: str = "cb") -> CrossReduce:
    """Reduce a DLNM to a one-dimensional summary.

    ``type="overall"`` gives the overall cumulative exposure-response,
    ``type="lag"`` the exposure-response at lag ``value``, and
    ``type="var"`` the lag-response at predictor value ``value``. The
    reduced coefficients/vcov can be used, e.g., in multivariate
    meta-analysis.
    """
    if not isinstance(basis, CrossBasis):
        raise TypeError("the first argument must be a CrossBasis object")
    if type not in ("overall", "var", "lag"):
        raise ValueError("'type' must be one of 'overall', 'var', 'lag'")
    if type != "overall":
        if value is None:
            raise ValueError("'value' must be provided for type 'var' or 'lag'")
        value = float(np.asarray(value).ravel()[0])
        if type == "lag" and not (basis.lag[0] <= value <= basis.lag[1]):
            raise ValueError("'value' of lag-specific effects must be within the lag range")
    else:
        value = None
    lag = basis.lag if lag is None else mklag(lag)
    if not np.array_equal(lag, basis.lag) and basis.arglag.get("fun") == "integer":
        raise ValueError("prediction for lag sub-period not allowed for type 'integer'")
    if not (0 < ci_level < 1):
        raise ValueError("'ci_level' must be numeric and between 0 and 1")

    coef, vcov, model_link, model_class = _resolve_coef(basis, model, coef, vcov, model_link, name, "cb")

    if at is not None and np.asarray(at).ndim == 2:
        raise ValueError("argument 'at' must be a vector")
    rng = basis.range
    at = mkat(at, from_, to, by, rng, lag, bylag)
    cen = mkcen(cen, basis, rng)

    ncol = basis.ncol
    if type == "overall":
        lagbasis = basis.basis_lag(seqlag(lag))
        M = np.kron(np.eye(ncol // lagbasis.shape[1]), np.ones((1, lagbasis.shape[0])) @ lagbasis)
        newbasis = basis.basis_var(at)
        if cen is not None:
            newbasis = newbasis - basis.basis_var([cen])
    elif type == "lag":
        lagbasis = basis.basis_lag([value])
        M = np.kron(np.eye(ncol // lagbasis.shape[1]), lagbasis)
        newbasis = basis.basis_var(at)
        if cen is not None:
            newbasis = newbasis - basis.basis_var([cen])
    else:
        varbasis = basis.basis_var([value])
        if cen is not None:
            varbasis = varbasis - basis.basis_var([cen])
        M = np.kron(varbasis, np.eye(ncol // varbasis.shape[1]))
        newbasis = basis.basis_lag(seqlag(lag, bylag))

    newcoef = M @ coef
    newvcov = M @ vcov @ M.T
    fit = newbasis @ newcoef
    se = _se(newbasis, newvcov)
    return CrossReduce(type=type, value=value, coef=newcoef, vcov=newvcov, basis=newbasis,
                       lag=lag, bylag=bylag, fit=fit, se=se,
                       predvar=None if type == "var" else at, cen=cen, ci_level=ci_level,
                       model_class=model_class, model_link=model_link)
