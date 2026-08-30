"""Core objects: ``OneBasis`` and ``CrossBasis``.

A *one-dimensional basis* transforms a vector ``x`` into a matrix of basis
variables through one of the functions in :mod:`dlnmpy.basis`.

A *cross-basis* is the bi-dimensional extension used by distributed lag
(non-)linear models. Given a basis for the predictor space (``argvar``) and a
basis for the lag space (``arglag``), the cross-basis columns are the
tensor product of the two. For a time series ``x`` the lagged occurrences are
built internally; alternatively ``x`` may already be a matrix of exposure
histories (rows = observations, columns = lags ``lag[0]..lag[1]``).

The algebra, following Gasparrini, Armstrong & Kenward (2010) and the R
implementation (``crossbasis.R``), is: with ``R`` the (n x vx) predictor
basis, ``Q_v`` the (n x L) matrix of lagged copies of column ``v`` of ``R``,
and ``C`` the (L x vl) lag basis, the cross-basis column for predictor basis
``v`` and lag basis ``l`` is ``Q_v @ C[:, l]``. Columns are ordered with the
lag index varying fastest: ``v1.l1, v1.l2, ..., v2.l1, ...``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .basis import PRED_ARGS, get_basis_function
from .lag import lag_matrix, mklag, seqlag

__all__ = ["OneBasis", "CrossBasis", "onebasis", "crossbasis", "normalise_args"]


# ----------------------------------------------------------------------------
# argument normalisation (accept R spellings)
# ----------------------------------------------------------------------------
_ALIASES = {
    "Boundary.knots": "boundary_knots",
    "Boundary_knots": "boundary_knots",
    "thr.value": "thr_value",
    "thr": "thr_value",
    "type": "fun",
}


def normalise_args(args: dict | None) -> dict:
    """Return a copy of ``args`` with R-style argument names mapped to the
    Python spellings used by :mod:`dlnmpy.basis`."""
    out: dict[str, Any] = {}
    for k, v in (args or {}).items():
        out[_ALIASES.get(k, k)] = v
    return out


def _apply_fun(fun, x, args: dict):
    f = get_basis_function(fun)
    basis, attrs = f(x, **args)
    basis = np.asarray(basis, dtype=float)
    if basis.ndim == 1:
        basis = basis[:, None]
    return basis, attrs


# ----------------------------------------------------------------------------
class _MatrixLike:
    """Mixin giving numpy-array behaviour to the basis containers."""

    matrix: np.ndarray

    def __array__(self, dtype=None, copy=None):
        return self.matrix if dtype is None else self.matrix.astype(dtype)

    @property
    def shape(self):
        return self.matrix.shape

    def __len__(self):
        return self.matrix.shape[0]

    def __getitem__(self, item):
        return self.matrix[item]

    @property
    def ncol(self) -> int:
        return self.matrix.shape[1]


@dataclass
class OneBasis(_MatrixLike):
    """A one-dimensional basis matrix with the information to rebuild it.

    Attributes
    ----------
    matrix : ndarray (n, k)
    fun : str or callable
    attrs : dict
        Attributes returned by the basis function (knots, thresholds, ...).
    range : (float, float)
        Range of the original ``x``.
    cen : float or None
        Centring value stored for later use in prediction (optional).
    """

    matrix: np.ndarray
    fun: Any
    attrs: dict = field(default_factory=dict)
    range: tuple = (np.nan, np.nan)
    cen: Any = None
    index: Any = None

    @property
    def colnames(self) -> list[str]:
        return [f"b{i + 1}" for i in range(self.ncol)]

    @property
    def fun_name(self) -> str:
        return self.fun if isinstance(self.fun, str) else getattr(self.fun, "__name__", "custom")

    def pred_args(self) -> dict:
        """Arguments to re-apply the same transformation to new data."""
        keys = PRED_ARGS.get(self.fun_name)
        if keys is None:  # custom callable: pass back everything it returned
            return dict(self.attrs)
        return {k: self.attrs[k] for k in keys if k in self.attrs}

    def transform(self, x) -> np.ndarray:
        """Apply the stored transformation to new values ``x``."""
        basis, _ = _apply_fun(self.fun, x, self.pred_args())
        return basis

    def to_dataframe(self, name: str = "b", index=None):
        """DataFrame of the basis columns.

        The index of the pandas object the basis was built from is carried
        over, so ``data.join(basis.to_dataframe(...))`` aligns even when
        ``data`` does not have a default RangeIndex. Pass ``index`` to
        override it.
        """
        import pandas as pd
        idx = index if index is not None else self.index
        return pd.DataFrame(self.matrix, columns=[f"{name}_{c}" for c in self.colnames], index=idx)

    def __repr__(self):
        return f"OneBasis(fun={self.fun_name!r}, shape={self.shape}, attrs={_short(self.attrs)})"


def onebasis(x, fun="ns", cen=None, **kwargs) -> OneBasis:
    """Generate a basis matrix for a vector ``x``.

    Parameters
    ----------
    x : array_like
    fun : str or callable
        Basis function name (``"ns"`` by default) or a callable
        ``f(x, **kwargs) -> (matrix, attrs)``.
    cen : float, optional
        Centring value kept for the prediction stage (not applied here).
    **kwargs
        Arguments passed to the basis function (R spellings accepted).
    """
    index = getattr(x, "index", None)
    x = np.asarray(x, dtype=float).ravel()
    rng = (float(np.nanmin(x)), float(np.nanmax(x)))
    args = normalise_args(kwargs)
    args.pop("cen", None)
    basis, attrs = _apply_fun(fun, x, args)
    return OneBasis(matrix=basis, fun=fun, attrs=attrs, range=rng, cen=cen, index=index)


# ----------------------------------------------------------------------------
@dataclass
class CrossBasis(_MatrixLike):
    """A cross-basis matrix for a distributed lag (non-)linear model.

    Attributes
    ----------
    matrix : ndarray (n, vx * vl)
    df : (int, int)
        Number of basis columns in the predictor and lag spaces.
    range : (float, float)
        Range of the predictor.
    lag : ndarray [lag0, lag1]
    argvar, arglag : dict
        Function name and attributes for each dimension, used to rebuild
        the bases at prediction time.
    group : int or None
        Number of groups if ``group`` was used.
    """

    matrix: np.ndarray
    df: tuple
    range: tuple
    lag: np.ndarray
    argvar: dict
    arglag: dict
    group: int | None = None
    index: Any = None

    @property
    def colnames(self) -> list[str]:
        vx, vl = self.df
        return [f"v{v + 1}.l{l + 1}" for v in range(vx) for l in range(vl)]

    def to_dataframe(self, name: str = "cb", index=None):
        """Return a DataFrame with columns ``{name}_v{i}_l{j}``, ready to be
        joined to a data frame and used in a statsmodels formula.

        The index of the pandas object the cross-basis was built from is
        carried over, so ``data.join(cb.to_dataframe("cb"))`` aligns even
        when ``data`` does not have a default RangeIndex (after a groupby or
        a filter, say). Pass ``index`` to override it.
        """
        import pandas as pd
        cols = [f"{name}_{c.replace('.', '_')}" for c in self.colnames]
        idx = index if index is not None else self.index
        return pd.DataFrame(self.matrix, columns=cols, index=idx)

    def basis_var(self, x) -> np.ndarray:
        """Predictor-space basis evaluated at ``x``."""
        a = dict(self.argvar)
        fun = a.pop("fun")
        a.pop("cen", None)
        return _apply_fun(fun, x, a)[0]

    def basis_lag(self, lags) -> np.ndarray:
        """Lag-space basis evaluated at ``lags``."""
        a = dict(self.arglag)
        fun = a.pop("fun")
        return _apply_fun(fun, lags, a)[0]

    def summary(self) -> str:
        vx, vl = self.df
        lines = ["CROSSBASIS FUNCTIONS",
                 f"observations: {self.shape[0]}",
                 f"range: {self.range[0]:g} to {self.range[1]:g}",
                 f"lag period: {self.lag[0]} {self.lag[1]}",
                 f"total df: {vx * vl}",
                 "",
                 "BASIS FOR VAR:"] + _fmt_args(self.argvar) + \
                ["", "BASIS FOR LAG:"] + _fmt_args(self.arglag)
        return "\n".join(lines)

    def __repr__(self):
        return (f"CrossBasis(shape={self.shape}, lag={tuple(int(v) for v in self.lag)}, "
                f"argvar={_short(self.argvar)}, arglag={_short(self.arglag)})")


def crossbasis(x, lag=None, argvar=None, arglag=None, group=None) -> CrossBasis:
    """Generate a cross-basis matrix.

    Parameters
    ----------
    x : array_like
        A time series vector, or a matrix of exposure histories with
        ``lag[1] - lag[0] + 1`` columns (see :func:`dlnmpy.exphist`).
    lag : int or (int, int)
        Maximum lag, or ``(min lag, max lag)``. Defaults to ``(0, ncol(x)-1)``
        for matrix input.
    argvar : dict
        Basis specification for the predictor space, e.g.
        ``{"fun": "ns", "df": 5}`` or ``{"fun": "thr", "thr_value": 25}``.
        Defaults to ``fun="ns"``.
    arglag : dict
        Basis specification for the lag space. Defaults to a single stratum
        over all lags (a moving sum) when empty or when ``lag[0] == lag[1]``.
        An intercept is included by default, as in R.
    group : array_like, optional
        Grouping factor for time series made of several consecutive series
        (lags are not computed across groups).
    """
    argvar = normalise_args(argvar)
    arglag = normalise_args(arglag)
    index = getattr(x, "index", None)
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    lag = np.array([0, x.shape[1] - 1]) if lag is None else mklag(lag)
    nlag = int(lag[1] - lag[0]) + 1
    if x.shape[1] not in (1, nlag):
        raise ValueError("ncol(x) must be 1 (time series) or equal to the lag period "
                         "(matrix of lagged occurrences)")

    # basis for the predictor space
    varfun = argvar.pop("fun", "ns")
    # R keeps a cen given inside argvar as attr(cb, "argvar")$cen and mkcen()
    # uses it as the default reference for every later crosspred/crossreduce
    varcen = argvar.pop("cen", None)
    basisvar, attrvar = _apply_fun(varfun, x.ravel(order="F"), argvar)

    # basis for the lag space
    if len(arglag) == 0 or nlag == 1:
        arglag = {"fun": "strata", "df": 1, "intercept": True}
    lagfun = arglag.pop("fun", "ns")
    if "intercept" not in arglag and _has_intercept(lagfun):
        arglag["intercept"] = True
    arglag.pop("cen", None)
    basislag, attrlag = _apply_fun(lagfun, seqlag(lag), arglag)

    if group is not None:
        group = np.asarray(group)
        if x.shape[1] > 1:
            raise ValueError("'group' allowed only for time series data")
        _, counts = np.unique(group, return_counts=True)
        if counts.min() <= lag[1] - lag[0]:
            raise ValueError("each group must have length > diff(lag)")

    # cross-basis computation
    n = x.shape[0]
    vx, vl = basisvar.shape[1], basislag.shape[1]
    cb = np.zeros((n, vx * vl))
    for v in range(vx):
        if x.shape[1] == 1:
            mat = lag_matrix(basisvar[:, v], seqlag(lag), group=group)
        else:
            mat = basisvar[:, v].reshape((n, nlag), order="F")
        with np.errstate(divide="ignore", over="ignore", invalid="ignore", under="ignore"):
            cb[:, vl * v: vl * (v + 1)] = mat @ basislag

    av = {"fun": varfun, **_select_pred_args(varfun, attrvar)}
    av["cen"] = varcen
    al = {"fun": lagfun, **_select_pred_args(lagfun, attrlag)}
    return CrossBasis(matrix=cb, df=(vx, vl),
                      range=(float(np.nanmin(x)), float(np.nanmax(x))),
                      lag=lag, argvar=av, arglag=al,
                      group=None if group is None else int(np.unique(group).size),
                      index=index)


def _select_pred_args(fun, attrs: dict) -> dict:
    name = fun if isinstance(fun, str) else getattr(fun, "__name__", "custom")
    keys = PRED_ARGS.get(name)
    if keys is None:
        return dict(attrs)
    return {k: attrs[k] for k in keys if k in attrs}


def _has_intercept(fun) -> bool:
    import inspect
    try:
        f = get_basis_function(fun)
        return "intercept" in inspect.signature(f).parameters
    except (ValueError, TypeError):
        return False


def _short(d: dict) -> dict:
    out = {}
    for k, v in d.items():
        if isinstance(v, np.ndarray):
            out[k] = f"array{v.shape}" if v.size > 6 else np.round(v, 4).tolist()
        else:
            out[k] = v
    return out


def _fmt_args(d: dict) -> list[str]:
    out = []
    for k, v in d.items():
        if v is None:
            continue
        if isinstance(v, np.ndarray):
            if v.ndim > 1:
                v = f"matrix {v.shape}"
            else:
                v = " ".join(f"{t:g}" for t in v)
        out.append(f"{k}: {v}")
    return out
