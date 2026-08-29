"""Basis functions for the predictor and lag spaces.

Each function takes a vector ``x`` and returns ``(matrix, attrs)``. ``attrs``
records everything needed to rebuild the *same* transformation on new values
at the prediction stage (knots, thresholds, scale ...). This mirrors how the R
package stores the arguments as attributes of the basis matrix and later
re-applies them inside ``crosspred``/``crossreduce``.

Available functions (names follow the R package):

- ``lin``      linear
- ``poly``     polynomial (scaled)
- ``strata``   indicator variables for intervals
- ``thr``      high/low/double threshold (hockey-stick) functions
- ``integer``  one indicator per integer value (unconstrained lag)
- ``ns``       natural cubic spline (``splines::ns``)
- ``bs``       B-spline (``splines::bs``)
- ``ps``       penalised B-spline (P-spline), with its difference penalty

``PRED_ARGS`` lists, for each function, which attributes are passed back when
the transformation is re-applied (R uses ``formals()`` for this).
"""

from __future__ import annotations

import numpy as np

from . import _splines
from ._rcompat import median, quantile7

__all__ = ["lin", "poly", "strata", "thr", "integer", "ns", "bs", "ps",
           "BASIS_FUNCTIONS", "PRED_ARGS", "get_basis_function"]


def _asvec(x) -> np.ndarray:
    return np.asarray(x, dtype=float).ravel()


# ----------------------------------------------------------------------------
def lin(x, intercept: bool = False):
    """Linear basis: ``x`` (optionally preceded by a column of ones)."""
    x = _asvec(x)
    basis = x[:, None]
    if intercept:
        basis = np.column_stack((np.ones(x.size), basis))
    return basis, {"intercept": bool(intercept)}


# ----------------------------------------------------------------------------
def poly(x, degree: int = 1, scale=None, intercept: bool = False):
    """Polynomial basis ``(x/scale)^k`` for ``k = 1..degree`` (``0..degree``
    with intercept). ``scale`` defaults to ``max(abs(x))``."""
    x = _asvec(x)
    degree = int(degree)
    if scale is None:
        scale = float(np.nanmax(np.abs(x)))
    powers = np.arange(1 - int(intercept), degree + 1)
    basis = np.power.outer(x / scale, powers)
    return basis, {"degree": degree, "scale": float(scale), "intercept": bool(intercept)}


# ----------------------------------------------------------------------------
def strata(x, df: int = 1, breaks=None, ref: int = 1, intercept: bool = False):
    """Indicator variables for strata of ``x`` defined by ``breaks``
    (left-closed intervals ``[a, b)``). Without ``breaks``, ``df`` strata are
    formed at equally spaced quantiles. ``ref`` (1-based) is the reference
    stratum dropped from the basis (``0`` keeps all, only with intercept)."""
    x = _asvec(x)
    rng = (np.nanmin(x), np.nanmax(x))
    intercept = bool(intercept)
    if breaks is not None:
        breaks = np.unique(np.asarray(breaks, dtype=float).ravel())
    elif df - intercept > 0:
        k = int(df) - intercept
        breaks = quantile7(x, np.arange(1, k + 1) / (k + 1))
    else:
        breaks = None
    df = (0 if breaks is None else breaks.size) + intercept
    # cut(x, c(range[1]-0.0001, breaks, range[2]+0.0001), right=FALSE)
    edges = np.concatenate(([rng[0] - 0.0001], [] if breaks is None else breaks, [rng[1] + 0.0001]))
    cat = np.searchsorted(edges, x, side="right") - 1  # 0-based stratum, NaN -> -1 or len
    nlev = edges.size - 1
    basis = np.zeros((x.size, nlev))
    valid = (cat >= 0) & (cat < nlev) & ~np.isnan(x)
    basis[np.nonzero(valid)[0], cat[valid]] = 1.0
    basis[np.isnan(x), :] = np.nan
    ref = int(ref)
    if ref not in range(0, basis.shape[1] + 1):
        raise ValueError("wrong value in 'ref' argument. See help('strata')")
    if not intercept and ref == 0:
        ref = 1
    if breaks is not None:
        if ref != 0:
            basis = np.delete(basis, ref - 1, axis=1)
        if intercept and ref != 0:
            basis = np.column_stack((np.ones(x.size), basis))
    attrs = {"df": int(df), "breaks": breaks, "ref": ref, "intercept": intercept}
    return basis, attrs


# ----------------------------------------------------------------------------
def thr(x, thr_value=None, side=None, intercept: bool = False):
    """Threshold (hockey-stick) functions.

    ``side`` is ``"h"`` (linear above the threshold, zero below), ``"l"``
    (linear below) or ``"d"`` (double threshold: two columns, linear below
    ``thr_value[0]`` and above ``thr_value[-1]``). The default is ``"d"`` when
    two thresholds are given and ``"h"`` otherwise; the default threshold is
    the median of ``x``."""
    x = _asvec(x)
    if thr_value is None:
        thr_value = np.array([median(x)])
    else:
        thr_value = np.sort(np.atleast_1d(np.asarray(thr_value, dtype=float)))
    if side is None:
        side = "d" if thr_value.size > 1 else "h"
    if side not in ("h", "l", "d"):
        raise ValueError("'side' must be one of 'h', 'l', 'd'")
    thr_value = thr_value[[0, -1]] if side == "d" else thr_value[[0]]
    if side == "h":
        basis = np.maximum(x - thr_value[0], 0)[:, None]
    elif side == "l":
        basis = (-np.minimum(x - thr_value[0], 0))[:, None]
    else:
        basis = np.column_stack((-np.minimum(x - thr_value[0], 0),
                                 np.maximum(x - thr_value[1], 0)))
    if intercept:
        basis = np.column_stack((np.ones(x.size), basis))
    # keep thr_value as scalar for single-threshold sides, as R does
    tv = thr_value if side == "d" else thr_value
    return basis, {"thr_value": tv, "side": side, "intercept": bool(intercept)}


# ----------------------------------------------------------------------------
def integer(x, values=None, intercept: bool = False):
    """One indicator column per (integer) value of ``x``; the first level is
    dropped unless ``intercept`` is True. Used for unconstrained DLMs."""
    x = _asvec(x)
    levels = np.sort(np.unique(x[~np.isnan(x)])) if values is None else np.asarray(values, dtype=float).ravel()
    basis = (x[:, None] == levels[None, :]).astype(float)
    # as in R's factor(x, levels): values not among the levels become NA
    basis[np.isnan(x) | ~np.isin(x, levels), :] = np.nan
    intercept = bool(intercept)
    if basis.shape[1] > 1:
        if not intercept:
            basis = basis[:, 1:]
    else:
        intercept = True
    return basis, {"values": levels, "intercept": intercept}


# ----------------------------------------------------------------------------
def ns(x, df=None, knots=None, intercept: bool = False, Boundary_knots=None, **kw):
    """Natural cubic spline (see :func:`dlnmpy._splines.ns`)."""
    bk = kw.pop("boundary_knots", Boundary_knots)
    if kw:
        raise TypeError(f"unexpected arguments for ns: {sorted(kw)}")
    return _splines.ns(x, df=df, knots=knots, intercept=intercept, boundary_knots=bk)


def bs(x, df=None, knots=None, degree: int = 3, intercept: bool = False,
       Boundary_knots=None, **kw):
    """B-spline (see :func:`dlnmpy._splines.bs`)."""
    bk = kw.pop("boundary_knots", Boundary_knots)
    if kw:
        raise TypeError(f"unexpected arguments for bs: {sorted(kw)}")
    return _splines.bs(x, df=df, knots=knots, degree=degree, intercept=intercept,
                       boundary_knots=bk)


# ----------------------------------------------------------------------------
def ps(x, df: int = 10, knots=None, degree: int = 3, intercept: bool = False,
       fx: bool = False, S=None, diff: int = 2):
    """P-spline basis (B-spline with equally spaced knots) and its
    difference penalty matrix ``S`` (port of ``dlnm::ps``)."""
    x = _asvec(x)
    rng = (np.nanmin(x), np.nanmax(x))
    nax = np.isnan(x)
    nas = bool(np.any(nax))
    if nas:
        x = x[~nax]
    degree = int(degree)
    if degree < 1:
        raise ValueError("'degree' must be integer >= 1")
    intercept = bool(intercept)
    if knots is None or np.asarray(knots).size == 2:
        nik = int(df) - degree + 2 - intercept
        if nik <= 1:
            raise ValueError("basis dimension too small for b-spline degree")
        span = rng[1] - rng[0]
        if knots is not None and np.asarray(knots).size == 2:
            kk = np.asarray(knots, dtype=float)
            xl = kk.min() - span * 0.001
            xu = kk.max() + span * 0.001
        else:
            xl = x.min() - span * 0.001
            xu = x.max() + span * 0.001
        dx = (xu - xl) / (nik - 1)
        knots = np.linspace(xl - dx * degree, xu + dx * degree, nik + 2 * degree)
    else:
        knots = np.asarray(knots, dtype=float).ravel()
        df = knots.size - degree - 2 + intercept
        if df - degree <= 1:
            raise ValueError("basis dimension too small for b-spline degree")
    basis = _splines.spline_design(knots, x, degree + 1, 0, outer_ok=True)
    if not intercept:
        basis = basis[:, 1:]
    if nas:
        nmat = np.full((nax.size, basis.shape[1]), np.nan)
        nmat[~nax, :] = basis
        basis = nmat
    diff = int(diff)
    if diff < 1:
        raise ValueError("'diff' must be an integer >=1")
    if fx:
        S = None
    elif S is None:
        D = np.diff(np.eye(basis.shape[1] + (not intercept)), n=diff, axis=0)
        S = D.T @ D
        S = (S + S.T) / 2
        if not intercept:
            S = S[1:, 1:]
    else:
        S = np.asarray(S, dtype=float)
        if S.shape != (basis.shape[1], basis.shape[1]):
            raise ValueError("dimensions of 'S' not compatible")
    attrs = {"df": int(df), "knots": knots, "degree": degree, "intercept": intercept,
             "fx": bool(fx), "S": S, "diff": diff}
    return basis, attrs


# ----------------------------------------------------------------------------
BASIS_FUNCTIONS = {
    "lin": lin, "poly": poly, "strata": strata, "thr": thr, "integer": integer,
    "ns": ns, "bs": bs, "ps": ps,
}

# Attributes re-used when the transformation is applied to new data
# (R: intersection of formals(fun) and attributes(basis)).
PRED_ARGS = {
    "lin": ["intercept"],
    "poly": ["degree", "scale", "intercept"],
    "strata": ["df", "breaks", "ref", "intercept"],
    "thr": ["thr_value", "side", "intercept"],
    "integer": ["values", "intercept"],
    "ns": ["knots", "intercept", "boundary_knots"],
    "bs": ["knots", "degree", "intercept", "boundary_knots"],
    "ps": ["df", "knots", "degree", "intercept", "fx", "S", "diff"],
}

# Functions that accept an 'intercept' argument (all of them here); kept for
# parity with the R logic that only adds intercept=TRUE to the lag basis when
# the function has such an argument.
HAS_INTERCEPT = set(BASIS_FUNCTIONS)


def get_basis_function(fun):
    """Resolve ``fun`` (a name or a callable) to a callable."""
    if callable(fun):
        return fun
    try:
        return BASIS_FUNCTIONS[fun]
    except KeyError:
        raise ValueError(f"unknown basis function '{fun}'. Choose from {sorted(BASIS_FUNCTIONS)}") from None
