"""Helpers to place knots for spline and strata bases."""

from __future__ import annotations

import numpy as np

from .lag import mklag

__all__ = ["logknots", "equalknots"]


def _nk_from_df(fun: str, df: int, degree: int, intercept: bool) -> int:
    if fun == "ns":
        return df - 1 - int(intercept)
    if fun == "bs":
        return df - degree - int(intercept)
    if fun == "strata":
        return df - int(intercept)
    raise ValueError("'fun' must be one of 'ns', 'bs', 'strata'")


def logknots(x, nk=None, fun: str = "ns", df: int = 1, degree: int = 3,
             intercept: bool = True) -> np.ndarray:
    """Knots at equally spaced values on the log scale of the lags.

    ``x`` is a lag range (length 1 or 2) or a vector whose range is used.
    The number of knots is ``nk`` or derived from ``fun``/``df``/``degree``/
    ``intercept`` for the lag basis (note ``intercept=True`` by default, as
    the lag basis in a cross-basis includes an intercept).
    """
    x = np.atleast_1d(np.asarray(x, dtype=float))
    rng = mklag(x).astype(float) if x.size < 3 else np.array([np.nanmin(x), np.nanmax(x)])
    if rng[1] - rng[0] == 0:
        raise ValueError("range must be >0")
    if nk is None:
        nk = _nk_from_df(fun, df, degree, intercept)
    if nk < 1:
        raise ValueError("choice of arguments defines no knots")
    span = rng[1] - rng[0]
    return rng[0] + np.exp(((1 + np.log(span)) / (nk + 1)) * np.arange(1, nk + 1) - 1)


def equalknots(x, nk=None, fun: str = "ns", df: int = 1, degree: int = 3,
               intercept: bool = False) -> np.ndarray:
    """Knots at equally spaced values along the range of ``x``."""
    x = np.asarray(x, dtype=float)
    rng = np.array([np.nanmin(x), np.nanmax(x)])
    if nk is None:
        nk = _nk_from_df(fun, df, degree, intercept)
    if nk < 1:
        raise ValueError("choice of arguments defines no knots")
    span = rng[1] - rng[0]
    return rng[0] + (span / (nk + 1)) * np.arange(1, nk + 1)
