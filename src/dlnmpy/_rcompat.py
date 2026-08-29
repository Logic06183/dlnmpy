"""Small numerical helpers that reproduce base-R behaviour exactly.

The R package ``dlnm`` leans on a few base-R functions whose precise
semantics matter for numerical equivalence: ``pretty()`` (used to choose the
default prediction grid and the automatic centring value), ``quantile()``
(type 7, used to place knots), ``cut()`` with ``right=FALSE`` (used by
``strata``) and ``seq()``. Each is reimplemented here from the R source so
that a port in any language can copy the same algorithm.
"""

from __future__ import annotations

import math
import sys

import numpy as np

__all__ = ["pretty", "quantile7", "seq", "median", "r_range"]

_DBL_EPSILON = sys.float_info.epsilon
_DBL_MIN = sys.float_info.min
_DBL_MAX = sys.float_info.max


def _r_pretty(lo: float, up: float, ndiv: int, min_n: int, shrink_sml: float,
              h: float, h5: float, f_min: float, eps_correction: int,
              return_bounds: bool):
    """Port of ``R_pretty`` in ``src/appl/pretty.c`` (R >= 4.2).

    Returns ``(lo, up, ndiv, unit)``.
    """
    rounding_eps = 1e-10
    lo_, up_ = lo, up
    dx = up_ - lo_

    if dx == 0 and up_ == 0:
        cell = 1.0
        i_small = True
    else:
        cell = max(abs(lo_), abs(up_))
        U = 1 + (1 / (1 + h) if h5 >= 1.5 * h + 0.5 else 1.5 / (1 + h5))
        U *= max(1, ndiv) * _DBL_EPSILON
        i_small = dx < cell * U * 3

    if i_small:
        if cell > 10:
            cell = 9 + cell / 10
        cell *= shrink_sml
        if min_n > 1:
            cell /= min_n
    else:
        cell = dx
        if math.isfinite(dx):
            if ndiv > 1:
                cell /= ndiv
        else:
            if ndiv >= 2:
                cell = up_ / ndiv - lo_ / ndiv

    subsmall = f_min * _DBL_MIN
    if subsmall == 0.0:
        subsmall = _DBL_MIN
    max_f = 1.25
    if cell < subsmall:
        cell = subsmall
    elif cell > _DBL_MAX / max_f:
        cell = _DBL_MAX / max_f

    base = 10.0 ** math.floor(math.log10(cell))
    unit = base
    U = 2 * base
    if U - cell < h * (cell - unit):
        unit = U
        U = 5 * base
        if U - cell < h5 * (cell - unit):
            unit = U
            U = 10 * base
            if U - cell < h * (cell - unit):
                unit = U

    ns = math.floor(lo_ / unit + rounding_eps)
    nu = math.ceil(up_ / unit - rounding_eps)

    if eps_correction and (eps_correction > 1 or not i_small):
        E_ = _DBL_EPSILON
        D_max = _DBL_MAX * (1.0 - math.ldexp(E_, -1))
        if lo_ < 0.0:
            lo *= (1 + E_)
        elif lo_ > 0:
            lo *= (1 - E_)
        else:
            lo = -min(unit, _DBL_MIN)
        if up_ < 0.0:
            up *= (1 - E_)
        elif up_ > 0.0:
            if up_ < D_max:
                up *= (1 + E_)
        else:
            up = min(unit, _DBL_MIN)

    while ns * unit > lo + rounding_eps * unit:
        ns -= 1
    while not math.isfinite(ns * unit):
        ns += 1
    while nu * unit < up - rounding_eps * unit:
        nu += 1
    while not math.isfinite(nu * unit):
        nu -= 1

    k = int(0.5 + nu - ns)
    if k < min_n:
        k = min_n - k
        if lo_ == 0.0 and ns == 0.0 and up_ != 0.0:
            nu += k
        elif up_ == 0.0 and nu == 0.0 and lo_ != 0.0:
            ns -= k
        elif ns >= 0.0:
            nu += k // 2
            ns -= k // 2 + k % 2
        else:
            ns -= k // 2
            nu += k // 2 + k % 2
        ndiv = min_n
    else:
        ndiv = k

    if return_bounds:
        if ns * unit < lo:
            lo = ns * unit
        if nu * unit > up:
            up = nu * unit
    else:
        lo, up = ns, nu
    return lo, up, ndiv, unit


def pretty(x, n: float = 5, min_n=None, shrink_sml: float = 0.75,
           high_u_bias: float = 1.5, u5_bias=None, eps_correct: int = 0,
           f_min: float = 2.0 ** -20) -> np.ndarray:
    """Equivalent of R's ``pretty.default``.

    Computes a sequence of about ``n + 1`` equally spaced "round" values
    covering the range of ``x``. Non-integer ``n`` is truncated, as R's
    ``asInteger`` does. ``min_n`` defaults to ``n %/% 3``.
    """
    x = np.asarray(x, dtype=float).ravel()
    x = x[np.isfinite(x)]
    if x.size == 0:
        return x
    if min_n is None:
        min_n = math.floor(n / 3)  # R: n %/% 3, evaluated before coercion
    if u5_bias is None:
        u5_bias = 0.5 + 1.5 * high_u_bias
    ndiv = int(n)  # asInteger truncates towards zero
    min_n = int(min_n)
    lo, up, ndiv, _unit = _r_pretty(float(x.min()), float(x.max()), ndiv, min_n,
                                    shrink_sml, high_u_bias, u5_bias, f_min,
                                    int(eps_correct), True)
    s = np.linspace(lo, up, ndiv + 1)
    if not eps_correct and ndiv:
        delta = (up - lo) / ndiv
        small = np.abs(s) < 1e-14 * delta
        s[small] = 0.0
    return s


def quantile7(x, probs) -> np.ndarray:
    """R's default ``quantile(x, probs, type=7)``, ignoring NaN."""
    x = np.asarray(x, dtype=float).ravel()
    x = np.sort(x[~np.isnan(x)])
    probs = np.atleast_1d(np.asarray(probs, dtype=float))
    n = x.size
    if n == 0:
        return np.full(probs.shape, np.nan)
    # Port of stats:::quantile.default for type 7 (a = b = 1).
    fuzz = 4 * _DBL_EPSILON
    nppm = 1 + probs * (n - 1)
    j = np.floor(nppm + fuzz).astype(int)
    h = nppm - j
    h[np.abs(h) < fuzz] = 0.0
    xp = np.concatenate(([x[0], x[0]], x, [x[-1], x[-1]]))  # pad both ends
    lo = xp[j + 1]  # R index j+2 (1-based) -> j+1 (0-based)
    hi = xp[j + 2]
    qs = lo.copy()
    qs[h == 1] = hi[h == 1]
    other = (0 < h) & (h < 1) & (lo != hi)
    qs[other] = ((1 - h) * lo + h * hi)[other]
    return qs


def median(x) -> float:
    """R's ``median`` (mean of the two middle order statistics, NaN removed)."""
    x = np.asarray(x, dtype=float).ravel()
    x = np.sort(x[~np.isnan(x)])
    n = x.size
    if n == 0:
        return float("nan")
    half = (n + 1) // 2
    if n % 2 == 1:
        return float(x[half - 1])
    return float((x[half - 1] + x[half]) / 2.0)


def seq(start: float, stop: float, by: float = 1.0) -> np.ndarray:
    """R's ``seq(from, to, by)`` (with its fuzz on the number of elements)."""
    if by == 0:
        raise ValueError("'by' must be non-zero")
    n = (stop - start) / by
    if n < -1e-10:
        raise ValueError("wrong sign in 'by' argument")
    nn = int(math.floor(n + 1e-10))
    out = start + np.arange(nn + 1, dtype=float) * by
    # R clips the last element so it never overshoots 'to'
    if by > 0:
        out[out > stop] = stop
    else:
        out[out < stop] = stop
    return out


def r_range(x) -> tuple[float, float]:
    x = np.asarray(x, dtype=float)
    return float(np.nanmin(x)), float(np.nanmax(x))
