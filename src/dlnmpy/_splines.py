"""B-spline machinery replicating R's ``splines`` package.

``spline_design`` is a port of ``splines::splineDesign`` (and the C routine
``spline_basis`` underneath it), ``ns`` and ``bs`` are ports of
``splines::ns`` and ``splines::bs``. The ports are deliberately literal so
that the matrices agree with R to machine precision, including the
linear/Taylor extrapolation beyond the boundary knots that R applies at
prediction time.
"""

from __future__ import annotations

import math
import warnings

import numpy as np

from ._rcompat import quantile7

__all__ = ["spline_design", "ns", "bs"]


# ----------------------------------------------------------------------------
# Low-level: port of spline_basis() from src/library/splines/src/splines.c
# ----------------------------------------------------------------------------

def _set_cursor(knots: np.ndarray, x: float, order: int) -> tuple[int, bool]:
    nk = knots.size
    curs = -1
    boundary = False
    for i in range(nk):
        if knots[i] >= x:
            curs = i
        if knots[i] > x:
            break
    if curs > nk - order:
        last_legit = nk - order
        if x == knots[last_legit]:
            boundary = True
            curs = last_legit
    return curs, boundary


def _basis_funcs(knots: np.ndarray, curs: int, x: float, order: int) -> np.ndarray:
    ordm1 = order - 1
    rdel = np.array([knots[curs + i] - x for i in range(ordm1)])
    ldel = np.array([x - knots[curs - (i + 1)] for i in range(ordm1)])
    b = np.zeros(order)
    b[0] = 1.0
    for j in range(1, ordm1 + 1):
        saved = 0.0
        for r in range(j):
            den = rdel[r] + ldel[j - 1 - r]
            if den != 0:
                term = b[r] / den
                b[r] = saved + rdel[r] * term
                saved = ldel[j - 1 - r] * term
            else:
                if r != 0 or rdel[r] != 0.0:
                    b[r] = saved
                saved = 0.0
        b[j] = saved
    return b


def _evaluate(knots: np.ndarray, curs: int, boundary: bool, a: np.ndarray,
              x: float, nder: int, order: int) -> float:
    ordm1 = order - 1
    if boundary and nder == ordm1:
        return 0.0
    a = a.copy()
    outer = ordm1
    ti = curs
    while nder > 0:
        nder -= 1
        for inner in range(outer):
            lpt = ti - outer + inner
            a[inner] = outer * (a[inner + 1] - a[inner]) / (knots[lpt + outer] - knots[lpt])
        outer -= 1
    rdel = np.array([knots[curs + i] - x for i in range(outer)])
    ldel = np.array([x - knots[curs - (i + 1)] for i in range(outer)])
    while outer > 0:
        outer -= 1
        for inner in range(outer + 1):
            lp = ldel[outer - inner]
            rp = rdel[inner]
            a[inner] = (a[inner + 1] * lp + a[inner] * rp) / (rp + lp)
    return a[0]


def _spline_basis(knots: np.ndarray, order: int, x: np.ndarray,
                  derivs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(values (nx, order), offsets (nx,))`` like R's C routine."""
    nx = x.size
    nk = knots.size
    val = np.empty((nx, order))
    offsets = np.empty(nx, dtype=int)
    nd = derivs.size
    for i in range(nx):
        xi = x[i]
        curs, boundary = _set_cursor(knots, xi, order)
        io = curs - order
        offsets[i] = io
        der_i = int(derivs[i % nd])
        if io < 0 or io > nk:
            val[i, :] = np.nan
        elif der_i > 0:
            if der_i >= order:
                raise ValueError(f"derivs = {der_i} >= ord = {order}, but should be in {{0,..,ord-1}}")
            for ii in range(order):
                a = np.zeros(order)
                a[ii] = 1.0
                val[i, ii] = _evaluate(knots, curs, boundary, a, xi, der_i, order)
        else:
            val[i, :] = _basis_funcs(knots, curs, xi, order)
    return val, offsets


def _basis_funcs_vectorised(knots: np.ndarray, curs: np.ndarray, x: np.ndarray,
                            order: int) -> np.ndarray:
    """Vectorised version of ``_basis_funcs`` for the derivative-free path."""
    ordm1 = order - 1
    nx = x.size
    idx = np.arange(ordm1)
    rdel = knots[curs[:, None] + idx[None, :]] - x[:, None]
    ldel = x[:, None] - knots[curs[:, None] - (idx[None, :] + 1)]
    b = np.zeros((nx, order))
    b[:, 0] = 1.0
    for j in range(1, ordm1 + 1):
        saved = np.zeros(nx)
        for r in range(j):
            den = rdel[:, r] + ldel[:, j - 1 - r]
            nz = den != 0
            term = np.zeros(nx)
            term[nz] = b[nz, r] / den[nz]
            new_b = saved + rdel[:, r] * term
            new_saved = ldel[:, j - 1 - r] * term
            # den == 0 branch
            if r != 0:
                new_b[~nz] = saved[~nz]
            else:
                keep = (~nz) & (rdel[:, r] != 0.0)
                new_b[keep] = saved[keep]
                new_b[(~nz) & ~keep] = b[(~nz) & ~keep, r]
            new_saved[~nz] = 0.0
            b[:, r] = new_b
            saved = new_saved
        b[:, j] = saved
    return b


def spline_design(knots, x, ord: int = 4, derivs=0, outer_ok: bool = False) -> np.ndarray:
    """Port of ``splines::splineDesign``.

    Parameters
    ----------
    knots : array_like
        Full knot sequence (boundary knots repeated ``ord`` times, as R does).
    x : array_like
        Points at which to evaluate the basis.
    ord : int
        Order of the spline (degree + 1); 4 for cubic.
    derivs : int or array_like
        Derivative order for each ``x`` (recycled).
    outer_ok : bool
        If True, ``x`` outside the inner knots gives rows of zeros instead of
        raising an error (mirrors ``outer.ok=TRUE``).
    """
    knots = np.sort(np.asarray(knots, dtype=float).ravel())
    x = np.asarray(x, dtype=float).ravel()
    nk = knots.size
    nx = x.size
    derivs = np.atleast_1d(np.asarray(derivs, dtype=int))
    if nk <= 0:
        raise ValueError("must have at least 'ord' knots")
    if ord > nk or ord < 1:
        raise ValueError("'ord' must be positive integer, at most the number of knots")
    if not outer_ok and nk < 2 * ord - 1:
        raise ValueError(f"need at least 2*ord -1 (={2 * ord - 1}) knots")
    degree = ord - 1
    need_outer = bool(np.any((x < knots[ord - 1]) | (knots[nk - degree - 1] < x)))
    in_x = np.ones(nx, dtype=bool)
    x_out = False
    xx = x
    if need_outer:
        if not outer_ok:
            raise ValueError(
                f"the 'x' data must be in the range {knots[ord - 1]:g} to "
                f"{knots[nk - degree - 1]:g} unless you set 'outer_ok=True'")
        in_x = (knots[0] <= x) & (x <= knots[-1])
        x_out = not bool(np.all(in_x))
        if x_out:
            xx = x[in_x]
        dkn = np.diff(knots)[::-1]
        first_pos = int(np.argmax(dkn > 0)) + 1 if np.any(dkn > 0) else 0
        n_right = max(0, ord - first_pos) if first_pos else 0
        knots = np.concatenate((np.repeat(knots[0], degree), knots, np.repeat(knots[-1], n_right)))
    # evaluate
    if np.all(derivs == 0):
        curs = np.empty(xx.size, dtype=int)
        for i, xi in enumerate(xx):
            curs[i], _ = _set_cursor(knots, xi, ord)
        io = curs - ord
        temp = np.full((xx.size, ord), np.nan)
        ok = (io >= 0) & (io <= knots.size)
        if np.any(ok):
            temp[ok] = _basis_funcs_vectorised(knots, curs[ok], xx[ok], ord)
        offsets = io
    else:
        temp, offsets = _spline_basis(knots, ord, xx, derivs)
    ncoef = nk - ord
    design = np.zeros((nx, ncoef))
    rows = np.nonzero(in_x)[0] if (need_outer and x_out) else np.arange(nx)
    ii = np.repeat(rows, ord)
    jj = (np.arange(ord)[None, :] + offsets[:, None]).ravel()  # 0-based column
    vals = temp.ravel()
    if need_outer:
        jj = jj - degree
    okj = (jj >= 0) & (jj < ncoef)
    design[ii[okj], jj[okj]] = vals[okj]
    return design


# ----------------------------------------------------------------------------
# ns() and bs(): ports of splines::ns and splines::bs
# ----------------------------------------------------------------------------

def _qr_qty(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """R's ``qr.qty(qr(A), B)`` = t(Q) %*% B with the full Q of A."""
    q, _ = np.linalg.qr(A, mode="complete")
    return q.T @ B


def ns(x, df=None, knots=None, intercept: bool = False, boundary_knots=None):
    """Natural cubic spline basis, port of ``splines::ns``.

    Returns ``(basis, attrs)`` where ``attrs`` holds ``degree``, ``knots``,
    ``Boundary.knots`` and ``intercept`` exactly as R stores them.
    """
    x = np.asarray(x, dtype=float).ravel()
    nax = np.isnan(x)
    nas = bool(np.any(nax))
    if nas:
        x = x[~nax]
    bk_given = boundary_knots is not None
    if bk_given:
        boundary_knots = np.sort(np.asarray(boundary_knots, dtype=float).ravel())
        ol = x < boundary_knots[0]
        or_ = x > boundary_knots[1]
        outside = ol | or_
    else:
        if x.size == 1:
            boundary_knots = x * np.array([7, 9]) / 8
        else:
            boundary_knots = np.array([np.min(x), np.max(x)])
        ol = or_ = outside = np.zeros(x.size, dtype=bool)
    mk_knots = df is not None and knots is None
    if mk_knots:
        n_iknots = int(df) - 1 - int(intercept)
        if n_iknots < 0:
            warnings.warn(f"'df' was too small; have used {1 + int(intercept)}", stacklevel=2)
            n_iknots = 0
        if n_iknots > 0:
            probs = np.linspace(0, 1, n_iknots + 2)[1:-1]
            knots = quantile7(x[~outside], probs)
        else:
            knots = np.array([])
    else:
        knots = np.array([]) if knots is None else np.asarray(knots, dtype=float).ravel()
        if not np.all(np.isfinite(knots)):
            raise ValueError("non-finite knots")
        n_iknots = knots.size
    if mk_knots and knots.size:
        knots = _shove_knots(knots, boundary_knots)
    a_knots = np.sort(np.concatenate((np.repeat(boundary_knots, 4), knots)))
    if np.any(outside):
        basis = np.zeros((x.size, n_iknots + 4))
        if np.any(ol):
            k_pivot = boundary_knots[0]
            xl = np.column_stack((np.ones(ol.sum()), x[ol] - k_pivot))
            tt = spline_design(a_knots, np.repeat(k_pivot, 2), 4, [0, 1])
            basis[ol, :] = xl @ tt
        if np.any(or_):
            k_pivot = boundary_knots[1]
            xr = np.column_stack((np.ones(or_.sum()), x[or_] - k_pivot))
            tt = spline_design(a_knots, np.repeat(k_pivot, 2), 4, [0, 1])
            basis[or_, :] = xr @ tt
        inside = ~outside
        if np.any(inside):
            basis[inside, :] = spline_design(a_knots, x[inside], 4)
    else:
        basis = spline_design(a_knots, x, 4)
    const = spline_design(a_knots, boundary_knots, 4, [2, 2])
    if not intercept:
        const = const[:, 1:]
        basis = basis[:, 1:]
    basis = _qr_qty(const.T, basis.T).T[:, 2:]
    n_col = basis.shape[1]
    if nas:
        nmat = np.full((nax.size, n_col), np.nan)
        nmat[~nax, :] = basis
        basis = nmat
    attrs = {"degree": 3, "knots": knots, "boundary_knots": boundary_knots,
             "intercept": bool(intercept)}
    return basis, attrs


def _shove_knots(knots: np.ndarray, boundary_knots: np.ndarray, strict: bool = True) -> np.ndarray:
    """R's "shoving 'interior' knots matching boundary knots to inside".

    When *every* interior knot coincides with a boundary knot, ``splines::ns``
    stops but ``splines::bs`` only warns and leaves the knots alone; ``strict``
    selects between the two. This is reached by any exposure with a floor
    (zero-inflated rainfall, a detection limit), where the lower quantiles tie.
    """
    knots = knots.copy()
    lr_eq = np.isin([knots.min(), knots.max()], boundary_knots)
    if np.any(lr_eq):
        warnings.warn("shoving 'interior' knots matching boundary knots to inside", stacklevel=2)
    for side, piv in ((0, boundary_knots[0]), (1, boundary_knots[1])):
        if not lr_eq[side]:
            continue
        where = "left" if side == 0 else "right"
        i = knots == piv
        if np.all(i):
            if strict:
                raise ValueError(f"all interior knots match {where} boundary knot")
            warnings.warn(f"all interior knots match {where} boundary knot", stacklevel=2)
            continue
        if side == 0:
            knots[i] = knots[i] + (knots[knots > piv].min() - piv) / 8
        else:
            knots[i] = knots[i] - (piv - knots[knots < piv].max()) / 8
    return knots


def bs(x, df=None, knots=None, degree: int = 3, intercept: bool = False,
       boundary_knots=None):
    """B-spline basis, port of ``splines::bs``.

    Returns ``(basis, attrs)`` with ``degree``, ``knots``, ``Boundary.knots``
    and ``intercept``.
    """
    degree = int(degree)
    ord = 1 + degree
    if ord <= 1:
        raise ValueError("'degree' must be integer >= 1")
    x = np.asarray(x, dtype=float).ravel()
    nax = np.isnan(x)
    nas = bool(np.any(nax))
    if nas:
        x = x[~nax]
    if boundary_knots is not None:
        boundary_knots = np.sort(np.asarray(boundary_knots, dtype=float).ravel())
        ol = x < boundary_knots[0]
        or_ = x > boundary_knots[1]
        outside = ol | or_
    else:
        boundary_knots = np.array([np.min(x), np.max(x)])
        ol = or_ = outside = np.zeros(x.size, dtype=bool)
    mk_knots = df is not None and knots is None
    if mk_knots:
        n_iknots = int(df) - ord + (1 - int(intercept))
        if n_iknots < 0:
            warnings.warn(f"'df' was too small; have used {ord - (1 - int(intercept))}", stacklevel=2)
            n_iknots = 0
        if n_iknots > 0:
            probs = np.linspace(0, 1, n_iknots + 2)[1:-1]
            knots = quantile7(x[~outside], probs)
        else:
            knots = np.array([])
    else:
        knots = np.array([]) if knots is None else np.asarray(knots, dtype=float).ravel()
        if not np.all(np.isfinite(knots)):
            raise ValueError("non-finite knots")
    if mk_knots and knots.size:
        knots = _shove_knots(knots, boundary_knots, strict=False)
    a_knots = np.sort(np.concatenate((np.repeat(boundary_knots, ord), knots)))
    if np.any(outside):
        derivs = np.arange(0, degree + 1)
        scalef = np.array([math.gamma(i) for i in range(1, ord + 1)])
        basis = np.zeros((x.size, a_knots.size - degree - 1))
        e = 1 / 4
        if np.any(ol):
            k_pivot = (1 - e) * boundary_knots[0] + e * a_knots[ord]
            xl = np.column_stack((np.ones(ol.sum()),
                                  np.power.outer(x[ol] - k_pivot, np.arange(1, degree + 1))))
            tt = spline_design(a_knots, np.repeat(k_pivot, ord), ord, derivs)
            basis[ol, :] = xl @ (tt / scalef[:, None])
        if np.any(or_):
            k_pivot = (1 - e) * boundary_knots[1] + e * a_knots[a_knots.size - ord - 1]
            xr = np.column_stack((np.ones(or_.sum()),
                                  np.power.outer(x[or_] - k_pivot, np.arange(1, degree + 1))))
            tt = spline_design(a_knots, np.repeat(k_pivot, ord), ord, derivs)
            basis[or_, :] = xr @ (tt / scalef[:, None])
        inside = ~outside
        if np.any(inside):
            basis[inside, :] = spline_design(a_knots, x[inside], ord)
    else:
        basis = spline_design(a_knots, x, ord)
    if not intercept:
        basis = basis[:, 1:]
    n_col = basis.shape[1]
    if nas:
        nmat = np.full((nax.size, n_col), np.nan)
        nmat[~nax, :] = basis
        basis = nmat
    attrs = {"degree": degree, "knots": knots, "boundary_knots": boundary_knots,
             "intercept": bool(intercept)}
    return basis, attrs
