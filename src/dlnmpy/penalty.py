"""Penalty matrices for penalised DLNMs (port of ``cbPen.R``).

For a cross-basis built with ``ps`` (P-spline) bases the smoothness
penalties on the predictor and lag dimensions are expanded to the
dimension of the cross-basis: ``S_var (x) I`` and ``I (x) S_lag`` (Kronecker
products), each rescaled by its largest eigenvalue. Additional lag penalties
(``add_slag``) can be supplied as vectors (diagonal) or matrices.

The result is a dict of penalty matrices with their ranks, ready to be used
with a penalised likelihood or a mixed-model fitter (mgcv's ``gam`` plays
that role in R; in Python see the notes in ``docs/penalized.md``).
"""

from __future__ import annotations

import numpy as np

from .core import CrossBasis, OneBasis

__all__ = ["cbpen", "findrank"]


def findrank(X: np.ndarray) -> int:
    ev = np.linalg.eigvalsh(X)
    return int(np.sum(ev > ev.max() * np.finfo(float).eps * 10))


def _rescale(X: np.ndarray) -> np.ndarray:
    return X / np.linalg.eigvalsh(X).max()


def cbpen(cb, sp=-1, add_slag=None) -> dict:
    """Penalty matrices for a cross-basis with penalised bases."""
    if isinstance(cb, OneBasis):
        df = (cb.ncol, 1)
        argvar = {"fun": cb.fun, **cb.attrs}
        arglag = {"fun": "strata"}
        one = True
    elif isinstance(cb, CrossBasis):
        df = cb.df
        argvar, arglag = cb.argvar, cb.arglag
        one = False
    else:
        raise TypeError("first argument must be a CrossBasis or OneBasis object")

    funs = (argvar.get("fun"), arglag.get("fun"))
    fx = (funs[0] not in ("ps", "cr") or bool(argvar.get("fx", False)),
          funs[1] not in ("ps", "cr") or bool(arglag.get("fx", False)))
    S = {}
    if not fx[0]:
        S["Svar"] = np.kron(np.asarray(argvar["S"]), np.eye(df[1]))
    if not fx[1]:
        S["Slag"] = np.kron(np.eye(df[0]), np.asarray(arglag["S"]))
    S = {k: _rescale(v) for k, v in S.items()}

    if add_slag is not None:
        if one:
            raise ValueError("penalties on lag not allowed for class 'OneBasis'")
        slist = add_slag if isinstance(add_slag, (list, tuple)) else [add_slag]
        for i, X in enumerate(slist):
            X = np.asarray(X, dtype=float)
            if X.ndim == 1:
                X = np.diag(X)
            if X.shape != (df[1], df[1]):
                raise ValueError("terms in add_slag with dimensions not consistent with basis for lag")
            X = X / np.linalg.eigvalsh(X).max()
            S[f"Slag{i + 2}"] = np.kron(np.eye(df[0]), X)

    if not S:
        raise ValueError("no penalization defined")
    rank = {k: findrank(v) for k, v in S.items()}
    sp = np.atleast_1d(np.asarray(sp, dtype=float))
    if sp.size == 1:
        sp = np.repeat(sp, len(S))
    if sp.size != len(S):
        raise ValueError("'sp' must be numeric and consistent with number of penalty terms")
    return {**S, "rank": rank, "sp": dict(zip(S.keys(), sp.tolist()))}
