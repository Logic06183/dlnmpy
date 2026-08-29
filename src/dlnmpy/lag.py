"""Lag utilities: lag ranges, lagged matrices and exposure histories."""

from __future__ import annotations

import numpy as np

from ._rcompat import seq

__all__ = ["mklag", "seqlag", "lag_matrix", "exphist"]


def mklag(lag) -> np.ndarray:
    """Normalise a lag specification to an integer pair ``[lag0, lag1]``.

    A scalar ``L >= 0`` means ``[0, L]``; a negative scalar ``L`` means
    ``[L, 0]``. A pair is returned as given (rounded).
    """
    lag = np.atleast_1d(np.asarray(lag, dtype=float)).ravel()
    if lag.size > 2 or lag.size == 0:
        raise ValueError("'lag' must be an integer vector of length 1 or 2")
    if lag.size == 1:
        lag = np.array([lag[0], 0.0]) if lag[0] < 0 else np.array([0.0, lag[0]])
    if lag[1] - lag[0] < 0:
        raise ValueError("lag[0] must be <= lag[1]")
    return np.rint(lag).astype(int)


def seqlag(lag, by: float = 1.0) -> np.ndarray:
    """Sequence of lag values from ``lag[0]`` to ``lag[1]`` in steps ``by``."""
    lag = np.asarray(lag)
    if by == 1:
        return np.arange(lag[0], lag[1] + 1, dtype=float)
    return seq(float(lag[0]), float(lag[1]), float(by))


def lag_matrix(v, k, group=None) -> np.ndarray:
    """Matrix of lagged copies of the series ``v`` (port of ``tsModel::Lag``).

    Column ``j`` holds ``v`` shifted by ``k[j]`` positions (positive lags
    shift forward in time, so the first ``k[j]`` rows are NaN). With
    ``group``, the series is broken at group boundaries and lagged within
    each group, as for seasonal or multi-city series stacked vertically.
    """
    v = np.asarray(v, dtype=float).ravel()
    k = np.atleast_1d(np.asarray(k)).astype(int)
    if np.max(np.abs(k)) >= v.size:
        raise ValueError("largest lag in 'k' must be less than 'length(v)'")

    def lag_f(x: np.ndarray) -> np.ndarray:
        n = x.size
        out = np.full((n, k.size), np.nan)
        for i, lg in enumerate(k):
            if lg > 0:
                if lg < n:
                    out[lg:, i] = x[: n - lg]
            elif lg < 0:
                if -lg < n:
                    out[: n + lg, i] = x[-lg:]
            else:
                out[:, i] = x
        return out

    if group is None:
        return lag_f(v)
    group = np.asarray(group)
    out = np.full((v.size, k.size), np.nan)
    for g in np.unique(group):
        idx = np.nonzero(group == g)[0]
        out[idx, :] = lag_f(v[idx])
    return out


def exphist(exp, times=None, lag=None, fill: float = 0.0) -> np.ndarray:
    """Build a matrix of exposure histories (port of ``dlnm::exphist``).

    Given an exposure profile ``exp`` observed at consecutive integer times
    ``1..len(exp)``, returns for each of ``times`` a row with the exposures
    experienced at lags ``lag[0]..lag[1]`` before it. Times outside the
    observed profile are padded with ``fill``.
    """
    exp = np.asarray(exp, dtype=float).ravel()
    lag = np.array([0, exp.size - 1]) if lag is None else mklag(lag)
    times = np.arange(1, exp.size + 1) if times is None else np.rint(np.atleast_1d(times)).astype(int)
    left = max(0, lag[1] + 1 - int(times.min()))
    right = max(0, int(times.max()) - exp.size - lag[0])
    ext = np.concatenate((np.full(left, fill), exp, np.full(right, fill)))
    rows = []
    for t in times:
        lo = t - lag[1] + left  # 1-based indices into ext
        hi = t - lag[0] + left
        rows.append(ext[lo - 1: hi][::-1])  # R: exp[seq(x[1], x[2])] with x = t - (lag - left)
    hist = np.vstack(rows)
    return hist
