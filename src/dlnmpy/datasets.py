"""Example datasets shipped with the R package, exported to CSV.

- ``chicagoNMMAPS``: daily mortality, weather and pollution in Chicago,
  1987-2000 (from the NMMAPS study; Samet et al. 2000). Columns: date,
  time, year, month, doy, dow, death, cvd, resp, temp, dptp, rhum, pm10, o3.
- ``drug``: simulated trial with time-varying drug exposure histories.
- ``nested``: simulated nested case-control study with occupational
  exposure histories.
"""

from __future__ import annotations

from importlib import resources

import pandas as pd

__all__ = ["chicago_nmmaps", "drug", "nested", "simulate_cities"]


def _read(name: str) -> pd.DataFrame:
    with resources.files("dlnmpy.data").joinpath(name).open("rb") as f:
        return pd.read_csv(f)


def chicago_nmmaps() -> pd.DataFrame:
    """Chicago NMMAPS daily time series (5114 days)."""
    df = _read("chicagoNMMAPS.csv")
    df["date"] = pd.to_datetime(df["date"])
    return df


def drug() -> pd.DataFrame:
    return _read("drug.csv")


def nested() -> pd.DataFrame:
    return _read("nested.csv")


def simulate_cities(n_cities: int = 12, n_days: int = 2000, seed: int = 1, cen: float = 18.0) -> pd.DataFrame:
    """Simulate daily temperature and counts for several locations with a
    known exposure-lag-response association, for examples and tests.

    The true overall cumulative curve is ``slope * (t - cen)^2`` with a
    1.6-fold steeper heat side, distributed over lags 0-10 with exponentially
    decaying weights; ``slope`` varies with each city's mean temperature so
    that a meta-regression on mean temperature explains the heterogeneity.
    Returns a DataFrame with columns ``city``, ``time``, ``tmean``, ``y``.
    """
    import numpy as np
    rng = np.random.default_rng(seed)
    w = np.exp(-np.arange(11) / 3)
    w /= w.sum()
    frames = []
    for k in range(1, n_cities + 1):
        shift = rng.normal(0, 3)
        t = np.arange(1, n_days + 1)
        tmean = 15 + shift + 9 * np.sin(2 * np.pi * t / 365.25) + rng.normal(0, 3, n_days)
        slope = 0.0025 * (1 + 0.05 * shift)

        def f(x):
            return slope * (x - cen) ** 2 * np.where(x > cen, 1.6, 1.0)

        Q = np.column_stack([np.r_[[np.nan] * lag, tmean[: n_days - lag]] for lag in range(11)])
        eta = np.nansum(f(Q) * w, axis=1)
        eta[np.isnan(Q).any(axis=1)] = np.nan
        mu = np.exp(np.log(30) + 0.1 * np.sin(2 * np.pi * t / 365.25 + 1) + eta)
        y = rng.poisson(np.where(np.isnan(mu), 30, mu)).astype(float)
        y[np.isnan(mu)] = np.nan
        frames.append(pd.DataFrame({"city": k, "time": t, "tmean": tmean, "y": y}))
    return pd.concat(frames, ignore_index=True)
