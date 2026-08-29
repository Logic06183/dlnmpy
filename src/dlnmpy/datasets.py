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

__all__ = ["chicago_nmmaps", "drug", "nested"]


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
