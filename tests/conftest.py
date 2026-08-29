import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

FIX = Path(__file__).parent / "fixtures"


@pytest.fixture(scope="session")
def cases():
    with open(FIX / "cases.json") as f:
        return json.load(f)


@pytest.fixture(scope="session")
def onebasis_fixtures():
    with open(FIX / "onebasis.json") as f:
        return json.load(f)


@pytest.fixture(scope="session")
def helpers():
    with open(FIX / "helpers.json") as f:
        return json.load(f)


@pytest.fixture(scope="session")
def chicago():
    import dlnmpy as dl
    return dl.datasets.chicago_nmmaps()


def rcsv(name: str) -> np.ndarray:
    return pd.read_csv(FIX / f"{name}.csv").to_numpy(dtype=float)


def assert_close(a, b, atol=1e-10, rtol=0.0, msg=""):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    assert a.shape == b.shape, f"{msg} shape {a.shape} != {b.shape}"
    assert np.array_equal(np.isnan(a), np.isnan(b)), f"{msg} NaN pattern differs"
    np.testing.assert_allclose(np.nan_to_num(a), np.nan_to_num(b), atol=atol, rtol=rtol, err_msg=msg)
