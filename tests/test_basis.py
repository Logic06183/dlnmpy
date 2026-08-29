"""Basis functions must reproduce R's dlnm/splines output to machine precision."""

import json

import numpy as np
import pytest

import dlnmpy as dl
from dlnmpy._rcompat import median, quantile7

from conftest import FIX, assert_close


def _spec_names():
    with open(FIX / "onebasis.json") as f:
        return list(json.load(f).keys())


@pytest.mark.parametrize("name", _spec_names())
def test_onebasis_matches_r(onebasis_fixtures, name):
    c = onebasis_fixtures[name]
    spec = dict(c["spec"])
    fun = spec.pop("fun")
    ob = dl.onebasis(np.array(c["x"], dtype=float), fun, **spec)
    assert_close(ob.matrix, c["basis"], atol=1e-12, msg=f"{name} basis")
    # transformation re-applied on new values (prediction stage), including
    # values outside the original range
    pred = ob.transform(np.array(c["xpred"], dtype=float))
    assert_close(pred, c["basis_pred"], atol=1e-12, msg=f"{name} prediction basis")


def test_onebasis_attributes(onebasis_fixtures):
    c = onebasis_fixtures["ns_df5"]
    ob = dl.onebasis(np.array(c["x"]), "ns", df=5)
    assert_close(ob.attrs["knots"], c["attributes"]["knots"], atol=1e-12)
    assert_close(ob.attrs["boundary_knots"], c["attributes"]["Boundary.knots"], atol=1e-12)
    c = onebasis_fixtures["strata_df3"]
    ob = dl.onebasis(np.array(c["x"]), "strata", df=3)
    assert_close(ob.attrs["breaks"], c["attributes"]["breaks"], atol=1e-12)
    c = onebasis_fixtures["ps_df10"]
    ob = dl.onebasis(np.array(c["x"]), "ps", df=10)
    assert_close(ob.attrs["S"], c["attributes"]["S"], atol=1e-12)
    assert_close(ob.attrs["knots"], c["attributes"]["knots"], atol=1e-10)


def test_ns_with_missing_values():
    with open(FIX / "onebasis_na.json") as f:
        c = json.load(f)
    x = np.array([np.nan if v is None else v for v in c["x"]], dtype=float)
    ob = dl.onebasis(x, "ns", df=4)
    ref = np.array([[np.nan if v is None else v for v in row] for row in c["basis"]], dtype=float)
    assert_close(ob.matrix, ref, atol=1e-12)


def test_r_style_argument_names(chicago):
    a = dl.onebasis(chicago.temp, "thr", **{"thr.value": 25})
    b = dl.onebasis(chicago.temp, "thr", thr_value=25)
    c = dl.onebasis(chicago.temp, "thr", thr=25)
    assert_close(a.matrix, b.matrix)
    assert_close(a.matrix, c.matrix)
    a = dl.onebasis(chicago.temp, "ns", knots=[0, 20], **{"Boundary.knots": [-20, 30]})
    b = dl.onebasis(chicago.temp, "ns", knots=[0, 20], boundary_knots=[-20, 30])
    assert_close(a.matrix, b.matrix)


def test_quantile_and_median():
    x = np.array([3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0])
    # R: quantile(x, c(.1,.25,.5,.9)) -> 1.0 1.75 3.5 6.9
    assert_close(quantile7(x, [0.1, 0.25, 0.5, 0.9]), [1.0, 1.75, 3.5, 6.9], atol=1e-12)
    assert median(x) == 3.5
    assert median([1, 2, np.nan, 4]) == 2.0


def test_custom_callable_basis(chicago):
    def mybasis(x, power=2):
        x = np.asarray(x, dtype=float)
        return np.column_stack((x, x ** power)), {"power": power}

    ob = dl.onebasis(chicago.temp, mybasis, power=3)
    assert ob.shape == (len(chicago), 2)
    assert_close(ob.transform([2.0]), [[2.0, 8.0]])
