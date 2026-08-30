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


def test_small_df_warns_like_r(chicago):
    with pytest.warns(UserWarning, match="too small"):
        dl.onebasis(chicago.temp, "ns", df=0)
    with pytest.warns(UserWarning, match="too small"):
        dl.onebasis(chicago.temp, "bs", df=2)


def test_bs_warns_where_ns_stops_on_tied_boundary_knots():
    """splines::ns stops when every interior knot lands on a boundary knot;
    splines::bs only warns. Any exposure with a floor (zero-inflated rainfall,
    a detection limit) reaches this."""
    import warnings

    x = np.concatenate([np.zeros(70), np.linspace(0.5, 50, 30)])
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        b = dl.onebasis(x, "bs", df=5)
    assert b.matrix.shape == (100, 5)
    assert_close(b.matrix[0], [0, 1, 0, 0, 0], atol=1e-12)
    assert any("boundary knot" in str(m.message) for m in w)

    # not all knots tie here, so both languages only shove and warn
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        assert dl.onebasis(x, "ns", df=5).matrix.shape == (100, 5)
    assert any("shoving" in str(m.message) for m in w)

    # all of them tie here: ns stops, bs does not
    xz = np.concatenate([np.zeros(90), np.linspace(0.5, 50, 10)])
    with pytest.raises(ValueError, match="interior knots match"):
        dl.onebasis(xz, "ns", df=5)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert dl.onebasis(xz, "bs", df=5).matrix.shape == (100, 5)


def test_strata_rejects_tied_breaks_and_is_na_outside_the_edges():
    """R's cut() raises on non-unique breaks. Without the check the basis has
    an empty stratum, hence a zero column and a singular design, silently."""
    x = np.concatenate([np.zeros(60), np.linspace(1, 50, 40)])
    with pytest.raises(ValueError, match="not unique"):
        dl.onebasis(x, "strata", df=4)
    with pytest.raises(ValueError, match="not unique"):
        dl.onebasis(np.full(50, 3.0), "strata", df=3)
    # cut() gives NA, not a zero row, outside the outer edges
    xh = np.arange(1, 9, dtype=float) * 1e12
    assert np.all(np.isnan(dl.onebasis(xh, "strata", df=3).matrix[-1]))
    # and an ordinary case is untouched
    b = dl.onebasis(np.linspace(0, 10, 100), "strata", df=3)
    assert b.matrix.shape == (100, 3) and np.linalg.matrix_rank(b.matrix) == 3


def test_thr_default_threshold_propagates_na_as_in_r():
    """dlnm's thr uses median(x) without na.rm, so a missing value makes the
    threshold NA. Dropping NaN would silently give a different threshold, and
    different results, from the same R code."""
    b = dl.onebasis(np.array([1, 2, np.nan, 4, 5, 6, 7.0]), "thr")
    assert np.isnan(b.attrs["thr_value"]).all()
    assert np.isnan(b.matrix).all()
    # no missing values: unchanged
    b2 = dl.onebasis(np.array([1, 2, 3, 4, 5.0]), "thr")
    assert_close(b2.attrs["thr_value"], [3.0], atol=0)
