"""attrdl / findmin / mmt / attr_table against Gasparrini's reference R functions."""

import json

import numpy as np
import pytest

import dlnmpy as dl
from dlnmpy.attribution import _coefsim_from_normals, attr_table, attrdl, findmin, mmt

from conftest import FIX, assert_close


@pytest.fixture(scope="module")
def setup(chicago):
    with open(FIX / "attribution.json") as f:
        r = json.load(f)
    cb = dl.crossbasis(chicago.temp, lag=21,
                       argvar={"fun": "bs", "degree": 2, "knots": dl.quantile7(chicago.temp, [0.10, 0.75, 0.90])},
                       arglag={"knots": dl.logknots(21, 3)})
    kw = dict(coef=r["coef"], vcov=r["vcov"], model_link="log")
    return r, cb, kw


def test_findmin_and_mmt(chicago, setup):
    r, cb, kw = setup
    q = r["q"]
    assert findmin(cb, from_=q[0], to=q[1], by=0.1, **{k: kw[k] for k in ("coef", "vcov")}) == r["mmt"]
    # simulated minima with the same normal draws as R
    # eigenvector signs are arbitrary, so R's normal draws cannot be replayed
    # exactly; check the simulated coefficients have the right covariance and
    # then use R's own coefficient draws for the exact comparison
    k = len(r["coef"])
    L = _coefsim_from_normals(np.zeros(k), np.array(r["vcov"]), np.eye(k))
    assert_close(L @ L.T, r["vcov"], atol=1e-14)
    coefsim = np.array(r["coefsim"])
    mins = findmin(cb, from_=q[0], to=q[1], by=0.1, sim=True, coefsim=coefsim, coef=r["coef"], vcov=r["vcov"])
    assert_close(mins, r["minsim"], atol=1e-9)
    res = mmt(cb, coef=r["coef"], vcov=r["vcov"], x=chicago.temp, coefsim=coefsim)
    assert res.mmt == r["mmt"]
    assert res.low <= res.mmt <= res.high
    assert 0 <= res.percentile <= 100
    assert "MMT" in repr(res)


def test_attrdl_point_estimates(chicago, setup):
    r, cb, kw = setup
    x, y, cen = chicago.temp, chicago.death, r["mmt"]
    assert np.isclose(attrdl(x, cb, y, cen=cen, **kw), r["af_back"], atol=1e-12)
    assert np.isclose(attrdl(x, cb, y, type="an", cen=cen, **kw), r["an_back"], atol=1e-8)
    assert np.isclose(attrdl(x, cb, y, dir="forw", cen=cen, **kw), r["af_forw"], atol=1e-12)
    assert np.isclose(attrdl(x, cb, y, type="an", dir="forw", cen=cen, **kw), r["an_forw"], atol=1e-8)
    assert_close(attrdl(x, cb, y, type="an", cen=cen, tot=False, **kw)[:200], r["an_back_daily"], atol=1e-9)
    assert_close(attrdl(x, cb, y, dir="forw", cen=cen, tot=False, **kw)[:200], r["af_forw_daily"], atol=1e-12)
    assert np.isclose(attrdl(x, cb, y, type="an", cen=cen, range=(-100, cen), **kw), r["an_cold"], atol=1e-8)
    assert np.isclose(attrdl(x, cb, y, type="an", cen=cen, range=(cen, 100), **kw), r["an_heat"], atol=1e-8)
    p975 = dl.quantile7(x, [0.975])[0]
    assert np.isclose(attrdl(x, cb, y, type="an", cen=cen, range=(p975, 100), **kw), r["an_extreme_heat"], atol=1e-8)
    # cold + heat is close to, but not exactly, the total (AF = 1 - exp(-sum) is
    # not additive); R gives the same
    assert abs(r["an_cold"] + r["an_heat"] - r["an_back"]) / r["an_back"] < 0.01


def test_attrdl_matrix_inputs_and_reduced(chicago, setup):
    r, cb, kw = setup
    x, y, cen = chicago.temp.to_numpy(), chicago.death.to_numpy(), r["mmt"]
    Q = dl.lag_matrix(x, np.arange(22))
    assert np.isclose(attrdl(Q, cb, y, type="an", cen=cen, **kw), r["an_back_matrix"], atol=1e-8)
    fut = dl.lag_matrix(y, -np.arange(22))
    assert np.isclose(attrdl(x, cb, fut, dir="forw", cen=cen, **kw), r["af_forw_casemat"], atol=1e-12)
    red = dl.crossreduce(cb, cen=cen, **kw)
    assert_close(red.coef, r["red_coef"], atol=1e-10)
    assert np.isclose(attrdl(x, cb, y, coef=red.coef, vcov=red.vcov, model_link="log", dir="forw", cen=cen),
                      r["af_reduced_forw"], atol=1e-12)
    with pytest.raises(ValueError):
        attrdl(x, cb, y, coef=red.coef, vcov=red.vcov, model_link="log", dir="back", cen=cen)
    with pytest.raises(ValueError):
        attrdl(x, cb, y, cen=cen, tot=False, sim=True, **kw)


def test_attrdl_simulation(chicago, setup):
    r, cb, kw = setup
    x, y, cen = chicago.temp, chicago.death, r["mmt"]
    coefsim = np.array(r["coefsim"])
    ansim = attrdl(x, cb, y, type="an", cen=cen, sim=True, coefsim=coefsim, **kw)
    assert_close(ansim, r["an_back_sim_fixed"], atol=1e-6)
    # R's own random draws: agree in distribution (same nsim=200, so loose)
    rsim = np.array(r["an_back_sim"])
    assert abs(np.mean(ansim) - np.mean(rsim)) < 3 * np.std(rsim) / np.sqrt(len(rsim)) * 3


def test_attr_table(chicago, setup):
    r, cb, kw = setup
    tab = attr_table(chicago.temp, cb, chicago.death, cen=r["mmt"], nsim=300, seed=1, **kw)
    assert list(tab.component) == ["total", "cold", "heat", "extreme cold", "moderate cold", "moderate heat", "extreme heat"]
    t = tab.set_index("component")
    assert np.isclose(t.loc["total", "an"], r["an_back"], atol=1e-6)
    assert abs(t.loc["cold", "an"] + t.loc["heat", "an"] - t.loc["total", "an"]) / t.loc["total", "an"] < 0.01
    assert np.isclose(t.loc["cold", "an"], r["an_cold"], atol=1e-6)
    assert (t.an_low <= t.an).all() and (t.an <= t.an_high).all()


def test_attrdl_from_model(chicago, setup):
    r, cb, kw = setup
    nst = dl.onebasis(chicago.time, "ns", df=98)
    X = chicago.join(cb.to_dataframe("cb")).join(nst.to_dataframe("nst"))
    fit = dl.fit_glm("death ~ " + " + ".join(list(cb.to_dataframe("cb").columns) + list(nst.to_dataframe("nst").columns)) + " + C(dow)", X)
    assert np.isclose(attrdl(chicago.temp, cb, chicago.death, fit, cen=r["mmt"], name="cb"), r["af_back"], atol=1e-8)
    assert findmin(cb, fit, from_=r["q"][0], to=r["q"][1], by=0.1, name="cb") == r["mmt"]
