"""Penalised DLNMs (REML/ML smoothing) versus mgcv::gam with paraPen."""

import json

import numpy as np
import pytest

import dlnmpy as dl
from dlnmpy.penalized import fit_pgam, fit_pglm

from conftest import FIX, assert_close

pytest.importorskip("statsmodels")


@pytest.fixture(scope="module")
def R():
    with open(FIX / "penalized.json") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def setup(chicago):
    cb7 = dl.crossbasis(chicago.temp, lag=10, argvar={"fun": "ps", "df": 8}, arglag={"fun": "ps", "df": 5})
    nst = dl.onebasis(chicago.time, "ns", df=98)
    X = chicago.join(cb7.to_dataframe("cb")).join(nst.to_dataframe("nst"))
    f = "death ~ " + " + ".join(list(cb7.to_dataframe("cb").columns) + list(nst.to_dataframe("nst").columns)) + " + C(dow)"
    return cb7, X, f


def _cb(m):
    return m.params.filter(like="cb_").to_numpy()


def _vp(m):
    return m.cov_params().filter(like="cb_").filter(like="cb_", axis=0).to_numpy()


def test_criterion_at_fixed_sp(R, setup):
    """The REML score at given smoothing parameters equals mgcv's gcv.ubre."""
    cb7, X, f = setup
    pen = dl.cbpen(cb7)
    for fam, key in [("poisson", "poisson_fixed"), ("quasipoisson", "quasi_fixed")]:
        r = R["ps"][key]
        m = fit_pgam(f, X, {"cb": pen}, family=fam, sp=r["spfix"])
        assert abs(m.reml - r["score"]) < 1e-5
        assert_close(_cb(m), r["coef"], atol=1e-8)
        assert_close(_vp(m), r["Vp"], atol=1e-9)
        assert abs(m.scale - r["scale"]) < 1e-6
        assert abs(m.edf_by("cb_") - r["edf_cb"]) < 1e-6
        assert abs(m.deviance - r["deviance"]) < 1e-6


@pytest.mark.parametrize("fam,key,meth", [("poisson", "poisson", "reml"), ("quasipoisson", "quasi", "reml"),
                                          ("poisson", "ml", "ml")])
def test_smoothing_parameter_selection(R, setup, fam, key, meth):
    cb7, X, f = setup
    r = R["ps"][key]
    m = fit_pgam(f, X, {"cb": dl.cbpen(cb7)}, family=fam, method=meth)
    assert m.converged
    assert abs(m.reml - r["score"]) < 1e-4
    np.testing.assert_allclose(m.sp, r["sp"], rtol=1e-3)
    assert_close(_cb(m), r["coef"], atol=2e-4)
    assert_close(_vp(m), r["Vp"], atol=2e-4)
    assert abs(m.scale - r["scale"]) < 1e-4
    assert abs(m.edf_by("cb_") - r["edf_cb"]) < 1e-2
    if key == "poisson":
        p = dl.crosspred(cb7, m, cen=20, by=1, name="cb")
        pr = R["ps"]["pred"]
        assert_close(p.predvar, pr["predvar"])
        assert_close(p.allfit, pr["allfit"], atol=1e-4)
        assert_close(p.allse, pr["allse"], atol=1e-4)
        assert p.model_link == "log"


def test_cr_crossbasis_and_added_lag_penalty(chicago, R, setup):
    nst = dl.onebasis(chicago.time, "ns", df=98)
    cb9 = dl.crossbasis(chicago.temp, lag=21, argvar={"fun": "cr", "df": 6}, arglag={"fun": "cr", "df": 5})
    X = chicago.join(cb9.to_dataframe("cb")).join(nst.to_dataframe("nst"))
    f = "death ~ " + " + ".join(list(cb9.to_dataframe("cb").columns) + list(nst.to_dataframe("nst").columns)) + " + C(dow)"
    r = R["cr"]["quasi"]
    m = fit_pgam(f, X, {"cb": dl.cbpen(cb9)}, family="quasipoisson")
    assert abs(m.reml - r["score"]) < 1e-4
    np.testing.assert_allclose(m.sp, r["sp"], rtol=1e-4)
    assert_close(_cb(m), r["coef"], atol=1e-6)
    assert abs(m.scale - r["scale"]) < 1e-5
    cb7, X7, f7 = setup
    pen = dl.cbpen(cb7, add_slag=np.ones(5))
    r = R["ps_addslag"]["poisson"]
    assert_close(pen["Slag2"], r and R["ps_addslag"]["pen"]["Slag2"], atol=1e-12)
    m = fit_pgam(f7, X7, {"cb": pen}, family="poisson")
    assert abs(m.reml - r["score"]) < 1e-4
    assert_close(_cb(m), r["coef"], atol=5e-3)
    assert abs(m.edf_by("cb_") - r["edf_cb"]) < 1e-2
    assert "smoothing parameters" in m.summary()


def test_fit_pglm_gaussian_and_errors():
    rng = np.random.default_rng(0)
    x = np.linspace(0, 1, 200)
    y = np.sin(2 * np.pi * x) + rng.normal(0, 0.3, 200)
    B, attrs = dl.ps(x, df=12)
    S = np.zeros((13, 13))
    S[1:, 1:] = attrs["S"]
    Xd = np.column_stack([np.ones(200), B])
    m = fit_pglm(y, Xd, [S], family="gaussian")
    assert m.link == "identity" and 2 < m.edf_total < 12
    fitted = Xd @ m.params.to_numpy()
    assert np.mean((fitted - np.sin(2 * np.pi * x)) ** 2) < 0.02
    with pytest.raises(ValueError):
        fit_pglm(y, Xd, [S[:5, :5]], family="gaussian")
    with pytest.raises(ValueError):
        fit_pglm(y, Xd, [S], family="gaussian", method="gcv")
