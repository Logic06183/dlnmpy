import json

import numpy as np
import pytest
from conftest import FIX, assert_close

import dlnmpy as dl
from dlnmpy.uncertainty import bootstrap, bootstrap_ci, empirical_ci, model_grid, qaic, simulate_pred


@pytest.fixture(scope="module")
def fitted(chicago, cases):
    cb = dl.crossbasis(chicago.temp, lag=21, argvar={"fun": "ns", "df": 4}, arglag={"fun": "ns", "df": 4})
    m = cases["ex5"]["model"]
    return cb, np.array(m["coef"]), np.array(m["vcov"])


def test_bootstrap_matches_delta_method(fitted):
    cb, coef, vcov = fitted
    pred = dl.crosspred(cb, coef=coef, vcov=vcov, model_link="log", cen=21, by=1)
    sim = simulate_pred(cb, coef, vcov, nsim=3000, seed=1, model_link="log", cen=21, by=1)
    assert sim["allfit"].shape == (3000, pred.predvar.size)
    # simulated mean and sd of the linear predictor agree with the delta method
    np.testing.assert_allclose(sim["allfit"].mean(axis=0), pred.allfit, atol=3 * pred.allse.max() / np.sqrt(3000) + 1e-3)
    np.testing.assert_allclose(sim["allfit"].std(axis=0), pred.allse, rtol=0.1)
    low, high = empirical_ci(np.exp(sim["allfit"]))
    assert np.all(low <= pred.allRRfit + 1e-12) and np.all(high >= pred.allRRfit - 1e-12)
    out = bootstrap_ci(pred, nsim=500, seed=2, basis=cb)
    assert out["low"].shape == pred.predvar.shape
    # a derived quantity: ratio of RR at 33 vs -20 (non-linear in the coefficients)
    i33, im20 = pred._var_index(33), pred._var_index(-20)
    ratio = bootstrap(lambda c: np.exp(dl.crosspred(cb, coef=c, vcov=vcov, model_link="log", cen=21, at=[-20, 33]).allfit[1]
                                       - dl.crosspred(cb, coef=c, vcov=vcov, model_link="log", cen=21, at=[-20, 33]).allfit[0]),
                      coef, vcov, nsim=200, seed=3)
    assert ratio.shape == (200,)
    assert np.isclose(np.median(ratio), pred.allRRfit[i33] / pred.allRRfit[im20], rtol=0.2)


def test_qaic_and_model_grid(chicago):
    pytest.importorskip("statsmodels")
    nst = dl.onebasis(chicago.time, "ns", df=98)

    def fit(spec):
        cb = dl.crossbasis(chicago.temp, lag=spec["lag"], argvar={"fun": "ns", "df": spec["df_var"]}, arglag={"fun": "ns", "df": 3})
        X = chicago.join(cb.to_dataframe("cb")).join(nst.to_dataframe("nst"))
        return dl.fit_glm("death ~ " + " + ".join(list(cb.to_dataframe("cb").columns) + list(nst.to_dataframe("nst").columns)) + " + C(dow)", X)

    specs = [{"lag": 7, "df_var": 3}, {"lag": 14, "df_var": 4}, {"lag": 21, "df_var": 4}]
    grid = model_grid(specs, fit)
    assert list(grid.columns[:2]) == ["lag", "df_var"]
    assert grid["delta"].iloc[0] == 0 and grid["criterion"].is_monotonic_increasing


def _fit_spec(chicago, nst, spec, family="quasipoisson"):
    cb = dl.crossbasis(chicago.temp, lag=spec["lag"], argvar={"fun": "ns", "df": spec["df_var"]},
                       arglag={"fun": "ns", "df": spec["df_lag"]})
    X = dl.design_matrix(chicago, ("cb", cb), ("nst", nst), intercept=False)
    return dl.fit_glm("death ~ " + " + ".join(X.columns) + " + C(dow)", chicago.join(X), family=family)


def test_qaic_matches_r(chicago):
    """QAIC against R (tools/make_fixtures_qaic.R).

    R's logLik() is NA for a quasipoisson glm, so the reference value is
    -2*sum(dpois(y, fitted, log=TRUE)) + 2*phi*k. statsmodels' results.llf
    is that log-likelihood divided by phi, so using it would scale the fit
    term by 1/phi -- see the note in dlnmpy.uncertainty.qaic.
    """
    pytest.importorskip("statsmodels")
    with open(FIX / "qaic.json") as f:
        r = json.load(f)
    nst = dl.onebasis(chicago.time, "ns", df=7 * 14)

    for spec in r["quasipoisson"]:
        res = _fit_spec(chicago, nst, spec)
        assert res.params.size == spec["rank"], f"parameter count differs from R for {spec}"
        assert_close(res.scale, spec["dispersion"], atol=1e-9, msg=f"dispersion {spec}")
        assert_close(qaic(res), spec["qaic"], atol=1e-6, msg=f"qaic {spec}")

    # with a known scale QAIC is the AIC, which R reports directly
    p = r["poisson"]
    res = _fit_spec(chicago, nst, p, family="poisson")
    assert_close(qaic(res), p["aic"], atol=1e-6, msg="poisson qaic == R AIC")


def test_model_grid_ranking_matches_r(chicago):
    """The ranking, not just the values: a criterion that scales the
    log-likelihood by the dispersion reverses this grid, because phi is
    largest for the most under-fitted specification."""
    pytest.importorskip("statsmodels")
    with open(FIX / "qaic.json") as f:
        r = json.load(f)
    nst = dl.onebasis(chicago.time, "ns", df=7 * 14)
    specs = [{k: s[k] for k in ("lag", "df_var", "df_lag")} for s in r["quasipoisson"]]
    grid = model_grid(specs, lambda s: _fit_spec(chicago, nst, s))

    r_order = [specs[i] for i in r["order"]]
    py_order = grid[["lag", "df_var", "df_lag"]].to_dict("records")
    assert py_order == r_order, f"ranking differs from R\nR : {r_order}\npy: {py_order}"
    assert_close(grid["criterion"].to_numpy(),
                 [r["quasipoisson"][i]["qaic"] for i in r["order"]], atol=1e-6, msg="ranked QAIC")
