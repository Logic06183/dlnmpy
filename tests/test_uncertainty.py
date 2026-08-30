import numpy as np
import pytest

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
    res = grid["results"].iloc[0]
    assert np.isclose(qaic(res), -2 * res.llf + 2 * res.scale * res.params.size)
