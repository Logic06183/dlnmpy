"""Multivariate meta-analysis versus the R package mixmeta (12 simulated locations)."""

import json

import numpy as np
import pandas as pd
import pytest

import dlnmpy as dl
from dlnmpy.meta import mixmeta, predict_reduced, stack_reduced, vech, xpnd

from conftest import FIX, assert_close


@pytest.fixture(scope="module")
def R():
    with open(FIX / "meta.json") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def inputs(R):
    y = np.array(R["y"])
    S = np.array(R["S_vech"])
    X = np.column_stack([np.ones(len(y)), R["meanT"]])
    return y, S, X


def _rows(d, key):
    return np.array([d[str(i + 1)][key] for i in range(len(d))])


def test_vech_roundtrip():
    M = np.array([[1.0, 2, 3], [2, 4, 5], [3, 5, 6]])
    assert_close(vech(M), [1, 2, 3, 4, 5, 6])
    assert_close(xpnd(vech(M)), M)


@pytest.mark.parametrize("key", ["reml_mm", "reml_mr", "ml_mm", "ml_mr", "fixed_mm", "fixed_mr"])
def test_mixmeta_matches_r(R, inputs, key):
    y, S, X = inputs
    meth, kind = key.split("_")
    f = R["fits"][key]
    mm = mixmeta(y, S, X=None if kind == "mm" else X, method=meth)
    tol = 1e-10 if meth == "fixed" else 1e-5
    assert_close(mm.coef_vec, f["coef"], atol=tol)
    assert_close(mm.vcov, f["vcov"], atol=tol)
    assert_close(mm.Psi, f["Psi"], atol=tol)
    # the optimum must be at least as good as R's (R stops at its own tolerance)
    assert mm.loglik >= f["logLik"] - 1e-9
    assert abs(mm.aic - f["AIC"]) < 1e-6 and abs(mm.bic - f["BIC"]) < 1e-6
    assert mm.df_residual == f["df_res"]
    q = mm.qtest()
    assert_close(q["Q"], f["Q"], atol=1e-9)
    assert_close(q["df"], f["Qdf"])
    assert_close(q["pvalue"], f["Qp"], atol=1e-9)
    assert_close(q["I2"], f["I2"], atol=1e-9)
    b = mm.blup(se=True)
    assert_close(b["blup"], _rows(f["blup"], "blup"), atol=1e-4)
    assert_close(b["se"], _rows(f["blup"], "se"), atol=1e-5)
    assert_close(b["pi_low"], _rows(f["blup"], "pi_lb"), atol=1e-4)
    assert_close(np.array(b["vcov"]), _rows(f["blup"], "vcov"), atol=1e-5)
    if kind == "mr":
        p = mm.predict(np.column_stack([np.ones(3), [10, 15, 20]]), se=True)
        assert_close(p["fit"], _rows(f["pred"], "fit"), atol=1e-6)
        assert_close(p["se"], _rows(f["pred"], "se"), atol=1e-6)
    assert "Cochran Q" in mm.summary()


def test_other_structures_and_univariate(R, inputs):
    y, S, X = inputs
    for key in ["reml_diag", "reml_id"]:
        f = R["fits"][key]
        mm = mixmeta(y, S, method="reml", bscov=key.split("_")[1])
        assert mm.loglik >= f["logLik"] - 1e-9  # R reported non-convergence for these
        assert_close(np.diag(mm.Psi) > 0, np.diag(np.array(f["Psi"])) > 0)
    f = R["fits"]["uni_reml"]
    mm = mixmeta(y[:, 0], S[:, 0], method="reml")
    assert_close(mm.coef_vec, np.atleast_1d(f["coef"]), atol=1e-8)
    assert_close(mm.Psi, np.atleast_2d(f["Psi"]), atol=1e-8)
    assert_close(mm.qtest()["I2"], np.atleast_1d(f["I2"]), atol=1e-9)


def test_two_stage_pipeline(R):
    """Stage 1 in Python (crossreduce per city), stage 2 pooled curve versus R."""
    sim = pd.read_csv(FIX / "meta_sim.csv")
    reds = []
    for k, d in sim.groupby("city"):
        d = d.reset_index(drop=True)
        cb = dl.crossbasis(d.tmean, lag=10, argvar={"fun": "ns", "knots": R["knots"], "boundary_knots": R["bk"]},
                           arglag={"fun": "ns", "df": 3})
        nt = dl.onebasis(d.time, "ns", df=18)
        Xd = d.join(cb.to_dataframe("cb")).join(nt.to_dataframe("nt"))
        fit = dl.fit_glm("y ~ " + " + ".join(list(cb.to_dataframe("cb").columns) + list(nt.to_dataframe("nt").columns)), Xd)
        reds.append(dl.crossreduce(cb, fit, cen=18, name="cb"))
        st = R["stage1"][k - 1]
        assert_close(reds[-1].coef, st["coef"], atol=1e-8)
        assert_close(reds[-1].vcov, st["vcov"], atol=1e-8)
    y, S = stack_reduced(reds)
    assert_close(y, R["y"], atol=1e-8)
    mm = mixmeta(y, S, method="reml")
    pooled = predict_reduced(cb, mm.coef_vec, mm.vcov, at=np.arange(-5, 41), cen=18)
    pc = R["pooled_curve"]
    assert_close(pooled.predvar, pc["at"])
    np.testing.assert_allclose(pooled.allRRfit, pc["allRRfit"], rtol=1e-5)
    np.testing.assert_allclose(pooled.allRRlow, pc["allRRlow"], rtol=1e-4)
    # BLUP curve for one city runs through the same helper
    b = mm.blup(se=True)
    city_curve = predict_reduced(cb, b["blup"][0], b["vcov"][0], at=np.arange(-5, 41), cen=18)
    assert city_curve.allRRfit.shape == (46,)
