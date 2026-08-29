"""End-to-end: fit with statsmodels and compare with R's glm + crosspred."""

import numpy as np
import pytest

import dlnmpy as dl

from conftest import assert_close

pytest.importorskip("statsmodels")


@pytest.fixture(scope="module")
def ex1_fit(chicago):
    cbpm = dl.crossbasis(chicago.pm10, lag=15, argvar={"fun": "lin"}, arglag={"fun": "poly", "degree": 4})
    cbtemp = dl.crossbasis(chicago.temp, lag=3, argvar={"df": 5}, arglag={"fun": "strata", "breaks": 1})
    nstime = dl.onebasis(chicago.time, "ns", df=7 * 14)
    X = dl.design_matrix(chicago, ("cbpm", cbpm), ("cbtemp", cbtemp), ("nstime", nstime), intercept=False)
    data = chicago.join(X)
    rhs = " + ".join(list(X.columns) + ["C(dow)"])
    fit = dl.fit_glm("death ~ " + rhs, data, family="quasipoisson")
    return cbpm, cbtemp, fit


def test_fit_matches_r_glm(ex1_fit, cases):
    cbpm, cbtemp, fit = ex1_fit
    m = cases["ex1"]["model_pm"]
    assert fit.aliased == ["nstime_b13", "nstime_b14"]  # R reports these as NA
    assert np.isclose(fit.scale, m["dispersion"], rtol=1e-6)
    assert_close(fit.params.filter(like="cbpm").to_numpy(), m["coef"], atol=1e-10)
    assert dl.get_link(fit) == "log"


def test_crosspred_from_model(ex1_fit, cases):
    cbpm, cbtemp, fit = ex1_fit
    p = dl.crosspred(cbpm, fit, at=np.arange(21), bylag=0.2, cumul=True, name="cbpm")
    r = cases["ex1"]["pred1_pm"]
    assert_close(p.allRRfit, r["allRRfit"], atol=1e-10)
    assert_close(p.allRRlow, r["allRRlow"], atol=1e-8)
    assert_close(p.allRRhigh, r["allRRhigh"], atol=1e-8)
    p2 = dl.crosspred(cbtemp, fit, by=1, cen=20, name="cbtemp")
    assert_close(p2.matfit, cases["ex1"]["pred1_temp"]["matfit"], atol=1e-10)
    assert_close(p2.matse, cases["ex1"]["pred1_temp"]["matse"], atol=1e-8)
    # wrong name -> informative error
    with pytest.raises(ValueError, match="coefficients matching"):
        dl.crosspred(cbpm, fit, name="nothere")


def test_extract_by_r_style_names():
    import pandas as pd

    class Fake:
        params = pd.Series([1.0, 2.0, 3.0, 4.0], index=["Intercept", "cb.v1.l1", "cb.v1.l2", "x"])

        def cov_params(self):
            return pd.DataFrame(np.eye(4), index=self.params.index, columns=self.params.index)

    coef, vcov = dl.extract_coef_vcov(Fake(), "cb", "cb", 2)
    assert_close(coef, [2.0, 3.0])
    assert vcov.shape == (2, 2)
