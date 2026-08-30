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


def test_basis_name_is_not_a_prefix_match(chicago):
    """A basis called 'cb' must not also match the columns of 'cb2'."""
    pytest.importorskip("statsmodels")
    cb = dl.crossbasis(chicago.temp, lag=5, argvar={"fun": "lin"}, arglag={"fun": "strata", "breaks": 2})
    cb2 = dl.crossbasis(chicago.pm10, lag=5, argvar={"fun": "lin"}, arglag={"fun": "strata", "breaks": 2})
    X = dl.design_matrix(chicago, ("cb", cb), ("cb2", cb2), intercept=False)
    fit = dl.fit_glm("death ~ " + " + ".join(X.columns) + " + C(dow)", chicago.join(X),
                     family="quasipoisson")
    for name, basis in (("cb", cb), ("cb2", cb2)):
        coef, vcov = dl.extract_coef_vcov(fit, name, "cb", basis.ncol)
        expected = fit.params[[c for c in fit.params.index if c.startswith(name + "_")]].to_numpy()
        assert_close(coef, expected, atol=0)
        assert vcov.shape == (basis.ncol, basis.ncol)
    # and prediction picks the right block
    assert dl.crosspred(cb2, fit, name="cb2", by=20).allfit.size > 0


def test_to_dataframe_keeps_the_input_index(chicago):
    """data.join(cb.to_dataframe(...)) must align when data has been filtered
    or grouped and so no longer has a default RangeIndex."""
    import pandas as pd

    sub = chicago[chicago.year >= 1996]          # index starts at 1826
    assert sub.index[0] != 0
    cb = dl.crossbasis(sub.temp, lag=10, argvar={"fun": "ns", "df": 4}, arglag={"fun": "ns", "df": 3})
    ob = dl.onebasis(sub.temp, "ns", df=4)
    assert list(cb.to_dataframe("cb").index) == list(sub.index)
    assert list(ob.to_dataframe("ob").index) == list(sub.index)

    joined = sub.join(cb.to_dataframe("cb"))
    # only the leading `lag` rows are missing, not every row
    assert int(joined.filter(like="cb_").isna().all(axis=1).sum()) == 10
    assert_close(joined.filter(like="cb_").to_numpy()[10:], cb.matrix[10:], atol=0)

    # a numpy input has no index, and an explicit index still wins
    plain = dl.crossbasis(sub.temp.to_numpy(), lag=10, argvar={"fun": "ns", "df": 4},
                          arglag={"fun": "ns", "df": 3})
    assert isinstance(plain.to_dataframe("cb").index, pd.RangeIndex)
    assert list(cb.to_dataframe("cb", index=sub.index[::-1]).index) == list(sub.index[::-1])


def test_fit_glm_offset_and_weights_reach_the_model():
    """offset/exposure/freq_weights/var_weights are model arguments; GLM.fit()
    accepts arbitrary keywords and would silently discard them."""
    pytest.importorskip("statsmodels")
    import pandas as pd
    import statsmodels.api as sm
    import statsmodels.formula.api as smf

    rng = np.random.default_rng(1)
    n = 200
    lpop = rng.normal(size=n)
    d = pd.DataFrame({"x": rng.normal(size=n), "lpop": lpop})
    d["y"] = rng.poisson(np.exp(0.5 + 0.3 * d.x + lpop))

    plain = dl.fit_glm("y ~ x", d, family="poisson")
    off = dl.fit_glm("y ~ x", d, family="poisson", offset=lpop)
    ref = smf.glm("y ~ x", d, family=sm.families.Poisson(), offset=lpop).fit()
    assert not np.allclose(off.params, plain.params), "offset was ignored"
    assert_close(off.params.to_numpy(), ref.params.to_numpy(), atol=1e-9)
    # exposure is offset(log(exposure))
    exp_ = dl.fit_glm("y ~ x", d, family="poisson", exposure=np.exp(lpop))
    assert_close(exp_.params.to_numpy(), off.params.to_numpy(), atol=1e-9)
    # weights reach the model too
    w = dl.fit_glm("y ~ x", d, family="poisson", freq_weights=np.full(n, 2.0))
    assert w.df_resid != plain.df_resid
    # and survive both the NaN drop and the aliased-column refit
    d2 = d.copy()
    d2.loc[5, "x"] = np.nan
    d2["x2"] = d2.x * 2
    m2 = dl.fit_glm("y ~ x + x2", d2, family="poisson", offset=lpop)
    ref2 = smf.glm("y ~ x", d2, family=sm.families.Poisson(), offset=lpop, missing="drop").fit()
    assert m2.nobs == ref2.nobs == n - 1
    assert_close(m2.params[["Intercept", "x"]].to_numpy(), ref2.params.to_numpy(), atol=1e-9)


def test_get_link_compound_names():
    """'log' is a substring of 'cloglog' and 'loglog'; matching it first
    reports those as a log link and makes crosspred exponentiate them."""
    sm = pytest.importorskip("statsmodels.api")

    class Fake:
        def __init__(self, fam):
            self.family = fam

    L = sm.families.links
    assert dl.get_link(Fake(sm.families.Binomial(link=L.CLogLog()))) == "cloglog"
    assert dl.get_link(Fake(sm.families.Binomial(link=L.LogLog()))) == "loglog"
    assert dl.get_link(Fake(sm.families.Poisson(link=L.Log()))) == "log"
    assert dl.get_link(Fake(sm.families.Binomial(link=L.Logit()))) == "logit"
    assert dl.get_link(Fake(sm.families.Binomial(link=L.Probit()))) == "probit"
