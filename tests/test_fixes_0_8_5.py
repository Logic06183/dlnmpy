"""Regression tests for the defects found by the statistical review for the
manuscript (unreleased 0.8.5)."""

import numpy as np
import pytest

import dlnmpy as dl
from dlnmpy.model import get_link


def _large_value_prediction():
    rng = np.random.default_rng(0)
    x = rng.uniform(1e5, 1.0002e5, 300)          # e.g. air pressure in Pa
    cb = dl.crossbasis(x, lag=3, argvar={"fun": "lin"}, arglag={"fun": "poly", "degree": 1})
    return dl.crosspred(cb, coef=[1e-3, 2e-3], vcov=np.eye(2) * 1e-8,
                        at=[100000, 100001, 100002], model_link="log")


def test_prediction_lookup_picks_the_exact_row_on_a_fine_grid_of_large_values():
    # np.isclose's default rtol=1e-5 matched 100002 to the row for 100001
    cp = _large_value_prediction()
    assert [cp._var_index(v) for v in (100000, 100001, 100002)] == [0, 1, 2]
    np.testing.assert_allclose(cp.slice_var(100002).fit, np.exp(cp.matfit[2]), rtol=1e-12)
    with pytest.raises(ValueError, match="must match values used for prediction"):
        cp._var_index(100001.5)


def test_prediction_lookup_still_tolerates_floating_point_noise():
    x = np.linspace(-5, 30, 400)
    cb = dl.crossbasis(x, lag=2, argvar={"fun": "lin"}, arglag={"fun": "integer"})
    cp = dl.crosspred(cb, coef=[0.01, 0.02, 0.03], vcov=np.eye(3) * 1e-6, from_=0, to=1,
                      by=0.1, model_link="log", cen=0)
    assert cp._var_index(0.1 * 3) == int(np.argmin(np.abs(cp.predvar - 0.3)))


def test_logc_link_is_not_reported_as_a_log_link():
    sm = pytest.importorskip("statsmodels.api")
    y = np.r_[np.zeros(5), np.ones(5)]
    X = np.column_stack([np.ones(10), np.arange(10.0)])
    logc = sm.GLM(y, X, family=sm.families.Binomial(link=sm.families.links.LogC()))
    assert get_link(logc) == "logc"
    for link, name in ((sm.families.links.Log(), "log"), (sm.families.links.CLogLog(), "cloglog"),
                       (sm.families.links.Logit(), "logit")):
        m = sm.GLM(y, X, family=sm.families.Binomial(link=link))
        assert get_link(m) == name


def _chicago_fit(family, lag=3, rows=slice(None)):
    pytest.importorskip("statsmodels")
    d = dl.datasets.chicago_nmmaps().iloc[rows].reset_index(drop=True)
    cb = dl.crossbasis(d.temp, lag=lag, argvar={"fun": "ns", "df": 3}, arglag={"fun": "integer"})
    X = dl.design_matrix(d, ("cb", cb), ("ns_time", dl.onebasis(d.time, "ns", df=4)), intercept=False)
    return d, cb, dl.fit_glm("death ~ " + " + ".join(X.columns), d.join(X), family=family)


def test_attr_table_refuses_a_model_without_a_log_or_logit_link():
    # attrdl refused a Gaussian model; attr_table returned an AF of -23826 for one
    d, cb, m = _chicago_fit("gaussian", rows=slice(0, 1500))
    with pytest.raises(ValueError, match="log or logit link"):
        dl.attr_table(d.temp, cb, d.death, m, cen=20, nsim=10, seed=1)


def test_model_grid_warns_when_the_fits_use_different_rows():
    # longer maximum lags drop more initial days, so the criteria are not comparable
    def fit(spec):
        return _chicago_fit("quasipoisson", lag=spec["lag"], rows=slice(0, 1500))[2]
    with pytest.warns(UserWarning, match="different numbers of observations"):
        grid = dl.model_grid([{"lag": 3}, {"lag": 10}], fit)
    assert len(grid) == 2


def test_fit_glm_dispersion_does_not_count_rows_with_zero_prior_weight():
    # R's glm.fit uses n.ok = n - sum(weights == 0) for the residual df
    pytest.importorskip("statsmodels")
    import pandas as pd
    rng = np.random.default_rng(2)
    n = 400
    d = pd.DataFrame({"x": rng.normal(size=n)})
    d["y"] = rng.poisson(np.exp(1 + 0.3 * d.x)).astype(float)
    w = np.ones(n)
    w[:100] = 0
    a = dl.fit_glm("y ~ x", d, family="quasipoisson", var_weights=w)
    b = dl.fit_glm("y ~ x", d.iloc[100:], family="quasipoisson")
    np.testing.assert_allclose(a.scale, b.scale, rtol=1e-8)
    np.testing.assert_allclose(np.asarray(a.bse), np.asarray(b.bse), rtol=1e-6)
    assert a.df_resid == b.df_resid == 298


def test_fit_glm_aliasing_uses_glm_fit_tolerance():
    # glm.fit tests aliasing at min(1e-7, epsilon / 1000) = 1e-11, not lm.fit's 1e-7
    pytest.importorskip("statsmodels")
    import pandas as pd
    rng = np.random.default_rng(5)
    n = 300
    d = pd.DataFrame({"x1": rng.normal(size=n), "x2": rng.normal(size=n)})
    d["x3"] = d.x1 + d.x2 + 1e-9 * rng.normal(size=n)       # nearly, not exactly, aliased
    d["x4"] = d.x1 + d.x2                                    # exactly aliased
    d["y"] = rng.poisson(np.exp(0.5 + 0.1 * d.x1)).astype(float)
    assert dl.fit_glm("y ~ x1 + x2 + x3", d, family="poisson").aliased == []
    assert dl.fit_glm("y ~ x1 + x2 + x4", d, family="poisson").aliased == ["x4"]


def test_fit_clogit_reports_non_convergence_on_separated_data():
    # coxph warns that a coefficient may be infinite; the step-halving path
    # used to stop at a non-stationary point and say nothing
    pytest.importorskip("statsmodels")
    import pandas as pd
    rng = np.random.default_rng(3)
    nset, nper = 50, 4
    grp = np.repeat(np.arange(nset), nper)
    x = rng.gamma(2, 2500, size=nset * nper)
    y = np.zeros(nset * nper)
    for s in range(nset):
        ix = np.arange(s * nper, (s + 1) * nper)
        y[ix[np.argmax(x[ix])]] = 1                          # the case always has the largest x
    X = pd.DataFrame({"x": x, "x2": rng.normal(size=x.size)})
    with pytest.warns(Warning, match="did not converge"):
        res = dl.fit_clogit(y, X, grp)
    assert res.converged is False
    rng = np.random.default_rng(1)
    g = np.repeat(np.arange(100), 4)
    yy = np.zeros(400)
    yy[::4] = 1
    ok = dl.fit_clogit(yy, pd.DataFrame({"x": rng.normal(size=400)}), g)
    assert ok.converged is True
