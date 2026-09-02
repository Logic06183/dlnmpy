"""Regression tests for the defects fixed in 0.6.0."""

import numpy as np
import pytest

import dlnmpy as dl


@pytest.fixture(scope="module")
def two_penalised_bases(chicago):
    pytest.importorskip("statsmodels")
    d = chicago.iloc[:1500]
    cb = dl.crossbasis(d.temp, lag=10, argvar={"fun": "ps", "df": 6}, arglag={"fun": "ps", "df": 5})
    cb_o3 = dl.crossbasis(d.o3, lag=2, argvar={"fun": "ps", "df": 5}, arglag={"fun": "strata", "df": 1})
    ns_time = dl.onebasis(d.time, "ns", df=12)
    X = dl.design_matrix(d, ("cb", cb), ("cb_o3", cb_o3), ("ns_time", ns_time), intercept=False)
    f = "death ~ " + " + ".join(X.columns) + " + C(dow)"
    return d.join(X), f, cb, cb_o3, X


def test_penalty_prefix_does_not_capture_a_longer_basis_name(two_penalised_bases):
    data, f, cb, cb_o3, X = two_penalised_bases
    # 'cb' used to match the 'cb_o3_*' columns too, so the penalty could not
    # be placed (a shape error) and edf_by('cb') silently summed both bases
    res = dl.fit_pgam(f, data, {"cb": dl.cbpen(cb), "cb_o3": dl.cbpen(cb_o3)}, family="quasipoisson")
    n_cb = sum(c.startswith("cb_v") for c in X.columns)
    n_o3 = sum(c.startswith("cb_o3_") for c in X.columns)
    assert n_cb == cb.ncol and n_o3 == cb_o3.ncol
    assert res.edf_by("cb") <= n_cb + 1e-9
    assert res.edf_by("cb_o3") <= n_o3 + 1e-9
    assert np.isclose(res.edf_by("cb") + res.edf_by("cb_o3"),
                      np.sum(res.edf[[i for i, c in enumerate(res.params.index) if c.startswith("cb")]]))
    # and the prediction machinery agrees on which columns are which
    pred = dl.crosspred(cb, res, cen=21, at=[-10, 30], name="cb")
    assert pred.coef.size == cb.ncol


def test_qaic_accepts_penalised_fits(two_penalised_bases):
    data, f, cb, cb_o3, X = two_penalised_bases
    res = dl.fit_pgam(f, data, {"cb": dl.cbpen(cb), "cb_o3": dl.cbpen(cb_o3)}, family="quasipoisson")
    q = dl.qaic(res)
    assert np.isfinite(q)
    # the parameter count is the effective df, so QAIC is below the unpenalised count version
    from scipy.special import gammaln
    y, mu = res.endog, res.fitted_values
    llf = float(np.sum(y * np.log(mu) - mu - gammaln(y + 1)))
    assert np.isclose(q, -2 * llf + 2 * res.scale * res.edf_total)
    assert q < -2 * llf + 2 * res.scale * res.params.size


def test_overlay_labels_must_match_curves(chicago, cases):
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    cb = dl.crossbasis(chicago.temp, lag=21, argvar={"fun": "ns", "df": 4}, arglag={"fun": "ns", "df": 4})
    m = cases["ex5"]["model"]
    p = dl.crosspred(cb, coef=m["coef"], vcov=m["vcov"], model_link=m["link"], cen=21, by=1)
    with pytest.raises(ValueError, match="labels"):
        dl.plot.overlay_slices(p, var=[-20, 0, 33], labels=["cold"])


def test_qtest_pvalue_is_not_rounded_to_zero():
    rng = np.random.default_rng(0)
    m, k = 30, 1
    y = rng.normal(0, 1, (m, k))  # large heterogeneity against small S: Q around 600
    S = np.full((m, k, k), 0.05)
    fit = dl.mixmeta(y, S, method="fixed")
    q = fit.qtest()
    assert q["pvalue"][0] > 0  # 1 - cdf would return exactly 0
    assert q["pvalue"][0] < 1e-30


def test_ps_too_small_for_prediction_fails_at_construction():
    # R builds ps(df=4, degree=3) and only fails inside crosspred(); here the
    # message comes from the basis function, where 'df' is being chosen
    x = np.linspace(0, 10, 50)
    with pytest.raises(ValueError, match="at least 5"):
        dl.ps(x, df=4)
    with pytest.raises(ValueError, match="at least 5"):
        dl.ps(x, df=4, intercept=True)
    basis, attrs = dl.ps(x, df=5, intercept=True)  # smallest that survives the round trip
    again, _ = dl.ps(x[:10], knots=attrs["knots"], degree=attrs["degree"], intercept=True)
    assert again.shape == (10, basis.shape[1])
