"""Regression tests for the defects fixed in 0.7.0 (third audit)."""

import warnings

import numpy as np
import pytest

import dlnmpy as dl
from dlnmpy.predict import CenteringWarning


@pytest.fixture(scope="module")
def chi(chicago):
    pytest.importorskip("statsmodels")
    d = chicago.iloc[:2000].copy()
    cb = dl.crossbasis(d.temp, lag=10, argvar={"fun": "ns", "df": 4}, arglag={"fun": "ns", "df": 3})
    ns_time = dl.onebasis(d.time, "ns", df=20)
    X = dl.design_matrix(d, ("cb", cb), ("ns_time", ns_time), intercept=False)
    f = "death ~ " + " + ".join(X.columns) + " + C(dow)"
    return d.join(X), f, cb


# --- fit_glm -----------------------------------------------------------------
def test_exposure_is_not_logged_twice_on_the_aliased_refit(chi):
    data, f, cb = chi
    pop = np.full(len(data), 2.7e6)
    ref = dl.fit_glm(f, data, family="quasipoisson", exposure=pop)
    data = data.assign(dup=data["cb_v1_l1"])   # an aliased column forces the refit path
    al = dl.fit_glm(f + " + dup", data, family="quasipoisson", exposure=pop)
    assert al.aliased == ["dup"]
    # the refit used to pass log(exposure) back as `exposure`, logging it
    # again and shifting every coefficient
    common = [c for c in ref.params.index]
    np.testing.assert_allclose(al.params[common].to_numpy(), ref.params.to_numpy(), rtol=1e-8)
    p_ref = dl.crosspred(cb, ref, cen=20, at=[0, 30])
    p_al = dl.crosspred(cb, al, cen=20, at=[0, 30])
    np.testing.assert_allclose(p_al.allRRfit, p_ref.allRRfit, rtol=1e-8)


# --- bootstrap_ci ------------------------------------------------------------
def test_bootstrap_ci_keeps_an_uncentred_prediction_uncentred(chi, cases):
    data, f, cb = chi
    m = cases["ex5"]["model"]
    cb21 = dl.crossbasis(data.temp, lag=21, argvar={"fun": "ns", "df": 4}, arglag={"fun": "ns", "df": 4})
    p = dl.crosspred(cb21, coef=m["coef"], vcov=m["vcov"], model_link=m["link"], at=[-20, 0, 20, 30], cen=False)
    assert p.cen is None
    b = dl.bootstrap_ci(p, nsim=200, seed=1, basis=cb21)
    # the draws used to be centred automatically while the fit stayed
    # uncentred, so the point estimate fell outside its own interval
    assert np.all(b["low"] <= b["fit"]) and np.all(b["fit"] <= b["high"])


# --- crossreduce -------------------------------------------------------------
def test_crossreduce_refuses_a_vector_value(chi, cases):
    data, f, cb = chi
    m = cases["ex5"]["model"]
    cb21 = dl.crossbasis(data.temp, lag=21, argvar={"fun": "ns", "df": 4}, arglag={"fun": "ns", "df": 4})
    with pytest.raises(ValueError, match="scalar"):
        dl.crossreduce(cb21, coef=m["coef"], vcov=m["vcov"], model_link=m["link"], type="lag", value=[0, 5], cen=20)
    with pytest.raises(ValueError, match="scalar"):
        dl.crossreduce(cb21, coef=m["coef"], vcov=m["vcov"], model_link=m["link"], type="var", value=[10, 20], cen=20)


# --- qaic ---------------------------------------------------------------------
def test_qaic_refuses_non_poisson_families(chi):
    data, f, cb = chi
    g = dl.fit_glm(f, data, family="gaussian")
    with pytest.raises(ValueError, match="Poisson"):
        dl.qaic(g)
    data = data.assign(high=(data.death > data.death.median()).astype(float))
    b = dl.fit_glm(f.replace("death ~", "high ~"), data, family="quasibinomial")
    with pytest.raises(ValueError, match="Poisson"):
        dl.qaic(b)
    assert np.isfinite(dl.qaic(dl.fit_glm(f, data, family="quasipoisson")))


# --- centring message ---------------------------------------------------------
def test_automatic_centring_warns_and_explicit_centring_does_not(chi, cases):
    data, f, cb = chi
    m = cases["ex5"]["model"]
    cb21 = dl.crossbasis(data.temp, lag=21, argvar={"fun": "ns", "df": 4}, arglag={"fun": "ns", "df": 4})
    with pytest.warns(CenteringWarning, match="Automatically set"):
        p = dl.crosspred(cb21, coef=m["coef"], vcov=m["vcov"], model_link=m["link"], at=[0, 30])
    assert p.cen is not None
    with warnings.catch_warnings():
        warnings.simplefilter("error", CenteringWarning)
        dl.crosspred(cb21, coef=m["coef"], vcov=m["vcov"], model_link=m["link"], at=[0, 30], cen=20)
        dl.crosspred(cb21, coef=m["coef"], vcov=m["vcov"], model_link=m["link"], at=[0, 30], cen=False)
        # a value stored in the basis is not "unspecified"
        cbc = dl.crossbasis(data.temp, lag=21, argvar={"fun": "ns", "df": 4, "cen": 20},
                            arglag={"fun": "ns", "df": 4})
        dl.crosspred(cbc, coef=m["coef"], vcov=m["vcov"], model_link=m["link"], at=[0, 30])


# --- basis-level validation ---------------------------------------------------
def test_cr_needs_three_knots_and_poly_a_positive_degree():
    x = np.linspace(0, 50, 60)
    with pytest.raises(ValueError, match="3"):
        dl.basis.cr(x, knots=[10, 30])
    with pytest.raises(ValueError, match="degree"):
        dl.crossbasis(x, lag=2, argvar={"fun": "poly", "degree": 0})


def test_group_with_missing_labels_is_refused():
    x = np.linspace(0, 50, 60)
    g = np.repeat([1.0, 2.0], 30)
    g[5] = np.nan
    with pytest.raises(ValueError, match="missing"):
        dl.crossbasis(x, lag=3, argvar={"fun": "lin"}, group=g)
    with pytest.raises(ValueError, match="length"):
        dl.crossbasis(x, lag=3, argvar={"fun": "lin"}, group=g[:-1])


def test_cbpen_refuses_a_zero_penalty():
    from dlnmpy.penalty import _rescale
    with pytest.raises(ValueError, match="positive eigenvalue"):
        _rescale(np.zeros((2, 2)))


# --- attribution --------------------------------------------------------------
def test_attrdl_lags_within_groups(cases):
    sim = dl.datasets.simulate_cities(3, 800, seed=2)
    sim["y"] = sim["y"].fillna(30.0)
    x, y, g = sim.tmean.to_numpy(), sim.y.to_numpy(), sim.city.to_numpy()
    cb = dl.crossbasis(x, lag=10, argvar={"fun": "ns", "df": 3}, arglag={"fun": "ns", "df": 3}, group=g)
    rng = np.random.default_rng(0)
    coef = rng.normal(0, 0.02, cb.ncol)
    vcov = np.eye(cb.ncol) * 1e-5
    with pytest.warns(UserWarning, match="group"):
        af_cross = dl.attrdl(x, cb, y, coef=coef, vcov=vcov, model_link="log", type="af", cen=20)
    af_group = dl.attrdl(x, cb, y, coef=coef, vcov=vcov, model_link="log", type="af", cen=20, group=g)
    # reference: mask the cases in the first `lag` rows of every group, the
    # rows where the cross-basis itself is NaN, and lag the whole series
    y_masked = y.astype(float).copy()
    y_masked[np.isnan(cb.matrix).any(axis=1)] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        af_ref = dl.attrdl(x, cb, y_masked, coef=coef, vcov=vcov, model_link="log", type="af", cen=20)
    assert af_group != af_cross
    assert np.isclose(af_group, af_ref)
    # the daily contributions are NaN exactly where the cross-basis is
    daily = dl.attrdl(x, cb, y, coef=coef, vcov=vcov, model_link="log", type="an", cen=20, group=g, tot=False)
    assert np.array_equal(np.isnan(daily), np.isnan(cb.matrix).any(axis=1))


def test_attr_table_components_never_straddle_the_reference(chi):
    data, f, cb = chi
    m = dl.fit_glm(f, data, family="quasipoisson")
    plo = np.asarray(dl._rcompat.quantile7(data.temp.to_numpy(), 0.025)).ravel()[0]
    cen = float(plo) - 3.0                                  # below the 2.5th percentile
    t = dl.attr_table(data.temp, cb, data.death, m, cen=cen, nsim=100, seed=0)
    rng = dict(zip(t.component, t.range))
    assert rng["extreme cold"] == (rng["total"][0], cen)   # cut at the reference
    assert rng["moderate cold"] == (cen, cen)              # empty
    assert rng["moderate heat"][0] == cen
    for lo, hi in t.range:
        assert lo <= hi
        assert not (lo < cen < hi) or (lo, hi) in (rng["total"],)


# --- fit_clogit ---------------------------------------------------------------
def test_fit_clogit_drops_rows_with_missing_values():
    pytest.importorskip("statsmodels")
    rng = np.random.default_rng(3)
    n = 400
    x = rng.normal(size=n)
    cb = dl.crossbasis(x, lag=2, argvar={"fun": "lin"}, arglag={"fun": "strata", "df": 1})
    X = cb.to_dataframe("cb")
    strata = np.repeat(np.arange(n // 4), 4)
    y = np.tile([1, 0, 0, 0], n // 4).astype(float)
    res = dl.fit_clogit(y, X, strata)               # first two rows are NaN by construction
    assert res.n_dropped == 2
    ref = dl.fit_clogit(y[2:], X.iloc[2:], strata[2:])
    np.testing.assert_allclose(np.asarray(res.params), np.asarray(ref.params), rtol=1e-10)


# --- meta / penalised ---------------------------------------------------------
def test_mixmeta_refuses_missing_outcomes():
    y = np.array([[0.1, 0.2], [0.3, np.nan], [0.2, 0.1]])
    S = np.array([np.eye(2) * 0.01] * 3)
    with pytest.raises(ValueError, match="missing"):
        dl.meta.mixmeta(y, S)


def test_fit_pglm_accepts_a_scalar_sp(chicago):
    pytest.importorskip("statsmodels")
    d = chicago.iloc[:800]
    ob = dl.onebasis(d.temp, "ps", df=8)
    X = np.column_stack([np.ones(len(d)), ob.matrix])
    S = np.zeros((X.shape[1],) * 2)
    S[1:, 1:] = ob.attrs["S"]
    a = dl.fit_pglm(d.death.to_numpy(), X, [S], family="quasipoisson", sp=1.0)
    b = dl.fit_pglm(d.death.to_numpy(), X, [S], family="quasipoisson", sp=[1.0])
    np.testing.assert_allclose(a.params, b.params)
    with pytest.raises(ValueError, match="one value per penalty"):
        dl.fit_pglm(d.death.to_numpy(), X, [S], family="quasipoisson", sp=[1.0, 2.0])


# --- workflow -----------------------------------------------------------------
def test_dlnm_reports_an_additive_effect_for_a_gaussian_model(chicago):
    pytest.importorskip("statsmodels")
    d = chicago.iloc[:1500]
    g = dl.dlnm(d, "death", "temp", lag=10, time="time", dow="dow", family="gaussian")
    t = g.rr_at([1, 99])
    assert "effect" in t.columns and "rr" not in t.columns
    assert np.all(np.abs(t.effect) < 100)      # deaths per day, not exp(deaths per day)
    assert "effect at 1st percentile" in g.summary()
    q = dl.dlnm(d, "death", "temp", lag=10, time="time", dow="dow")
    assert "rr" in q.rr_at([1]).columns and "RR at 1st percentile" in q.summary()


def test_dlnm_group_gets_intercepts_and_per_group_seasonality():
    pytest.importorskip("statsmodels")
    sim = dl.datasets.simulate_cities(3, 1500, seed=2)
    sim["y"] = sim["y"].fillna(30.0)
    fit = dl.dlnm(sim, "y", "tmean", lag=10, time="time", df_per_year=4, group="city")
    names = list(fit.model.params.index)
    assert any(n.startswith("C(city)") for n in names)
    assert any("ns_time_g0" in n for n in names) and any("ns_time_g2" in n for n in names)
    assert fit.group is not None
    # attribution lags within groups (no warning about crossing boundaries)
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        t = fit.attributable(nsim=50, seed=1)
    assert np.isfinite(t.an).all()


def test_dlnm_mmt_is_cached_and_figure_has_the_interval_band(chicago):
    pytest.importorskip("statsmodels")
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    d = chicago.iloc[:1500]
    fit = dl.dlnm(d, "death", "temp", lag=10, time="time", dow="dow")
    a = fit.mmt(nsim=200, seed=4)
    assert fit.mmt(nsim=200, seed=4) is a
    assert fit.mmt(nsim=200, seed=5) is not a
    fresh = dl.dlnm(d, "death", "temp", lag=10, time="time", dow="dow")
    ax = fresh.figure()          # first call used to draw without the MMT band
    # axvspan returns a Polygon before matplotlib 3.9 and a Rectangle after
    assert any(type(p).__name__ in ("Rectangle", "Polygon") for p in ax.patches)
    import matplotlib.pyplot as plt
    plt.close("all")
