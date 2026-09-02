"""The one-call workflow must give exactly what the explicit pipeline gives."""

import numpy as np
import pytest

import dlnmpy as dl

pytest.importorskip("statsmodels")


@pytest.fixture(scope="module")
def fits(chicago):
    d = chicago
    argvar = {"fun": "bs", "degree": 2, "knots": dl.percentile_knots(d.temp, [10, 75, 90])}
    arglag = {"fun": "ns", "knots": dl.logknots(21, 3)}
    fit = dl.dlnm(d, outcome="death", exposure="temp", lag=21, argvar=argvar, arglag=arglag,
                  time="time", df_per_year=7, dow="dow")
    # the same by hand
    cb = dl.crossbasis(d.temp, lag=21, argvar=dict(argvar), arglag=dict(arglag))
    ns_time = dl.onebasis(d.time, "ns", df=7 * 14)
    X = dl.design_matrix(d, ("cb", cb), ("ns_time", ns_time), intercept=False)
    model = dl.fit_glm("death ~ " + " + ".join(X.columns) + " + C(dow)", d.join(X), family="quasipoisson")
    return fit, cb, model


def test_dlnm_matches_the_explicit_pipeline(fits, chicago):
    fit, cb, model = fits
    np.testing.assert_allclose(fit.cb.matrix[30:40], cb.matrix[30:40], atol=1e-12)
    np.testing.assert_allclose(fit.coef, dl.crosspred(cb, model, cen=21, name="cb").coef, atol=1e-12)
    np.testing.assert_allclose(fit.vcov, dl.crosspred(cb, model, cen=21, name="cb").vcov, atol=1e-12)
    p1 = fit.predict(cen=21, by=1)
    p2 = dl.crosspred(cb, model, cen=21, by=1, name="cb")
    np.testing.assert_allclose(p1.allRRfit, p2.allRRfit, atol=1e-12)
    np.testing.assert_allclose(p1.matse, p2.matse, atol=1e-12)


def test_mmt_is_the_default_reference(fits, chicago):
    fit, cb, model = fits
    m = fit.mmt(nsim=300, seed=1)
    ref = dl.mmt(cb, model, x=chicago.temp, nsim=300, seed=1, name="cb")
    assert m.mmt == ref.mmt and m.low == ref.low and m.high == ref.high
    p = fit.predict(by=1)
    assert p.cen == m.mmt
    rr = fit.rr_at([1, 99])
    assert list(rr.columns) == ["percentile", "temp", "rr", "low", "high"]
    np.testing.assert_allclose(rr.temp.to_numpy(), dl.quantile7(chicago.temp, [0.01, 0.99]))
    q = dl.crosspred(cb, model, at=rr.temp.to_numpy(), cen=m.mmt, name="cb")
    np.testing.assert_allclose(rr.rr.to_numpy(), q.allRRfit, atol=1e-12)
    # the RR at the reference itself is 1
    assert np.isclose(fit.predict(at=[m.mmt]).allRRfit[0], 1.0)


def test_attributable_and_percentiles(fits, chicago):
    fit, cb, model = fits
    fit.mmt(nsim=200, seed=1)
    tab = fit.attributable(nsim=200, seed=2)
    ref = dl.attr_table(chicago.temp, cb, chicago.death, model, cen=fit._mmt.mmt, nsim=200, seed=2, name="cb")
    np.testing.assert_allclose(tab.af.to_numpy(), ref.af.to_numpy(), atol=1e-12)
    np.testing.assert_allclose(tab.an_low.to_numpy(), ref.an_low.to_numpy(), atol=1e-8)
    assert list(tab.component) == ["total", "cold", "heat", "extreme cold", "moderate cold", "moderate heat", "extreme heat"]
    pct = fit.percentile(fit.quantile([1, 50, 99]))
    assert np.all(np.abs(pct - [1, 50, 99]) < 1.5)  # ties in the temperature series
    assert np.isclose(dl.percentile_of([1, 2, 3, 4], 2)[0], 50)


def test_to_dataframe_long_format(fits):
    fit, cb, model = fits
    p = fit.predict(cen=21, at=[-20, 0, 30], bylag=0.5, cumul=True)
    df = p.to_dataframe()
    assert df.shape == (3 * len(p.predlag), 6)
    row = df[(df["var"] == 30) & (df["lag"] == 2.5)].iloc[0]
    assert np.isclose(row.fit, p.slice_var(30).set_index("lag").loc[2.5, "fit"])
    dfc = p.to_dataframe(cumul=True)
    assert dfc.shape == (3 * 22, 6)
    assert np.isclose(dfc[(dfc["var"] == -20) & (dfc["lag"] == 21)].fit.iloc[0], p.allRRfit[0])
    raw = p.to_dataframe(exp=False)
    assert np.isclose(raw.fit.iloc[0], p.matfit[0, 0])


def test_defaults_penalised_and_offset(chicago):
    d = chicago.iloc[:1200].copy()
    d["pop"] = 1e6 + 10 * np.arange(len(d))
    d["logpop"] = np.log(d["pop"])
    f = dl.dlnm(d, "death", "temp", lag=3, time="time", offset="logpop", controls=["pm10"])
    assert "pm10" in f.formula and f.cb.arglag["fun"] == "ns"
    assert f.model.model.offset is not None
    # ps bases go through the penalised fitter by default
    g = dl.dlnm(d, "death", "temp", lag=5, argvar={"fun": "ps", "df": 6}, arglag={"fun": "ps", "df": 5}, time="time")
    assert type(g.model).__name__ == "PenalizedGLMResults"
    assert np.isfinite(g.qaic())
    assert "lag 0-5" in g.summary(nsim=50)


def test_new_figures(fits, chicago):
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    fit, cb, model = fits
    fit.mmt(nsim=100, seed=1)
    ax = fit.figure(xlab="Temperature")
    assert hasattr(ax, "hist_ax")
    assert ax.get_xlabel() == "" and ax.hist_ax.get_xlabel() == "Temperature"
    fig, a = plt.subplots()
    ax2 = dl.plot.plot_overall_risk(fit.predict(by=1), x=chicago.temp, mmt=fit._mmt, ax=a)
    assert ax2 is a
    tab = fit.attributable(nsim=100, seed=1)
    ax3 = dl.plot.plot_attributable(tab, components=["cold", "heat"])
    assert [t.get_text() for t in ax3.get_yticklabels()] == ["Cold", "Heat"]
    ax4 = dl.plot.plot_attributable(tab, kind="an")
    assert "number" in ax4.get_xlabel()
    fit.plot("3d", by=2)
    with dl.plot.style():
        assert mpl.rcParams["axes.spines.top"] is False
    plt.close("all")


# --- time= given as dates rather than a day index ---------------------------
# A date column is the obvious thing to pass and used to be silently fatal:
# datetime64 cast to float is nanoseconds, so fourteen years became 1.2e12
# "years" and the seasonal spline asked for 8e12 df, killing the process.


def test_datetime_time_column_matches_day_index(chicago):
    """time='date' must give exactly what the equivalent day index gives."""
    by_index = dl.dlnm(chicago, outcome="death", exposure="temp", lag=10,
                       time="time", df_per_year=7, dow="dow")
    by_date = dl.dlnm(chicago, outcome="death", exposure="temp", lag=10,
                      time="date", df_per_year=7, dow="dow")

    at = np.arange(-20, 31, 2.0)
    a = by_index.predict(at=at, cen=20.0)
    b = by_date.predict(at=at, cen=20.0)
    assert np.allclose(a.allRRfit, b.allRRfit, rtol=1e-10, atol=1e-12)
    assert np.allclose(a.allse, b.allse, rtol=1e-10, atol=1e-12)


def test_time_index_converts_datetimes_to_days(chicago):
    from dlnmpy.workflow import _time_index, _years

    days = _time_index(chicago["date"])
    assert days[0] == 0.0
    assert np.isclose(_years(days), _years(_time_index(chicago["time"])), atol=1e-9)
    assert np.isclose(_years(days), 14.0, atol=0.01)


def test_object_dtype_dates_are_understood(chicago):
    from dlnmpy.workflow import _time_index

    as_objects = chicago["date"].dt.date.astype(object)
    assert np.allclose(_time_index(as_objects), _time_index(chicago["date"]))


def test_absurd_df_per_year_raises_rather_than_exhausting_memory(chicago):
    """The guard that catches any unit error, not just this one."""
    with pytest.raises(ValueError, match="degrees of freedom"):
        dl.dlnm(chicago.head(200), outcome="death", exposure="temp", lag=5,
                time="time", df_per_year=400)
