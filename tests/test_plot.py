import pytest

import dlnmpy as dl

mpl = pytest.importorskip("matplotlib")
mpl.use("Agg")


@pytest.fixture(scope="module")
def pred(chicago, cases):
    cb = dl.crossbasis(chicago.temp, lag=21, argvar={"fun": "ns", "df": 4}, arglag={"fun": "ns", "df": 4})
    m = cases["ex5"]["model"]
    p = dl.crosspred(cb, coef=m["coef"], vcov=m["vcov"], model_link=m["link"], cen=21, by=1, cumul=True)
    rd = dl.crossreduce(cb, coef=m["coef"], vcov=m["vcov"], model_link=m["link"], cen=21, by=1)
    return p, rd


def test_plot_types(pred):
    import matplotlib.pyplot as plt
    p, rd = pred
    assert p.plot("overall", xlab="Temperature", ylab="RR") is not None
    plt.close("all")
    assert p.plot("slices", var=33) is not None
    plt.close("all")
    axes = p.plot("slices", var=[-20, 33], lag=[0, 5])
    assert len(axes) == 4
    plt.close("all")
    assert p.plot("slices", var=33, cumul=True, ci="bars") is not None
    plt.close("all")
    ax = p.plot("contour", xlab="Temperature", ylab="Lag (days)", zlab="RR")
    assert ax.get_ylabel() == "Lag (days)"
    plt.close("all")
    ax = p.plot("3d", zlab="RR")
    assert ax.get_zlabel() == "RR" and ax.get_ylabel() == "Lag"
    plt.close("all")
    assert rd.plot(ci="lines") is not None
    plt.close("all")
    with pytest.raises(ValueError):
        p.plot("slices")


def test_overlay_summary_and_themes(pred):
    import matplotlib.pyplot as plt
    from dlnmpy.plot import overlay_slices, set_theme, summary_figure
    p, rd = pred
    for theme in ["journal", "colour"]:
        set_theme(theme)
        ax = overlay_slices(p, var=[-20, 0, 33], ci="area")
        assert ax.get_legend() is not None and len(ax.get_lines()) >= 3
        plt.close("all")
        ax = overlay_slices(p, lag=[0, 5], legend=False)
        assert ax.get_legend() is None
        plt.close("all")
        fig, axes = summary_figure(p, var=[-20, 33], lag=[0, 5], xlab="Temperature")
        assert axes.shape == (2, 2)
        plt.close(fig)
    with pytest.raises(ValueError):
        set_theme("neon")
    with pytest.raises(ValueError):
        overlay_slices(p, var=[1], lag=[2])
    set_theme("journal")


def test_cumulative_slices_use_integer_lags(chicago, cases):
    """The cumulative outcomes have one column per integer lag. Resolving a
    lag against pred.predlag silently plots a different lag when bylag != 1."""
    import numpy as np

    from dlnmpy.plot import plot_crosspred

    cb = dl.crossbasis(chicago.pm10, lag=15, argvar={"fun": "lin"},
                       arglag={"fun": "poly", "degree": 4})
    coef = np.asarray(cases["ex1"]["pred1_pm"]["coef"], float)
    vcov = np.asarray(cases["ex1"]["pred1_pm"]["vcov"], float)
    p = dl.crosspred(cb, coef=coef, vcov=vcov, model_link="log",
                     at=np.arange(21), bylag=0.2, cumul=True)
    assert p.predlag.size != np.shape(p.cumfit)[1]  # 76 vs 16: the trap
    i10 = int(np.argmin(np.abs(p.predvar - 10)))
    truth = np.exp(np.asarray(p.cumfit, float)[i10])
    for lg in (0, 1, 2, 3, 4, 10, 15):
        ax = plot_crosspred(p, "slices", lag=lg, cumul=True)
        line = ax.get_lines()[0]
        got = float(np.interp(10, line.get_xdata(), line.get_ydata()))
        mpl.pyplot.close("all")
        assert abs(got - truth[lg]) < 1e-9, f"lag {lg}: plotted {got}, expected {truth[lg]}"
    with pytest.raises(ValueError, match="integer lags"):
        plot_crosspred(p, "slices", lag=2.5, cumul=True)


def test_overlay_legend_title_and_single_value(pred):
    """The legend indexes the curves, not the x axis; and summary_figure must
    survive a single value, for which overlay_slices draws no legend."""
    from dlnmpy.plot import overlay_slices, summary_figure

    p, _ = pred
    ax = overlay_slices(p, var=[-10.0, 0.0, 20.0], xlab="Lag since exposure")
    assert ax.get_legend().get_title().get_text() == "Var"
    mpl.pyplot.close("all")
    summary_figure(p, var=[20.0])             # used to raise AttributeError
    mpl.pyplot.close("all")
