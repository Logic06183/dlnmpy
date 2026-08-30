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
