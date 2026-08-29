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
    assert p.plot("contour") is not None
    plt.close("all")
    assert p.plot("3d") is not None
    plt.close("all")
    assert rd.plot(ci="lines") is not None
    plt.close("all")
    with pytest.raises(ValueError):
        p.plot("slices")
