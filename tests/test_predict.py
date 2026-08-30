"""crosspred / crossreduce versus R, using R's own coefficients so that the
comparison isolates the DLNM algebra from the regression fitter."""

import numpy as np
import pytest

import dlnmpy as dl
from dlnmpy.predict import mkat, mkcen

from conftest import assert_close


def _pred(cb, m, **kw):
    return dl.crosspred(cb, coef=m["coef"], vcov=m["vcov"], model_link=m["link"], **kw)


def _check_pred(p, r, atol=1e-12):
    assert_close(p.predvar, r["predvar"], msg="predvar")
    assert_close(p.matfit, r["matfit"], atol=atol, msg="matfit")
    assert_close(p.matse, r["matse"], atol=atol, msg="matse")
    assert_close(p.allfit, r["allfit"], atol=atol, msg="allfit")
    assert_close(p.allse, r["allse"], atol=atol, msg="allse")
    assert (p.cen is None and r.get("cen") is None) or np.isclose(p.cen, r["cen"])
    if "cumfit" in r:
        assert_close(p.cumfit, r["cumfit"], atol=atol, msg="cumfit")
        assert_close(p.cumse, r["cumse"], atol=atol, msg="cumse")
    if "allRRfit" in r:
        assert_close(p.allRRfit, r["allRRfit"], atol=1e-10, msg="allRRfit")
        assert_close(p.allRRlow, r["allRRlow"], atol=1e-10, msg="allRRlow")
        assert_close(p.allRRhigh, r["allRRhigh"], atol=1e-10, msg="allRRhigh")
        assert_close(p.matRRfit, r["matRRfit"], atol=1e-10, msg="matRRfit")
    else:
        assert_close(p.alllow, r["alllow"], atol=1e-10, msg="alllow")
        assert_close(p.allhigh, r["allhigh"], atol=1e-10, msg="allhigh")


def test_example1_dlm(chicago, cases):
    cb = dl.crossbasis(chicago.pm10, lag=15, argvar={"fun": "lin"}, arglag={"fun": "poly", "degree": 4})
    p = _pred(cb, cases["ex1"]["model_pm"], at=np.arange(21), bylag=0.2, cumul=True)
    _check_pred(p, cases["ex1"]["pred1_pm"])
    assert p.predlag.size == 76
    # vignette: RR for a 10-unit increase in PM10
    assert np.isclose(p.allRRfit[10], cases["ex1"]["pred1_pm"]["allRRfit"][10])
    cb = dl.crossbasis(chicago.temp, lag=3, argvar={"df": 5}, arglag={"fun": "strata", "breaks": 1})
    p = _pred(cb, cases["ex1"]["model_temp"], by=1, cen=20)
    _check_pred(p, cases["ex1"]["pred1_temp"])


def test_example2_thresholds(chicago, cases):
    seas = chicago[chicago.month.isin([6, 7, 8, 9])]
    cb = dl.crossbasis(seas.o3, lag=5, argvar={"fun": "thr", "thr": 40.3}, arglag={"fun": "integer"}, group=seas.year)
    p = _pred(cb, cases["ex2"]["model_o3"], at=list(range(66)) + [40.3, 50.3])
    _check_pred(p, cases["ex2"]["pred2_o3"])
    with pytest.raises(ValueError):
        _pred(cb, cases["ex2"]["model_o3"], bylag=0.5)
    cb = dl.crossbasis(seas.temp, lag=10, argvar={"fun": "thr", "thr": [15, 25]},
                       arglag={"fun": "strata", "breaks": [2, 6]}, group=seas.year)
    p = _pred(cb, cases["ex2"]["model_temp"], by=1)
    _check_pred(p, cases["ex2"]["pred2_temp"])


def test_example3_dlnm(chicago, cases):
    vk = dl.equalknots(chicago.temp, fun="bs", df=5, degree=2)
    lk = dl.logknots(30, 3)
    cb = dl.crossbasis(chicago.temp, lag=30, argvar={"fun": "bs", "knots": vk}, arglag={"knots": lk})
    m = cases["ex3"]["model_temp"]
    _check_pred(_pred(cb, m, cen=21, by=1), cases["ex3"]["pred3_temp"])
    p = _pred(cb, m, cen=21, at=[-20, 0, 27, 33], lag=(0, 10), bylag=0.5)
    _check_pred(p, cases["ex3"]["pred3_temp_sub"])
    with pytest.raises(ValueError):
        _pred(cb, m, cen=21, lag=(0, 10), cumul=True)
    cb = dl.crossbasis(chicago.pm10, lag=1, argvar={"fun": "lin"}, arglag={"fun": "strata"})
    _check_pred(_pred(cb, cases["ex3"]["model_pm"], at=np.arange(31)), cases["ex3"]["pred3_pm"])


def _check_red(rd, r):
    assert_close(rd.coef, r["coef"], atol=1e-12, msg="coef")
    assert_close(rd.vcov, r["vcov"], atol=1e-14, msg="vcov")
    assert_close(rd.fit, r["fit"], atol=1e-12, msg="fit")
    assert_close(rd.se, r["se"], atol=1e-12, msg="se")
    assert_close(rd.basis, r["basis"], atol=1e-12, msg="basis")
    if r.get("RRfit") is not None:
        assert_close(rd.RRlow, r["RRlow"], atol=1e-10)


def test_example4_reduce(chicago, cases):
    lk = dl.logknots(30, 3)
    cb = dl.crossbasis(chicago.temp, lag=30, argvar={"fun": "thr", "thr": [10, 25]}, arglag={"knots": lk})
    m = cases["ex4"]["model"]
    p = _pred(cb, m, by=1)
    _check_pred(p, cases["ex4"]["pred4"])
    kw = dict(coef=m["coef"], vcov=m["vcov"], model_link=m["link"])
    redall = dl.crossreduce(cb, **kw)
    _check_red(redall, cases["ex4"]["redall"])
    assert redall.coef.size == 2
    redlag = dl.crossreduce(cb, type="lag", value=5, **kw)
    _check_red(redlag, cases["ex4"]["redlag"])
    redvar = dl.crossreduce(cb, type="var", value=33, **kw)
    _check_red(redvar, cases["ex4"]["redvar"])
    assert redvar.coef.size == 5
    _check_red(dl.crossreduce(cb, type="var", value=33, bylag=0.25, **kw), cases["ex4"]["redvar_by"])
    # reduced overall association equals the full prediction
    assert_close(redall.fit, p.allfit, atol=1e-12)
    assert_close(redall.se, p.allse, atol=1e-12)
    with pytest.raises(ValueError):
        dl.crossreduce(cb, type="lag", **kw)
    with pytest.raises(ValueError):
        dl.crossreduce(cb, type="lag", value=40, **kw)


def test_example5_centred_reduce_and_gaussian(chicago, cases):
    cb = dl.crossbasis(chicago.temp, lag=21, argvar={"fun": "ns", "df": 4}, arglag={"fun": "ns", "df": 4})
    m = cases["ex5"]["model"]
    p = _pred(cb, m, cen=21, by=1, cumul=True)
    _check_pred(p, cases["ex5"]["pred5"])
    kw = dict(coef=m["coef"], vcov=m["vcov"], model_link=m["link"])
    _check_red(dl.crossreduce(cb, cen=21, by=1, **kw), cases["ex5"]["red5all"])
    _check_red(dl.crossreduce(cb, type="var", value=33, cen=21, **kw), cases["ex5"]["red5var"])
    _check_red(dl.crossreduce(cb, type="lag", value=3, cen=21, by=1, **kw), cases["ex5"]["red5lag"])
    mg = cases["ex5"]["model_gauss"]
    pg = _pred(cb, mg, cen=21, by=1)
    _check_pred(pg, cases["ex5"]["pred5g"])
    assert not pg.is_exp


def test_onebasis_prediction(chicago, cases):
    ob = dl.onebasis(chicago.temp, "ns", df=4)
    r = cases["ex7"]
    p = dl.crosspred(ob, coef=r["ob_coef"], vcov=r["ob_vcov"], model_link="log", cen=21, by=1)
    _check_pred(p, r["predob"])
    assert tuple(p.lag) == (0, 0)


def test_mkat_and_centring(chicago, cases):
    import json
    from conftest import FIX
    with open(FIX / "mkat.json") as f:
        mk = json.load(f)
    cb = dl.crossbasis(chicago.temp, lag=3, argvar={"df": 3}, arglag={"df": 2})
    assert_close(mkat(None, None, None, None, cb.range, cb.lag, 1), mk["default"])
    assert_close(mkat(None, None, None, 1, cb.range, cb.lag, 1), mk["by1"])
    assert_close(mkat(None, None, None, 2.5, cb.range, cb.lag, 1), mk["by2_5"])
    assert_close(mkat(None, -10, 30, None, cb.range, cb.lag, 1), mk["from_to"])
    assert_close(mkat(None, -10.5, 30, 3, cb.range, cb.lag, 1), mk["from_to_by"])
    assert_close(mkat([30, 0, 0, -20], None, None, None, cb.range, cb.lag, 1), mk["at_unsorted"])
    assert mkcen(None, cb, cb.range) == mk["auto_cen"]
    assert mkcen(False, cb, cb.range) is None
    assert mkcen(10, cb, cb.range) == 10
    cbl = dl.crossbasis(chicago.temp, lag=3, argvar={"fun": "lin"})
    assert mkcen(None, cbl, cbl.range) is None
    assert mkcen(True, cbl, cbl.range) is None
    assert mkcen(20, cbl, cbl.range) == 20


def test_exposure_history_matrix_at(chicago, cases):
    cb = dl.crossbasis(chicago.temp, lag=3, argvar={"df": 3}, arglag={"df": 2})
    coef = np.linspace(-0.01, 0.01, cb.ncol)
    vcov = np.eye(cb.ncol) * 1e-5
    at = np.array([[20, 21, 22, 23], [30, 30, 30, 30]], dtype=float)
    p = dl.crosspred(cb, coef=coef, vcov=vcov, model_link="log", at=at, cen=20)
    assert p.matfit.shape == (2, 4)
    # the row with a constant history at 30 equals the vector prediction at 30
    q = dl.crosspred(cb, coef=coef, vcov=vcov, model_link="log", at=[30], cen=20)
    assert_close(p.allfit[1], q.allfit[0], atol=1e-12)


def test_dataframe_helpers(chicago, cases):
    cb = dl.crossbasis(chicago.temp, lag=21, argvar={"fun": "ns", "df": 4}, arglag={"fun": "ns", "df": 4})
    p = _pred(cb, cases["ex5"]["model"], cen=21, by=1, cumul=True)
    ov = p.overall()
    assert list(ov.columns) == ["var", "fit", "low", "high"] and len(ov) == p.predvar.size
    sv = p.slice_var(33)
    assert len(sv) == 22
    assert len(p.slice_var(33, cumul=True)) == 22
    sl = p.slice_lag(0)
    assert_close(sl["fit"], np.exp(p.matfit[:, 0]))
    assert "centered at: 21" in p.summary()
    rd = dl.crossreduce(cb, cen=21, coef=cases["ex5"]["model"]["coef"], vcov=cases["ex5"]["model"]["vcov"], model_link="log")
    assert list(rd.to_dataframe().columns) == ["var", "fit", "low", "high"]


def test_argvar_cen_is_kept_by_crossbasis(chicago):
    """R stores a cen given inside argvar on the cross-basis and mkcen() uses
    it as the default reference; dropping it silently re-centres the curve."""
    x = chicago.temp.to_numpy()[:400]
    cb = dl.crossbasis(x, lag=5, argvar={"fun": "ns", "df": 3, "cen": 21},
                       arglag={"fun": "ns", "df": 2})
    assert cb.argvar["cen"] == 21
    coef = np.arange(1, cb.ncol + 1) / 50
    vcov = np.eye(cb.ncol) / 100
    p = dl.crosspred(cb, coef=coef, vcov=vcov, model_link="log", at=[-10, 0, 10, 20, 30])
    assert p.cen == 21
    # R: crosspred(cb, coef=coef, vcov=V, model.link="log", at=c(-10,0,10,20,30))$allRRfit
    assert_close(p.allRRfit, [0.852108974969, 0.856612999448, 0.889354799627,
                              0.985363679309, 1.165800623334], atol=1e-11)
    # an explicit cen still wins, and a basis without a stored cen is unchanged
    assert dl.crosspred(cb, coef=coef, vcov=vcov, model_link="log", at=[0], cen=10).cen == 10
    cb2 = dl.crossbasis(x, lag=5, argvar={"fun": "ns", "df": 3}, arglag={"fun": "ns", "df": 2})
    assert dl.crosspred(cb2, coef=coef, vcov=vcov, model_link="log", at=[0]).cen == 5.0


def test_numpy_boolean_cen(chicago):
    """np.bool_ is not a subclass of bool, so a numpy boolean used to fall
    through mkcen's branches and be coerced to 1.0 or 0.0."""
    x = chicago.temp.to_numpy()[:400]
    cb = dl.crossbasis(x, lag=5, argvar={"fun": "ns", "df": 3}, arglag={"fun": "ns", "df": 2})
    coef = np.arange(1, cb.ncol + 1) / 50
    vcov = np.eye(cb.ncol) / 100
    kw = {"coef": coef, "vcov": vcov, "model_link": "log", "at": [0]}
    for truthy in (True, np.True_):
        assert dl.crosspred(cb, cen=truthy, **kw).cen == 5.0
    for falsy in (False, np.False_):
        assert dl.crosspred(cb, cen=falsy, **kw).cen is None
