"""Cross-basis matrices, lag utilities and knot helpers versus R."""

import numpy as np
import pytest

import dlnmpy as dl

from conftest import assert_close, rcsv


def test_example1(chicago):
    cb = dl.crossbasis(chicago.pm10, lag=15, argvar={"fun": "lin"}, arglag={"fun": "poly", "degree": 4})
    assert_close(cb.matrix, rcsv("cb1_pm"), atol=1e-11)
    assert cb.df == (1, 5) and tuple(cb.lag) == (0, 15)
    cb = dl.crossbasis(chicago.temp, lag=3, argvar={"df": 5}, arglag={"fun": "strata", "breaks": 1})
    assert_close(cb.matrix, rcsv("cb1_temp"), atol=1e-12)
    assert cb.argvar["fun"] == "ns" and cb.arglag["fun"] == "strata"
    assert cb.colnames[:3] == ["v1.l1", "v1.l2", "v2.l1"]


def test_example2_group(chicago):
    seas = chicago[chicago.month.isin([6, 7, 8, 9])]
    cb = dl.crossbasis(seas.o3, lag=5, argvar={"fun": "thr", "thr": 40.3},
                       arglag={"fun": "integer"}, group=seas.year)
    assert_close(cb.matrix, rcsv("cb2_o3"), atol=1e-12)
    assert cb.group == 14
    cb = dl.crossbasis(seas.temp, lag=10, argvar={"fun": "thr", "thr": [15, 25]},
                       arglag={"fun": "strata", "breaks": [2, 6]}, group=seas.year)
    assert_close(cb.matrix, rcsv("cb2_temp"), atol=1e-11)


def test_example3_splines_and_knots(chicago, cases):
    vk = dl.equalknots(chicago.temp, fun="bs", df=5, degree=2)
    lk = dl.logknots(30, 3)
    assert_close(vk, cases["ex3"]["varknots"], atol=1e-12)
    assert_close(lk, cases["ex3"]["lagknots"], atol=1e-12)
    cb = dl.crossbasis(chicago.temp, lag=30, argvar={"fun": "bs", "knots": vk}, arglag={"knots": lk})
    assert_close(cb.matrix, rcsv("cb3_temp"), atol=1e-12)
    cb = dl.crossbasis(chicago.pm10, lag=1, argvar={"fun": "lin"}, arglag={"fun": "strata"})
    assert_close(cb.matrix, rcsv("cb3_pm"), atol=1e-12)


def test_example4_5(chicago, cases):
    lk = dl.logknots(30, 3)
    cb = dl.crossbasis(chicago.temp, lag=30, argvar={"fun": "thr", "thr": [10, 25]}, arglag={"knots": lk})
    assert_close(cb.matrix, rcsv("cb4"), atol=1e-11)
    cb = dl.crossbasis(chicago.temp, lag=21, argvar={"fun": "ns", "df": 4}, arglag={"fun": "ns", "df": 4})
    assert_close(cb.matrix, rcsv("cb5"), atol=1e-12)


def test_matrix_input_and_exphist(chicago, cases):
    Q = dl.exphist(chicago.temp.to_numpy()[:200], lag=(2, 8))
    assert_close(Q, rcsv("exphist_Q"))
    cb = dl.crossbasis(Q, lag=(2, 8), argvar={"fun": "ns", "df": 3}, arglag={"fun": "poly", "degree": 2})
    assert_close(cb.matrix, rcsv("cb6"), atol=1e-12)
    assert_close(dl.exphist([1, 2, 3, 4, 5], times=[2, 5, 7], lag=3), cases["ex6"]["Q2"])
    assert_close(dl.exphist([1, 2, 3, 4, 5], lag=(1, 3), fill=-1), cases["ex6"]["Q3"])


def test_ps_and_penalty(chicago, cases):
    cb = dl.crossbasis(chicago.temp, lag=10, argvar={"fun": "ps", "df": 8}, arglag={"fun": "ps", "df": 5})
    assert_close(cb.matrix, rcsv("cb7"), atol=1e-12)
    pen = dl.cbpen(cb)
    assert_close(pen["Svar"], cases["ex6"]["pen7"]["Svar"], atol=1e-12)
    assert_close(pen["Slag"], cases["ex6"]["pen7"]["Slag"], atol=1e-12)
    assert list(pen["rank"].values()) == cases["ex6"]["pen7"]["rank"]


def test_lag_zero(chicago):
    cb = dl.crossbasis(chicago.temp, lag=0, argvar={"fun": "bs", "df": 5})
    assert_close(cb.matrix, rcsv("cb8"), atol=1e-12)
    assert cb.df == (5, 1)


def test_lag_matrix():
    m = dl.lag_matrix([1, 2, 3, 4, 5], [0, 1, 2])
    assert_close(m[:, 0], [1, 2, 3, 4, 5])
    assert_close(m[:, 1], [np.nan, 1, 2, 3, 4])
    assert_close(m[:, 2], [np.nan, np.nan, 1, 2, 3])
    g = dl.lag_matrix([1, 2, 3, 4, 5, 6], [1], group=[1, 1, 1, 2, 2, 2])
    assert_close(g[:, 0], [np.nan, 1, 2, np.nan, 4, 5])
    assert_close(dl.lag_matrix([1, 2, 3], [-1])[:, 0], [2, 3, np.nan])


def test_mklag_and_seqlag():
    assert tuple(dl.mklag(5)) == (0, 5)
    assert tuple(dl.mklag(-3)) == (-3, 0)
    assert tuple(dl.mklag((2, 8))) == (2, 8)
    with pytest.raises(ValueError):
        dl.mklag((5, 2))
    assert_close(dl.seqlag((0, 2), 0.5), [0, 0.5, 1, 1.5, 2])


def test_knot_helpers(chicago, helpers):
    assert_close(dl.logknots(30, df=4), helpers["logknots_30_df4"], atol=1e-12)
    assert_close(dl.logknots(30, fun="bs", df=5), np.atleast_1d(helpers["logknots_30_df5_bs"]), atol=1e-12)
    assert_close(dl.logknots([3, 21], 2), helpers["logknots_range"], atol=1e-12)
    assert_close(dl.logknots(np.arange(16), df=3), np.atleast_1d(helpers["logknots_vec"]), atol=1e-12)
    assert_close(dl.logknots(10, fun="strata", df=3), helpers["logknots_strata_df3"], atol=1e-12)
    assert_close(dl.equalknots(chicago.temp, df=4), helpers["equalknots_temp_ns"], atol=1e-12)
    assert_close(dl.equalknots(chicago.temp, nk=3), helpers["equalknots_nk"], atol=1e-12)


def test_pretty(helpers):
    for c in helpers["pretty_cases"]:
        a = c["args"]
        res = dl.pretty(a[:2], n=a[2]) if len(a) == 3 else dl.pretty(a)
        assert_close(res, c["result"], atol=1e-12, msg=f"pretty{a}")


def test_dataframe_and_summary(chicago):
    cb = dl.crossbasis(chicago.temp, lag=3, argvar={"df": 3}, arglag={"df": 2})
    df = cb.to_dataframe("cb")
    assert list(df.columns)[:2] == ["cb_v1_l1", "cb_v1_l2"]
    assert "BASIS FOR LAG" in cb.summary()
    assert np.asarray(cb).shape == cb.shape


def test_input_validation(chicago):
    with pytest.raises(ValueError):
        dl.crossbasis(np.ones((10, 3)), lag=5)
    with pytest.raises(ValueError):
        dl.crossbasis(chicago.temp[:20], lag=10, group=np.repeat([1, 2], 10))
