"""Offline checks against the R results of the published reproductions.

The full reproductions (examples/*_vs_r.py, about an hour in all) download
the authors' data. The stages below need only what is stored in
tests/reference/ and tests/fixtures/, so they run with the ordinary test
suite: the second stage of Gasparrini et al. (2015) and of Gasparrini and
Armstrong (2013) from R's own first-stage reductions, R's simulated
coefficients of Vicedo-Cabrera et al. (2019) from R's normal deviates and
eigenvectors, and every component of attr_table with its intervals.
"""

import json
from pathlib import Path

import numpy as np
import pytest

import dlnmpy as dl
from dlnmpy.attribution import attr_table, attrdl, simulate_coef

from conftest import FIX

REF = Path(__file__).parent / "reference"
mixmeta = pytest.importorskip("mixmetapy").mixmeta


def _ref(name):
    return json.loads((REF / f"{name}.json").read_text())


def test_lancet2015_meta_regression_and_blups_from_r_reductions():
    r = _ref("lancet2015")
    fs, st = r["first_stage"]["converged"], r["settings"]
    rr = r["second_stage"]["converged_mixmeta_reltol14"]            # R at its optimum
    Y, S = np.array(fs["coef"]), [np.array(v) for v in fs["vcov"]]
    X = np.column_stack([np.ones(len(Y)), st["avgtmean"], st["rangetmean"]])
    m = mixmeta(Y, S, X, method="reml")
    # 30 parameters from 10 regions on a flat likelihood: agreement is limited
    # by the optimisers, not by the computation
    np.testing.assert_allclose(m.coef, np.array(rr["coef"]), atol=1e-5)
    np.testing.assert_allclose(m.Psi, np.array(rr["Psi"]), atol=1e-6)
    assert abs(m.loglik - rr["logLik"]) < 1e-7
    np.testing.assert_allclose(m.blup(se=True)["blup"], np.array(rr["blup"]), atol=1e-5)


@pytest.mark.parametrize("key", ["all", "hot", "cold"])
def test_bmcmrm2013_pooling_and_meta_regression_from_r_reductions(key):
    r = _ref("bmcmrm2013")
    red, r14 = r["converged"]["first_stage"]["reductions"][key], r["reltol14"]
    Y, S = np.array(red["coef"]), [np.array(v) for v in red["vcov"]]
    m = mixmeta(Y, S, method="reml")
    p = r14["pooled"][key]
    np.testing.assert_allclose(np.ravel(m.coef), np.ravel(p["coef"]), atol=1e-8)
    np.testing.assert_allclose(m.Psi, np.array(p["Psi"]), atol=1e-8)
    assert abs(m.loglik - p["logLik"]) < 1e-8
    q = m.qtest()
    np.testing.assert_allclose([np.atleast_1d(q[k])[0] for k in ("Q", "df", "I2")],
                               [np.atleast_1d(p["qtest"][k])[0] for k in ("Q", "df", "I2")], rtol=1e-8)
    lat = np.column_stack([np.ones(len(Y)), r["settings"]["lat"]])
    mr = mixmeta(Y, S, lat, method="reml")
    np.testing.assert_allclose(np.ravel(mr.coef), np.ravel(r14["metareg"][key]["coef"]), atol=1e-5)
    assert abs(mr.loglik - r14["metareg"][key]["logLik"]) < 1e-8


def test_vicedo2019_simulated_coefficients_replay_r_draws():
    r = _ref("vicedo2019")
    coef = np.asarray(r["red"]["coef"], float)
    k = coef.size
    vcov = np.asarray(r["red"]["vcov"], float).reshape(k, k)
    nsim = r["nsim"]
    for g, z in r["normals"].items():
        Z = np.asarray(z, float).reshape((nsim, k), order="F")
        E = np.asarray(r["eigvecs"][g], float).reshape((k, k), order="F")
        sims = simulate_coef(coef, vcov, normals=Z, eigvecs=E)
        head = np.asarray(r["coefsim_head"][g]).reshape((10, k), order="F")
        np.testing.assert_allclose(sims[:, :10].T, head, rtol=1e-10, atol=1e-13)


@pytest.fixture(scope="module")
def chicago_attr(chicago):
    r = json.loads((FIX / "attribution.json").read_text())
    cb = dl.crossbasis(chicago.temp, lag=21,
                       argvar={"fun": "bs", "degree": 2,
                               "knots": dl.quantile7(chicago.temp, [0.10, 0.75, 0.90])},
                       arglag={"knots": dl.logknots(21, 3)})
    return r, cb, dict(coef=r["coef"], vcov=r["vcov"], model_link="log")


def test_attr_table_components_match_attrdl_r(chicago, chicago_attr):
    # attrdl.R values for the total, cold, heat and extreme heat
    r, cb, kw = chicago_attr
    t = attr_table(chicago.temp, cb, chicago.death, cen=r["mmt"], nsim=50, seed=1, **kw).set_index("component")
    for comp, key in (("total", "an_back"), ("cold", "an_cold"), ("heat", "an_heat"),
                      ("extreme heat", "an_extreme_heat")):
        assert np.isclose(t.loc[comp, "an"], r[key], rtol=1e-10), comp
    np.testing.assert_allclose(t.loc["total", "af"], r["af_back"], rtol=1e-10)


def test_attr_table_intervals_are_the_quantiles_of_shared_draws(chicago, chicago_attr):
    # every component's interval comes from one set of draws, the same as
    # attrdl called with those draws (which is checked against R's own draws
    # in test_attribution.py::test_attrdl_simulation)
    r, cb, kw = chicago_attr
    nsim, seed, cen = 300, 11, r["mmt"]
    t = attr_table(chicago.temp, cb, chicago.death, cen=cen, nsim=nsim, seed=seed, **kw)
    draws = simulate_coef(np.asarray(r["coef"]), np.asarray(r["vcov"]), nsim, seed)
    for _, row in t.iterrows():
        sims = attrdl(chicago.temp, cb, chicago.death, type="an", cen=cen, range=row["range"],
                      sim=True, coefsim=draws, **kw)
        np.testing.assert_allclose([row.an_low, row.an_high], np.quantile(sims, [0.025, 0.975]), rtol=1e-12)
