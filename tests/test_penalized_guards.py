"""Targeted tests for three penalised-fitting defects that the ledger
(docs/defects.md) listed without one: D16, D32 and D38."""

import numpy as np
import pandas as pd
import pytest

import dlnmpy as dl
import dlnmpy.penalized as pen

pytest.importorskip("statsmodels")


def _data(seed=0, n=600):
    rng = np.random.default_rng(seed)
    x = 15 + 8 * np.sin(np.arange(n) / 58) + rng.normal(0, 2, n)
    cb = dl.crossbasis(x, lag=5, argvar={"fun": "ps", "df": 6}, arglag={"fun": "ps", "df": 5})
    pop = rng.uniform(800, 1200, n)
    eta = np.log(pop / 1000) + 2 + 0.002 * (x - 15) ** 2
    d = pd.DataFrame({"y": rng.poisson(np.exp(eta)).astype(float), "logpop": np.log(pop / 1000)})
    return d.join(cb.to_dataframe("cb")), cb


def test_d16_full_length_offset_follows_the_rows_dropped_for_missing_values():
    d, cb = _data()
    d.loc[[40, 41, 300], "y"] = np.nan                    # missing outcomes, beyond the lag rows
    f = "y ~ " + " + ".join(c for c in d.columns if c.startswith("cb"))
    full = dl.fit_pgam(f, d, penalties={"cb": dl.cbpen(cb)}, family="poisson",
                       sp=[1.0, 1.0], offset=d.logpop.to_numpy())
    ok = d.dropna()
    sub = dl.fit_pgam(f, ok, penalties={"cb": dl.cbpen(cb)}, family="poisson",
                      sp=[1.0, 1.0], offset=ok.logpop.to_numpy())
    np.testing.assert_allclose(np.asarray(full.params), np.asarray(sub.params), rtol=1e-10, atol=1e-12)


def test_d32_converged_is_false_when_the_search_fails(monkeypatch):
    d, cb = _data(1)
    f = "y ~ " + " + ".join(c for c in d.columns if c.startswith("cb"))
    real = pen.optimize.minimize

    def stalls(fun, x0, method=None, **kw):
        if method == "Nelder-Mead":                       # stop far from the optimum, unsuccessfully
            x = np.asarray(x0, float) + 3.0
            return pen.optimize.OptimizeResult(x=x, fun=fun(x), success=False, nit=1)
        return real(fun, x0, method=method, **kw)

    monkeypatch.setattr(pen.optimize, "minimize", stalls)
    with pytest.warns(pen.ConvergenceWarning, match="did not converge"):
        res = dl.fit_pgam(f, d, penalties={"cb": dl.cbpen(cb)}, family="poisson", method="reml")
    assert res.converged is False


def test_d38_a_singular_trial_point_does_not_abort_the_search(monkeypatch):
    d, cb = _data(2)
    f = "y ~ " + " + ".join(c for c in d.columns if c.startswith("cb"))
    ref = dl.fit_pgam(f, d, penalties={"cb": dl.cbpen(cb)}, family="poisson", method="reml")
    real, calls = pen._pirls, {"n": 0}

    def sometimes_singular(*a, **k):
        calls["n"] += 1
        if calls["n"] == 12:                              # one trial point inside the search
            raise np.linalg.LinAlgError("Singular matrix")
        return real(*a, **k)

    monkeypatch.setattr(pen, "_pirls", sometimes_singular)
    res = dl.fit_pgam(f, d, penalties={"cb": dl.cbpen(cb)}, family="poisson", method="reml")
    assert calls["n"] > 12 and np.all(np.isfinite(np.asarray(res.params)))
    np.testing.assert_allclose(res.reml, ref.reml, rtol=1e-6)
