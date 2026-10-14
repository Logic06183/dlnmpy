"""Regression tests for the defects found while committing the four
published reproductions (release 0.8.4)."""

import numpy as np
import pytest

import dlnmpy as dl
from dlnmpy.attribution import simulate_coef
from dlnmpy.core import normalise_args


def test_int_is_partial_match_for_intercept():
    # the authors' scripts write `int = TRUE`; R's partial matching reads it as intercept
    assert normalise_args({"int": True}) == {"intercept": True}
    assert normalise_args({"int": False, "Bound": [0, 1]}) == {"intercept": False,
                                                              "boundary_knots": [0, 1]}
    x = np.linspace(0, 10, 50)
    a = dl.onebasis(x, "ns", df=3, int=True)
    b = dl.onebasis(x, "ns", df=3, intercept=True)
    np.testing.assert_array_equal(np.asarray(a), np.asarray(b))


def _rand_vcov(k, seed):
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((k, k))
    return A @ A.T + k * np.eye(k)


def test_simulate_coef_sign_does_not_depend_on_eigensolver():
    # eigenvectors are defined up to sign; the draws must not depend on which
    # sign the LAPACK routine happens to return
    import scipy.linalg as sl
    V = _rand_vcov(6, 1)
    coef = np.arange(6.0)
    Z = np.random.default_rng(2).standard_normal((50, 6))
    ref = simulate_coef(coef, V, normals=Z)
    for driver in ("evr", "evd", "ev"):
        w, U = sl.eigh(V, driver=driver)
        o = np.argsort(w)[::-1]
        U = U[:, o] * np.where(np.random.default_rng(3).random(6) < 0.5, -1, 1)   # arbitrary signs
        big = np.argmax(np.abs(U), axis=0)
        U = U * np.where(U[big, np.arange(6)] < 0, -1.0, 1.0)
        np.testing.assert_allclose(simulate_coef(coef, V, normals=Z, eigvecs=U), ref, atol=1e-10)


def test_simulate_coef_eigvecs_reproduces_mvrnorm_construction():
    # with R's eigenvectors, coef + U diag(sqrt(w)) Z' exactly (MASS::mvrnorm)
    V = _rand_vcov(5, 4)
    coef = np.ones(5)
    w, U = np.linalg.eigh(V)
    o = np.argsort(w)[::-1]
    w, U = w[o], U[:, o]
    U[:, [0, 3]] *= -1                                  # signs as some other solver chose them
    Z = np.random.default_rng(5).standard_normal((20, 5))
    expected = coef[:, None] + U @ (np.sqrt(w)[:, None] * Z.T)
    np.testing.assert_allclose(simulate_coef(coef, V, normals=Z, eigvecs=U), expected, atol=1e-12)
    with pytest.raises(ValueError):
        simulate_coef(coef, V, normals=Z, eigvecs=U[:, :4])


def test_fit_clogit_keeps_names_on_step_halving_path(monkeypatch):
    # force the fallback that runs when statsmodels' Newton diverges to NaN
    # (as on the radon cross-basis of Gasparrini 2014), and check the result
    # still carries the column names crosspred(name=...) needs
    pytest.importorskip("statsmodels")
    import pandas as pd
    from statsmodels.discrete.conditional_models import ConditionalLogit
    orig = ConditionalLogit.fit

    def diverging_fit(self, *a, **k):
        res = orig(self, *a, **k)
        res._results.params = np.full_like(np.asarray(res.params, dtype=float), np.nan)
        return res

    monkeypatch.setattr(ConditionalLogit, "fit", diverging_fit)
    rng = np.random.default_rng(0)
    n_sets, m = 60, 5
    g = np.repeat(np.arange(n_sets), m)
    X = pd.DataFrame({"cb_v1_l1": rng.standard_normal(n_sets * m),
                      "z": rng.standard_normal(n_sets * m)})
    y = np.zeros(n_sets * m)
    for s in range(n_sets):
        idx = np.where(g == s)[0]
        y[idx[np.argmax(X.cb_v1_l1.to_numpy()[idx] + rng.standard_normal(m))]] = 1
    res = dl.fit_clogit(y, X, g)
    assert np.all(np.isfinite(np.asarray(res.params)))
    assert list(res.params.index) == ["cb_v1_l1", "z"]
