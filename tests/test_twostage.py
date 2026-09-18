"""Two-stage designs: the bridge between dlnmpy and mixmetapy.

The meta-analysis itself moved to the package mixmetapy in 0.8.0 and is
tested there against R's mixmeta. What stays here is the boundary: dlnmpy
must import and run without mixmetapy, the helpers must work on plain
arrays, the deprecated ``dlnmpy.meta`` must either forward to mixmetapy or
say how to install it, and the full pipeline must still match R when both
packages are present.
"""

import json
import subprocess
import sys
import textwrap

import numpy as np
import pandas as pd
import pytest

import dlnmpy as dl

from conftest import FIX, assert_close


@pytest.fixture(scope="module")
def R():
    with open(FIX / "meta.json") as f:
        return json.load(f)


# a child interpreter in which 'import mixmetapy' fails, as if it were not installed
_BLOCK = textwrap.dedent("""
    import sys
    class _Block:
        def find_spec(self, name, path=None, target=None):
            if name == "mixmetapy" or name.startswith("mixmetapy."):
                raise ModuleNotFoundError("No module named 'mixmetapy'", name=name)
    sys.meta_path.insert(0, _Block())
""")


def _run_without_mixmetapy(code):
    return subprocess.run([sys.executable, "-c", _BLOCK + textwrap.dedent(code)],
                          capture_output=True, text=True, timeout=120)


def test_single_location_dlnm_runs_without_mixmetapy():
    r = _run_without_mixmetapy("""
        import sys
        import numpy as np
        import dlnmpy as dl
        x = np.sin(np.arange(400) / 20) * 10 + 15
        cb = dl.crossbasis(x, lag=5, argvar={"fun": "ns", "df": 3}, arglag={"fun": "ns", "df": 3})
        coef = np.linspace(-0.01, 0.01, cb.matrix.shape[1])
        p = dl.crosspred(cb, coef=coef, vcov=np.eye(coef.size) * 1e-4, model_link="log", cen=15, by=1)
        red = dl.crossreduce(cb, coef=coef, vcov=np.eye(coef.size) * 1e-4, model_link="log", cen=15)
        assert "mixmetapy" not in sys.modules
        print("ok", p.allRRfit.shape, red.coef.shape)
    """)
    assert r.returncode == 0, r.stderr
    assert r.stdout.startswith("ok")


def test_meta_says_how_to_install_when_mixmetapy_is_missing():
    r = _run_without_mixmetapy("""
        import dlnmpy
        try:
            import dlnmpy.meta
        except ImportError as e:
            print(e)
        else:
            print("imported")
    """)
    assert r.returncode == 0, r.stderr
    assert "pip install dlnmpy[twostage]" in r.stdout
    r = _run_without_mixmetapy("""
        import dlnmpy as dl
        try:
            dl.mixmeta
        except ImportError as e:
            print(e)
    """)
    assert "pip install dlnmpy[twostage]" in r.stdout, r.stderr


def test_meta_forwards_to_mixmetapy_with_a_warning():
    mixmetapy = pytest.importorskip("mixmetapy")
    import dlnmpy.meta
    with pytest.warns(DeprecationWarning, match="use mixmetapy.mixmeta"):
        f = dlnmpy.meta.mixmeta
    assert f is mixmetapy.mixmeta
    with pytest.warns(DeprecationWarning, match="use mixmetapy.MixMeta"):
        assert dl.MixMeta is mixmetapy.MixMeta
    with pytest.warns(DeprecationWarning, match="use dlnmpy.stack_reduced"):
        assert dlnmpy.meta.stack_reduced is dl.stack_reduced


def test_predict_reduced_takes_plain_arrays():
    """A curve rebuilt by predict_reduced from a reduction's own coefficients
    must equal the reduction's curve: no meta-analysis object involved."""
    x = np.sin(np.arange(600) / 30) * 12 + 16
    cb = dl.crossbasis(x, lag=7, argvar={"fun": "ns", "df": 4}, arglag={"fun": "ns", "df": 3})
    rng = np.random.default_rng(3)
    coef = rng.normal(0, 0.02, cb.matrix.shape[1])
    A = rng.normal(0, 0.01, (coef.size, coef.size))
    vcov = A @ A.T + np.eye(coef.size) * 1e-5
    red = dl.crossreduce(cb, coef=coef, vcov=vcov, model_link="log", cen=16)
    rebuilt = dl.predict_reduced(cb, np.asarray(red.coef), np.asarray(red.vcov), at=red.predvar, cen=16)
    assert isinstance(rebuilt, dl.CrossPred)
    assert_close(rebuilt.allRRfit, red.RRfit, atol=1e-12)
    assert_close(rebuilt.allRRlow, red.RRlow, atol=1e-12)


def test_two_stage_pipeline(R):
    """Stage 1 in Python (crossreduce per city), stage 2 pooled curve versus R."""
    pytest.importorskip("statsmodels")
    mixmeta = pytest.importorskip("mixmetapy").mixmeta
    sim = pd.read_csv(FIX / "meta_sim.csv")
    reds = []
    for k, d in sim.groupby("city"):
        d = d.reset_index(drop=True)
        cb = dl.crossbasis(d.tmean, lag=10, argvar={"fun": "ns", "knots": R["knots"], "boundary_knots": R["bk"]},
                           arglag={"fun": "ns", "df": 3})
        nt = dl.onebasis(d.time, "ns", df=18)
        Xd = d.join(cb.to_dataframe("cb")).join(nt.to_dataframe("nt"))
        fit = dl.fit_glm("y ~ " + " + ".join(list(cb.to_dataframe("cb").columns) + list(nt.to_dataframe("nt").columns)), Xd)
        reds.append(dl.crossreduce(cb, fit, cen=18, name="cb"))
        st = R["stage1"][k - 1]
        assert_close(reds[-1].coef, st["coef"], atol=1e-8)
        assert_close(reds[-1].vcov, st["vcov"], atol=1e-8)
    y, S = dl.stack_reduced(reds)
    assert_close(y, R["y"], atol=1e-8)
    mm = mixmeta(y, S, method="reml")
    pooled = dl.predict_reduced(cb, mm.coef_vec, mm.vcov, at=np.arange(-5, 41), cen=18)
    pc = R["pooled_curve"]
    assert_close(pooled.predvar, pc["at"])
    np.testing.assert_allclose(pooled.allRRfit, pc["allRRfit"], rtol=1e-5)
    np.testing.assert_allclose(pooled.allRRlow, pc["allRRlow"], rtol=1e-4)
    # BLUP curve for one city runs through the same helper
    b = mm.blup(se=True)
    city_curve = dl.predict_reduced(cb, b["blup"][0], b["vcov"][0], at=np.arange(-5, 41), cen=18)
    assert city_curve.allRRfit.shape == (46,)
