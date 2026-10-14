"""Reproduce Gasparrini & Leone (2014), attributable risk from DLNMs.

    Gasparrini A, Leone M. Attributable risk from distributed lag models.
    BMC Medical Research Methodology 2014;14:55. doi:10.1186/1471-2288-14-55

A line-by-line translation of the authors' code and data,
https://github.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata at commit
68a8d444f85c673868bdd1511619dac93af6f563 (01.runmodel.R, 02.results.R and
the attrdl() calls of 03.graphs.R), on London 1993-2006: temperature and
all-cause mortality, quadratic B-spline of temperature with two equally
spaced knots, natural spline of lag 0-25 with three log-spaced knots,
centred at 20 C.

Every cell of Tables 1 and 2 for London is recomputed (overall, cold, heat,
mild and extreme cold and heat; backward and forward; attributable number
and fraction) with 95% empirical intervals from 5,000 draws, and checked
against tests/reference/leone2014.json, written by
tools/reproduce_leone2014.R from the authors' own scripts. The intervals
replay R's standard normal deviates, so they are compared without Monte
Carlo noise. They are replayed twice: through simulate_coef(normals=..., eigvecs=...),
and with this model's eigenvectors given the signs R's eigen() chose (an
eigenvector is defined only up to sign, and numpy and R pick differently
for some columns, which sends the same deviates to different draws).

The data (london.csv, 650 kB) are downloaded once into examples/ from the
pinned commit and checked against a hardcoded SHA-256.

Tolerances, |python - R| <= atol + rtol * |R| (values in TOL below). Two
references:

(a) The authors' fit as published. glm() stops at its default epsilon of
    1e-8 after 3 iterations, 1.5e-8 (relative) short of convergence in the
    coefficients, and summary.glm() computes the dispersion with the
    working weights of the previous iteration, 6.5e-6 above the Pearson
    statistic of the final fit. statsmodels iterates further. Tolerances
    follow from that: rtol 1e-7 on coefficients, 1e-9 absolute on
    fractions, 1e-8 relative on numbers and relative risks, 2e-5 on
    anything scaled by the dispersion (covariance) and 1e-6 / 2e-5 on what
    scales with its square root (confidence limits / eCI bounds).
(b) The same model refitted in R with epsilon = 1e-15 (not in the authors'
    code). This separates porting error from R's stopping rule; tolerances
    there are 1e-13 absolute on fractions and 1e-10 to 1e-12 relative
    elsewhere.
"""

import hashlib
import io
import json
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

import dlnmpy as dl

SHA = "68a8d444f85c673868bdd1511619dac93af6f563"
DATA_URL = ("https://raw.githubusercontent.com/gasparrini/"
            f"2014_gasparrini_BMCmrm_Rcodedata/{SHA}/london.csv")
DATA_SHA256 = "62387cba84d5dcf0443fc048b3df0f3b589e6c16d525277b7222e79ec75f64d1"
CACHE = Path(__file__).with_name("london.csv")
REFERENCE = Path(__file__).resolve().parents[1] / "tests" / "reference" / "leone2014.json"

# London rows of Tables 1 and 2 as printed in the paper (b-AF, f-AF, %), for
# the record only: the intervals there come from other random draws
PAPER = {("overall", "back"): "13.59 (10.04-17.09)", ("overall", "forw"): "13.41 (9.72-16.87)",
         ("cold", "back"): "12.95 (9.32-16.38)", ("cold", "forw"): "12.84 (9.38-16.33)",
         ("heat", "back"): "0.66 (0.52-0.80)", ("heat", "forw"): "0.57 (0.45-0.68)",
         ("extreme_cold", "back"): "0.55 (0.45-0.64)", ("extreme_cold", "forw"): "0.47 (0.40-0.53)",
         ("mild_cold", "back"): "12.48 (8.86-15.88)", ("mild_cold", "forw"): "12.38 (8.98-15.78)",
         ("mild_heat", "back"): "0.31 (0.23-0.38)", ("mild_heat", "forw"): "0.29 (0.22-0.35)",
         ("extreme_heat", "back"): "0.36 (0.29-0.43)", ("extreme_heat", "forw"): "0.28 (0.23-0.33)"}


# |python - R| <= atol + rtol * |R|; see the module docstring for the reasons
TOL = dict(
    # (a) against the authors' fit, which glm() stops 1.5e-8 short of convergence
    coef=1e-7,      # rtol, cross-basis coefficients
    vcov=2e-5,      # rtol, and atol relative to max|vcov|: dispersion and covariance
    af=1e-9,        # atol, attributable fractions (0.003-0.14)
    an=1e-8,        # rtol, attributable numbers (2,400-115,000)
    af_ci=2e-5,     # rtol, eCI bounds of fractions (scale with sqrt(dispersion))
    an_ci=2e-5,     # rtol, eCI bounds of numbers
    daily=1e-7,     # atol, daily attributable deaths (up to ~60)
    rr=1e-8,        # rtol, relative risks
    rr_ci=1e-6,     # rtol, their 95% limits (scale with sqrt(dispersion))
    # (b) against the same model run to convergence in R
    c_coef=1e-10, c_vcov=1e-9, c_af=1e-13, c_an=1e-12, c_rr=1e-11)


def load() -> pd.DataFrame:
    if not CACHE.exists():
        print(f"downloading {DATA_URL} ...")
        with urllib.request.urlopen(DATA_URL) as r:
            raw = r.read()
        if hashlib.sha256(raw).hexdigest() != DATA_SHA256:
            raise RuntimeError("london.csv does not match the pinned SHA-256")
        CACHE.write_bytes(raw)
    raw = CACHE.read_bytes()
    if hashlib.sha256(raw).hexdigest() != DATA_SHA256:
        raise RuntimeError(f"{CACHE} does not match the pinned SHA-256; delete it to re-download")
    return pd.read_csv(io.BytesIO(raw))


class Checker:
    """Collects |python - R| for every compared scalar."""

    def __init__(self):
        self.n = 0
        self.fail = 0
        self.maxabs = 0.0
        self.groups = {}

    def __call__(self, group, py, r, atol=0.0, rtol=0.0):
        py = np.asarray(py, dtype=float).ravel()
        r = np.asarray([np.nan if v is None else v for v in np.ravel(np.asarray(r, dtype=object))],
                       dtype=float)
        if py.shape != r.shape:
            raise ValueError(f"{group}: shape {py.shape} vs R {r.shape}")
        nan_ok = np.isnan(py) == np.isnan(r)
        both = ~np.isnan(py) & ~np.isnan(r)
        d = np.abs(py[both] - r[both])
        bad = int((~nan_ok).sum() + (d > atol + rtol * np.abs(r[both])).sum())
        worst = float(d.max()) if d.size else 0.0
        rel = float((d / np.maximum(np.abs(r[both]), 1e-300)).max()) if d.size else 0.0
        self.n += py.size
        self.fail += bad
        self.maxabs = max(self.maxabs, worst)
        g = self.groups.setdefault(group, [0, 0.0, 0.0, 0])
        g[0] += py.size
        g[1] = max(g[1], worst)
        g[2] = max(g[2], rel)
        g[3] += bad


def main() -> int:
    t0 = time.perf_counter()
    ref = json.loads(REFERENCE.read_text())
    lndn = load()
    tmean = lndn.tmean.to_numpy(float)
    death = lndn.death.to_numpy(float)

    # ---- 01.runmodel.R ------------------------------------------------------
    vk = dl.equalknots(tmean, fun="bs", df=4, degree=2)
    maxlag = 25
    lk = dl.logknots(maxlag, 3)
    cen = 20
    cb = dl.crossbasis(tmean, lag=maxlag, argvar={"fun": "bs", "degree": 2, "knots": vk},
                       arglag={"knots": lk})
    perc = dl.quantile7(tmean, [0.01, 0.99])

    seas = dl.onebasis(lndn.time.to_numpy(float), "ns", df=10 * 14)
    X = dl.design_matrix(lndn, ("cb", cb), ("seas", seas), intercept=False)
    model = dl.fit_glm("death ~ " + " + ".join(X.columns) + " + C(dow)", lndn.join(X),
                       family="quasipoisson")

    rng_ = np.linspace(np.nanmin(tmean), np.nanmax(tmean), 25)
    cp1 = dl.crosspred(cb, model, at=rng_, cen=cen, name="cb")
    cp2 = dl.crosspred(cb, model, at=perc, bylag=0.2, cen=cen, name="cb")
    cp3 = dl.crosspred(cb, model, at=[25, 27, 29], bylag=0.2, cen=cen, name="cb")

    # ---- 02.results.R and Tables 1-2 ---------------------------------------
    names = [c for c in model.params.index if c.startswith("cb")]
    coef = model.params[names].to_numpy()
    vcov = model.cov_params().loc[names, names].to_numpy()
    s = ref["settings"]
    Z = np.asarray(ref["normals"], dtype=float)           # (nsim, k), R's rnorm draws

    # Two replays of R's draws. (1) simulate_coef(normals=Z), the documented
    # route. (2) The same computation, coef + U diag(sqrt(w)) Z', with this
    # model's own eigenvectors U given the signs of R's eigen(): an
    # eigenvector is only defined up to sign, numpy's eigh and R's eigen
    # (LAPACK dsyevd vs dsyevr) choose differently for some of them, and a
    # flipped column sends the same deviates to different draws.
    E = np.asarray(ref["model"]["eigen_vectors"], dtype=float).reshape(vcov.shape, order="F")
    coefsim_api = dl.simulate_coef(coef, vcov, normals=Z, eigvecs=E)
    w, U = np.linalg.eigh(vcov)
    order = np.argsort(w)[::-1]
    w, U = np.maximum(w[order], 0), U[:, order]
    flip = np.sign(np.sum(U * np.asarray(ref["model"]["eigen_vectors"]), axis=0))
    coefsim = coef[:, None] + (U * flip) @ (np.sqrt(w)[:, None] * Z.T)

    ranges = {"overall": None, "cold": (-100, cen), "heat": (cen, 100),
              "extreme_cold": (-100, perc[0]), "mild_cold": (perc[0], cen),
              "mild_heat": (cen, perc[1]), "extreme_heat": (perc[1], 100)}

    def eci(sims):
        return np.quantile(sims, [0.025, 0.975])           # R's quantile type 7

    # attrdl() returns the number as af * den with den = sum(cases), for the
    # point estimate and for every draw alike (as attrdl.R does), so one call
    # per quantity gives both; R's calls with type="an" are the same product.
    # Both sets of draws go through a single call. This halves the runtime,
    # which is dominated by rebuilding the 5,114 x 26 lagged basis per call.
    den = np.nansum(death)
    nsim = Z.shape[0]
    attr = {}
    for comp, rg in ranges.items():
        for d in ("back", "forw"):
            kw = dict(dir=d, cen=cen, range=rg, name="cb")
            a = attr[comp, d] = {}
            a["af"] = dl.attrdl(tmean, cb, death, model, type="af", **kw)
            a["an"] = a["af"] * den
            afsim = dl.attrdl(tmean, cb, death, model, type="af", sim=True,
                              coefsim=np.hstack([coefsim, coefsim_api]), **kw)
            for tag, sims in (("", afsim[:nsim]), ("_api", afsim[nsim:])):
                a["af_ci" + tag] = eci(sims)
                a["an_ci" + tag] = eci(sims * den)

    daily = {}
    for comp, rg in (("cold", (-100, cen)), ("heat", (cen, 100))):
        for d in ("back", "forw"):
            daily[f"{comp}_{d}"] = dl.attrdl(tmean, cb, death, model, tot=False, type="an",
                                             dir=d, cen=cen, range=rg, name="cb")

    # ---- compare with R -----------------------------------------------------
    chk = Checker()
    chk("settings: knots, percentiles", np.r_[vk, lk, perc],
        np.r_[s["varknots"], s["lagknots"], s["perc"]], atol=1e-12)
    # (a) against the authors' fit as published (glm defaults)
    m = ref["model"]
    vscale = np.abs(np.asarray(m["vcov"])).max()
    chk("cb coefficients", coef, m["coef"], rtol=TOL["coef"])
    chk("dispersion, vcov", np.r_[model.scale, vcov.ravel()], np.r_[m["dispersion"], np.ravel(m["vcov"])],
        rtol=TOL["vcov"], atol=TOL["vcov"] * vscale)
    for row in ref["attributable"]:
        a = attr[row["component"], row["dir"]]
        chk("attributable fractions", a["af"], row["af"], atol=TOL["af"])
        chk("attributable numbers", a["an"], row["an"], rtol=TOL["an"])
        chk("eCI, R draws (sign-aligned)", a["af_ci"], row["af_ci"], rtol=TOL["af_ci"])
        chk("eCI, R draws (sign-aligned)", a["an_ci"], row["an_ci"], rtol=TOL["an_ci"])
        chk("eCI, simulate_coef(normals=, eigvecs=)", a["af_ci_api"], row["af_ci"], rtol=TOL["af_ci"])
        chk("eCI, simulate_coef(normals=, eigvecs=)", a["an_ci_api"], row["an_ci"], rtol=TOL["an_ci"])
    for k, v in daily.items():
        chk("daily attributable numbers", v, ref["daily"][k], atol=TOL["daily"])
    for key, cp in (("cp1", cp1), ("cp2", cp2)):
        p = ref["predictions"][key]
        chk("relative risks", cp.allRRfit, p["allRRfit"], rtol=TOL["rr"])
        chk("RR confidence limits", np.r_[cp.allRRlow, cp.allRRhigh],
            np.r_[p["allRRlow"], p["allRRhigh"]], rtol=TOL["rr_ci"])
    p = ref["predictions"]["cp3"]
    chk("relative risks", cp3.matRRfit.ravel(), np.ravel(p["matRRfit"]), rtol=TOL["rr"])
    chk("RR confidence limits", np.r_[cp3.matRRlow.ravel(), cp3.matRRhigh.ravel()],
        np.r_[np.ravel(p["matRRlow"]), np.ravel(p["matRRhigh"])], rtol=TOL["rr_ci"])

    # (b) against the same model run to convergence in R (epsilon = 1e-15)
    mc = ref["model_converged"]
    chk("converged: coefficients", coef, mc["coef"], rtol=TOL["c_coef"])
    chk("converged: dispersion, vcov", np.r_[model.scale, vcov.ravel()],
        np.r_[mc["dispersion"], np.ravel(mc["vcov"])], rtol=TOL["c_vcov"], atol=TOL["c_vcov"] * vscale)
    for row in mc["attributable"]:
        a = attr[row["component"], row["dir"]]
        chk("converged: AF", a["af"], row["af"], atol=TOL["c_af"])
        chk("converged: AN", a["an"], row["an"], rtol=TOL["c_an"])
    chk("converged: RR and limits", np.r_[cp2.allRRfit, cp2.allRRlow, cp2.allRRhigh],
        np.r_[mc["cp2"]["allRRfit"], mc["cp2"]["allRRlow"], mc["cp2"]["allRRhigh"]], rtol=TOL["c_rr"])

    # ---- report -------------------------------------------------------------
    v = ref["versions"]
    print(f"Gasparrini & Leone 2014, London 1993-2006: {int(s['ndeaths']):,} deaths, "
          f"centred at {cen} C (percentile {100 * s['cen_percentile']:.1f}), "
          f"cut-offs {perc[0]:.1f} / {perc[1]:.1f} C")
    print(f"R reference: {v['R']}, dlnm {v['dlnm']}, authors' commit {v['authors_commit'][:7]}")
    print(f"overall RR at 1st / 99th percentile: "
          f"{cp2.allRRfit[0]:.3f} ({cp2.allRRlow[0]:.3f}-{cp2.allRRhigh[0]:.3f}) / "
          f"{cp2.allRRfit[1]:.3f} ({cp2.allRRlow[1]:.3f}-{cp2.allRRhigh[1]:.3f})\n")
    print(f"  {'component':13s} {'dir':4s}  {'AF % (95% eCI), python':24s} {'AN':>8s}   paper")
    for (comp, d), a in attr.items():
        lo, hi = 100 * a["af_ci"]
        print(f"  {comp:13s} {d:4s}  {100 * a['af']:5.2f} ({lo:5.2f}-{hi:5.2f})        "
              f"{a['an']:8.0f}   {PAPER[comp, d]}")
    print("\n  (paper intervals used different random draws; point estimates agree to 2 dp)\n")
    print(f"  {'compared with R':30s} {'n':>6s} {'max abs diff':>13s} {'max rel diff':>13s}  fail")
    for g, (n, mx, rel, bad) in chk.groups.items():
        print(f"  {g:30s} {n:6d} {mx:13.2e} {rel:13.2e}  {bad:4d}")
    if chk.groups["eCI, simulate_coef(normals=, eigvecs=)"][3]:
        print(f"\n  simulate_coef(normals=R draws) does not reproduce R's draws: numpy's eigh "
              f"returned {int((flip < 0).sum())} of {flip.size} eigenvectors with the opposite "
              f"sign to R's eigen(), so the same deviates map to different coefficient draws.")
    print(f"\n  runtime {time.perf_counter() - t0:.1f} s")
    print(f"leone2014: {chk.n} comparisons with R, max abs difference {chk.maxabs:.2e}, "
          f"{chk.fail} failures")
    return 1 if chk.fail else 0


if __name__ == "__main__":
    sys.exit(main())
