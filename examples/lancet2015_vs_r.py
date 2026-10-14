"""Stage-by-stage comparison of examples/lancet_2015.py with the authors' R code.

    Gasparrini A, Guo Y, Hashizume M, et al. Mortality risk attributable to
    high and low ambient temperature: a multicountry observational study.
    The Lancet 2015;386(9991):369-375. doi:10.1016/S0140-6736(14)62114-0

Recomputes the five stages of
https://github.com/gasparrini/2015_gasparrini_Lancet_Rcodedata (commit
7c3e03c425529d2f210b1659a01bbc582448b9c0) on the 10 regions of England and
Wales with dlnmpy and mixmetapy, as examples/lancet_2015.py does, and checks
every intermediate against tests/reference/lancet2015.json, written by
tools/reproduce_lancet2015.R from the authors' scripts. The empirical
intervals replay R's standard normal deviates and R's eigenvectors through
simulate_coef(normals=, eigvecs=), so they are compared without Monte Carlo
noise. examples/lancet_2015.py is left as it is; it draws its own.

The data (regEngWales.csv, 7 MB) are downloaded once into examples/ from the
pinned commit and checked against a hardcoded SHA-256.

Three R references, |python - R| <= atol + rtol * |R| (values in TOL below):

(a) "published": the authors' code as written. glm() stops at its default
    epsilon = 1e-8 and computes the dispersion from the working weights of
    the previous IRLS iteration, so the reduced covariances differ by up to
    2e-6 (relative); the second stage is mvmeta(). Everything downstream
    inherits that.
(b) "converged": the same first stage refitted in R with
    glm.control(epsilon = 1e-14), the second stage with mixmeta() (whose
    REML log-likelihood constant mixmetapy follows) at its default control,
    and stages 3-4 rerun from those BLUPs. The first stage then agrees to
    ~1e-13.
(c) "reltol14": the second stage of (b) with mixmeta()'s optimiser
    tolerance at 1e-14 in place of its default 1.5e-8, i.e. at the REML
    optimum.

The meta-regression is the weak point. It estimates 15 fixed coefficients
and 15 parameters of Psi from 10 regions (the authors' script says "clearly
overparameterized"), so the REML surface is flat: a log-likelihood 1e-8
short of the optimum leaves the parameters ~1e-5 away. That is where R's
default stopping rule leaves both mvmeta() and mixmeta() (they land 5e-6
apart on identical inputs; at reltol 1e-14 they agree to 6e-8), while
mixmetapy stops within 1e-9 of the optimum. Against (a) and (b) the
meta-regression parameters, the BLUPs and what is computed from them are
therefore held to 5e-4 relative, against (c) to 5e-5. The minimum-mortality
percentile is an integer and is required to be exact.
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
from mixmetapy import mixmeta  # pip install dlnmpy[twostage]

SHA = "7c3e03c425529d2f210b1659a01bbc582448b9c0"
DATA_URL = ("https://raw.githubusercontent.com/gasparrini/"
            f"2015_gasparrini_Lancet_Rcodedata/{SHA}/regEngWales.csv")
DATA_SHA256 = "6f7b1471e641ed3b403b4d912388e57c1007fb91c3329b525fb1f689e31f29c4"
CACHE = Path(__file__).with_name("regEngWales.csv")
REFERENCE = Path(__file__).resolve().parents[1] / "tests" / "reference" / "lancet2015.json"

VARFUN, VARDEGREE, VARPER = "bs", 2, [0.10, 0.75, 0.90]
LAG, LAGNK, DFSEAS = 21, 3, 8

# |python - R| <= atol + rtol * |R|; reasons in the module docstring
TOL = dict(
    exact=1e-12,      # rtol, settings: knots, means, ranges, death totals
    # (a) against the authors' code as published
    p_coef=1e-8,      # atol, reduced coefficients: glm() stops at a relative
                      # deviance change of 1e-8; IRLS converges quadratically, so
                      # the coefficients are left far closer than that
    p_vcov=1e-5,      # rtol on max|vcov|: R's dispersion uses the working weights
                      # of the previous IRLS iteration (as in leone2014)
    p_meta=5e-4,      # rtol on max|.|, meta-regression coef, vcov, Psi vs mvmeta
    p_blup=5e-4,      # rtol on max|.|, BLUPs and their vcov
    p_mmt=1e-9,       # atol, minimum-mortality temperature (a quantile of the data)
    p_an=5e-4,        # rtol, attributable deaths; eCI bounds relative to the point
    p_af=5e-4,        # rtol, attributable fractions (%); eCI bounds likewise
    # (b) against the same pipeline run to convergence in R, with mixmeta
    c_coef=1e-10,     # rtol on max|.|, reduced coefficients
    c_vcov=1e-9,      # rtol on max|.|, reduced covariance
    c_loglik=1e-6,    # atol, REML log-likelihood (52.39)
    c_meta=5e-4,      # rtol on max|.|, meta-regression coef, vcov, Psi: R's default
                      # optimiser tolerance on a flat surface (see the docstring)
    c_blup=5e-4,      # rtol on max|.|, BLUPs and their vcov
    c_mmt=1e-9,       # atol, minimum-mortality temperature
    c_an=5e-4,        # rtol, attributable deaths; eCI bounds relative to the point
    c_af=5e-4,        # rtol, attributable fractions (%); eCI bounds likewise
    # (c) second stage only, against mixmeta with reltol = 1e-14 (the optimum)
    t_loglik=1e-8,    # atol
    t_meta=5e-5)      # rtol on max|.|, coef, vcov, Psi, BLUPs


def load() -> pd.DataFrame:
    if not CACHE.exists():
        print(f"downloading {DATA_URL} ...")
        with urllib.request.urlopen(DATA_URL) as r:          # noqa: S310
            raw = r.read()
        if hashlib.sha256(raw).hexdigest() != DATA_SHA256:
            raise RuntimeError("regEngWales.csv does not match the pinned SHA-256")
        CACHE.write_bytes(raw)
    raw = CACHE.read_bytes()
    if hashlib.sha256(raw).hexdigest() != DATA_SHA256:
        raise RuntimeError(f"{CACHE} does not match the pinned SHA-256; delete it to re-download")
    d = pd.read_csv(io.BytesIO(raw), index_col=0)
    d["date"] = pd.to_datetime(d["date"])
    return d


class Checker:
    """Collects |python - R| for every compared scalar."""

    def __init__(self):
        self.n = 0
        self.fail = 0
        self.maxabs = 0.0
        self.groups = {}

    def __call__(self, group, py, r, atol=0.0, rtol=0.0, scale=None):
        """``scale`` (scalar or array) replaces |R| in the relative bound: the
        largest entry of a matrix, or the point estimate for an interval
        bound, which can sit near zero (heat in the North East: -8.7 deaths)."""
        py = np.asarray(py, dtype=float).ravel()
        r = np.asarray([np.nan if v is None else v for v in np.ravel(np.asarray(r, dtype=object))],
                       dtype=float)
        if py.shape != r.shape:
            raise ValueError(f"{group}: shape {py.shape} vs R {r.shape}")
        nan_ok = np.isnan(py) == np.isnan(r)
        both = ~np.isnan(py) & ~np.isnan(r)
        d = np.abs(py[both] - r[both])
        ref = (np.abs(r[both]) if scale is None
               else np.broadcast_to(np.abs(np.asarray(scale, dtype=float)).ravel(), r.shape)[both])
        bad = int((~nan_ok).sum() + (d > atol + rtol * ref).sum())
        worst = float(d.max()) if d.size else 0.0
        rel = float((d / np.maximum(ref, 1e-300)).max()) if d.size else 0.0
        self.n += py.size
        self.fail += bad
        self.maxabs = max(self.maxabs, worst)
        g = self.groups.setdefault(group, [0, 0.0, 0.0, 0])
        g[0] += py.size
        g[1] = max(g[1], worst)
        g[2] = max(g[2], rel)
        g[3] += bad

    def block(self, group, py, r, rtol):
        """A matrix compared relative to its largest entry."""
        r = np.asarray(r, dtype=float)
        self(group, py, r, rtol=rtol, scale=float(np.abs(r).max()))


def main() -> int:
    t0 = time.perf_counter()
    ref = json.loads(REFERENCE.read_text())
    s = ref["settings"]
    d = load()

    # ---- 00.prepdata.R: regions ordered by name -----------------------------
    regions = list(dict.fromkeys(d.regnames))
    citynames = ["North East", "North West", "Yorkshire & Humber", "East Midlands",
                 "West Midlands", "East", "London", "South East", "South West", "Wales"]
    cities = (pd.DataFrame({"city": regions, "cityname": citynames})
              .sort_values("cityname").reset_index(drop=True))
    if list(cities.city) != s["regions"]:
        raise RuntimeError("region order differs from R's")
    dlist = {c: d[d.regnames == c] for c in cities.city}
    m = len(cities)

    # ---- 01.firststage.R ----------------------------------------------------
    coefs, vcovs, keep = [], [], {}
    for city in cities.city:
        data = dlist[city]
        tmean = data.tmean.to_numpy(float)
        knots = dl.quantile7(tmean, VARPER)
        cb = dl.crossbasis(data.tmean, lag=LAG,
                           argvar={"fun": VARFUN, "knots": knots, "degree": VARDEGREE},
                           arglag={"knots": dl.logknots(LAG, LAGNK)})
        datenum = (data.date - pd.Timestamp("1970-01-01")).dt.days
        seas = dl.onebasis(datenum, "ns", df=DFSEAS * data.year.nunique())
        X = dl.design_matrix(data, ("cb", cb), ("seas", seas), intercept=False)
        model = dl.fit_glm("death ~ " + " + ".join(X.columns) + " + C(dow)",
                           data.join(X), family="quasipoisson")
        red = dl.crossreduce(cb, model, cen=float(tmean.mean()), name="cb")
        coefs.append(np.ravel(red.coef))
        vcovs.append(np.asarray(red.vcov, float))
        keep[city] = (cb, tmean, knots, data.death.to_numpy(float))
    t1 = time.perf_counter()

    # ---- 02.secondstage.R: meta-regression, BLUPs, minimum mortality -------
    Y, S = np.vstack(coefs), np.stack(vcovs)
    avgtmean = np.array([keep[c][1].mean() for c in cities.city])
    rangetmean = np.array([np.ptp(keep[c][1]) for c in cities.city])
    mv = mixmeta(Y, S, np.column_stack([np.ones(m), avgtmean, rangetmean]),
                 method="reml", bscov="unstr")
    bl = mv.blup(se=True)
    blups, blupvcov = np.asarray(bl["blup"], float), [np.asarray(v, float) for v in bl["vcov"]]
    t2 = time.perf_counter()

    minperc, mintemp = [], []
    for i, city in enumerate(cities.city):
        _, tmean, knots, _ = keep[city]
        predvar = dl.quantile7(tmean, np.arange(1, 100) / 100)
        bvar = dl.onebasis(predvar, VARFUN, knots=knots, degree=VARDEGREE,
                           **{"Boundary.knots": [tmean.min(), tmean.max()]})
        p = int(np.argmin(bvar.matrix @ blups[i])) + 1
        minperc.append(p)
        mintemp.append(float(dl.quantile7(tmean, [p / 100])[0]))
    minperc, mintemp = np.array(minperc), np.array(mintemp)

    # ---- 03.attr.R: attributable deaths, forward, eCIs from R's draws -------
    Z = np.asarray(ref["normals"], dtype=float)            # (region, comp, nsim, k)
    nsim = Z.shape[2]
    totdeath = np.array([keep[c][3].sum() for c in cities.city])
    an = np.zeros((m, 3))
    sims = {ch: np.zeros((m, 3, nsim)) for ch in ("published", "converged")}
    for i, city in enumerate(cities.city):
        cb, tmean, _, cases = keep[city]
        cen = mintemp[i]
        E = {ch: np.asarray(ref["attributable"][ch]["eigen_vectors"][i], float)
             for ch in sims}
        for j, rng in enumerate([None, (-100, cen), (cen, 100)]):
            kw = dict(coef=blups[i], vcov=blupvcov[i], type="an", dir="forw",
                      cen=cen, range=rng, name="cb")
            an[i, j] = dl.attrdl(tmean, cb, cases, **kw)
            # both chains' replays in one call: same deviates, R's eigenvectors
            # of that chain's BLUP covariance (an eigenvector's sign is arbitrary)
            cs = [dl.simulate_coef(blups[i], blupvcov[i], normals=Z[i, j], eigvecs=E[ch])
                  for ch in sims]
            out = dl.attrdl(tmean, cb, cases, sim=True, coefsim=np.hstack(cs), **kw)
            for c, ch in enumerate(sims):
                sims[ch][i, j] = out[c * nsim:(c + 1) * nsim]
    t3 = time.perf_counter()

    # ---- 04.tables.R ---------------------------------------------------------
    def table(sim):
        def q(a, p):
            return np.quantile(a, p, axis=-1)        # R's quantile type 7
        tot = sim.sum(axis=0)
        return dict(an=an, anlow=q(sim, 0.025), anhigh=q(sim, 0.975),
                    antot=an.sum(axis=0), antotlow=q(tot, 0.025), antothigh=q(tot, 0.975),
                    af=an / totdeath[:, None] * 100, aflow=q(sim, 0.025) / totdeath[:, None] * 100,
                    afhigh=q(sim, 0.975) / totdeath[:, None] * 100,
                    aftot=an.sum(axis=0) / totdeath.sum() * 100,
                    aftotlow=q(tot, 0.025) / totdeath.sum() * 100,
                    aftothigh=q(tot, 0.975) / totdeath.sum() * 100)
    tabs = {ch: table(sims[ch]) for ch in sims}

    # ---- compare with R ------------------------------------------------------
    chk = Checker()
    chk("settings", np.r_[np.ravel([keep[c][2] for c in cities.city]), avgtmean, rangetmean,
                          totdeath],
        np.r_[np.ravel(s["knots"]), s["avgtmean"], s["rangetmean"], s["totdeath"]],
        rtol=TOL["exact"])

    for ch, P in (("published", "p_"), ("converged", "c_")):
        fs = ref["first_stage"][ch]
        lab = "(a)" if ch == "published" else "(b)"
        Rc, Rv = np.asarray(fs["coef"], float), np.asarray(fs["vcov"], float)
        if ch == "published":
            chk(f"{lab} 1 reduced coef", Y, Rc, atol=TOL["p_coef"])
        else:
            chk.block(f"{lab} 1 reduced coef", Y, Rc, TOL["c_coef"])
        for i in range(m):
            chk.block(f"{lab} 1 reduced vcov", S[i], Rv[i], TOL[P + "vcov"])

        ss = ref["second_stage"][f"{ch}_{'mvmeta' if ch == 'published' else 'mixmeta'}"]
        meta = TOL[P + "meta"]
        if ch == "converged":
            chk(f"{lab} 2 meta-regression logLik", mv.loglik, ss["logLik"], atol=TOL["c_loglik"])
        chk.block(f"{lab} 2 meta-regression coef", mv.coef, ss["coef"], meta)
        chk.block(f"{lab} 2 meta-regression vcov", mv.vcov, ss["vcov"], meta)
        chk.block(f"{lab} 2 Psi", mv.Psi, ss["Psi"], meta)
        chk.block(f"{lab} 2 BLUPs", blups, ss["blup"], TOL[P + "blup"])
        for i in range(m):
            chk.block(f"{lab} 2 BLUP vcov", blupvcov[i], ss["blupvcov"][i], TOL[P + "blup"])

        at = ref["attributable"][ch]
        chk(f"{lab} 3 MMT percentile", minperc, at["minperc"], atol=0)
        chk(f"{lab} 3 MMT percentile", np.median(minperc), at["minperccountry"], atol=0)
        chk(f"{lab} 3 MMT (C)", mintemp, at["mintemp"], atol=TOL[P + "mmt"])
        t = tabs[ch]
        for key in ("an", "antot"):
            chk(f"{lab} 4 AN", t[key], at[key], rtol=TOL[P + "an"])
        for key, pt in (("anlow", "an"), ("anhigh", "an"), ("antotlow", "antot"),
                        ("antothigh", "antot")):
            chk(f"{lab} 4 AN eCI (R draws)", t[key], at[key], rtol=TOL[P + "an"], scale=at[pt])
        for key in ("af", "aftot"):
            chk(f"{lab} 4 AF %", t[key], at[key], rtol=TOL[P + "af"])
        for key, pt in (("aflow", "af"), ("afhigh", "af"), ("aftotlow", "aftot"),
                        ("aftothigh", "aftot")):
            chk(f"{lab} 4 AF % eCI (R draws)", t[key], at[key], rtol=TOL[P + "af"], scale=at[pt])

    # (c) the second stage against R's mixmeta run to its optimum
    ss = ref["second_stage"]["converged_mixmeta_reltol14"]
    chk("(c) 2 meta-regression logLik", mv.loglik, ss["logLik"], atol=TOL["t_loglik"])
    chk.block("(c) 2 meta-regression coef", mv.coef, ss["coef"], TOL["t_meta"])
    chk.block("(c) 2 meta-regression vcov", mv.vcov, ss["vcov"], TOL["t_meta"])
    chk.block("(c) 2 Psi", mv.Psi, ss["Psi"], TOL["t_meta"])
    chk.block("(c) 2 BLUPs", blups, ss["blup"], TOL["t_meta"])

    # ---- report ---------------------------------------------------------------
    v = ref["versions"]
    tp, tc = ref["attributable"]["published"], ref["attributable"]["converged"]
    print(f"Gasparrini et al. 2015, England & Wales {d.year.min()}-{d.year.max()}: "
          f"{int(totdeath.sum()):,} deaths, {m} regions")
    print(f"R reference: {v['R']}, dlnm {v['dlnm']}, mvmeta {v['mvmeta']}, "
          f"mixmeta {v['mixmeta']}, authors' commit {v['authors_commit'][:7]}, "
          f"seed {s['seed']}, nsim {s['nsim']}\n")
    print(f"  {'':22s} {'python':>22s} {'R published':>22s} {'R converged':>22s}")
    print(f"  {'MMT percentile':22s} {np.median(minperc):22.1f} {tp['minperccountry']:22.1f} "
          f"{tc['minperccountry']:22.1f}")
    for j, nm in enumerate(("AF total %", "AF cold %", "AF heat %")):
        def f(t):
            return f"{t['aftot'][j]:.3f} ({t['aftotlow'][j]:.3f}-{t['aftothigh'][j]:.3f})"
        print(f"  {nm:22s} {f(tabs['converged']):>22s} {f(tp):>22s} {f(tc):>22s}")
    print("  (python's interval uses the converged chain's eigenvectors; the published"
          " chain's\n   BLUP covariances give R's eigen() different signs, hence other draws)\n")

    rr = ref["second_stage"]
    mvc, mxc = rr["converged_mvmeta"], rr["converged_mixmeta"]
    print("  R against R, same converged inputs: mvmeta vs mixmeta max abs difference "
          f"coef {np.abs(np.subtract(mvc['coef'], mxc['coef'])).max():.1e}, "
          f"Psi {np.abs(np.subtract(mvc['Psi'], mxc['Psi'])).max():.1e}, "
          f"BLUPs {np.abs(np.subtract(mvc['blup'], mxc['blup'])).max():.1e}; "
          f"logLik constant {rr['loglik_constant_mvmeta_minus_mixmeta']:.4f}")
    print(f"  logLik: python {mv.loglik:.10f}, R mixmeta {mxc['logLik']:.10f} (default control), "
          f"{rr['converged_mixmeta_reltol14']['logLik']:.10f} (reltol 1e-14)")
    print("  (a) authors' code as published (glm epsilon 1e-8, mvmeta); (b) glm epsilon 1e-14,\n"
          "  mixmeta at its default control; (c) as (b), mixmeta with reltol 1e-14\n")
    print(f"  {'compared with R':32s} {'n':>6s} {'max abs diff':>13s} {'max rel diff':>13s}  fail")
    for g, (n, mx, rel, bad) in chk.groups.items():
        print(f"  {g:32s} {n:6d} {mx:13.2e} {rel:13.2e}  {bad:4d}")
    print(f"\n  runtime {time.perf_counter() - t0:.1f} s (first stage {t1 - t0:.1f} s, "
          f"meta-regression {t2 - t1:.1f} s, attributable {t3 - t2:.1f} s)")
    print(f"lancet2015_vs_r: {chk.n} comparisons with R, max abs difference {chk.maxabs:.2e}, "
          f"{chk.fail} failures")
    return 1 if chk.fail else 0


if __name__ == "__main__":
    sys.exit(main())
