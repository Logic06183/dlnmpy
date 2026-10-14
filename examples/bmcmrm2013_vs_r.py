"""Stage-by-stage comparison of examples/bmcmrm_2013.py with the authors' R code.

    Gasparrini A, Armstrong B. Reducing and meta-analysing estimates from
    distributed lag non-linear models. BMC Medical Research Methodology
    2013;13:1. doi:10.1186/1471-2288-13-1

Recomputes the analysis of
https://github.com/gasparrini/2013_gasparrini_BMCmrm_Rcodedata (commit
d542c9314901b00d89dae6d5000a599b7700d6da) on the 10 regions of England and
Wales with dlnmpy and mixmetapy, as examples/bmcmrm_2013.py does, and checks
every intermediate against tests/reference/bmcmrm2013.json, written by
tools/reproduce_bmcmrm2013.R from the authors' scripts: the five reductions
in each region, QAIC, the five pooled fits with their Q tests and I2, the
fixed-effects fit, the three meta-regressions on latitude with their Wald
tests and predictions, and every pooled curve. examples/bmcmrm_2013.py is
left as it is.

The data (regEngWales.csv, 7 MB) are downloaded once into examples/ from the
pinned commit and checked against a hardcoded SHA-256.

Three R references, |python - R| <= atol + rtol * |R| (values in TOL below):

(a) "published": the authors' code as written. glm() stops at its default
    epsilon = 1e-8 and reports a dispersion from the working weights of the
    previous IRLS iteration; the second stage is mvmeta(), whose REML
    log-likelihood differs from mixmeta's by a constant that depends only on
    the dimensions (added back before comparing).
(b) "converged": 03.firststage.R rerun in R with glm.control(epsilon =
    1e-14), the second stage with mixmeta() at its default control.
(c) "reltol14": as (b), with mixmeta()'s optimiser tolerance at 1e-14, so
    R sits at the REML optimum (second stage only).

The pooled fits estimate 4 or 5 coefficients and 10 or 15 parameters of Psi
from 10 regions, and the meta-regressions twice the coefficients. The REML
surface is flat, so a log-likelihood that agrees to 1e-8 leaves the
parameters agreeing only to ~1e-5, depending on where each optimiser stops.
R's default stopping rule (relative change below 1.5e-8) is the larger part
of that: against (c) the agreement is tighter. Tolerances on Psi, on the
pooled coefficients and on what is predicted from them reflect that, and
are relative to the largest entry of each matrix, since Psi has entries
near zero.
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
from scipy.stats import chi2

import dlnmpy as dl
from mixmetapy import mixmeta  # pip install dlnmpy[twostage]

SHA = "d542c9314901b00d89dae6d5000a599b7700d6da"
DATA_URL = ("https://raw.githubusercontent.com/gasparrini/"
            f"2013_gasparrini_BMCmrm_Rcodedata/{SHA}/regEngWales.csv")
DATA_SHA256 = "6f7b1471e641ed3b403b4d912388e57c1007fb91c3329b525fb1f689e31f29c4"
CACHE = Path(__file__).with_name("regEngWales.csv")
REFERENCE = Path(__file__).resolve().parents[1] / "tests" / "reference" / "bmcmrm2013.json"
LAT = [54.84815, 53.58832, 53.72352, 52.85539, 52.53304,
       52.03734, 51.50583, 51.24213, 51.05361, 52.02615]
KEYS = ("all", "all2", "all3", "hot", "cold")

# |python - R| <= atol + rtol * |R|; "block" tolerances are relative to the
# largest entry of the matrix. Reasons in the module docstring.
TOL = dict(
    exact=1e-12,       # rtol, knots, boundaries, grids, latitude quartiles
    # (a) against the authors' code as published
    p_red=1e-7,        # block rtol, reduced coefficients (glm epsilon 1e-8)
    p_redv=2e-5,       # block rtol, reduced covariance: R's dispersion comes from
                       # the working weights of the previous IRLS iteration
    p_qaic=5e-7,       # rtol, QAIC (~40,000 per region; its phi term inherits that)
    p_meta=5e-3,       # block rtol, pooled and meta-regression coef, vcov, Psi, and
                       # Wald p-values: mvmeta's default stopping rule on a flat
                       # surface, on inputs that differ by the above
    p_loglik=1e-4,     # atol, log-likelihoods (mvmeta's constant added back); the
                       # log-determinants of 50 within-region covariances inherit
                       # their 4e-6 relative difference
    p_q=5e-4,          # rtol, Cochran Q and p-values; I2 to 100 * this in points
    p_pred=1e-4,       # rtol, pooled RR and CI limits, latitude predictions
    # (b) against the converged first stage and mixmeta at its default control
    c_red=1e-10,       # block rtol, reduced coefficients
    c_redv=1e-8,       # block rtol, reduced covariance
    c_qaic=1e-10,      # rtol, QAIC
    c_meta=5e-3,       # block rtol, pooled and meta-regression parameters: mixmeta's
                       # default stopping rule (the hot fit stops 5e-6 short in
                       # logLik and 2e-3 of max|Psi| away; see (c))
    c_loglik=1e-5,     # atol, log-likelihoods, same reason
    c_q=1e-6,          # rtol, Cochran Q (a fixed-effects quantity, no optimiser)
    c_pred=5e-5,       # rtol, pooled RR and CI limits, latitude predictions
    # (c) second stage only, against mixmeta at reltol 1e-14: both at the optimum
    t_meta=1e-5,       # block rtol, coef, vcov, Psi
    t_loglik=1e-8,     # atol, log-likelihoods
    fixed=1e-9)        # block rtol, fixed-effects fit (closed form, no optimiser)


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
    return pd.read_csv(io.BytesIO(raw), index_col=0)


class Checker:
    """Collects |python - R| for every compared scalar."""

    def __init__(self):
        self.n = 0
        self.fail = 0
        self.maxabs = 0.0
        self.groups = {}

    def __call__(self, group, py, r, atol=0.0, rtol=0.0, scale=None):
        """``scale`` replaces |R| in the relative bound (the largest entry of a
        matrix, for matrices with entries near zero)."""
        py = np.asarray(py, dtype=float).ravel()
        r = np.asarray([np.nan if v is None else v for v in np.ravel(np.asarray(r, dtype=object))],
                       dtype=float)
        if py.shape != r.shape:
            raise ValueError(f"{group}: shape {py.shape} vs R {r.shape}")
        nan_ok = np.isnan(py) == np.isnan(r)
        both = ~np.isnan(py) & ~np.isnan(r)
        d = np.abs(py[both] - r[both])
        ref = np.abs(r[both]) if scale is None else np.full(d.shape, float(scale))
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
        r = np.asarray(r, dtype=float)
        self(group, py, r, rtol=rtol, scale=float(np.abs(r).max()) or 1.0)


def qtest(fit):
    q = fit.qtest()
    return dict(Q=q["Q"], df=q["df"], pvalue=q["pvalue"],
                I2=float((q["Q"][0] - q["df"][0]) / q["Q"][0] * 100))


def wald(fit):
    """R's fwald(): joint test of the latitude coefficients."""
    k = fit.k
    b = fit.coef_vec[-k:]
    V = fit.vcov[-k:, -k:]
    return float(chi2.sf(b @ np.linalg.solve(V, b), k))


def main() -> int:
    t0 = time.perf_counter()
    ref = json.loads(REFERENCE.read_text())
    s = ref["settings"]
    d = load()
    regions = list(dict.fromkeys(d.regnames))
    if regions != s["regions"]:
        raise RuntimeError("region order differs from R's")
    dlist = {r: d[d.regnames == r] for r in regions}

    # ---- 03.firststage.R ----------------------------------------------------
    ranges = np.array([[x.tmean.min(), x.tmean.max()] for x in dlist.values()])
    bound = ranges.mean(axis=0)
    varknots = dl.equalknots(bound, fun="bs", degree=2, df=4)
    lagknots = dl.logknots(21, df=5, intercept=True)
    argvar = {"fun": "bs", "degree": 2, "boundary_knots": bound, "knots": varknots}
    arglag = {"fun": "ns", "knots": lagknots}
    arglag_const = {"fun": "strata", "df": 1}

    red = {k: ([], []) for k in KEYS}
    qaic = {m: [] for m in ("m1", "m2", "m3")}
    for region in regions:
        sub = dlist[region]
        tmean = sub.tmean.to_numpy(float)
        cb = dl.crossbasis(tmean, lag=(0, 21), argvar=argvar, arglag=arglag)
        cb2 = dl.crossbasis(tmean, lag=(0, 3), argvar=argvar, arglag=arglag_const)
        cb3 = dl.crossbasis(tmean, lag=(0, 21), argvar=argvar, arglag=arglag_const)
        cb2.matrix[:21, :] = np.nan               # R: cb2[0:21,] <- NA
        seas = dl.onebasis(sub.time.to_numpy(float), "ns", df=10 * 14)
        fits = {}
        for name, basis in (("cb", cb), ("cb2", cb2), ("cb3", cb3)):
            X = dl.design_matrix(sub, (name, basis), ("seas", seas), intercept=False)
            fits[name] = dl.fit_glm("death ~ " + " + ".join(X.columns) + " + C(dow)",
                                    sub.join(X), family="quasipoisson")
        for m, name in zip(qaic, ("cb", "cb2", "cb3")):
            qaic[m].append(dl.uncertainty.qaic(fits[name]))
        for key, basis, name, kw in (
                ("all", cb, "cb", {}), ("all2", cb2, "cb2", {}), ("all3", cb3, "cb3", {}),
                ("hot", cb, "cb", {"type": "var", "value": 22}),
                ("cold", cb, "cb", {"type": "var", "value": 0})):
            cr = dl.crossreduce(basis, fits[name], cen=17, name=name, **kw)
            red[key][0].append(np.ravel(cr.coef))
            red[key][1].append(np.asarray(cr.vcov, float))
        realised_arglag = dict(cb.arglag)
    Y = {k: np.vstack(red[k][0]) for k in KEYS}
    S = {k: np.stack(red[k][1]) for k in KEYS}
    t1 = time.perf_counter()

    # ---- 04.secondstage.R, 05.mainresults.R, 06.metareg.R --------------------
    pooled = {k: mixmeta(Y[k], S[k], method="reml") for k in KEYS}
    fixed_hot = mixmeta(Y["hot"], S["hot"], method="fixed")
    latX = np.column_stack([np.ones(len(regions)), LAT])
    metareg = {k: mixmeta(Y[k], S[k], latX, method="reml") for k in ("all", "hot", "cold")}
    t2 = time.perf_counter()

    a = dict(argvar)
    bvar = dl.onebasis(np.arange(bound[0], bound[1] + 1e-9, 0.1), a.pop("fun"), **a)
    al = dict(realised_arglag)
    lagfun = al.pop("fun")
    al.pop("cen", None)
    xlag = np.arange(0, 211) / 10
    blag = dl.onebasis(xlag, lagfun, **al)

    def cpvar(coef, vcov, **kw):
        return dl.crosspred(bvar, coef=np.ravel(coef), vcov=vcov, model_link="log", **kw)

    def cplag(coef, vcov, **kw):
        return dl.crosspred(blag, coef=np.ravel(coef), vcov=vcov, model_link="log", **kw)

    grid = dict(by=0.1, from_=bound[0], to=bound[1], cen=17)
    curves = {
        "cpall": cpvar(pooled["all"].coef, pooled["all"].vcov, **grid),
        "cpall2": cpvar(pooled["all2"].coef, pooled["all2"].vcov, **grid),
        "cpall3": cpvar(pooled["all3"].coef, pooled["all3"].vcov, **grid),
        "cphot": cplag(pooled["hot"].coef, pooled["hot"].vcov, at=xlag, cen=17),
        "cpcold": cplag(pooled["cold"].coef, pooled["cold"].vcov, at=xlag, cen=17),
        "cphot2": cplag(fixed_hot.coef, fixed_hot.vcov, at=xlag)}
    val = np.round(dl.quantile7(LAT, [0.25, 0.75]), 1)
    latpred = {}
    for k in ("all", "hot", "cold"):
        p = metareg[k].predict(np.column_stack([np.ones(2), val]), se=True)
        latpred[k] = [(p["fit"][j], p["vcov"][j]) for j in range(2)]
        for j, q in enumerate(("25", "75")):
            fit, V = latpred[k][j]
            curves[f"cp{k}lat{q}"] = (cpvar(fit, V, by=0.2, cen=17) if k == "all"
                                      else cplag(fit, V, by=0.2))
    cp = curves["cpall"]
    mmt = float(cp.predvar[np.argmin(cp.allRRfit)])
    i22, i0 = int(np.argmin(np.abs(cp.predvar - 22))), int(np.argmin(np.abs(cp.predvar)))
    rr22 = [cp.allRRfit[i22], cp.allRRlow[i22], cp.allRRhigh[i22]]
    rr0 = [cp.allRRfit[i0], cp.allRRlow[i0], cp.allRRhigh[i0]]

    # ---- compare with R -------------------------------------------------------
    chk = Checker()
    chk("settings", np.r_[bound, varknots, lagknots, val],
        np.r_[s["bound"], s["varknots"], s["lagknots"],
              ref["published"]["second_stage"]["latpred"]["val"]], rtol=TOL["exact"])
    const = ref["loglik_constant_mvmeta_minus_mixmeta"]

    def cmp_fit(lab, py, r, tol_meta, tol_ll, tol_q, ll_const=0.0, q=True):
        chk.block(f"{lab} coef", py.coef, r["coef"], tol_meta)
        chk.block(f"{lab} vcov", py.vcov, r["vcov"], tol_meta)
        if py.method != "fixed":
            chk.block(f"{lab} Psi", py.Psi, r["Psi"], tol_meta)
        chk(f"{lab} logLik", py.loglik + ll_const, r["logLik"], atol=tol_ll)
        if q:
            qt = qtest(py)
            chk(f"{lab} Q test", np.r_[qt["Q"], qt["df"]], np.r_[r["qtest"]["Q"], r["qtest"]["df"]],
                rtol=tol_q)
            chk(f"{lab} Q test", qt["pvalue"], r["qtest"]["pvalue"], rtol=tol_q, atol=1e-15)
            chk(f"{lab} Q test", qt["I2"], r["qtest"]["I2"], atol=tol_q * 100)

    for ch, P, lab in (("published", "p_", "(a)"), ("converged", "c_", "(b)")):
        fs, ss = ref[ch]["first_stage"], ref[ch]["second_stage"]
        for k in KEYS:
            chk.block(f"{lab} 1 reduced coef", Y[k], fs["reductions"][k]["coef"], TOL[P + "red"])
            for i in range(len(regions)):
                chk.block(f"{lab} 1 reduced vcov", S[k][i], fs["reductions"][k]["vcov"][i],
                          TOL[P + "redv"])
        for m in qaic:
            chk(f"{lab} 1 QAIC per region", qaic[m], fs["qaic"][m], rtol=TOL[P + "qaic"])
        chk(f"{lab} 1 QAIC totals", [sum(qaic[m]) for m in qaic], fs["qaic_total"],
            rtol=TOL[P + "qaic"])

        for k in KEYS:
            cmp_fit(f"{lab} 2 pooled", pooled[k], ss["pooled"][k], TOL[P + "meta"],
                    TOL[P + "loglik"], TOL[P + "q"],
                    ll_const=const["pooled"][k] if ch == "published" else 0.0)
        cmp_fit(f"{lab} 2 fixed effects", fixed_hot, ss["fixed_hot"], TOL["fixed"]
                if ch == "converged" else TOL["p_meta"], TOL[P + "loglik"], TOL[P + "q"],
                ll_const=const["fixed_hot"] if ch == "published" else 0.0)
        for k in ("all", "hot", "cold"):
            r = ss["metareg"][k]
            cmp_fit(f"{lab} 2 meta-regression", metareg[k], r, TOL[P + "meta"],
                    TOL[P + "loglik"], TOL[P + "q"],
                    ll_const=const["metareg"][k] if ch == "published" else 0.0)
        chk(f"{lab} 2 Wald test, latitude", [wald(metareg[k]) for k in ("all", "hot", "cold")],
            ss["wald_lat"], rtol=TOL[P + "meta"])
        for k in ("all", "hot", "cold"):
            for j in range(2):
                fit, V = latpred[k][j]
                rp = ss["latpred"][k][j]
                chk.block(f"{lab} 2 latitude predictions", fit, rp["fit"], TOL[P + "pred"])
                chk.block(f"{lab} 2 latitude predictions", V, rp["vcov"], TOL[P + "meta"])
        for name, c in curves.items():
            rc = ss["curves"][name]
            chk(f"{lab} 3 curves: grid", c.predvar, rc["predvar"], rtol=TOL["exact"])
            chk(f"{lab} 3 curves: RR", c.allRRfit, rc["allRRfit"], rtol=TOL[P + "pred"])
            chk(f"{lab} 3 curves: RR limits", np.r_[c.allRRlow, c.allRRhigh],
                np.r_[rc["allRRlow"], rc["allRRhigh"]], rtol=TOL[P + "pred"])
        chk(f"{lab} 3 MMT, RR at 22 and 0 C", mmt, ss["mmt"], atol=1e-12)
        chk(f"{lab} 3 MMT, RR at 22 and 0 C", np.r_[rr22, rr0], np.r_[ss["rr22"], ss["rr0"]],
            rtol=TOL[P + "pred"])

    tt = ref["reltol14"]
    for k in KEYS:
        cmp_fit("(c) 2 pooled", pooled[k], tt["pooled"][k], TOL["t_meta"], TOL["t_loglik"], 0, q=False)
    for k in ("all", "hot", "cold"):
        cmp_fit("(c) 2 meta-regression", metareg[k], tt["metareg"][k], TOL["t_meta"],
                TOL["t_loglik"], 0, q=False)

    # ---- report ----------------------------------------------------------------
    v = ref["versions"]
    pub = ref["published"]["second_stage"]
    print(f"Gasparrini & Armstrong 2013, England & Wales {d.year.min()}-{d.year.max()}, "
          f"{len(regions)} regions")
    print(f"R reference: {v['R']}, dlnm {v['dlnm']}, mvmeta {v['mvmeta']}, "
          f"mixmeta {v['mixmeta']}, authors' commit {v['authors_commit'][:7]}\n")
    def f3(x):
        return f"{x[0]:.3f} ({x[1]:.3f}-{x[2]:.3f})"
    print(f"  {'':26s} {'python':>22s} {'R published':>22s}")
    print(f"  {'MMT (C)':26s} {mmt:22.1f} {pub['mmt']:22.1f}")
    print(f"  {'pooled RR at 22 C':26s} {f3(rr22):>22s} {f3(pub['rr22']):>22s}")
    print(f"  {'pooled RR at 0 C':26s} {f3(rr0):>22s} {f3(pub['rr0']):>22s}")
    i2 = "/".join(f"{qtest(pooled[k])['I2']:.1f}" for k in ("all", "hot", "cold"))
    i2r = "/".join(f"{pub['pooled'][k]['qtest']['I2']:.1f}" for k in ("all", "hot", "cold"))
    print(f"  {'I2 overall / 22 C / 0 C':26s} {i2:>22s} {i2r:>22s}")
    for j, m in enumerate(qaic):
        print(f"  {'QAIC ' + m:26s} {sum(qaic[m]):22.1f} "
              f"{ref['published']['first_stage']['qaic_total'][j]:22.1f}")
    print("\n  logLik of the pooled fits: python / R mixmeta default / R mixmeta reltol 1e-14")
    for k in KEYS:
        print(f"    {k:5s} {pooled[k].loglik:.10f} / "
              f"{ref['converged']['second_stage']['pooled'][k]['logLik']:.10f} / "
              f"{tt['pooled'][k]['logLik']:.10f}")
    for k in ("all", "hot", "cold"):
        print(f"    {k + '~lat':9s} {metareg[k].loglik:.10f} / "
              f"{ref['converged']['second_stage']['metareg'][k]['logLik']:.10f} / "
              f"{tt['metareg'][k]['logLik']:.10f}")
    print("\n  (a) authors' code as published (glm epsilon 1e-8, mvmeta); (b) glm epsilon 1e-14,\n"
          "  mixmeta at its default control; (c) as (b), mixmeta with reltol 1e-14\n")
    print(f"  {'compared with R':34s} {'n':>6s} {'max abs diff':>13s} {'max rel diff':>13s}  fail")
    for g, (n, mx, rel, bad) in chk.groups.items():
        print(f"  {g:34s} {n:6d} {mx:13.2e} {rel:13.2e}  {bad:4d}")
    print(f"\n  runtime {time.perf_counter() - t0:.1f} s (first stage {t1 - t0:.1f} s, "
          f"second stage {t2 - t1:.1f} s)")
    print(f"bmcmrm2013_vs_r: {chk.n} comparisons with R, max abs difference {chk.maxabs:.2e}, "
          f"{chk.fail} failures")
    return 1 if chk.fail else 0


if __name__ == "__main__":
    sys.exit(main())
