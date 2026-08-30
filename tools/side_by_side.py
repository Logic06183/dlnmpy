"""Python half of the side-by-side validation (see side_by_side.R).

Each analysis is run end to end in Python, fitting the model with
statsmodels, and compared with R's results (coefficients, dispersion, and
every crosspred quantity). Run from the repository root:

    python tools/side_by_side.py
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

import dlnmpy as dl

HERE = Path(__file__).resolve().parents[1] / "tests" / "side_by_side"
R = json.load(open(HERE / "r_results.json"))
rows = []


def report(section, name, a, b, tol):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if a.shape != b.shape:
        rows.append((section, name, f"SHAPE {a.shape} vs {b.shape}", "FAIL"))
        return
    d = float(np.nanmax(np.abs(a - b))) if a.size else 0.0
    # scale the tolerance by the magnitude of the R value. A flat absolute
    # bound asks for ~10 significant figures from quantities as large as the
    # nested vcov (entries up to 35, built from a finite-difference Hessian),
    # which failed on arm64/Accelerate and passed on x86_64 -- the comparison
    # was platform-dependent rather than wrong. Read as a relative tolerance
    # this is *stricter* for the many small-magnitude quantities.
    ok = bool(np.all(np.abs(a - b) <= 1e-10 + tol * np.abs(b))) if a.size else True
    rows.append((section, name, f"{d:.2e}", "ok" if ok else "FAIL"))


def compare_pred(section, p, r, tol=1e-8):
    report(section, "predvar", p.predvar, r["predvar"], 1e-12)
    report(section, "coef", p.coef, r["coef"], tol)
    report(section, "vcov", p.vcov, r["vcov"], tol)
    report(section, "matfit", p.matfit, r["matfit"], tol)
    report(section, "matse", p.matse, r["matse"], tol)
    report(section, "allfit", p.allfit, r["allfit"], tol)
    report(section, "allse", p.allse, r["allse"], tol)
    report(section, "alllow", p.alllow, r["alllow"], tol)
    report(section, "allhigh", p.allhigh, r["allhigh"], tol)
    if r.get("cumfit") is not None:
        report(section, "cumfit", p.cumfit, r["cumfit"], tol)
    assert p.ci_level == r["ci.level"], "ci level mismatch"
    assert p.model_link == r["link"], f"link {p.model_link} vs {r['link']}"


# A. drug trial (Gaussian OLS, exposure-history matrix) -----------------------
drug = dl.datasets.drug()
Qdrug = pd.read_csv(HERE / "Qdrug.csv").to_numpy(float)
cbdrug = dl.crossbasis(Qdrug, lag=27, argvar={"fun": "lin"}, arglag={"fun": "ns", "knots": [9, 18]})
X = drug.join(cbdrug.to_dataframe("cb"))
mdrug = dl.fit_glm("out ~ " + " + ".join(cbdrug.to_dataframe("cb").columns) + " + C(sex)", X, family="gaussian")
pdrug = dl.crosspred(cbdrug, mdrug, at=np.arange(21) * 5, ci_level=0.90, name="cb")
# (the sex factor is parameterised with a different reference level in patsy, so compare the cross-basis part)
report("drug", "cb coefs", mdrug.params.filter(like="cb_").to_numpy(), np.array(R["drug"]["coef_all"])[1:1 + cbdrug.ncol], 1e-8)
report("drug", "residual sd", np.sqrt(mdrug.scale), R["drug"]["sigma"], 1e-8)
compare_pred("drug", pdrug, R["drug"]["pred"])

# B. nested case-control (conditional logistic regression) --------------------
nested = dl.datasets.nested()
Qnest = pd.read_csv(HERE / "Qnest.csv").to_numpy(float)
# rebuild Qnest in Python too and check it against R's
Qpy = np.vstack([dl.exphist(np.repeat(np.r_[0, 0, 0, row[4:14]], 5), times=[row[3]], lag=(3, 40))
                 for row in nested[["id", "case", "riskset", "age"] + [f"exp{a}" for a in range(15, 65, 5)]].to_numpy(float)])
report("nested", "exphist matrix", Qpy, Qnest, 1e-12)
cbnest = dl.crossbasis(Qnest, lag=(3, 40), argvar={"fun": "bs", "degree": 2, "df": 3},
                       arglag={"fun": "ns", "knots": [10, 30], "intercept": False})
# NB: statsmodels' default optimiser (BFGS) stops on this flat likelihood with
# coefficients 0.025 away from R's although the log-likelihood agrees to 1e-5;
# fit_clogit uses Newton-Raphson with a tight tolerance, which matches R.
mnest = dl.fit_clogit(nested["case"], cbnest.to_dataframe("cb"), groups=nested["riskset"])
report("nested", "clogit coefs", mnest.params.to_numpy(), R["nested"]["coef_all"], 1e-8)
pnest = dl.crosspred(cbnest, mnest, cen=0, at=np.arange(21) * 5, name="cb")
compare_pred("nested", pnest, R["nested"]["pred"])

# C. Chicago logistic regression, thr var + integer lag -----------------------
chi = dl.datasets.chicago_nmmaps()
chi["high"] = (chi["death"] > chi["death"].median()).astype(int)
cbC = dl.crossbasis(chi.temp, lag=7, argvar={"fun": "thr", "thr": [5, 25]}, arglag={"fun": "integer"})
nst = dl.onebasis(chi.time, "ns", df=7 * 14)
XC = chi.join(cbC.to_dataframe("cb")).join(nst.to_dataframe("nst"))
mC = dl.fit_glm("high ~ " + " + ".join(list(cbC.to_dataframe("cb").columns) + list(nst.to_dataframe("nst").columns)) + " + C(dow)", XC, family="binomial")
report("chicago_logit", "cb coefs", mC.params.filter(like="cb_").to_numpy(), np.array(R["chicago_logit"]["coef_all"])[1:1 + cbC.ncol], 1e-8)
pC = dl.crosspred(cbC, mC, at=[-20, -10, 0, 10, 20, 30, 33], cumul=True, name="cb")
compare_pred("chicago_logit", pC, R["chicago_logit"]["pred"])

# D. Chicago quasi-Poisson, poly var + strata lag, matrix 'at', 80% CI --------
cbD = dl.crossbasis(chi.temp, lag=10, argvar={"fun": "poly", "degree": 3}, arglag={"fun": "strata", "breaks": [1, 4]})
XD = chi.join(cbD.to_dataframe("cb")).join(nst.to_dataframe("nst"))
mD = dl.fit_glm("death ~ " + " + ".join(list(cbD.to_dataframe("cb").columns) + list(nst.to_dataframe("nst").columns)) + " + C(dow)", XD, family="quasipoisson")
report("chicago_poly", "dispersion", mD.scale, R["chicago_poly"]["dispersion"], 1e-8)
pD = dl.crosspred(cbD, mD, cen=20, by=2, ci_level=0.80, name="cb")
compare_pred("chicago_poly", pD, R["chicago_poly"]["pred"])
atmat = dl.exphist(chi.temp.to_numpy(), times=[200, 2000, 4000], lag=10)
report("chicago_poly", "exphist at-matrix", atmat, R["chicago_poly"]["atmat"], 1e-12)
pDm = dl.crosspred(cbD, mD, cen=20, at=atmat, name="cb")
rm = R["chicago_poly"]["pred_mat"]
report("chicago_poly", "matrix-at matfit", pDm.matfit, rm["matfit"], 1e-8)
report("chicago_poly", "matrix-at allfit", pDm.allfit, rm["allfit"], 1e-8)
report("chicago_poly", "matrix-at allse", pDm.allse, rm["allse"], 1e-8)
rD = dl.crossreduce(cbD, mD, cen=20, by=2, ci_level=0.8, name="cb")
report("chicago_poly", "reduce coef", rD.coef, R["chicago_poly"]["red"]["coef"], 1e-8)
report("chicago_poly", "reduce vcov", rD.vcov, R["chicago_poly"]["red"]["vcov"], 1e-8)
report("chicago_poly", "reduce se", rD.se, R["chicago_poly"]["red"]["se"], 1e-8)

# E. simulated multi-city data with known truth ---------------------------------
sim = pd.read_csv(HERE / "sim.csv")
truth = np.array(R["sim"]["truth_all"])
for k, st in enumerate(R["sim"]["stage1"], start=1):
    d = sim[sim.city == k].reset_index(drop=True)
    cb = dl.crossbasis(d.tmean, lag=10, argvar={"fun": "ns", "knots": [10, 18, 25]}, arglag={"fun": "ns", "df": 3})
    nt = dl.onebasis(d.time, "ns", df=8 * 3)
    Xs = d.join(cb.to_dataframe("cb")).join(nt.to_dataframe("nt"))
    m = dl.fit_glm("y ~ " + " + ".join(list(cb.to_dataframe("cb").columns) + list(nt.to_dataframe("nt").columns)), Xs, family="quasipoisson")
    p = dl.crosspred(cb, m, cen=18, at=np.arange(0, 31, 2), name="cb")
    r = dl.crossreduce(cb, m, cen=18, at=np.arange(0, 31, 2), name="cb")
    sec = f"sim city {k}"
    report(sec, "dispersion", m.scale, st["dispersion"], 1e-8)
    compare_pred(sec, p, st["pred"])
    report(sec, "reduce coef", r.coef, st["red_coef"], 1e-8)
    report(sec, "reduce vcov", r.vcov, st["red_vcov"], 1e-8)
    # sanity against the simulated truth: 95% CI should cover the true curve at most points
    cover = np.mean((p.alllow <= truth) & (truth <= p.allhigh))
    rows.append((sec, "truth inside 95% CI", f"{cover:.0%} of grid points", "ok" if cover >= 0.8 else "CHECK"))

# ------------------------------------------------------------------------------
w = max(len(r[1]) for r in rows)
print(f"{'analysis':16s} {'quantity':{w}s} {'max |Py - R|':>16s}  status")
for s, n, d, ok in rows:
    print(f"{s:16s} {n:{w}s} {d:>16s}  {ok}")
fails = [r for r in rows if r[3] == "FAIL"]
print(f"\n{len(rows)} comparisons, {len(fails)} failures")
raise SystemExit(1 if fails else 0)
