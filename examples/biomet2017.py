"""Reproduce the penalised DLNMs of Gasparrini et al. (2017) Biometrics against mgcv.

    Gasparrini A, Scheipl F, Armstrong B, Kenward MG. A penalized framework
    for distributed lag non-linear models. Biometrics 2017;73(3):938-948.
    doi:10.1111/biom.12645

The authors' code and data are public at
https://github.com/gasparrini/2017_gasparrini_Biomet_Rcodedata. This script
translates the penalised ("external" method) models of their two examples:

    ex1_gam1      London temperature, ps x ps cross-basis, default penalties
    ex1_gam2      as gam1 with the doubly varying lag penalties (eq. 8b, 7a)
    ex1_gam3      double-threshold exposure-response, penalised lag
    ex2_gam_ubre  uranium miners, cr x cr with a ridge lag penalty (eq. 7b),
                  Poisson, mgcv's default UBRE criterion, as published
    ex2_gam_reml  the same model with REML (not in the authors' code)

and compares every fit with the reference values written by
``tools/reproduce_biomet2017.R`` to ``tests/reference/biomet2017.json``.

Two kinds of comparison are made.

(a) Counted. dlnmpy fitted at mgcv's smoothing parameters (``sp=``): the
    penalty matrices and ranks from ``cbpen``, the REML score, coefficients,
    standard errors, the cross-basis covariance block, edf, scale and the
    overall cumulative curve from ``crosspred``. These are deterministic given
    the smoothing parameters and should agree tightly.

(b) Information. dlnmpy choosing its own smoothing parameters by REML. The
    optimisers differ (mgcv: Newton on analytic derivatives; dlnmpy: BFGS and
    Nelder-Mead on numerical ones), and on a flat criterion they can stop at
    different smoothing parameters with almost the same score. This is
    reported (score difference, sp ratio, largest difference in the predicted
    RR) and only counted as a failure when dlnmpy's optimum is *worse* than
    mgcv's by more than SCORE_TOL. UBRE has no counterpart in dlnmpy, so the
    published ex2 model has comparison (a) only.

Runtime (2 cores): about 2 minutes by default, (a) for all five fits and (b)
for the three London models; ``--quick`` runs (a) only, in about 1 minute.
Most of that is building the miners' cr cross-basis: dlnmpy evaluates ``cr``
one value at a time in Python, and Qx holds 3 million. ``--all`` adds (b) for
the miners REML model, which took 26 minutes (28 in all): the Nelder-Mead stage of
``fit_pglm`` exhausts its 2000 iterations there, because the criterion,
evaluated from warm-started P-IRLS, varies by ~5e-9 between calls at the same
smoothing parameters, more than its 1e-9 stopping tolerance on the score.
The reference JSON takes about 20 minutes to regenerate in R.

Data are downloaded once from the pinned commit and checked against
hard-coded SHA-256 digests; set DLNMPY_CACHE to choose where they are kept.
The script exits non-zero if any comparison fails. Its last line is
``biomet2017: N comparisons with R, max abs difference X, F failures``.

Deviations from the authors' code: the internal method (``s(Q, L, bs="cb")``),
the unpenalised GLMs, the Q-AIC knot search, the simulation study and the
plots are not reproduced. The authors' ``int=F`` is written ``intercept=False``
(dlnmpy has accepted ``int=`` as R's partial match since 0.8.4).
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
import urllib.request
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import dlnmpy as dl

COMMIT = "531c5c8bec3d477a32e5fd66e7f37f1e4fdead1d"
BASE = f"https://raw.githubusercontent.com/gasparrini/2017_gasparrini_Biomet_Rcodedata/{COMMIT}/"
DATA = {
    "london.csv": ("example1/london.csv",
                   "51b1753f307df245ee8678cc7fcbf9e76226e84498e055be44ed4769edceae64"),
    "uminers.csv": ("example2/uminers.csv",
                    "478e5a56db6aa0100024a44e739f4f6904dfecfbb0ff06b050a96dca79a4df7f"),
}
CACHE = Path(os.environ.get("DLNMPY_CACHE", Path.home() / ".cache" / "dlnmpy")) / f"biomet2017-{COMMIT[:7]}"
REFERENCE = Path(__file__).resolve().parent.parent / "tests" / "reference" / "biomet2017.json"

# ---- tolerances --------------------------------------------------------------
# (a) at matched smoothing parameters. mgcv stops its inner P-IRLS on a relative
# change in penalised deviance of 1e-7 (gam.control epsilon), so its
# coefficients are not exact either; these bounds leave room for that and
# nothing more. Differences are absolute, except coefficients and standard
# errors, which are scaled by max(1, |R value|): in the miners model the
# intercept, bs(age) and year columns are nearly collinear (condition number of
# the penalised information ~1e13), their standard errors are ~40, and dlnmpy
# itself moves them by ~7e-8 when sp changes by 1e-12 relative.
TOL = {
    "penalty": 1e-12,     # cbpen vs cbPen matrices (entries O(1) after rescaling)
    "rank": 0,
    "score": 1e-6,        # REML score (O(1e3)), Laplace criterion at fixed sp
    "coef": 1e-6,         # scaled
    "se": 1e-7,           # scaled
    "vcov": 1e-9,
    "edf": 1e-6,
    "scale": 1e-8,
    "allfit": 1e-6,       # overall cumulative log-RR
    "allse": 1e-7,
    "data": 1e-8,         # ex2 exposure histories (column sums), relative
}
# (b) dlnmpy's own selection fails only if its REML score is higher (worse)
# than mgcv's by more than this, in log-likelihood units. 1e-3 is far below
# anything that would matter for inference.
SCORE_TOL = 1e-3


# ---- data --------------------------------------------------------------------
def fetch(name: str) -> Path:
    rel, digest = DATA[name]
    path = CACHE / name
    if not path.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        print(f"downloading {BASE + rel}")
        with urllib.request.urlopen(BASE + rel, timeout=60) as r:   # noqa: S310
            raw = r.read()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise RuntimeError(f"{name}: SHA-256 mismatch, refusing to use the download")
        path.write_bytes(raw)
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise RuntimeError(f"{path}: SHA-256 mismatch; delete it and rerun")
    return path


# ---- bookkeeping -------------------------------------------------------------
class Tally:
    def __init__(self):
        self.n, self.fail, self.maxdiff, self.rows = 0, 0, 0.0, []

    def check(self, model: str, what: str, py, r, tol: float, mode: str = "abs"):
        """mode "abs": max |py - r|; "scaled": max |py - r| / max(1, |r|), for
        quantities whose size varies by orders of magnitude; "rel": relative to
        the largest |r| (data checks, not counted in the max abs difference)."""
        py, r = np.asarray(py, float), np.asarray(r, float)
        if py.shape != r.shape:
            raise ValueError(f"{model} {what}: shape {py.shape} vs R {r.shape}")
        raw = np.abs(py - r)
        absd = float(raw.max()) if raw.size else 0.0
        if mode == "abs":
            d = absd
        elif mode == "scaled":
            d = float(np.max(raw / np.maximum(1.0, np.abs(r)))) if raw.size else 0.0
        else:
            d = absd / max(float(np.max(np.abs(r))), 1e-300)
        ok = d <= tol
        self.n += 1
        self.fail += not ok
        if mode != "rel":
            self.maxdiff = max(self.maxdiff, absd)
        self.rows.append((model, f"{what}{' (rel)' if mode == 'rel' else ' (scaled)' if mode == 'scaled' else ''}",
                          py.size, d, tol, ok))
        return d


def r_mat(x) -> np.ndarray:
    return np.atleast_2d(np.asarray(x, float))


def dummies(s: pd.Series) -> np.ndarray:
    """Treatment contrasts with alphabetical levels, as R's factor()."""
    lev = sorted(s.unique())
    return np.column_stack([(s == lv).to_numpy(float) for lv in lev[1:]])


def fit_compare(tally: Tally, info: list, label: str, ref: dict, y, X, cb, cbidx, pen,
                family: str, cen, select: bool):
    """Fit at mgcv's sp (counted), and optionally with dlnmpy's own sp (info)."""
    names = [k for k in pen if k not in ("rank", "sp")]
    # penalties and ranks
    if names != list(np.atleast_1d(ref["penalty"]["names"])):
        raise ValueError(f"{label}: penalty names {names} vs R {ref['penalty']['names']}")
    for k in names:
        tally.check(label, f"penalty {k}", pen[k], r_mat(ref["penalty"]["mats"][k]), TOL["penalty"])
    tally.check(label, "ranks", [pen["rank"][k] for k in names], np.atleast_1d(ref["penalty"]["rank"]), TOL["rank"])

    p = X.shape[1]
    Slist = []
    for k in names:
        full = np.zeros((p, p))
        full[np.ix_(cbidx, cbidx)] = pen[k]
        Slist.append(full)
    # At fixed sp the coefficients do not depend on the criterion, so the
    # published UBRE model is refitted with method="reml" at mgcv's UBRE sp.
    method = "reml"
    sp_r = np.atleast_1d(np.asarray(ref["sp"], float))
    predvar = np.atleast_1d(np.asarray(ref["pred"]["predvar"], float))

    def predict(fit):
        b = fit.params.to_numpy()[cbidx]
        V = fit.vcov[np.ix_(cbidx, cbidx)]
        return dl.crosspred(cb, coef=b, vcov=V, model_link="log", at=predvar, cen=cen)

    t = time.perf_counter()
    fit = dl.fit_pglm(y, X, Slist, family=family, method=method, sp=sp_r)
    t_a = time.perf_counter() - t
    if ref["method"] == "REML":
        tally.check(label, "REML score", fit.reml, ref["score"], TOL["score"])
    tally.check(label, "coefficients", fit.params, ref["coef"], TOL["coef"], mode="scaled")
    tally.check(label, "standard errors", fit.bse, ref["se"], TOL["se"], mode="scaled")
    tally.check(label, "vcov (cross-basis)", fit.vcov[np.ix_(cbidx, cbidx)], r_mat(ref["vcov_cb"]), TOL["vcov"])
    tally.check(label, "edf", fit.edf, ref["edf"], TOL["edf"])
    if family.startswith("quasi"):
        tally.check(label, "scale", fit.scale, ref["scale"], TOL["scale"])
    pr = predict(fit)
    tally.check(label, "crosspred allfit", np.ravel(pr.allfit), ref["pred"]["allfit"], TOL["allfit"])
    tally.check(label, "crosspred allse", np.ravel(pr.allse), ref["pred"]["allse"], TOL["allse"])
    print(f"  {label}: (a) fit at mgcv's sp in {t_a:.1f}s, edf(cb) {fit.edf[cbidx].sum():.4f} "
          f"vs R {ref['edf_cb']:.4f}")

    if not (select and label in SELECT):
        return
    t = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        own = dl.fit_pglm(y, X, Slist, family=family, method="reml")
    t_b = time.perf_counter() - t
    pro = predict(own)
    dscore = own.reml - ref["score"]
    ratio = own.sp / sp_r
    drr = float(np.max(np.abs(np.exp(np.ravel(pro.allfit)) - np.exp(np.asarray(ref["pred"]["allfit"])))))
    ok = dscore <= SCORE_TOL
    tally.n += 1
    tally.fail += not ok
    tally.rows.append((label, "own sp: score not worse than mgcv", 1, max(dscore, 0.0), SCORE_TOL, ok))
    info.append(dict(model=label, seconds=t_b, converged=own.converged, sp_py=own.sp, sp_r=sp_r,
                     ratio=ratio, dscore=dscore, drr=drr, edf_py=float(own.edf[cbidx].sum()),
                     edf_r=ref["edf_cb"], ok=ok,
                     warnings=sorted({f"{w.category.__name__}: {str(w.message)[:80]}" for w in caught})))


# ---- example 1: London -------------------------------------------------------
def example1(ref: dict, tally: Tally, info: list):
    london = pd.read_csv(fetch("london.csv"))
    ns_time = np.asarray(dl.onebasis(london.time.to_numpy(float), "ns", df=10 * 14))
    dow = dummies(london.dow)

    def design(cb):
        B = np.asarray(cb)
        X = np.column_stack([np.ones(len(london)), B, ns_time, dow])
        keep = ~np.isnan(X).any(axis=1)               # R drops the first 25 days
        return london.death.to_numpy(float)[keep], X[keep], np.arange(1, 1 + B.shape[1])

    # GAM WITH DEFAULT PENALTIES
    cbgam1 = dl.crossbasis(london.tmean, lag=25, argvar={"fun": "ps", "df": 9},
                           arglag={"fun": "ps", "df": 10})
    y, X, idx = design(cbgam1)
    fit_compare(tally, info, "ex1_gam1", ref["ex1_gam1"], y, X, cbgam1, idx, dl.cbpen(cbgam1),
                "quasipoisson", 20, True)

    # GAM WITH DOUBLY VARYING PENALTY ON THE LAG
    cbgam2 = dl.crossbasis(london.tmean, lag=25, argvar={"fun": "ps", "df": 9},
                           arglag={"fun": "ps", "df": 10, "fx": True})
    C = np.asarray(dl.onebasis(np.arange(26.0), "ps", df=10, intercept=True))
    D = np.diff(np.eye(25 + 1), n=2, axis=0)
    P = np.diag(np.arange(0, 25 - 1, dtype=float) ** 2)
    Slag1 = C.T @ D.T @ P @ D @ C
    Slag2 = np.diag(np.repeat([0.0, 1.0], [6, 4]))
    tally.check("ex1_gam2", "addSlag Slag1", Slag1, r_mat(ref["ex1_gam2"]["addslag"]["Slag1"]),
                TOL["penalty"] * max(1.0, np.abs(Slag1).max()))
    y, X, idx = design(cbgam2)
    fit_compare(tally, info, "ex1_gam2", ref["ex1_gam2"], y, X, cbgam2, idx,
                dl.cbpen(cbgam2, add_slag=[Slag1, Slag2]), "quasipoisson", 20, True)

    # GAM WITH ONE DIMENSION UNPENALIZED (DOUBLE THRESHOLD; no centring, as in R)
    cbgam3 = dl.crossbasis(london.tmean, lag=25, argvar={"fun": "thr", "thr": [17, 21]},
                           arglag={"fun": "ps", "df": 10})
    y, X, idx = design(cbgam3)
    fit_compare(tally, info, "ex1_gam3", ref["ex1_gam3"], y, X, cbgam3, idx,
                dl.cbpen(cbgam3, add_slag=[Slag1, Slag2]), "quasipoisson", None, True)


# ---- example 2: uranium miners -----------------------------------------------
def fexpfull(start, end, exp):
    """ex2_01.prep.R: yearly exposure profile at ages 0..99 from 5-year totals."""
    age = np.arange(19) * 5
    start = max(round(start - 0.5), 6)
    end = min(round(end + 0.5), 84)
    age[((age - start) < 0) & ((start - age) < 5)] = start
    age[((age - end) > 0) & ((age - end) < 5)] = end
    d = np.diff(age)
    with np.errstate(divide="ignore", invalid="ignore"):
        per = np.asarray(exp, float) / d
    return np.r_[np.repeat(per, d), np.zeros(10)]


def prepare2():
    u = pd.read_csv(fetch("uminers.csv"), index_col=0).sort_values("record").reset_index(drop=True)
    rcols = [f"rexp{5 * k}" for k in range(1, 19)]
    scols = [f"sexp{5 * k}" for k in range(1, 19)]
    R, S = u[rcols].to_numpy(float), u[scols].to_numpy(float)
    cols = {k: u[k].to_numpy(float) for k in ("record", "ind", "agest", "ageexit", "byr",
                                               "rdnstar", "rendage", "smkstar", "sendage")}
    ids, outs, ages, Qx, Qz = [], [], [], [], []
    for i in range(len(u)):
        c = {k: v[i] for k, v in cols.items()}
        age = np.arange(round(c["agest"]), round(c["ageexit"]) + 1)   # R's round: half to even
        ids.append(np.full(age.size, c["record"]))
        outs.append(np.r_[np.zeros(age.size - 1), c["ind"]])
        ages.append(age)
        Qx.append(dl.exphist(fexpfull(c["rdnstar"], c["rendage"], R[i]), age, lag=(2, 40)))
        Qz.append(dl.exphist(fexpfull(c["smkstar"], c["sendage"], S[i]), age, lag=(2, 40)))
    age = np.concatenate(ages)
    main = pd.DataFrame({"id": np.concatenate(ids), "out": np.concatenate(outs), "age": age,
                         "year": 1000 + np.repeat(cols["byr"], [a.size for a in ages]) + age})
    return main, np.vstack(Qx), np.vstack(Qz)


def example2(ref: dict, tally: Tally, info: list):
    t = time.perf_counter()
    main, Qx, Qz = prepare2()
    print(f"  ex2 data: {len(main)} person-years, prepared in {time.perf_counter() - t:.1f}s")
    rd = ref["ex2_data"]
    tally.check("ex2_data", "rows, events, sum(age), sum(year)",
                [len(main), main.out.sum(), main.age.sum(), main.year.sum()],
                [rd["nrow"], rd["sum_out"], rd["sum_age"], rd["sum_year"]], 0)
    tally.check("ex2_data", "column sums of Qx", Qx.sum(0), rd["colsum_Qx"], TOL["data"], mode="rel")
    tally.check("ex2_data", "column sums of Qz", Qz.sum(0), rd["colsum_Qz"], TOL["data"], mode="rel")
    tally.check("ex2_data", "sum of squares Qx, Qz", [np.sum(Qx ** 2), np.sum(Qz ** 2)],
                [rd["sumsq_Qx"], rd["sumsq_Qz"]], TOL["data"], mode="rel")

    cbz = dl.crossbasis(Qz, lag=(2, 40), argvar={"fun": "ns", "knots": 2.5},
                        arglag={"fun": "strata", "breaks": 20})
    vk = np.r_[Qx.min(), dl.logknots(Qx, nk=9), Qx.max()]
    tally.check("ex2_spec", "knots of f(x)", vk, ref["ex2_spec"]["vk"], 1e-14, mode="rel")
    cbgam = dl.crossbasis(Qx, lag=(2, 40), argvar={"fun": "cr", "knots": vk},
                          arglag={"fun": "cr", "intercept": False})
    C = np.asarray(dl.onebasis(np.arange(2.0, 41.0), "cr", df=10, intercept=False))
    P = np.diag(np.repeat([0.0, 1.0], [30 - 2 + 1, 40 - 30]))
    Slag2 = C.T @ P @ C
    tally.check("ex2_spec", "addSlag Slag2", Slag2, r_mat(ref["ex2_spec"]["Slag2"]), TOL["penalty"])
    pen = dl.cbpen(cbgam, add_slag=Slag2)
    age_bs = np.asarray(dl.onebasis(main.age.to_numpy(float), "bs", df=5))
    B = np.asarray(cbgam)
    X = np.column_stack([np.ones(len(main)), B, np.asarray(cbz), age_bs, main.year.to_numpy(float)])
    y = main.out.to_numpy(float)
    idx = np.arange(1, 1 + B.shape[1])
    if X.shape[1] != len(ref["ex2_gam_reml"]["coef"]):
        raise ValueError(f"ex2 design has {X.shape[1]} columns, R has {len(ref['ex2_gam_reml']['coef'])}")
    fit_compare(tally, info, "ex2_gam_ubre", ref["ex2_gam_ubre"], y, X, cbgam, idx, pen, "poisson", 0, False)
    fit_compare(tally, info, "ex2_gam_reml", ref["ex2_gam_reml"], y, X, cbgam, idx, pen, "poisson", 0, True)


# ---- main --------------------------------------------------------------------
SELECT: set = set()          # models that get comparison (b); set in main()


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--quick", action="store_true",
                   help="comparison (a) only, about 1 minute")
    g.add_argument("--all", action="store_true",
                   help="also run (b) on the miners REML model (slow, see the module docstring)")
    args = ap.parse_args(argv)
    SELECT.clear()
    if not args.quick:
        SELECT.update({"ex1_gam1", "ex1_gam2", "ex1_gam3"})
    if args.all:
        SELECT.add("ex2_gam_reml")
    t0 = time.perf_counter()
    ref = json.loads(REFERENCE.read_text())
    v = ref["versions"]
    print(f"Gasparrini et al. (2017) Biometrics, authors' commit {v['authors_commit'][:7]}")
    print(f"reference: {v['R']}, dlnm {v['dlnm']}, mgcv {v['mgcv']}; dlnmpy {dl.__version__}")
    for name, (_, digest) in DATA.items():
        if v["data_sha256"][name] != digest:
            raise RuntimeError(f"reference JSON was built from a different {name}")
    tally, info = Tally(), []
    example1(ref, tally, info)
    example2(ref, tally, info)

    print("\n(a) dlnmpy at mgcv's smoothing parameters")
    print(f"  {'model':<14}{'quantity':<36}{'n':>6}{'max diff':>11}{'tol':>9}")
    for model, what, n, d, tol, ok in tally.rows:
        if what.startswith("own sp"):
            continue
        print(f"  {model:<14}{what:<36}{n:>6}{d:>11.2e}{tol:>9.0e}{'' if ok else '  FAIL'}")

    print(f"\n(b) dlnmpy choosing its own smoothing parameters (REML); "
          f"fails only if its score is worse than mgcv's by > {SCORE_TOL:g}")
    for r in info:
        sp = ", ".join(f"{a:.4g}/{b:.4g}" for a, b in zip(r["sp_py"], r["sp_r"]))
        print(f"  {r['model']}: {r['seconds']:.1f}s, converged={r['converged']}")
        print(f"    sp dlnmpy/mgcv: {sp}  (ratio {', '.join(f'{x:.3g}' for x in r['ratio'])})")
        print(f"    score dlnmpy - mgcv: {r['dscore']:+.3e}   edf(cb) {r['edf_py']:.3f} vs {r['edf_r']:.3f}"
              f"   max |RR diff| overall curve: {r['drr']:.2e}{'' if r['ok'] else '  FAIL'}")
        for w in r["warnings"]:
            print(f"    warning: {w}")

    print(f"\nruntime {time.perf_counter() - t0:.0f}s")
    print(f"biomet2017: {tally.n} comparisons with R, max abs difference {tally.maxdiff:.1e}, "
          f"{tally.fail} failures")
    return 1 if tally.fail else 0


if __name__ == "__main__":
    sys.exit(main())
