"""Reproduce Gasparrini (2014), the exposure-lag-response paper, and check it against R.

    Gasparrini A. Modeling exposure-lag-response associations with distributed
    lag non-linear models. Statistics in Medicine 2014;33(5):881-899.
    doi:10.1002/sim.5963

The paper takes DLNMs out of time series. Its example is the Colorado Plateau
uranium miners cohort: 3,347 men, 258 lung cancer deaths, and yearly radon
exposure (working level months, WLM/year) over a working life. The risk of
death at a given age depends on the whole exposure history at lags 2-40
years, so each row of the design is a 39-year history, and the model is a
Cox model on risk sets sampled from the cohort (case + up to 100 controls).

This is a translation of the author's scripts 01.prep.R to 06.additional.R
(https://github.com/gasparrini/2014_gasparrini_StatMed_Rcodedata), step by
step:

    01  risk sets from survSplit, sampled with set.seed(13041975), yearly
        exposure profiles from 5-year totals, exposure histories by exphist
    02  Table I, descriptive statistics
    03  Table II models 1-4 (DLMs, linear in radon), Figure 1, LR tests
    04  Table II models 5-8 (DLNMs), Figures 2-3, LR tests, peak lag
    05  Table III and Figure 5, predictions for exposure histories
    06  a user-defined right-constrained lag spline, the subset with
        max < 300 WLM/year, a log transform, and the knot position

The Cox fits use dl.fit_clogit with strata = risk sets. The two are the same
likelihood here: after survSplit each row (start, stop] contains exactly one
failure time, so it belongs to one risk set only, and there are no tied
failure times, so Efron, Breslow and exact all agree.

R's risk-set sampling is reproduced exactly by re-implementing R's
Mersenne-Twister seeding and its sample() with the "Rejection" index sampler
(R >= 3.6), so the Python data set is the authors' one row for row.

Each quantity is compared with tests/reference/statmed2014.json, written by
tools/reproduce_statmed2014.R from the authors' own scripts. The simulation
study (07.simcall.R, 08.simres.R) is not reproduced.

Tolerances (an element fails if |python - R| > tol * scale):
    exact     integer equality: counts, sampled risk sets, df, peak lags
    1e-12     data preparation and basis matrices, scale max(|R|, 1); sums
              over 25,906 rows and Table I use 1e-11 for summation order
    1e-6      anything downstream of a fitted model, on a statistical scale:
              coefficients and fits in units of R's standard error, vcov
              entries in units of sqrt(V_ii V_jj), standard errors, the
              log-likelihood, AIC, BIC and LR statistics relative to R's value.
              fit_clogit's covariance is a finite-difference Hessian, so
              agreement past ~1e-8 is not owed.

The fits are compared with R's coxph refitted to convergence (eps = 1e-14).
coxph's default stops when the log-likelihood changes by less than 1e-9,
which for model 8 leaves R's own coefficients 1.6e-6 standard errors short of
the maximum; the distance to the fits as the authors' scripts leave them is
printed but not counted. The crosspred results are the authors' own objects,
from those default fits, and agree within the 1e-6 tolerance regardless.

Needs the cohort data, ~0.6 MB, downloaded once from the pinned commit.

Slow: about 10 minutes, almost all of it in the 15 Cox fits. statsmodels'
ConditionalLogit, which dl.fit_clogit wraps, evaluates the score by a
pure-Python recursion over every risk set (~200 calls per set of 101 rows),
and fit_clogit's covariance differentiates that score numerically, so each
fit takes 30-90 s. R's coxph takes well under a second.
"""

import hashlib
import json
import math
import sys
import time
import urllib.request
import warnings
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

import dlnmpy as dl

SHA = "3ae639c53023b74fd1c332aa0bae3ea83b9dc594"
DATA_URL = ("https://raw.githubusercontent.com/gasparrini/"
            f"2014_gasparrini_StatMed_Rcodedata/{SHA}/uminers.csv")
DATA_SHA256 = "356cceab063a264d537bb14e75f2302197a43aabc2cdbe4a11667e40c280ce29"
CACHE = Path(__file__).with_name("uminers.csv")
REFERENCE = Path(__file__).resolve().parents[1] / "tests" / "reference" / "statmed2014.json"

TOL_DATA, TOL_SUM, TOL_FIT = 1e-12, 1e-11, 1e-6


def load() -> pd.DataFrame:
    if not CACHE.exists():
        print(f"downloading {DATA_URL} ...")
        with urllib.request.urlopen(DATA_URL) as r:
            raw = r.read()
        if hashlib.sha256(raw).hexdigest() != DATA_SHA256:
            raise RuntimeError("uminers.csv does not match the pinned SHA-256")
        CACHE.write_bytes(raw)
    raw = CACHE.read_bytes()
    if hashlib.sha256(raw).hexdigest() != DATA_SHA256:
        raise RuntimeError(f"{CACHE} does not match the pinned SHA-256; delete it and rerun")
    return pd.read_csv(CACHE, index_col=0)


# ---------------------------------------------------------------------------
# R's set.seed() and sample(), so the sampled risk sets are the authors' own
# ---------------------------------------------------------------------------
class RRandom:
    """Mersenne-Twister as seeded by R's set.seed(), with R's sample()."""

    def __init__(self, seed: int):
        s = seed & 0xFFFFFFFF
        for _ in range(50):                      # initial scrambling (RNG.c)
            s = (69069 * s + 1) & 0xFFFFFFFF
        words = []
        for _ in range(625):
            s = (69069 * s + 1) & 0xFFFFFFFF
            words.append(s)
        # words[0] is mti, which FixupSeeds sets to 624 (regenerate on first use)
        self.bg = np.random.MT19937()
        self.bg.state = {"bit_generator": "MT19937",
                         "state": {"key": np.array(words[1:], dtype=np.uint32), "pos": 624}}

    def unif_rand(self) -> float:
        v = int(self.bg.random_raw()) * 2.3283064365386963e-10
        half = 0.5 * 2.328306437080797e-10      # fixup(): keep inside (0, 1)
        return half if v <= 0 else (1 - half if 1 - v <= 0 else v)

    def _rbits(self, bits: int) -> int:
        v, n = 0, 0
        while n <= bits:
            v = 65536 * v + int(math.floor(self.unif_rand() * 65536))
            n += 16
        return v & ((1 << bits) - 1)

    def unif_index(self, dn: int) -> int:
        """R_unif_index with sample.kind = "Rejection"."""
        if dn <= 0:
            return 0
        bits = int(math.ceil(math.log2(dn)))
        while True:
            dv = self._rbits(bits)
            if dn > dv:
                return dv

    def sample(self, x, size: int) -> list:
        """R's sample(x, size) without replacement, including its habit of
        treating a single number n >= 1 as 1:n."""
        x = list(x)
        pop = list(range(1, int(x[0]) + 1)) if len(x) == 1 and x[0] >= 1 else x
        n = len(pop)
        if size > n:
            raise ValueError("cannot take a sample larger than the population")
        idx, out = list(range(n)), []
        for _ in range(size):
            j = self.unif_index(n)
            out.append(pop[idx[j]])
            idx[j] = idx[n - 1]
            n -= 1
        return out


# ---------------------------------------------------------------------------
# functions.R
# ---------------------------------------------------------------------------
def fexpfull(start, end, exp):
    """Yearly exposure at ages 0-99 from the 5-year totals (functions.R)."""
    age = np.arange(19) * 5.0
    start = max(np.round(start - 0.5), 6)       # np.round and R round both go half-to-even
    end = min(np.round(end + 0.5), 84)
    age[((age - start) < 0) & ((start - age) < 5)] = start
    age[((age - end) > 0) & ((age - end) < 5)] = end
    width = np.diff(age)
    return np.concatenate((np.repeat(np.asarray(exp, float) / width, width.astype(int)),
                           np.zeros(10)))


def bsrc(x, df=None, knots=None, degree=3, intercept=False, boundary_knots=None, red=2):
    """06.additional.R: a B-spline with its last `red` columns dropped, which
    constrains the lag curve to zero at the far end."""
    basis, attrs = dl.bs(x, df=df, knots=knots, degree=degree, intercept=intercept,
                         Boundary_knots=boundary_knots)
    if basis.shape[1] <= red:
        raise ValueError("dimension of the basis not compatible with 'red'")
    attrs = dict(attrs, red=red)
    return basis[:, :basis.shape[1] - red], attrs


def mylog(x):
    """06.additional.R: log(x + 1)."""
    return np.log(np.asarray(x, float) + 1), {}


# ---------------------------------------------------------------------------
# 01.prep.R
# ---------------------------------------------------------------------------
def prepare(uminers: pd.DataFrame) -> dict:
    uminers = uminers.sort_values("record", kind="stable").reset_index(drop=True)
    main = uminers.iloc[:, :9]

    # survSplit at every failure time, keeping only the rows that end at one:
    # each subject contributes a row to every risk set it is in
    ftime = np.sort(main.ageexit[main.ind == 1].unique())
    ftime_set = set(ftime.tolist())
    rows = []
    for rec, ind, agest, ageexit in main[["record", "ind", "agest", "ageexit"]].itertuples(index=False):
        cuts = ftime[(ftime > agest) & (ftime < ageexit)]
        starts = np.concatenate(([agest], cuts))
        stops = np.concatenate((cuts, [ageexit]))
        for a, b in zip(starts, stops):
            if b in ftime_set:
                rows.append((rec, a, b, int(ind == 1 and b == ageexit)))
    mainspl = pd.DataFrame(rows, columns=["record", "agest", "ageexit", "ind"])
    mainspl = mainspl.sort_values(["record", "ageexit"], kind="stable").reset_index(drop=True)
    mainspl["riskset"] = np.searchsorted(ftime, mainspl.ageexit) + 1

    # the case and (max) 100 controls from each risk set
    rng = RRandom(13041975)
    subset = []
    for _, x in mainspl.groupby("riskset", sort=True).groups.items():
        x = np.asarray(x) + 1                    # 1-based positions, as in R
        ind = mainspl.ind.to_numpy()[x - 1]
        subset += list(x[ind == 1]) + rng.sample(x[ind == 0], min(len(x) - 1, 100))
    mainsub = mainspl.iloc[np.asarray(subset) - 1]
    mainsub = mainsub.sort_values(["record", "ageexit"], kind="stable").reset_index(drop=True)

    # yearly exposure profiles and the matrices of exposure histories, lag 2-40
    cols = list(uminers.columns)
    rcols, scols = cols[16:34], cols[34:52]
    rexp = np.vstack([fexpfull(r.rdnstar, r.rendage, r[rcols]) for _, r in uminers.iterrows()])
    sexp = np.vstack([fexpfull(r.smkstar, r.sendage, r[scols]) for _, r in uminers.iterrows()])
    rec = mainsub.record.to_numpy(int)
    t = mainsub.ageexit.to_numpy(float)
    Qx = np.vstack([dl.exphist(rexp[r - 1], t[i], (2, 40)) for i, r in enumerate(rec)])
    Qz = np.vstack([dl.exphist(sexp[r - 1], t[i], (2, 40)) for i, r in enumerate(rec)])

    # calendar time: birth date on the 15th of the birth month, in years since 1970
    byear = mainsub.record.map(dict(zip(uminers.record, uminers.byr)))
    bmon = mainsub.record.map(dict(zip(uminers.record, uminers.bmon)))
    days = np.array([(date(int(f"1{y}"), int(m), 15) - date(1970, 1, 1)).days
                     for y, m in zip(byear, bmon)], float)
    caltime = days / 365.25 + t

    cbz = dl.crossbasis(Qz, lag=(2, 40), argvar={"fun": "ns", "knots": 2.5},
                        arglag={"fun": "strata", "breaks": 20})
    return dict(uminers=uminers, main=main, ftime=ftime, mainspl=mainspl,
                mainsub=mainsub, rexp=rexp, sexp=sexp, Qx=Qx, Qz=Qz, caltime=caltime, cbz=cbz)


# the 8 models of Table II (01.prep.R). R writes int=TRUE, partial matching
# for intercept=TRUE; dlnmpy needs the full name.
LAGLIST = [{"fun": "strata", "df": 1, "intercept": True},
           {"fun": "strata", "breaks": [10, 20, 30], "intercept": True},
           {"fun": "bs", "degree": 2, "knots": 13.3, "intercept": True},
           {"fun": "bs", "degree": 2, "knots": 13.3, "intercept": False}] * 2
_VAR = [{"fun": "lin"}, {"fun": "strata", "breaks": [26.7, 60.2, 122.2]},
        {"fun": "bs", "degree": 2, "knots": 60.2}]
VARLIST = [_VAR[i] for i in (0, 0, 0, 0, 2, 1, 2, 2)]


def coxfit(cbx, d, rows=None):
    """coxph(Surv(agest, ageexit, ind) ~ cbx + cbz + caltime, mainsub[rows, ])"""
    rows = slice(None) if rows is None else rows
    sub = d["mainsub"][rows]
    X = pd.concat([cbx.to_dataframe("cbx", index=sub.index),
                   pd.DataFrame(d["cbz"].matrix[rows], index=sub.index,
                                columns=d["cbz"].to_dataframe("cbz").columns),
                   pd.Series(d["caltime"][rows], index=sub.index, name="caltime")], axis=1)
    with warnings.catch_warnings():
        # in the subset analysis (max < 300 WLM/year) 57 risk sets lose their
        # case; statsmodels says it drops them, coxph drops them silently
        warnings.filterwarnings("ignore", message=r"Dropped \d+ groups")
        m = dl.fit_clogit(sub.ind, X, groups=sub.riskset)
    k = len(m.params)
    m.aic_cox = -2 * m.llf + 2 * k                   # AIC.coxph in functions.R
    m.bic_cox = -2 * m.llf + k * math.log(int(sub.ind.sum()))
    return m


def pred(cbx, m, **kw):
    if hasattr(m.params, "index"):
        return dl.crosspred(cbx, m, name="cbx", model_link="log", **kw)
    # fit_clogit's step-halving fallback (taken for the cumulative-exposure
    # models) returns unnamed parameters, so crosspred cannot find the
    # cross-basis by name; the cbx columns are the first ones in the design
    k = cbx.ncol
    return dl.crosspred(cbx, coef=np.asarray(m.params)[:k],
                        vcov=np.asarray(m.cov_params())[:k, :k], model_link="log", **kw)


def lrtest(m1, m2):
    stat = 2 * abs(m1.llf - m2.llf)
    df = abs(len(m1.params) - len(m2.params))
    return stat, df, 1 - chi2.cdf(stat, df)


def table1(d) -> np.ndarray:
    """Table I (02.descriptive.R)."""
    u, Qx, sexp = d["uminers"], d["Qx"], d["sexp"]
    def q(v):                                   # quantile() is type 7 in both
        return [np.median(v), v.min(), np.quantile(v, .25), np.quantile(v, .75), v.max()]

    def fvtab(v, ind):
        v = np.asarray(v, float)
        return q(v) + q(v[ind])

    def fmtab(m, ind):                          # positive entries only
        return q(m[m > 0]) + q(m[ind][m[ind] > 0])

    vcase = (u.ind == 1).to_numpy()
    mcase = (d["mainsub"].ind == 1).to_numpy()
    vsmk = (u.totsmk > 0).to_numpy()
    t = np.full((18, 10), np.nan)
    t[0, [0, 1, 5, 6]] = [len(u), 100, vcase.sum(), 100]
    dead = (u.status == 1).to_numpy()
    t[1, [0, 1, 5, 6]] = [dead.sum(), dead.sum() / t[0, 0] * 100, dead[vcase].sum(),
                          dead[vcase].sum() / t[0, 5] * 100]
    t[2, [0, 1, 5, 6]] = [vsmk.sum(), vsmk.sum() / t[0, 0] * 100, vsmk[vcase].sum(),
                          vsmk[vcase].sum() / t[0, 5] * 100]
    t[3] = fvtab(u.agest, vcase)
    t[4] = fvtab(u.ageexit - u.agest, vcase)
    t[6] = fvtab(u.rendage - u.rdnstar, vcase)
    t[7] = fvtab(u.totrdn, vcase)
    t[9] = fmtab(Qx, mcase)
    for row, (a, b) in zip(range(10, 14), ((0, 8), (8, 18), (18, 28), (28, 39))):
        t[row] = fmtab(Qx[:, a:b], mcase)
    us = u[vsmk]
    t[15] = fvtab(us.sendage - us.smkstar, vcase[vsmk])
    t[16] = fvtab(us.totsmk, vcase[vsmk])
    t[17] = fmtab(sexp[vsmk], vcase[vsmk])
    return t


# ---------------------------------------------------------------------------
# comparison with R
# ---------------------------------------------------------------------------
def _arr(x) -> np.ndarray:
    """R's JSON (null for NA) as a float array."""
    if isinstance(x, list):
        return np.array([_arr(v) for v in x], float) if x else np.array([], float)
    return np.asarray(np.nan if x is None else x, float)


class Checker:
    def __init__(self):
        self.groups: dict[str, dict] = {}

    def __call__(self, group, label, py, r, tol, scale=None):
        py = np.asarray(py, float)
        r = _arr(r)
        if r.ndim == 0 and py.size == 1:         # jsonlite unboxes length-1 vectors
            r = r.reshape(py.shape)
        g = self.groups.setdefault(group, {"n": 0, "abs": 0.0, "err": 0.0, "fail": 0,
                                           "tol": tol, "labels": []})
        g["tol"] = max(g["tol"], tol)
        if py.shape != r.shape:
            g["n"] += max(r.size, 1)
            g["fail"] += max(r.size, 1)
            g["labels"].append(f"{label}: shape {py.shape} vs R {r.shape}")
            return
        both_nan = np.isnan(py) & np.isnan(r)
        diff = np.where(both_nan, 0.0, np.abs(py - r))
        diff = np.where(np.isnan(diff), np.inf, diff)
        if tol == 0:
            err = diff
            bad = diff > 0
        else:
            s = np.maximum(np.abs(r), 1.0) if scale is None else np.asarray(scale, float)
            s = np.where(both_nan | ~np.isfinite(s), 1.0, s)
            s = np.maximum(s, 1e-8 * np.nanmax(s) if s.size else 1.0)
            err = diff / s
            bad = err > tol
        g["n"] += r.size
        if r.size:
            g["abs"] = max(g["abs"], float(np.max(diff)))
            g["err"] = max(g["err"], float(np.max(err)))
        nbad = int(bad.sum())
        g["fail"] += nbad
        if nbad:
            g["labels"].append(f"{label}: {nbad} of {r.size} beyond tolerance, "
                               f"max scaled error {float(np.max(err)):.2e}")

    def model(self, group, m, r):
        """Compare with R's fit refitted to convergence (see the R script); the
        distance to the fit as the authors' scripts leave it is only reported."""
        pub = r["published"]
        se = np.asarray(pub["se"])
        V = np.asarray(pub["vcov"])
        self.published = max(getattr(self, "published", 0.0),
                              float(np.max(np.abs(np.asarray(m.params) - pub["coef"]) / se)),
                              float(np.max(np.abs(np.asarray(m.cov_params()) - V)
                                           / np.sqrt(np.outer(np.diag(V), np.diag(V))))))
        r = r["converged"]
        se = np.asarray(r["se"])
        V = np.asarray(r["vcov"])
        self(group, "coef", m.params, r["coef"], TOL_FIT, se)
        self(group, "vcov", m.cov_params(), V, TOL_FIT, np.sqrt(np.outer(np.diag(V), np.diag(V))))
        self(group, "se", m.bse, se, TOL_FIT, se)
        self(group, "loglik", m.llf, r["loglik"], TOL_FIT, abs(r["loglik"]))
        self(group, "AIC", m.aic_cox, r["AIC"], TOL_FIT, abs(r["AIC"]))
        self(group, "BIC", m.bic_cox, r["BIC"], TOL_FIT, abs(r["BIC"]))

    def pred(self, group, p, r):
        self(group, "predvar", p.predvar, r["predvar"], TOL_DATA)
        se = np.asarray(r["matse"])
        self(group, "matfit", p.matfit, r["matfit"], TOL_FIT, se)
        self(group, "matse", p.matse, se, TOL_FIT, se)
        ase = np.asarray(r["allse"])
        self(group, "allfit", p.allfit, r["allfit"], TOL_FIT, ase)
        self(group, "allse", p.allse, ase, TOL_FIT, ase)
        for k in ("allRRfit", "allRRlow", "allRRhigh"):
            self(group, k, getattr(p, k), r[k], TOL_FIT, np.abs(r[k]))

    def basis(self, group, cb, r, rows):
        self(group, "dim", cb.shape, r["dim"], 0)
        self(group, "df", cb.df, r["df"], 0)
        self(group, "range", cb.range, r["range"], TOL_DATA)
        self(group, "colsums", np.asarray(cb).sum(axis=0), r["colsums"], TOL_SUM)
        self(group, "rows", np.asarray(cb)[rows], r["rows"], TOL_DATA)


def main() -> int:
    t0 = time.time()
    R = json.loads(REFERENCE.read_text())
    if R["versions"]["authors_commit"] != SHA:
        raise RuntimeError("reference was written from a different commit of the authors' code")
    check = Checker()

    # ---- 01 data preparation ----------------------------------------------
    u = load()
    d = prepare(u)
    rp = R["prep"]
    rows = np.asarray(rp["rows"], int) - 1
    ms = d["mainsub"]
    c = "data: risk sets, histories"
    check(c, "counts", [len(d["uminers"]), len(d["main"]), d["main"].ind.sum(), d["ftime"].size,
                        len(d["mainspl"]), len(ms), ms.ind.sum()],
          [rp["n_uminers"], rp["n_main"], rp["n_cases"], rp["n_ftime"], rp["n_mainspl"],
           rp["n_mainsub"], rp["n_mainsub_cases"]], 0)
    check(c, "ftime", d["ftime"], rp["ftime"], TOL_DATA)
    check(c, "sampled records", ms.record, rp["mainsub_record"], 0)
    check(c, "sampled risk sets", ms.riskset, rp["mainsub_riskset"], 0)
    check(c, "sum agest", ms.agest.sum(), rp["mainsub_agest_sum"], TOL_SUM)
    for k in ("rexp", "sexp"):
        check(c, f"{k} row sums", d[k].sum(axis=1), rp[f"{k}_rowsums"], TOL_DATA)
        check(c, f"{k} col sums", d[k].sum(axis=0), rp[f"{k}_colsums"], TOL_SUM)
    for k in ("Qx", "Qz"):
        check(c, f"{k} col sums", d[k].sum(axis=0), rp[f"{k}_colsums"], TOL_SUM)
        check(c, f"{k} rows", d[k][rows], rp[f"{k}_rows"], TOL_DATA)
    check(c, "caltime sum", d["caltime"].sum(), rp["caltime_sum"], TOL_SUM)
    check(c, "caltime rows", d["caltime"][rows], rp["caltime_rows"], TOL_DATA)
    check.basis("bases: cross-bases", d["cbz"], rp["cbz"], rows)

    # ---- 02 Table I -------------------------------------------------------
    t1 = table1(d)
    check("Table I", "table1", t1, R["table1"], TOL_SUM)

    # ---- 03/04 Table II ---------------------------------------------------
    Qx = d["Qx"]
    cbx, models, preds = {}, {}, {}
    for i in range(1, 9):
        cbx[i] = dl.crossbasis(Qx, lag=(2, 40), argvar=VARLIST[i - 1], arglag=LAGLIST[i - 1])
        models[i] = coxfit(cbx[i], d)
        preds[i] = pred(cbx[i], models[i], at=np.arange(26) * 10.0,
                        **({} if i <= 4 else {"cen": 0}))
        check.basis("bases: cross-bases", cbx[i], R["cbs"][f"cbx{i}"], rows)
        check.model("Table II: models 1-8", models[i], R["models"][f"model{i}"])
        check.pred("Figures 1-3: crosspred 1-8", preds[i], R["preds"][f"pred{i}"])
    rt2 = R["table2"]
    check("Table II: models 1-8", "df", [cbx[i].ncol for i in range(1, 9)], rt2["df"], 0)

    def lag1012(p):
        j = [list(p.predlag).index(v) for v in (10, 11, 12)]
        i = list(p.predvar).index(100)
        return {k: getattr(p, f"matRR{k}")[i, j] for k in ("fit", "low", "high")}

    for name, i in (("pred1", 1), ("pred4", 4), ("pred8", 8)):
        got = lag1012(preds[i])
        for k in ("fit", "low", "high"):
            check("HR 100 WLM at lags 10-12", f"{name} {k}", got[k], R["lag1012"][name][k],
                  TOL_FIT, np.abs(R["lag1012"][name][k]))

    lr = {"m3_m1": (3, 1), "m3_m4": (3, 4), "m8_m4": (8, 4), "m8_m5": (8, 5), "m8_m7": (8, 7)}
    lrres = {}
    for key, (a, b) in lr.items():
        stat, df, p = lrres[key] = lrtest(models[a], models[b])
        r = R["lrtests"][key]
        check("LR tests", f"{key} stat", stat, r["stat"], TOL_FIT, abs(r["stat"]))
        check("LR tests", f"{key} df", df, r["df"], 0)
        # 1 - pchisq() carries an absolute rounding error of ~1e-16 in both
        # languages, so small p-values are compared on an absolute floor
        check("LR tests", f"{key} p", p, r["pvalue"], TOL_FIT, max(abs(r["pvalue"]), 1e-6))

    peak = preds[8].predlag[np.argmax(preds[8].matRRfit, axis=1)]
    check("peak lag, model 8", "peak", peak, R["peak8"], 0)

    # ---- 05 Table III and Figure 5 ----------------------------------------
    hist = np.vstack([dl.exphist(np.full(10, 20.0), times=10, lag=(2, 40)),
                      dl.exphist(np.full(10, 100.0), times=10, lag=(2, 40)),
                      dl.exphist(np.full(20, 20.0), times=20, lag=(2, 40)),
                      dl.exphist(np.repeat([20.0, 0], [10, 10]), times=20, lag=(2, 40)),
                      dl.exphist(np.repeat([20.0, 0], [10, 30]), times=40, lag=(2, 40))])
    rpred = R["predictions"]
    check("Table III, Figure 5", "histories", hist, rpred["hist"], TOL_DATA)
    tab3 = {}
    for i in (1, 4, 5, 8):
        tab3[i] = pred(cbx[i], models[i], at=hist, cen=0)
        r = rpred["tab3"][f"model{i}"]
        for k in ("allfit", "allse", "allRRfit", "allRRlow", "allRRhigh"):
            scale = np.asarray(r["allse"]) if k in ("allfit", "allse") else np.abs(r[k])
            check("Table III, Figure 5", f"model{i} {k}", getattr(tab3[i], k), r[k], TOL_FIT, scale)
    profile = np.repeat([20.0, 0], [15, 45])
    dynhist = dl.exphist(profile, lag=(2, 40))
    check("Table III, Figure 5", "dynhist", dynhist, rpred["dynhist"], TOL_DATA)
    predhist = {}
    for i in (4, 5, 8):
        predhist[i] = pred(cbx[i], models[i], at=dynhist, cen=0)
        r = rpred[f"predhist{i}"]
        for k in ("allfit", "allse", "allRRfit", "allRRlow", "allRRhigh"):
            scale = np.asarray(r["allse"]) if k in ("allfit", "allse") else np.abs(r[k])
            check("Table III, Figure 5", f"predhist{i} {k}", getattr(predhist[i], k), r[k],
                  TOL_FIT, scale)

    # ---- 06 additional analyses -------------------------------------------
    ra = R["additional"]
    at = np.arange(26) * 10.0
    g = "06: constrained, subset, log, knots"
    cbxrc = dl.crossbasis(Qx, lag=(2, 40), argvar=VARLIST[7],
                          arglag={"fun": bsrc, "degree": 2, "knots": [13.3, 26.7],
                                  "intercept": False, "red": 2})
    modelrc = coxfit(cbxrc, d)
    check.basis("bases: cross-bases", cbxrc, ra["rc"]["cb"], rows)
    check.model(g, modelrc, ra["rc"]["model"])
    check.pred(g, pred(cbxrc, modelrc, at=at, cen=0), ra["rc"]["pred"])

    sub = Qx.max(axis=1) < 300
    check(g, "subset size", [sub.sum(), ms.record[sub].nunique()],
          [ra["sub"]["n"], ra["sub"]["subjects"]], 0)
    sub_models = {}
    for i in (4, 8):
        cb = dl.crossbasis(Qx[sub], lag=(2, 40), argvar=VARLIST[i - 1], arglag=LAGLIST[i - 1])
        sub_models[i] = m = coxfit(cb, d, rows=sub)
        check.model(g, m, ra["sub"][f"model{i}"])
        check.pred(g, pred(cb, m, at=at, cen=0), ra["sub"][f"pred{i}"])

    cbxlog = dl.crossbasis(Qx, lag=(2, 40), argvar={"fun": mylog}, arglag=LAGLIST[3])
    modellog = coxfit(cbxlog, d)
    check.basis("bases: cross-bases", cbxlog, ra["log"]["cb"], rows)
    check.model(g, modellog, ra["log"]["model"])
    check.pred(g, pred(cbxlog, modellog, at=at, cen=0), ra["log"]["pred"])

    knot_models = {}
    for key, knots in (("8b", [20]), ("8c", [26.6]), ("8d", [13.3, 26.6])):
        cb = dl.crossbasis(Qx, lag=(2, 40), argvar=VARLIST[7],
                           arglag={"fun": "bs", "df": 3, "degree": 2, "knots": knots,
                                   "intercept": False})
        knot_models[key] = m = coxfit(cb, d)
        check.model(g, m, ra["knots"][f"model{key}"])
        check.pred(g, pred(cb, m, at=at, cen=0), ra["knots"][f"pred{key}"])

    # ---- reductions of model 8 (not in the paper) -------------------------
    for key, kw in (("overall", {}), ("var100", {"type": "var", "value": 100}),
                    ("lag15", {"type": "lag", "value": 15})):
        red = dl.crossreduce(cbx[8], models[8], cen=0, name="cbx", model_link="log", **kw)
        r = R["reductions"][key]
        V = np.asarray(r["vcov"])
        se = np.sqrt(np.diag(V))
        check("crossreduce, model 8", f"{key} coef", np.ravel(red.coef), r["coef"], TOL_FIT, se)
        check("crossreduce, model 8", f"{key} vcov", red.vcov, V, TOL_FIT, np.sqrt(np.outer(se, se)))
        check("crossreduce, model 8", f"{key} fit", red.fit, r["fit"], TOL_FIT, np.asarray(r["se"]))
        check("crossreduce, model 8", f"{key} se", red.se, r["se"], TOL_FIT, np.asarray(r["se"]))

    elapsed = time.time() - t0

    # ---- summary ----------------------------------------------------------
    print(f"\nColorado Plateau uranium miners: {len(d['uminers'])} men, "
          f"{int(d['main'].ind.sum())} lung cancer deaths")
    print(f"  {len(ms)} rows in {d['ftime'].size} sampled risk sets "
          f"(case + up to 100 controls, R seed 13041975)")
    print("\n  Table II   f(x)    w(lag)  df     AIC     BIC")
    for i in range(1, 9):
        print(f"    model {i}  {VARLIST[i - 1]['fun']:6s}  {LAGLIST[i - 1]['fun']:6s} "
              f"{cbx[i].ncol:3d}  {models[i].aic_cox:7.1f} {models[i].bic_cox:7.1f}")
    print("\n  LR tests: " + ", ".join(
        f"{k.replace('_', ' vs ').replace('m', 'M')} {s:.2f} on {df} df (p={p:.3f})"
        for k, (s, df, p) in lrres.items()))
    for i in (4, 8):
        h = lag1012(preds[i])
        print(f"  HR for 100 WLM/year at lag 11, model {i}: "
              f"{h['fit'][1]:.4f} ({h['low'][1]:.4f}-{h['high'][1]:.4f})")
    print(f"  peak lag for 100 WLM/year, model 8: {peak[10]:.0f} years")
    print("\n  Table III, HR (95% CI)   " + "".join(f"model {i:<14d}" for i in (1, 4, 5, 8)))
    for s in range(5):
        print(f"    scenario {s + 1}            " + "".join(
            f"{tab3[i].allRRfit[s]:.2f} ({tab3[i].allRRlow[s]:.2f}-{tab3[i].allRRhigh[s]:.2f})  "
            for i in (1, 4, 5, 8)))
    print("  Figure 5, HR at year 20 after 15 years at 20 WLM/year: " + ", ".join(
        f"model {i} {predhist[i].allRRfit[19]:.2f}" for i in (4, 5, 8)))
    print(f"  AIC: right-constrained {modelrc.aic_cox:.1f}, log {modellog.aic_cox:.1f}, "
          + ", ".join(f"knots {k} {m.aic_cox:.1f}" for k, m in knot_models.items()))

    print(f"\n  comparison with R ({R['versions']['R']}, dlnm {R['versions']['dlnm']}, "
          f"survival {R['versions']['survival']})")
    m2 = models[2]
    path = "Newton" if hasattr(m2.params, "index") else "step-halving Newton fallback"
    for which in ("published", "converged"):
        r2 = R["models"]["model2"][which]
        print(f"    Table II model 2 ({path}) vs coxph {which:9s}: coef max abs diff "
              f"{np.max(np.abs(np.asarray(m2.params) - r2['coef'])):.1e}, SE max abs diff "
              f"{np.max(np.abs(np.asarray(m2.bse) - r2['se'])):.1e}")
    print(f"    all fits vs coxph as published (eps=1e-9, not counted): max "
          f"{check.published:.1e} SE; R's own refit at eps=1e-14 moves by the same amount")
    print(f"    {'group':40s} {'n':>6s} {'max abs diff':>13s} {'max scaled':>11s} {'tol':>7s}")
    total, worst, failures = 0, 0.0, 0
    for name, gr in check.groups.items():
        tol = "exact" if gr["tol"] == 0 else f"{gr['tol']:.0e}"
        status = "ok" if gr["fail"] == 0 else f"FAIL ({gr['fail']})"
        print(f"    {name:40s} {gr['n']:6d} {gr['abs']:13.2e} {gr['err']:11.2e} {tol:>7s}  {status}")
        for lab in gr["labels"]:
            print(f"      {lab}")
        total += gr["n"]
        worst = max(worst, gr["abs"])
        failures += gr["fail"]
    print(f"  (tol is the loosest in the group: fits and coefficients in units of R's"
          f" standard error, otherwise relative; {elapsed:.0f} s)")
    print(f"statmed2014: {total} comparisons with R, max abs difference {worst:.1e}, "
          f"{failures} failures")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
