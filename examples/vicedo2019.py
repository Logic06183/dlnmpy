"""Reproduce Vicedo-Cabrera, Sera & Gasparrini (2019) Epidemiology.

    Vicedo-Cabrera AM, Sera F, Gasparrini A. Hands-on tutorial on a modeling
    framework for projections of climate change impacts on health.
    Epidemiology 2019;30(3):321-329. doi:10.1097/EDE.0000000000000982

A translation of the authors' scripts at
https://github.com/gasparrini/2019_vicedo-cabrera_Epidem_Rcodedata (commit
313606cd), section by section:

    01  DLNM of daily deaths on temperature in London 1990-2012, and the MMT
    02  projected mortality series (daily mean deaths by day of year)
    03  bias correction of 5 GCMs x 2 RCPs (the authors' fhempel.r, ported here)
    04  extrapolation of the curve to the projected temperatures
    05  projected attributable deaths by decade, cold/heat, GCM and RCP
    06  ensemble estimates with Monte Carlo empirical intervals
    07  age-specific curves and the adaptation scenario

Every step is compared with the R output written by
tools/reproduce_vicedo2019.R (tests/reference/vicedo2019.json). The Monte
Carlo intervals replay R's own normal deviates through
``simulate_coef(normals=Z)``, so they are compared exactly, not up to
simulation error. The data (~9.5 MB) are downloaded once from the pinned
commit into examples/vicedo2019_data/ and checked against SHA-256 digests.

Run from anywhere: python examples/vicedo2019.py
"""

import hashlib
import json
import os
import sys
import time
import urllib.request
from pathlib import Path

# The fits are small (8,258 x 214); a threaded BLAS only adds overhead, and
# on a busy machine slowed this script four-fold. Override by setting these.
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import dlnmpy as dl  # noqa: E402

COMMIT = "313606cd0fbeec3cd52e72cd4fd5956faece9c48"
BASE = ("https://raw.githubusercontent.com/gasparrini/"
        f"2019_vicedo-cabrera_Epidem_Rcodedata/{COMMIT}/")
SHA256 = {
    "lndn_obs.csv": "9133393737710d5a666b41da41526fc5f9c9a153178047624190ebb6199f6e52",
    "lndn_rcp4p5.csv": "e7cb0d7bbd1660930eb7dd8adac33f8a2a745968859b6082f2c6931f14801430",
    "lndn_rcp8p5.csv": "d36bc4e26f1ad6eec882e18bf9e3813f6639e736a90b098a6ea072bf88ccac25",
}
HERE = Path(__file__).resolve().parent
DATA = HERE / "vicedo2019_data"
REF = HERE.parent / "tests" / "reference" / "vicedo2019.json"


def fetch(name: str) -> Path:
    path = DATA / name
    if not path.exists():
        DATA.mkdir(exist_ok=True)
        print(f"downloading {name} ...")
        with urllib.request.urlopen(BASE + name) as r:      # noqa: S310
            blob = r.read()
        tmp = path.with_suffix(".part")
        tmp.write_bytes(blob)
        tmp.rename(path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != SHA256[name]:
        raise SystemExit(f"{path}: SHA-256 {digest} does not match the pinned file; "
                         "delete it to download again")
    return path


# ----------------------------------------------------------------------------
# fhempel.r: the ISI-MIP trend-preserving bias correction (Hempel et al. 2013)
def fhempel(obs: pd.DataFrame, mod: pd.DataFrame, add=True, mult=True, output="series"):
    """Line-by-line port of the authors' ``fhempel``. ``obs`` is (date, series),
    ``mod`` is date plus one or more modelled series."""
    odate, mdate = obs.iloc[:, 0], mod.iloc[:, 0]
    out, corr = {}, {}
    for col in mod.columns[1:]:
        # period with no missing values in either series
        common = odate[odate.isin(set(mdate))]
        oi = np.flatnonzero(odate.isin(set(common)))
        mi = np.flatnonzero(mdate.isin(set(common)))
        notna = ~(np.isnan(obs.iloc[oi, 1].to_numpy(float))
                  | np.isnan(mod[col].iloc[mi].to_numpy(float)))
        oi, mi = oi[notna], mi[notna]
        if oi.size == 0:
            raise ValueError("no common non-missing days in 'obs' and 'mod'")
        month = odate.iloc[oi].dt.month.to_numpy()
        if np.unique(month).size != 12:
            raise ValueError("some months are not represented in the overlapping period")
        monthyear = (odate.iloc[oi].dt.year * 100 + odate.iloc[oi].dt.month).to_numpy()
        yo = obs.iloc[oi, 1].to_numpy(float)
        ym = mod[col].iloc[mi].to_numpy(float)

        # additive correction: difference in monthly means
        months = np.arange(1, 13)
        C = np.array([np.nanmean(yo[month == m]) - np.nanmean(ym[month == m]) for m in months])
        if not add:
            C[:] = 0
        # residuals from month-year means, multiplicative correction by month:
        # slope through the origin of sorted observed on sorted modelled residuals
        resobs = yo - pd.Series(yo).groupby(monthyear).transform("mean").to_numpy()
        resmod = ym - pd.Series(ym).groupby(monthyear).transform("mean").to_numpy()
        B = np.empty(12)
        for m in months:
            x, y = np.sort(resmod[month == m]), np.sort(resobs[month == m])
            B[m - 1] = np.linalg.lstsq(x[:, None], y, rcond=None)[0][0]
        if not mult:
            B[:] = 1
        if output == "correction":
            corr[col] = {"add": C, "mult": B}
            continue

        # interpolate the corrections across month boundaries (Hempel p. 228)
        day, mon = mdate.dt.day.to_numpy(), mdate.dt.month.to_numpy()
        my = (mdate.dt.year * 100 + mdate.dt.month).to_numpy()
        x = mod[col].to_numpy(float)
        nday = pd.Series(x).groupby(my).transform("size").to_numpy()
        d = (day - 1) / (nday - 1) - 0.5
        dm, d0, dp = 0.5 * (np.abs(d) - d), 1 - np.abs(d), 0.5 * (np.abs(d) + d)
        Cday = dm * np.roll(C, 1)[mon - 1] + d0 * C[mon - 1] + dp * np.roll(C, -1)[mon - 1]
        Bday = dm * np.roll(B, 1)[mon - 1] + d0 * B[mon - 1] + dp * np.roll(B, -1)[mon - 1]
        myavg = pd.Series(x).groupby(my).transform("mean").to_numpy()
        out[col] = myavg + Cday + (x - myavg) * Bday
    if output == "correction":
        return corr
    res = pd.DataFrame(out, index=mod.index)
    res.insert(0, mod.columns[0], mdate)
    return res


# ----------------------------------------------------------------------------
class Compare:
    """Collects |python - R| <= atol + rtol * |R| checks by group."""

    def __init__(self):
        self.rows, self.n, self.fail, self.maxabs = [], 0, 0, 0.0

    def __call__(self, name, py, r, rtol, atol=0.0):
        py = np.asarray(py, float).ravel()
        r = np.asarray(r, float).ravel()
        if py.shape != r.shape:
            print(f"  {name}: shape {py.shape} vs R {r.shape}")
            self.n += 1
            self.fail += 1
            self.rows.append((name, 1, np.inf, np.inf, rtol, atol, 1))
            return
        diff = np.abs(py - r)
        bad = int(np.sum(~(diff <= atol + rtol * np.abs(r))))
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(r != 0, diff / np.abs(r), diff)
        md, mr = float(diff.max()), float(rel.max())
        self.rows.append((name, r.size, md, mr, rtol, atol, bad))
        self.n += r.size
        self.fail += bad
        self.maxabs = max(self.maxabs, md)

    def report(self):
        print(f"\n{'quantity':38s} {'n':>6s} {'max abs':>9s} {'max rel':>9s} "
              f"{'rtol':>7s} {'atol':>7s} {'fail':>5s}")
        for name, n, md, mr, rtol, atol, bad in self.rows:
            print(f"{name:38s} {n:6d} {md:9.1e} {mr:9.1e} {rtol:7.0e} "
                  f"{atol:7.0e} {bad:5d}")


def rarr(d):
    """An R array exported by tools/reproduce_vicedo2019.R, as numpy."""
    return np.asarray(d["values"], float).reshape(d["dim"], order="F")


# ----------------------------------------------------------------------------
def main() -> int:
    t0 = time.time()
    ref = json.loads(REF.read_text())
    if ref["versions"]["authors_commit"] != COMMIT:
        raise SystemExit("reference JSON was made from a different commit of the authors' code")
    cmp = Compare()
    # Tolerances. Anything that passes through a GLM fit is held to 1e-8
    # relative: R stops IRLS on a relative deviance change of 1e-8 and
    # dlnmpy on 1e-10, so the two fits are the same optimum to roughly that
    # precision (observed: about 1e-10 and below). Pure arithmetic on data
    # (calibration, mortality series, draws) is held to 1e-12.
    FIT, ARITH = 1e-8, 1e-12

    # ---- 01 estimation of the exposure-response association ----------------
    obs = pd.read_csv(fetch("lndn_obs.csv"))
    obs["date"] = pd.to_datetime(obs["date"], format="%d/%m/%Y")
    tmean = obs.tmean.to_numpy(float)

    argvar = {"fun": "ns", "knots": dl.quantile7(tmean, [0.10, 0.75, 0.90]),
              "Bound": [tmean.min(), tmean.max()]}
    argvar1 = {"fun": "lin"}
    argvar2 = {"fun": "thr", "thr.value": dl.quantile7(tmean, [0.10, 0.90]), "side": "d"}
    maxlag = 21
    arglag = {"fun": "ns", "knots": dl.logknots(maxlag, nk=3)}

    cb = dl.crossbasis(obs.tmean, lag=maxlag, argvar=argvar, arglag=arglag)
    cb1 = dl.crossbasis(obs.tmean, lag=maxlag, argvar=argvar1, arglag=arglag)
    cb2 = dl.crossbasis(obs.tmean, lag=maxlag, argvar=argvar2, arglag=arglag)

    # all ~ cb + dow + ns(date, df=round(8*length(date)/365.25)), quasipoisson
    seas_df = int(round(8 * len(obs) / 365.25))
    datenum = (obs.date - pd.Timestamp("1970-01-01")).dt.days
    seas = dl.onebasis(datenum, "ns", df=seas_df)

    def fit(outcome, basis):
        X = dl.design_matrix(obs, ("cb", basis), ("seas", seas), intercept=False)
        return dl.fit_glm(f"Q('{outcome}') ~ " + " + ".join(X.columns) + " + C(dow)",
                          obs.join(X), family="quasipoisson")

    def cbpart(model):
        names = [c for c in model.params.index if c.startswith("cb_")]
        return (model.params[names].to_numpy(),
                model.cov_params().loc[names, names].to_numpy())

    m, m1, m2 = fit("all", cb), fit("all", cb1), fit("all", cb2)
    s1 = ref["stage1"]
    cmp("seasonal spline df", seas_df, s1["seas_df"], 0)
    cmp("fitted rows", m.nobs, s1["nobs"], 0)
    for key, model in (("ns", m), ("lin", m1), ("thr", m2)):
        b, V = cbpart(model)
        cmp(f"01 cross-basis coef ({key})", b, s1[key]["coef"], FIT)
        # covariance entries near zero are held to 1e-8 of the largest entry
        Vr = np.asarray(s1[key]["vcov"], float)
        cmp(f"01 cross-basis vcov ({key})", V, Vr, FIT, atol=FIT * np.abs(Vr).max())
    cmp("01 dispersion", m.scale, s1["dispersion"], FIT)
    cmp("01 deviance", m.deviance, s1["deviance"], FIT)

    # provisional centring, then the MMT as the minimum of the overall curve
    varcen = 19
    cp = dl.crosspred(cb, m, cen=varcen, by=0.1)
    cen = float(cp.predvar[np.argmin(cp.allRRfit)])
    cmp("01 prediction grid (cen=19)", cp.predvar, ref["cp_cen19"]["predvar"], ARITH)
    cmp("01 overall RR (cen=19)", cp.allRRfit, ref["cp_cen19"]["allRRfit"], FIT)
    # the MMT is a point of the 0.1 grid: the same point, up to how each
    # language accumulates seq(by=0.1)
    cmp("01 MMT", cen, ref["mmt"], ARITH)
    cmp("01 MMT via findmin()", dl.findmin(cb, m, by=0.1), ref["mmt"], ARITH)

    pred = dl.crosspred(cb, m, cen=cen, by=1)
    pred1 = dl.crosspred(cb1, m1, cen=cen, by=0.1)
    pred2 = dl.crosspred(cb2, m2, cen=cen, by=0.1)
    cmp("01 overall RR at MMT (ns)", pred.allRRfit, ref["pred"]["allRRfit"], FIT)
    cmp("01 overall RR 95% CI (ns)", np.r_[pred.allRRlow, pred.allRRhigh],
        np.r_[ref["pred"]["allRRlow"], ref["pred"]["allRRhigh"]], FIT)
    cmp("01 overall RR (linear)", pred1.allRRfit, ref["pred_lin"]["allRRfit"], FIT)
    cmp("01 overall RR (double threshold)", pred2.allRRfit, ref["pred_thr"]["allRRfit"], FIT)

    # ---- 02 projected temperature and mortality series ---------------------
    rcp4p5 = pd.read_csv(fetch("lndn_rcp4p5.csv"), parse_dates=["date"])
    rcp8p5 = pd.read_csv(fetch("lndn_rcp8p5.csv"), parse_dates=["date"])

    # mean deaths by day of year (format %j), first 365 days, gaps filled with
    # the mean of the neighbouring days, then repeated along the projection
    deathdoy = obs.groupby(obs.date.dt.dayofyear)["all"].mean().to_numpy()
    deathdoy = np.r_[deathdoy, np.full(max(0, 365 - deathdoy.size), np.nan)][:365]
    while np.isnan(deathdoy).any():
        isna = np.isnan(deathdoy)
        nb = np.column_stack([np.r_[deathdoy[1:], np.nan], np.r_[np.nan, deathdoy[:-1]]])
        deathdoy[isna] = np.nanmean(nb, axis=1)[isna]
    deathproj = np.resize(deathdoy, len(rcp8p5))
    cmp("02 deaths by day of year", deathdoy, ref["deathdoy"], ARITH)

    # ---- 03 bias correction --------------------------------------------------
    o2 = obs[["date", "tmean"]]
    cal = {"rcp4p5": fhempel(o2, rcp4p5), "rcp8p5": fhempel(o2, rcp8p5)}
    rc = ref["calibration"]
    for sc, raw in (("rcp4p5", rcp4p5), ("rcp8p5", rcp8p5)):
        corr = fhempel(o2, raw, output="correction")
        cmp(f"03 correction add+mult ({sc})",
            np.concatenate([np.r_[corr[g]["add"], corr[g]["mult"]] for g in raw.columns[1:]]),
            np.concatenate([np.r_[rc["correction"][sc][g]["add"], rc["correction"][sc][g]["mult"]]
                            for g in raw.columns[1:]]), ARITH)
        ser = rc["series"][sc]
        vals = cal[sc].iloc[:, 1:].to_numpy()
        rows = np.asarray(ser["rows"]) - 1
        cmp(f"03 calibrated series sample ({sc})", vals[rows].T,
            [ser["sample"][g] for g in raw.columns[1:]], ARITH)
        cmp(f"03 calibrated series sums ({sc})",
            np.r_[vals.sum(0), (vals ** 2).sum(0)], np.r_[ser["colsum"], ser["colsumsq"]], ARITH)
    rcp4p5, rcp8p5 = cal["rcp4p5"], cal["rcp8p5"]

    # ---- 04-06 extrapolation, projection, uncertainty ------------------------
    red = dl.crossreduce(cb, m, cen=cen)
    coef, vcov = np.ravel(red.coef), np.asarray(red.vcov, float)
    cmp("04 reduced coef", coef, ref["red"]["coef"], FIT)
    Vr = np.asarray(ref["red"]["vcov"], float)
    cmp("04 reduced vcov", vcov, Vr, FIT, atol=FIT * np.abs(Vr).max())
    deathperiod = deathdoy.sum() * 10
    cmp("05 deaths per decade", deathperiod, ref["deathperiod"], ARITH)

    histperiod = "1971-1989"
    projperiod = [f"{d}0-{str(d)[2]}9" for d in range(199, 210)]
    levels = [histperiod] + projperiod
    seqperiod = np.r_[np.zeros(365 * 19, int), np.repeat(np.arange(1, 12), 365 * 10)]
    temprange = ["tot", "cold", "heat"]
    gcm = {"GFDL-ESM2M": "tmean_gfdl", "HadGEM2-ES": "tmean_hadgem2",
           "IPSL-CM5A-LR": "tmean_ipsl", "MIROC-ESM-CHEM": "tmean_miroc",
           "NorESM1-M": "tmean_noresm1"}
    rcp = {"RCP4.5": rcp4p5, "RCP8.5": rcp8p5}
    nsim = int(ref["nsim"])
    nper = len(levels)
    # ansim[period, range, abs/rel, gcm, rcp, est + nsim]
    ansim = np.full((nper, 3, 2, len(gcm), len(rcp), nsim + 1), np.nan)

    def period_sums(an, keep):
        """tapply(an[keep], factor(seqperiod[keep]), sum), for 1 or many columns."""
        if not np.all(np.isin(np.arange(nper), seqperiod[keep])):
            # R's factor() would drop the empty decade and misalign the result
            raise RuntimeError("a decade has no days on one side of the MMT")
        an = an.reshape(len(seqperiod), -1)
        # the decades are contiguous blocks of rows
        return np.add.reduceat(np.where(keep[:, None], an, 0.0), starts, axis=0)

    Zs = {g: np.asarray(z, float).reshape((nsim, coef.size), order="F")
          for g, z in ref["normals"].items()}
    allrows = np.ones(len(seqperiod), bool)
    starts = np.flatnonzero(np.r_[True, np.diff(seqperiod) != 0])
    for i, (rname, tmeanproj) in enumerate(rcp.items()):
        for j, (gname, col) in enumerate(gcm.items()):
            x = tmeanproj[col].to_numpy(float)
            # (4) extrapolation: basis of the projected series centred on the MMT
            bvar = dl.onebasis(x, **argvar)
            cenvec = dl.onebasis([cen], **argvar)
            bvarcen = bvar.matrix - cenvec.matrix
            indheat = x > cen
            # (5) daily attributable deaths, summed by decade and range
            an = (1 - np.exp(-bvarcen @ coef)) * deathproj
            ansim[:, 0, 0, j, i, 0] = period_sums(an, allrows)[:, 0]
            ansim[:, 1, 0, j, i, 0] = period_sums(an, ~indheat)[:, 0]
            ansim[:, 2, 0, j, i, 0] = period_sums(an, indheat)[:, 0]
            # (6) uncertainty: R's draws, mvrnorm(nsim, coef, vcov) after
            # set.seed(13041975 + j), replayed from its normal deviates
            E = np.asarray(ref["eigvecs"][gname], dtype=float).reshape((coef.size, coef.size), order="F")
            coefsim = dl.simulate_coef(coef, vcov, normals=Zs[gname], eigvecs=E)   # (k, nsim)
            if i == 0:
                # inherits the (fitted) coef and vcov, hence the FIT tolerance
                cmp(f"06 first 10 draws ({gname})", coefsim[:, :10].T,
                    np.asarray(ref["coefsim_head"][gname]).reshape((10, coef.size), order="F"),
                    FIT)
            ansims = (1 - np.exp(-bvarcen @ coefsim)) * deathproj[:, None]
            ansim[:, 0, 0, j, i, 1:] = period_sums(ansims, allrows)
            ansim[:, 1, 0, j, i, 1:] = period_sums(ansims, ~indheat)
            ansim[:, 2, 0, j, i, 1:] = period_sums(ansims, indheat)

    # AN relative to the current decade
    cur = levels.index("2010-19")
    ansim[:, :, 1] = ansim[:, :, 0] - ansim[cur:cur + 1, :, 0]

    cmp("05 AN by decade/range/GCM/RCP", ansim[:, :, 0, :, :, 0], rarr(ref["gcm_est"]), FIT)
    gci = np.quantile(ansim[:, :, 0, :, :, 1:], [0.025, 0.975], axis=-1)
    cmp("06 per-GCM 95% eCI of AN", gci, rarr(ref["gcm_ci"]), FIT)

    # ensemble: mean over GCMs; intervals pooled over GCMs and simulations
    shape = (nper, 3, len(temprange), len(rcp))           # period, est/ci, range, rcp
    anabs, anrel = np.full(shape, np.nan), np.full(shape, np.nan)
    for out, k in ((anabs, 0), (anrel, 1)):
        for i in range(len(rcp)):
            a = ansim[:, :, k, :, i]
            out[:, 0, :, i] = a[..., 0].mean(axis=2)
            pooled = a[..., 1:].reshape(nper, 3, -1)
            out[:, 1, :, i] = np.quantile(pooled, 0.025, axis=2)
            out[:, 2, :, i] = np.quantile(pooled, 0.975, axis=2)
    afabs = anabs / deathperiod * 100
    afrel = anrel / deathperiod * 100
    for nm, a in (("anabs", anabs), ("anrel", anrel), ("afabs", afabs), ("afrel", afrel)):
        R = rarr(ref[nm])
        cmp(f"06 ensemble {nm} est", a[:, 0], R[:, 0], FIT)
        # rel AN in 2010-19 is 0 by construction in both; compare with an atol
        cmp(f"06 ensemble {nm} 95% eCI", a[:, 1:], R[:, 1:], FIT,
            atol=1e-9 if nm in ("anrel", "afrel") else 0.0)

    # ---- Figure 5: one GCM (NorESM1-M), RCP8.5, 2010-19 vs 2090-99 ------------
    dec = rcp8p5.date.dt.year.astype(str).str[:3]
    tmean2010 = rcp8p5.loc[dec == "201", "tmean_noresm1"].to_numpy()
    tmean2090 = rcp8p5.loc[dec == "209", "tmean_noresm1"].to_numpy()
    bvar = dl.onebasis(np.r_[tmean2010, tmean2090], **argvar)
    pred_exp = dl.crosspred(bvar, coef=coef, vcov=vcov, model_link="log", cen=cen, by=0.1)
    death2010 = deathproj[(dec == "201").to_numpy()]
    death2090 = deathproj[(dec == "209").to_numpy()]
    cenvec = dl.onebasis([cen], **argvar).matrix
    anhist = (1 - np.exp(-(dl.onebasis(tmean2010, **argvar).matrix - cenvec) @ coef)) * death2010
    anproj = (1 - np.exp(-(dl.onebasis(tmean2090, **argvar).matrix - cenvec) @ coef)) * death2090
    f5 = ref["fig5"]
    cmp("Fig 5 extrapolated curve grid", pred_exp.predvar, f5["predvar"], ARITH)
    cmp("Fig 5 extrapolated RR and CI",
        np.r_[pred_exp.allRRfit, pred_exp.allRRlow, pred_exp.allRRhigh],
        np.r_[f5["allRRfit"], f5["allRRlow"], f5["allRRhigh"]], FIT)
    cmp("Fig 5 AN 2010-19 and 2090-99", [anhist.sum(), anproj.sum()],
        [f5["anhist"], f5["anproj"]], FIT)

    # ---- 07 age-specific curves and adaptation -------------------------------
    agemmt, agecoef = [], []
    for outcome in ("all_0_64", "all_65_74", "all_75_84", "all_85plus"):
        ma = fit(outcome, cb)
        cpa = dl.crosspred(cb, ma, cen=varcen, by=0.1)
        agemmt.append(float(cpa.predvar[np.argmin(cpa.allRRfit)]))
        agecoef.append(cbpart(ma)[0])
    cmp("07 age-specific cross-basis coef", np.concatenate(agecoef),
        np.concatenate(ref["age"]["coef"]), FIT)
    cmp("07 age-specific MMT", agemmt, ref["age"]["mmt"], ARITH)
    ind = pred.predvar >= cen
    cmp("07 heat curve with 30% adaptation", np.exp(pred.allfit[ind] * 0.7), ref["adapt"], FIT)

    # ---- summary ---------------------------------------------------------------
    cmp.report()
    print(f"\nLondon, {obs.date.dt.year.min()}-{obs.date.dt.year.max()}: "
          f"MMT {cen:.1f} C; age groups {', '.join(f'{v:.1f}' for v in agemmt)} C")
    print("change in attributable fraction (%), 2010-19 to 2090-99, ensemble of 5 GCMs:")
    # the paper reports RCP8.5 only (section "Projection and Quantification")
    published = {("RCP8.5", "heat"): "[paper: 4.0 (0.7, 6.8)]",
                 ("RCP8.5", "cold"): "[paper: -3.3 (-4.3, -1.9)]"}
    for i, rname in enumerate(rcp):
        for rg in ("heat", "cold"):
            e, lo, hi = afrel[levels.index("2090-99"), :, temprange.index(rg), i]
            print(f"  {rname} {rg:4s}  {e:6.2f}  ({lo:.2f} to {hi:.2f})  "
                  f"{published.get((rname, rg), '')}")
    v = ref["versions"]
    print(f"R reference: {v['R']}, dlnm {v['dlnm']}, MASS {v['MASS']}; "
          f"authors' commit {v['authors_commit'][:8]}; dlnmpy {dl.__version__}; "
          f"{time.time() - t0:.0f} s")
    print(f"largest relative difference {max(r[3] for r in cmp.rows):.1e} "
          "(see the table: the largest are on near-zero covariance entries)")
    print(f"vicedo2019: {cmp.n} comparisons with R, max abs difference "
          f"{cmp.maxabs:.1e}, {cmp.fail} failures")
    return 1 if cmp.fail else 0


if __name__ == "__main__":
    sys.exit(main())
