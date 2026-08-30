"""Reproduce Gasparrini & Armstrong (2013), the DLNM meta-analysis paper.

    Gasparrini A, Armstrong B. Reducing and meta-analysing estimates from
    distributed lag non-linear models. BMC Medical Research Methodology
    2013;13:1.

The author's code and data are public at
https://github.com/gasparrini/2013_gasparrini_BMCmrm_Rcodedata; this is a
translation of it, on the same 10 regions of England and Wales used by
examples/lancet_2015.py but with a different question: how to reduce a
bi-dimensional exposure-lag-response surface to a set of one-dimensional
summaries, and how to pool those across locations.

It exercises the parts examples/lancet_2015.py does not:

    * crossreduce in two of its three forms -- the overall cumulative curve
      and the predictor-specific (type="var") lag-response at a given
      temperature
    * three competing lag specifications compared by QAIC
    * multivariate meta-analysis of each reduction, random and fixed effects
    * Cochran's Q and I-squared
    * meta-regression on latitude
    * prediction from the pooled coefficients through a onebasis

Needs the region data, ~7 MB, downloaded once from the repository above.
Reference values from the R code are in docs/validation.md.
"""

import io
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

import dlnmpy as dl

DATA_URL = ("https://raw.githubusercontent.com/gasparrini/"
            "2013_gasparrini_BMCmrm_Rcodedata/master/regEngWales.csv")
CACHE = Path(__file__).with_name("regEngWales.csv")
LAT = [54.84815, 53.58832, 53.72352, 52.85539, 52.53304,
       52.03734, 51.50583, 51.24213, 51.05361, 52.02615]


def load() -> pd.DataFrame:
    if not CACHE.exists():
        print(f"downloading {DATA_URL} ...")
        with urllib.request.urlopen(DATA_URL) as r:
            CACHE.write_bytes(r.read())
    return pd.read_csv(io.BytesIO(CACHE.read_bytes()), index_col=0)


def main() -> None:
    d = load()
    regions = list(dict.fromkeys(d.regnames))
    dlist = {r: d[d.regnames == r] for r in regions}

    # the same basis everywhere, so the reduced coefficients are comparable
    ranges = np.array([[s.tmean.min(), s.tmean.max()] for s in dlist.values()])
    bound = ranges.mean(axis=0)
    argvar = {"fun": "bs", "degree": 2, "boundary_knots": bound,
              "knots": dl.equalknots(bound, fun="bs", degree=2, df=4)}
    arglag = {"fun": "ns", "knots": dl.logknots(21, df=5, intercept=True)}
    arglag_const = {"fun": "strata", "df": 1}     # a constant over the lag window

    red = {k: ([], []) for k in ("all", "all2", "all3", "hot", "cold")}
    qaic = dict.fromkeys(("m1", "m2", "m3"), 0.0)
    for region in regions:
        sub = dlist[region]
        tmean = sub.tmean.to_numpy(float)
        cb = dl.crossbasis(tmean, lag=(0, 21), argvar=argvar, arglag=arglag)
        cb2 = dl.crossbasis(tmean, lag=(0, 3), argvar=argvar, arglag=arglag_const)
        cb3 = dl.crossbasis(tmean, lag=(0, 21), argvar=argvar, arglag=arglag_const)
        cb2.matrix[:21, :] = np.nan               # compare the models on the same days
        seas = dl.onebasis(sub.time.to_numpy(float), "ns", df=10 * 14)

        fits = {}
        for name, basis in (("cb", cb), ("cb2", cb2), ("cb3", cb3)):
            X = dl.design_matrix(sub, (name, basis), ("seas", seas), intercept=False)
            fits[name] = dl.fit_glm("death ~ " + " + ".join(X.columns) + " + C(dow)",
                                    sub.join(X), family="quasipoisson")
        for key, val in zip(("m1", "m2", "m3"), ("cb", "cb2", "cb3")):
            qaic[key] += dl.uncertainty.qaic(fits[val])

        for key, basis, name, kw in (
                ("all", cb, "cb", {}), ("all2", cb2, "cb2", {}), ("all3", cb3, "cb3", {}),
                ("hot", cb, "cb", {"type": "var", "value": 22}),
                ("cold", cb, "cb", {"type": "var", "value": 0})):
            cr = dl.crossreduce(basis, fits[name], cen=17, name=name, **kw)
            red[key][0].append(np.ravel(cr.coef))
            red[key][1].append(np.asarray(cr.vcov, float))
        realised_arglag = dict(cb.arglag)   # R: attr(cb, "arglag"), with the intercept

    # ---- second stage ------------------------------------------------------
    pooled = {k: dl.mixmeta(np.vstack(y), np.stack(S), method="reml")
              for k, (y, S) in red.items()}

    a = dict(argvar)
    bvar = dl.onebasis(np.arange(bound[0], bound[1] + 1e-9, 0.1), a.pop("fun"), **a)
    cp = dl.crosspred(bvar, coef=np.ravel(pooled["all"].coef), vcov=pooled["all"].vcov,
                      model_link="log", from_=bound[0], to=bound[1], by=0.1, cen=17)

    mmt = float(cp.predvar[np.argmin(cp.allRRfit)])

    def at(t):
        return int(np.argmin(np.abs(cp.predvar - t)))

    def i2(m):
        q = m.qtest()
        Q, df = np.atleast_1d(q["Q"])[0], np.atleast_1d(q["df"])[0]
        return float((Q - df) / Q * 100)

    print(f"\nEngland & Wales, {len(regions)} regions, {d.year.min()}-{d.year.max()}")
    print(f"  minimum-mortality temperature        {mmt:.1f} C            [paper: 17.1]")
    for t, ref in ((22, "1.101 (1.079-1.124)"), (0, "1.309 (1.245-1.376)")):
        i = at(t)
        print(f"  pooled RR at {t:2d} C                    "
              f"{cp.allRRfit[i]:.3f} ({cp.allRRlow[i]:.3f}-{cp.allRRhigh[i]:.3f})"
              f"   [paper: {ref}]")
    print(f"  I2  overall {i2(pooled['all']):.1f}%   at 22 C {i2(pooled['hot']):.1f}%"
          f"   at 0 C {i2(pooled['cold']):.1f}%      [paper: 63.8 / 16.4 / 63.5]")
    print("\n  QAIC, summed over regions (lower is better):")
    for key, label in (("m1", "B-spline of lag 0-21"), ("m2", "constant of lag 0-3"),
                       ("m3", "constant of lag 0-21")):
        print(f"    {label:24s} {qaic[key]:12.1f}")
    print(f"    -> the flexible lag model wins by "
          f"{min(qaic['m2'], qaic['m3']) - qaic['m1']:.0f} QAIC units")

    lat = np.column_stack([np.ones(len(regions)), LAT])
    mr = dl.mixmeta(np.vstack(red["all"][0]), np.stack(red["all"][1]), lat, method="reml")
    print(f"\n  meta-regression on latitude: {mr.coef.shape[0]} predictors x "
          f"{mr.coef.shape[1]} outcomes, logLik {mr.loglik:.4f}")

    # the predictor-specific reductions are lag-response curves at one temperature
    al = realised_arglag
    lagfun = al.pop("fun")
    al.pop("cen", None)
    blag = dl.onebasis(np.arange(0, 211) / 10, lagfun, **al)
    for key, t in (("hot", 22), ("cold", 0)):
        p = dl.crosspred(blag, coef=np.ravel(pooled[key].coef), vcov=pooled[key].vcov,
                         model_link="log", at=np.arange(0, 211) / 10, cen=17)
        peak = int(np.argmax(p.allRRfit))
        print(f"  pooled lag-response at {t:2d} C: peak RR {p.allRRfit[peak]:.4f} "
              f"at lag {p.predvar[peak]:.1f}")


if __name__ == "__main__":
    main()
