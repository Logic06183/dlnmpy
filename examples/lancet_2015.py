"""Reproduce Gasparrini et al. (2015) Lancet, England & Wales subset.

    Gasparrini A, Guo Y, Hashizume M, et al. Mortality risk attributable to
    high and low ambient temperature: a multicountry observational study.
    The Lancet 2015;386(9991):369-375.

The author's own code and data are public at
https://github.com/gasparrini/2015_gasparrini_Lancet_Rcodedata; this script
is a line-by-line translation of its five stages, using the 10 regions of
England and Wales (1993-2006, 7,573,716 deaths):

    1. a DLNM per region, reduced to the overall cumulative exposure-response
    2. multivariate meta-regression of the reduced coefficients, and BLUPs
    3. the minimum-mortality temperature per region from its BLUP
    4. attributable deaths, forward perspective, with empirical CIs
    5. the country-level table

Running it needs the region data, ~7 MB, downloaded once from the repository
above (it is not redistributed here). Reference values produced by the R
code, and the agreement of every intermediate, are in docs/validation.md.
"""

import io
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

import dlnmpy as dl

DATA_URL = ("https://raw.githubusercontent.com/gasparrini/"
            "2015_gasparrini_Lancet_Rcodedata/master/regEngWales.csv")
CACHE = Path(__file__).with_name("regEngWales.csv")

# specification of the published analysis
VARFUN, VARDEGREE, VARPER = "bs", 2, [0.10, 0.75, 0.90]
LAG, LAGNK, DFSEAS, NSIM = 21, 3, 8, 1000


def load() -> pd.DataFrame:
    if not CACHE.exists():
        print(f"downloading {DATA_URL} ...")
        with urllib.request.urlopen(DATA_URL) as r:          # noqa: S310
            CACHE.write_bytes(r.read())
    d = pd.read_csv(io.BytesIO(CACHE.read_bytes()), index_col=0)
    d["date"] = pd.to_datetime(d["date"])
    return d


def main() -> None:
    d = load()
    regions = list(dict.fromkeys(d.regnames))
    citynames = ["North East", "North West", "Yorkshire & Humber", "East Midlands",
                 "West Midlands", "East", "London", "South East", "South West", "Wales"]
    cities = (pd.DataFrame({"city": regions, "cityname": citynames})
              .sort_values("cityname").reset_index(drop=True))
    dlist = {c: d[d.regnames == c] for c in cities.city}

    # ---- stage 1: region-specific DLNMs, reduced to overall cumulative ------
    coefs, vcovs, keep, meta_x = [], [], {}, []
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
        keep[city] = (cb, tmean, knots)
        meta_x.append((tmean.mean(), tmean.max() - tmean.min()))

    # ---- stage 2: meta-regression on mean and range of temperature ---------
    Y, S = np.vstack(coefs), np.stack(vcovs)
    avgtmean = np.array([m[0] for m in meta_x])
    rangetmean = np.array([m[1] for m in meta_x])
    mv = dl.mixmeta(Y, S, np.column_stack([np.ones(len(Y)), avgtmean, rangetmean]),
                    method="reml", bscov="unstr")
    blups = np.asarray(mv.blup(se=True)["blup"], float)
    blupvcov = mv.blup(se=True)["vcov"]

    # ---- stage 3: minimum-mortality temperature from each BLUP -------------
    minperc, mintemp = [], []
    for i, city in enumerate(cities.city):
        _, tmean, knots = keep[city]
        predvar = dl.quantile7(tmean, np.arange(1, 100) / 100)
        bvar = dl.onebasis(predvar, VARFUN, knots=knots, degree=VARDEGREE,
                           **{"Boundary.knots": [tmean.min(), tmean.max()]})
        p = int(np.argmin(bvar.matrix @ blups[i])) + 1
        minperc.append(p)
        mintemp.append(float(dl.quantile7(tmean, [p / 100])[0]))

    # ---- stage 4: attributable deaths, forward perspective -----------------
    an = np.zeros((len(cities), 3))
    ansim = np.zeros((len(cities), 3, NSIM))
    totdeath = np.zeros(len(cities))
    for i, city in enumerate(cities.city):
        data = dlist[city]
        cb, tmean, _ = keep[city]
        cases = data.death.to_numpy(float)
        cen = mintemp[i]
        for j, rng in enumerate([None, (-100, cen), (cen, 100)]):
            kw = dict(coef=blups[i], vcov=blupvcov[i], type="an", dir="forw",
                      cen=cen, range=rng, name="cb")
            an[i, j] = dl.attrdl(tmean, cb, cases, **kw)
            ansim[i, j] = dl.attrdl(tmean, cb, cases, sim=True, nsim=NSIM,
                                    seed=1000 + i, **kw)
        totdeath[i] = cases.sum()

    # ---- stage 5: the country table ---------------------------------------
    tot = totdeath.sum()
    af = an.sum(axis=0) / tot * 100
    sim_tot = ansim.sum(axis=0)
    aflow = np.quantile(sim_tot, 0.025, axis=1) / tot * 100
    afhigh = np.quantile(sim_tot, 0.975, axis=1) / tot * 100

    print(f"\nEngland & Wales, {d.year.min()}-{d.year.max()}, {int(tot):,} deaths")
    print(f"minimum-mortality percentile (median over regions): {np.median(minperc):.1f}"
          "     [published code: 89.5]")
    print("\nattributable fraction (%)        estimate      95% empirical CI"
          "        [published code]")
    # the point estimates are deterministic and match R exactly; the interval is
    # simulated, so it differs from R's in the third digit by Monte Carlo error
    for j, nm, ref in ((0, "total", "8.94 (8.42-9.43)"), (1, "cold", "8.63 (8.13-9.12)"),
                       (2, "heat", "0.31 (0.27-0.34)")):
        print(f"  {nm:6s}                         {af[j]:5.2f}        "
              f"{aflow[j]:5.2f} to {afhigh[j]:5.2f}          [{ref}]")

    print("\nper region:")
    out = pd.DataFrame({"region": cities.cityname, "deaths": totdeath.astype(int),
                        "MMT_pct": minperc, "MMT_C": np.round(mintemp, 1),
                        "AF_total_%": np.round(an[:, 0] / totdeath * 100, 2),
                        "AF_cold_%": np.round(an[:, 1] / totdeath * 100, 2),
                        "AF_heat_%": np.round(an[:, 2] / totdeath * 100, 2)})
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
