"""Penalised DLNM (Gasparrini et al. 2017): P-spline cross-basis with smoothing
parameters selected by REML, the Python counterpart of
gam(death ~ cb + ..., family=quasipoisson, paraPen=list(cb=cbPen(cb)), method="REML").
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import dlnmpy as dl  # noqa: E402

dl.plot.set_style()
chicago = dl.datasets.chicago_nmmaps()

# cubic regression spline for temperature, P-spline for lags: generous df,
# the penalties decide the smoothness (see the note in docs/penalized.md on
# why "ps" on a predictor with sparse tails can over-fit the extremes)
cb = dl.crossbasis(chicago.temp, lag=21, argvar={"fun": "cr", "df": 8}, arglag={"fun": "ps", "df": 6})
pen = dl.cbpen(cb)                       # Svar, Slag, ranks
ns_time = dl.onebasis(chicago.time, "ns", df=7 * 14)
X = dl.design_matrix(chicago, ("cb", cb), ("ns_time", ns_time), intercept=False)
fit = dl.fit_pgam("death ~ " + " + ".join(X.columns) + " + C(dow)", chicago.join(X),
                  penalties={"cb": pen}, family="quasipoisson", method="reml")
print(fit.summary())
print(f"edf of the cross-basis: {fit.edf_by('cb_'):.1f} of {cb.ncol} coefficients\n")

pred = dl.crosspred(cb, fit, cen=21, by=1, name="cb")
fig = plt.figure(figsize=(10, 3.8))
ax = fig.add_subplot(121)
pred.plot("overall", ax=ax, xlab="Temperature (C)", ylab="RR", title="Overall cumulative (REML-penalised)")
ax = fig.add_subplot(122)
pred.plot("contour", ax=ax, xlab="Temperature (C)", ylab="Lag (days)", zlab="RR", title="Exposure-lag-response surface")
fig.tight_layout()
out = Path(__file__).parent / "figures" / "penalized.png"
out.parent.mkdir(exist_ok=True)
fig.savefig(out, dpi=130)
print("figure written to", out)
