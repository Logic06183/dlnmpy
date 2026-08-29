"""Two-stage multi-location analysis: location-specific DLNMs, reduced to
the exposure-response dimension, pooled by multivariate meta-analysis
(Gasparrini, Armstrong & Kenward 2012; Gasparrini et al. 2015 Lancet).

Uses 12 locations simulated by dlnmpy.datasets.simulate_cities (no external
data needed).
"""

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import dlnmpy as dl  # noqa: E402

sim = dl.datasets.simulate_cities(n_cities=12, n_days=2000, seed=1)
knots, bk, cen = [8, 15, 22], [-5, 40], 18

# ---- stage 1: one DLNM per location, reduced to the predictor space ----------
reduced, cbs, mean_temp = [], [], []
for city, d in sim.groupby("city"):
    d = d.reset_index(drop=True)
    # same knots everywhere so that the reduced coefficients are comparable
    cb = dl.crossbasis(d.tmean, lag=10, argvar={"fun": "ns", "knots": knots, "boundary_knots": bk},
                       arglag={"fun": "ns", "df": 3})
    ns_time = dl.onebasis(d.time, "ns", df=18)
    X = dl.design_matrix(d, ("cb", cb), ("ns_time", ns_time), intercept=False)
    fit = dl.fit_glm("y ~ " + " + ".join(X.columns), d.join(X), family="quasipoisson")
    reduced.append(dl.crossreduce(cb, fit, cen=cen, name="cb"))
    cbs.append(cb)
    mean_temp.append(d.tmean.mean())

# ---- stage 2: multivariate meta-analysis and meta-regression -----------------
y, S = dl.stack_reduced(reduced)
mm = dl.mixmeta(y, S, method="reml")
print(mm.summary(), "\n")

Xmeta = np.column_stack([np.ones(len(y)), mean_temp])
mr = dl.mixmeta(y, S, X=Xmeta, method="reml")
q0, q1 = mm.qtest(), mr.qtest()
print(f"meta-regression on mean temperature: I2 {q0['I2'][0]:.1f}% -> {q1['I2'][0]:.1f}%, "
      f"AIC {mm.aic:.2f} -> {mr.aic:.2f}\n")

# ---- pooled curve, BLUPs and plot ---------------------------------------------
grid = np.arange(-5, 41)
pooled = dl.predict_reduced(cbs[0], mm.coef_vec, mm.vcov, at=grid, cen=cen)
blups = mm.blup(se=True)

fig, ax = plt.subplots(figsize=(6.5, 4.2))
for i in range(len(reduced)):
    first = dl.predict_reduced(cbs[i], reduced[i].coef, reduced[i].vcov, at=grid, cen=cen)
    ax.plot(grid, first.allRRfit, color="0.75", lw=0.8, label="first stage" if i == 0 else None)
    bl = dl.predict_reduced(cbs[i], blups["blup"][i], blups["vcov"][i], at=grid, cen=cen)
    ax.plot(grid, bl.allRRfit, color="C0", lw=0.8, alpha=0.7, label="BLUP" if i == 0 else None)
ax.fill_between(grid, pooled.allRRlow, pooled.allRRhigh, color="C3", alpha=0.15)
ax.plot(grid, pooled.allRRfit, color="C3", lw=2, label="pooled (REML)")
ax.axhline(1, color="0.3", lw=0.8)
ax.set_xlabel("Temperature"); ax.set_ylabel("RR (ref. 18)"); ax.set_ylim(0.8, 4)
ax.legend(frameon=False)
ax.set_title("12 locations: first-stage, BLUP and pooled exposure-response")
fig.tight_layout()
out = Path(__file__).parent / "figures" / "two_stage.png"
out.parent.mkdir(exist_ok=True)
fig.savefig(out, dpi=130)
print("figure written to", out)
