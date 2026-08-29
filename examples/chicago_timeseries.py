"""Reproduction of the examples in the R vignette "Distributed lag linear
and non-linear models for time series data" (dlnmTS) with dlnmpy.

Run:  python examples/chicago_timeseries.py
Figures are written to examples/figures/.
"""

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import dlnmpy as dl  # noqa: E402

out = Path(__file__).parent / "figures"
out.mkdir(exist_ok=True)

chicago = dl.datasets.chicago_nmmaps()


def fit(data, *terms, extra="C(dow)"):
    """Fit death ~ terms + extra with a quasi-Poisson GLM."""
    X = dl.design_matrix(data, *terms, intercept=False)
    rhs = " + ".join(list(X.columns) + [extra])
    return dl.fit_glm("death ~ " + rhs, data.join(X), family="quasipoisson")


# ---------------------------------------------------------------------------
# Example 1: a simple DLM (PM10, linear; temperature, ns with 5 df)
# ---------------------------------------------------------------------------
cb1_pm = dl.crossbasis(chicago.pm10, lag=15, argvar={"fun": "lin"}, arglag={"fun": "poly", "degree": 4})
cb1_temp = dl.crossbasis(chicago.temp, lag=3, argvar={"df": 5}, arglag={"fun": "strata", "breaks": 1})
print(cb1_pm.summary(), "\n")
ns_time = dl.onebasis(chicago.time, "ns", df=7 * 14)

model1 = fit(chicago, ("cb1_pm", cb1_pm), ("cb1_temp", cb1_temp), ("ns_time", ns_time))
pred1_pm = dl.crosspred(cb1_pm, model1, at=np.arange(21), bylag=0.2, cumul=True, name="cb1_pm")

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
pred1_pm.plot("slices", var=10, ax=axes[0], ylab="RR", title="Lag-response curve for a 10-unit increase in PM10")
pred1_pm.plot("slices", var=10, ax=axes[1], cumul=True, ylab="Cumulative RR",
              title="Incremental cumulative effects")
fig.tight_layout()
fig.savefig(out / "example1_slices.png", dpi=130)

i = pred1_pm._var_index(10)
print(f"Example 1: overall RR for 10-unit PM10 increase = {pred1_pm.allRRfit[i]:.4f} "
      f"({pred1_pm.allRRlow[i]:.4f}, {pred1_pm.allRRhigh[i]:.4f})\n")

# ---------------------------------------------------------------------------
# Example 2: seasonal analysis (summer only, series grouped by year)
# ---------------------------------------------------------------------------
seas = chicago[chicago.month.isin([6, 7, 8, 9])].reset_index(drop=True)
cb2_o3 = dl.crossbasis(seas.o3, lag=5, argvar={"fun": "thr", "thr": 40.3}, arglag={"fun": "integer"}, group=seas.year)
cb2_temp = dl.crossbasis(seas.temp, lag=10, argvar={"fun": "thr", "thr": [15, 25]},
                         arglag={"fun": "strata", "breaks": [2, 6]}, group=seas.year)
ns_doy = dl.onebasis(seas.doy, "ns", df=4)
ns_time2 = dl.onebasis(seas.time, "ns", df=3)
model2 = fit(seas, ("cb2_o3", cb2_o3), ("cb2_temp", cb2_temp), ("ns_doy", ns_doy), ("ns_time", ns_time2))
pred2_o3 = dl.crosspred(cb2_o3, model2, at=list(range(66)) + [40.3, 50.3], name="cb2_o3")

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
pred2_o3.plot("slices", var=50.3, ci="bars", ci_level=0.80, ax=axes[0], ylab="RR",
              title="Lag-response, 10-unit increase above threshold (80% CI)")
pred2_o3.plot("overall", ci="lines", ax=axes[1], xlab="Ozone", ylab="RR",
              title="Overall cumulative association for 5 lags")
fig.tight_layout()
fig.savefig(out / "example2_ozone.png", dpi=130)
i = pred2_o3._var_index(50.3)
print(f"Example 2: overall RR for ozone 50.3 vs threshold = {pred2_o3.allRRfit[i]:.4f} "
      f"({pred2_o3.allRRlow[i]:.4f}, {pred2_o3.allRRhigh[i]:.4f})\n")

# ---------------------------------------------------------------------------
# Example 3: a bi-dimensional DLNM (temperature, lag 0-30)
# ---------------------------------------------------------------------------
cb3_pm = dl.crossbasis(chicago.pm10, lag=1, argvar={"fun": "lin"}, arglag={"fun": "strata"})
varknots = dl.equalknots(chicago.temp, fun="bs", df=5, degree=2)
lagknots = dl.logknots(30, 3)
cb3_temp = dl.crossbasis(chicago.temp, lag=30, argvar={"fun": "bs", "knots": varknots}, arglag={"knots": lagknots})
model3 = fit(chicago, ("cb3_pm", cb3_pm), ("cb3_temp", cb3_temp), ("ns_time", ns_time))
pred3_temp = dl.crosspred(cb3_temp, model3, cen=21, by=1, name="cb3_temp")

fig = plt.figure(figsize=(11, 4.2))
ax3d = fig.add_subplot(121, projection="3d")
pred3_temp.plot("3d", ax=ax3d, xlab="Temperature", ylab="RR", theta=200, phi=40, title="3D graph of temperature effect")
ax = fig.add_subplot(122)
pred3_temp.plot("contour", ax=ax, xlab="Temperature", ylab="RR", title="Contour plot")
fig.tight_layout()
fig.savefig(out / "example3_surface.png", dpi=130)

fig, ax = plt.subplots(figsize=(6, 4))
for v, col in zip([-20, 0, 27, 33], ["k", "C0", "C2", "C3"]):
    s = pred3_temp.slice_var(v)
    ax.plot(s.lag, s.fit, color=col, label=f"Temperature = {v}")
ax.axhline(1, color="0.5", lw=0.8)
ax.set_xlabel("Lag"); ax.set_ylabel("RR"); ax.set_ylim(0.95, 1.25); ax.legend()
ax.set_title("Lag-response curves for different temperatures, ref. 21C")
fig.tight_layout()
fig.savefig(out / "example3_slices.png", dpi=130)

axes = pred3_temp.plot("slices", var=[-20, 33], lag=[0, 5])
axes[0].figure.tight_layout()
axes[0].figure.savefig(out / "example3_slices_grid.png", dpi=130)

# ---------------------------------------------------------------------------
# Example 4: reducing a DLNM to one-dimensional summaries
# ---------------------------------------------------------------------------
cb4 = dl.crossbasis(chicago.temp, lag=30, argvar={"fun": "thr", "thr": [10, 25]}, arglag={"knots": lagknots})
model4 = fit(chicago, ("cb4", cb4), ("ns_time", ns_time))
pred4 = dl.crosspred(cb4, model4, by=1, name="cb4")
redall = dl.crossreduce(cb4, model4, name="cb4")
redlag = dl.crossreduce(cb4, model4, type="lag", value=5, name="cb4")
redvar = dl.crossreduce(cb4, model4, type="var", value=33, name="cb4")
print("Example 4: parameters -", "full:", pred4.coef.size, "| overall:", redall.coef.size,
      "| lag-specific:", redlag.coef.size, "| var-specific:", redvar.coef.size)

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
pred4.plot("overall", ax=axes[0], xlab="Temperature", ylab="RR", title="Overall cumulative association")
r = redall.to_dataframe()
axes[0].plot(r["var"], r.fit, color="C0", ls="--", label="Reduced")
axes[0].plot([], [], color="C3", label="Original"); axes[0].legend()
pred4.plot("slices", var=33, ax=axes[1], ylab="RR", title="Predictor-specific association at 33C")
r = redvar.to_dataframe()
axes[1].plot(r["lag"], r.fit, color="C0", ls="--")
fig.tight_layout()
fig.savefig(out / "example4_reduce.png", dpi=130)
print("max |reduced - full| overall:", float(np.max(np.abs(redall.fit - pred4.allfit))))
print(f"\nFigures written to {out}")
