"""The figures in the README: one call from the Chicago data to the
temperature-mortality results, then the standard figures in both themes.

Run:  python examples/gallery.py
Figures are written to examples/figures/.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import dlnmpy as dl  # noqa: E402

out = Path(__file__).parent / "figures"
out.mkdir(exist_ok=True)

chicago = dl.datasets.chicago_nmmaps()

fit = dl.dlnm(chicago, outcome="death", exposure="temp", lag=21,
              argvar={"fun": "bs", "degree": 2, "knots": dl.percentile_knots(chicago.temp, [10, 75, 90])},
              arglag={"fun": "ns", "knots": dl.logknots(21, 3)},
              time="time", df_per_year=7, dow="dow")
print(fit)
print(fit.mmt(nsim=2000, seed=1))
print(fit.rr_at([1, 2.5, 97.5, 99]).round(3).to_string(index=False))
tab = fit.attributable(nsim=2000, seed=1)
print(tab.drop(columns="range").round(4).to_string(index=False))
print(fit.summary())

for theme in ("colour", "journal"):
    dl.plot.set_theme(theme)
    ax = fit.figure(xlab="Temperature (°C)")
    ax.figure.savefig(out / f"overall_risk_{theme}.png", dpi=160)
    ax = dl.plot.plot_attributable(tab)
    ax.figure.savefig(out / f"attributable_{theme}.png", dpi=160)
    ax = fit.plot("3d", by=1, xlab="Temperature (°C)", zlab="RR")
    ax.figure.savefig(out / f"surface_{theme}.png", dpi=160)
    fig, _ = fit.summary_figure(xlab="Temperature (°C)")
    fig.savefig(out / f"summary_{theme}.png", dpi=160)
    plt.close("all")

print(f"\nFigures written to {out}")
