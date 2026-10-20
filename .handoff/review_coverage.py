"""Coverage of dlnmpy's overall cumulative RR intervals under a known DLNM.

Simulates daily counts from a known temperature-lag-response surface (Chicago
temperatures, Poisson counts), fits dl.dlnm() with the
basis that contains the truth, and records how often the 95% interval of the
overall cumulative RR covers the true value at three temperatures.

    python .handoff/review_coverage.py [nrep]
"""
import sys
import warnings

import numpy as np

import dlnmpy as dl

warnings.simplefilter("ignore")
nrep = int(sys.argv[1]) if len(sys.argv) > 1 else 200
chi = dl.datasets.chicago_nmmaps().iloc[:3000].reset_index(drop=True)
x = chi.temp.to_numpy()
lag = 7
argvar = {"fun": "ns", "knots": [5.0, 20.0], "boundary_knots": [-26.0, 33.0]}
arglag = {"fun": "ns", "df": 3}
cb = dl.crossbasis(x, lag=lag, argvar=argvar, arglag=arglag)

# a true surface inside the model space: coefficients chosen by least squares so
# that the overall cumulative curve is 0.0006 (t - 20)^2
grid = np.linspace(-20, 30, 51)
D = np.column_stack([dl.crosspred(cb, coef=np.eye(cb.ncol)[j], vcov=np.eye(cb.ncol), model_link="log",
                                  cen=20, at=grid).allfit for j in range(cb.ncol)])
eta_true = np.linalg.lstsq(D, 0.0006 * (grid - 20) ** 2, rcond=None)[0]
at = np.array([-10.0, 0.0, 28.0])
true_rr = dl.crosspred(cb, coef=eta_true, vcov=np.eye(cb.ncol), model_link="log", cen=20, at=at).allRRfit

lin = np.asarray(cb.matrix) @ eta_true
base = np.log(30) + 0.15 * np.sin(2 * np.pi * chi.time.to_numpy() / 365.25)
ok = np.isfinite(lin)
cover = np.zeros((nrep, at.size), bool)
rng = np.random.default_rng(2026)
for r in range(nrep):
    y = np.full(len(chi), np.nan)
    y[ok] = rng.poisson(np.exp(base[ok] + lin[ok]))
    d = chi.assign(sim=y)
    fit = dl.dlnm(d, outcome="sim", exposure="temp", lag=lag, argvar=argvar, arglag=arglag,
                  time="time", df_per_year=4, family="quasipoisson")
    p = fit.predict(at=at, cen=20)
    cover[r] = (p.allRRlow <= true_rr) & (true_rr <= p.allRRhigh)
cov = cover.mean(axis=0)
se = np.sqrt(cov * (1 - cov) / nrep)
for t, c, s in zip(at, cov, se):
    print(f"temperature {t:5.1f} C: true RR {true_rr[list(at).index(t)]:.3f}, coverage {c:.3f} (MC SE {s:.3f})")
print(f"nrep {nrep}")
