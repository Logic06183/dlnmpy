"""Temperature-attributable mortality in Chicago: the standard workflow.

Fit a DLNM, find the minimum mortality temperature with its uncertainty,
then compute deaths attributable to cold and heat with empirical intervals
(Gasparrini & Leone 2014; Tobías et al. 2017; Gasparrini et al. 2015).
"""

import numpy as np

import dlnmpy as dl

chicago = dl.datasets.chicago_nmmaps()

# cross-basis: quadratic B-spline for temperature (knots at the 10th, 75th and
# 90th percentiles), natural cubic spline for lags 0-21 with knots on the log scale
cb = dl.crossbasis(chicago.temp, lag=21,
                   argvar={"fun": "bs", "degree": 2, "knots": dl.quantile7(chicago.temp, [0.10, 0.75, 0.90])},
                   arglag={"knots": dl.logknots(21, 3)})
ns_time = dl.onebasis(chicago.time, "ns", df=7 * 14)
X = dl.design_matrix(chicago, ("cb", cb), ("ns_time", ns_time), intercept=False)
model = dl.fit_glm("death ~ " + " + ".join(X.columns) + " + C(dow)", chicago.join(X), family="quasipoisson")

# minimum mortality temperature, searched between the 1st and 99th percentiles
res = dl.mmt(cb, model, x=chicago.temp, nsim=2000, seed=1, name="cb")
print(res)

# attributable deaths and fractions, reference = MMT
tab = dl.attr_table(chicago.temp, cb, chicago.death, model, cen=res.mmt, nsim=2000, seed=1, name="cb")
with np.printoptions(precision=1):
    print(tab.assign(af=lambda d: 100 * d.af, af_low=lambda d: 100 * d.af_low, af_high=lambda d: 100 * d.af_high)
          .round({"an": 0, "an_low": 0, "an_high": 0, "af": 2, "af_low": 2, "af_high": 2})
          .to_string(index=False))

# the same from the forward perspective (exposure on day t contributing to the next 21 days)
af_forw = dl.attrdl(chicago.temp, cb, chicago.death, model, dir="forw", cen=res.mmt, name="cb")
print(f"\nforward-perspective total AF: {100 * af_forw:.2f}%")
