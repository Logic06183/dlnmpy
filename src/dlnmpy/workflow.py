"""One call from a daily data frame to the quantities a paper reports.

The R workflow for a temperature-mortality analysis is forty lines spread
over ``crossbasis``, ``glm``, ``crosspred``, ``findmin.R``, ``attrdl.R`` and
a plotting block, and every study re-types them. :func:`dlnm` collapses
that into one function whose result knows its own data, so the minimum
mortality temperature, the relative risks at chosen percentiles, the
attributable fractions and the standard figure are each one method call::

    fit = dl.dlnm(chicago, outcome="death", exposure="temp", lag=21,
                  argvar={"fun": "ns", "knots": dl.percentile_knots(chicago.temp, [10, 75, 90])},
                  arglag={"fun": "ns", "knots": dl.logknots(21, 3)},
                  time="time", df_per_year=7, dow="dow")
    fit.mmt()             # minimum mortality temperature with its interval and percentile
    fit.rr_at([1, 99])    # RR (95% CI) at the 1st and 99th percentiles, relative to the MMT
    fit.attributable()    # attributable numbers and fractions: total, cold, heat, extreme, moderate
    fit.figure()          # overall curve with the MMT, percentile axis and the exposure distribution

Everything is built from the same validated primitives (``crossbasis``,
``fit_glm``, ``crosspred``, ``mmt``, ``attr_table``), which remain available
for anything the wrapper does not cover.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from ._rcompat import quantile7
from .attribution import MMTResult, attr_table, mmt as _mmt
from .core import CrossBasis, crossbasis, onebasis
from .model import design_matrix, extract_coef_vcov, fit_glm
from .predict import CrossPred, CrossReduce, crosspred, crossreduce

__all__ = ["dlnm", "DLNM", "percentile_knots", "percentile_of"]


# ----------------------------------------------------------------------------
def percentile_knots(x, percentiles=(10, 75, 90)) -> np.ndarray:
    """Knots at the given percentiles of ``x`` (R's ``quantile(type=7)``,
    NaN ignored), the usual placement for the exposure-response basis."""
    return quantile7(np.asarray(x, dtype=float), np.asarray(percentiles, dtype=float) / 100)


def percentile_of(x, values) -> np.ndarray:
    """Empirical percentile (0-100) of each of ``values`` within ``x``."""
    xx = np.asarray(x, dtype=float)
    xx = np.sort(xx[~np.isnan(xx)])
    v = np.atleast_1d(np.asarray(values, dtype=float))
    return 100 * np.searchsorted(xx, v, side="right") / xx.size


def _years(t: np.ndarray) -> float:
    return (np.nanmax(t) - np.nanmin(t) + 1) / 365.25


# ----------------------------------------------------------------------------
@dataclass
class DLNM:
    """A fitted distributed lag non-linear model bound to its data.

    Attributes
    ----------
    cb : CrossBasis
    model : fitted results (statsmodels GLM or :class:`PenalizedGLMResults`)
    x, y : ndarray
        Exposure and outcome series in data order.
    formula : str
        The formula that was fitted.
    """

    cb: CrossBasis
    model: Any
    x: np.ndarray
    y: np.ndarray
    formula: str
    name: str = "cb"
    exposure: str = "x"
    outcome: str = "y"
    family: str = "quasipoisson"
    data: pd.DataFrame | None = field(default=None, repr=False)
    _mmt: MMTResult | None = field(default=None, repr=False)

    # --- basics ---------------------------------------------------------------
    @property
    def lag(self) -> tuple:
        return tuple(int(v) for v in self.cb.lag)

    @property
    def coef(self) -> np.ndarray:
        """Coefficients of the cross-basis columns."""
        return extract_coef_vcov(self.model, self.name, "cb", self.cb.ncol)[0]

    @property
    def vcov(self) -> np.ndarray:
        """Covariance matrix of the cross-basis coefficients."""
        return extract_coef_vcov(self.model, self.name, "cb", self.cb.ncol)[1]

    def quantile(self, percentiles) -> np.ndarray:
        """Exposure values at the given percentiles (0-100)."""
        return quantile7(self.x, np.atleast_1d(np.asarray(percentiles, dtype=float)) / 100)

    def percentile(self, values) -> np.ndarray:
        """Percentile (0-100) of exposure values."""
        return percentile_of(self.x, values)

    # --- centring -------------------------------------------------------------
    def mmt(self, percentiles=(1, 99), ci_level: float = 0.95, nsim: int = 5000, seed=None,
            by: float = 0.1) -> MMTResult:
        """Minimum mortality (minimum risk) exposure, searched between the
        given percentiles, with a simulation interval. Cached: later calls
        with default arguments return the same object, and it is the default
        reference for :meth:`predict`, :meth:`rr_at` and :meth:`attributable`."""
        res = _mmt(self.cb, self.model, x=self.x, percentiles=percentiles, ci_level=ci_level,
                   nsim=nsim, seed=seed, by=by, name=self.name)
        self._mmt = res
        return res

    def _cen(self, cen):
        if cen is not None:
            return cen
        if self._mmt is None:
            self.mmt(seed=0)
        return self._mmt.mmt

    # --- predictions ----------------------------------------------------------
    def predict(self, at=None, percentiles=None, cen=None, by=None, from_=None, to=None,
                lag=None, bylag: float = 1, cumul: bool = False, ci_level: float = 0.95) -> CrossPred:
        """Exposure-lag-response predictions (:func:`dlnmpy.crosspred`).

        ``percentiles`` gives the grid as percentiles of the exposure instead
        of values. ``cen`` defaults to the MMT (computed on first use);
        ``cen=False`` leaves the prediction uncentred.
        """
        if percentiles is not None:
            at = self.quantile(percentiles)
        if cen is None:
            cen = self._cen(None)
        return crosspred(self.cb, self.model, at=at, from_=from_, to=to, by=by, lag=lag, bylag=bylag,
                         cen=cen, cumul=cumul, ci_level=ci_level, name=self.name)

    def reduce(self, type: str = "overall", value=None, cen=None, **kwargs) -> CrossReduce:
        """One-dimensional reduction (:func:`dlnmpy.crossreduce`), e.g. the
        overall cumulative curve for a second-stage meta-analysis."""
        return crossreduce(self.cb, self.model, type=type, value=value, cen=self._cen(cen),
                           name=self.name, **kwargs)

    def rr_at(self, percentiles=(1, 2.5, 10, 90, 97.5, 99), cen=None, ci_level: float = 0.95) -> pd.DataFrame:
        """Overall cumulative relative risk at exposure percentiles, relative
        to ``cen`` (default the MMT): the numbers a results table reports."""
        pct = np.atleast_1d(np.asarray(percentiles, dtype=float))
        at = self.quantile(pct)
        p = self.predict(at=at, cen=cen, ci_level=ci_level)
        rows = []
        for q, v in zip(pct, at):
            i = p._var_index(v)
            rows.append({"percentile": q, self.exposure: v, "rr": p.allRRfit[i],
                         "low": p.allRRlow[i], "high": p.allRRhigh[i]})
        return pd.DataFrame(rows)

    def attributable(self, cen=None, dir: str = "back", extreme_percentiles=(2.5, 97.5),
                     ci_level: float = 0.95, nsim: int = 5000, seed=None) -> pd.DataFrame:
        """Attributable numbers and fractions for total, cold and heat, and
        their extreme and moderate parts (:func:`dlnmpy.attr_table`), with
        ``cen`` (default the MMT) as the counterfactual."""
        return attr_table(self.x, self.cb, self.y, self.model, cen=self._cen(cen), dir=dir,
                          extreme_percentiles=extreme_percentiles, ci_level=ci_level, nsim=nsim,
                          seed=seed, name=self.name)

    def qaic(self) -> float:
        from .uncertainty import qaic
        return qaic(self.model)

    # --- figures --------------------------------------------------------------
    def figure(self, cen=None, by=None, percentiles=(1, 10, 50, 90, 99), hist: bool = True,
               xlab=None, ylab="RR", title=None, ax=None, **kwargs):
        """The overall cumulative exposure-response curve with the MMT, a
        percentile axis and the exposure distribution (see
        :func:`dlnmpy.plot.plot_overall_risk`)."""
        from .plot import plot_overall_risk
        m = self._mmt if cen is None else None
        pred = self.predict(cen=cen, by=by)
        return plot_overall_risk(pred, x=self.x, mmt=m if m is not None else pred.cen, percentiles=percentiles,
                                 hist=hist, xlab=xlab or self.exposure, ylab=ylab, title=title, ax=ax, **kwargs)

    def plot(self, ptype=None, cen=None, by=None, **kwargs):
        """Any of the ``CrossPred`` plots (``overall``, ``slices``, ``contour``, ``3d``)."""
        return self.predict(cen=cen, by=by).plot(ptype, **kwargs)

    def summary_figure(self, cen=None, by=None, **kwargs):
        from .plot import summary_figure
        kwargs.setdefault("xlab", self.exposure)
        return summary_figure(self.predict(cen=cen, by=by), **kwargs)

    # --- report ---------------------------------------------------------------
    def summary(self, nsim: int = 2000, seed=0) -> str:
        m = self._mmt or self.mmt(nsim=nsim, seed=seed)
        rr = self.rr_at([1, 99])
        lines = [f"DLNM: {self.outcome} ~ {self.exposure}, lag {self.lag[0]}-{self.lag[1]}, {self.family}",
                 f"observations: {self.x.size} (fitted: {self._nobs()})   cross-basis df: {self.cb.df[0]} x {self.cb.df[1]}",
                 f"minimum-risk {self.exposure}: {m.mmt:.1f} ({int(100 * m.ci_level)}% CI {m.low:.1f} to {m.high:.1f}), "
                 f"percentile {m.percentile:.1f}",
                 f"RR at 1st percentile ({rr.iloc[0, 1]:.1f}): {rr.rr[0]:.3f} ({rr.low[0]:.3f}, {rr.high[0]:.3f})",
                 f"RR at 99th percentile ({rr.iloc[1, 1]:.1f}): {rr.rr[1]:.3f} ({rr.low[1]:.3f}, {rr.high[1]:.3f})"]
        try:
            lines.append(f"QAIC: {self.qaic():.1f}")
        except Exception:
            pass
        return "\n".join(lines)

    def _nobs(self) -> int:
        n = getattr(self.model, "nobs", None)
        return int(n) if n is not None else int(np.sum(~np.isnan(self.cb.matrix).any(axis=1)))

    def __repr__(self):
        return (f"DLNM({self.outcome} ~ {self.exposure}, lag={self.lag}, df={self.cb.df}, "
                f"family={self.family!r})")


# ----------------------------------------------------------------------------
def dlnm(data: pd.DataFrame, outcome: str, exposure: str, lag, argvar=None, arglag=None,
         time=None, df_per_year: float = 7, dow=None, controls=(), family: str = "quasipoisson",
         offset=None, group=None, penalised: bool | None = None, method: str = "reml",
         name: str = "cb", **fit_kwargs) -> DLNM:
    """Fit a DLNM from a data frame in one call.

    Parameters
    ----------
    data : DataFrame
        Daily (or otherwise regularly spaced) series, one row per unit of time.
    outcome, exposure : str
        Column names of the outcome counts and the exposure.
    lag : int or (int, int)
    argvar, arglag : dict
        Basis specifications as for :func:`dlnmpy.crossbasis`. Defaults: a
        natural cubic spline with knots at the 10th, 75th and 90th percentiles
        of the exposure; a natural cubic spline of the lags with knots
        equally spaced on the log scale (``logknots(lag, 3)``), or a single
        stratum for ``lag <= 1``.
    time : str, optional
        Column with the time index (days). A natural cubic spline with
        ``df_per_year`` df per year is added to control for season and trend.
    dow : str, optional
        Column with the day of the week, added as a factor.
    controls : sequence of str
        Extra formula terms on columns of ``data`` (``"C(holiday)"``, ``"pm10"``).
    family : str
        ``"quasipoisson"`` (default), ``"poisson"``, ``"gaussian"``, ``"binomial"``.
    offset : str, optional
        Column with a log-offset (e.g. log population).
    group : str, optional
        Column identifying separate series (lags are not computed across groups).
    penalised : bool, optional
        Fit by penalised likelihood with REML/ML smoothing
        (:func:`dlnmpy.fit_pgam`). Defaults to True when either basis is
        ``ps`` or ``cr``.
    name : str
        Prefix of the cross-basis columns in the design.
    **fit_kwargs
        Passed to :func:`dlnmpy.fit_glm` (or ``fit_pgam``).
    """
    from .knots import logknots
    from .lag import mklag

    x = np.asarray(data[exposure], dtype=float)
    y = np.asarray(data[outcome], dtype=float)
    lag_ = mklag(lag)
    if argvar is None:
        argvar = {"fun": "ns", "knots": percentile_knots(x, (10, 75, 90))}
    if arglag is None:
        arglag = {"fun": "ns", "knots": logknots(lag_, 3)} if lag_[1] - lag_[0] > 1 else {}
    cb = crossbasis(data[exposure], lag=lag_, argvar=dict(argvar), arglag=dict(arglag),
                    group=None if group is None else data[group])

    terms = [(name, cb)]
    rhs = []
    if time is not None:
        t = np.asarray(data[time], dtype=float)
        df_time = int(round(df_per_year * _years(t)))
        terms.append((f"ns_{time}", onebasis(data[time], "ns", df=max(df_time, 1))))
    if dow is not None:
        rhs.append(f"C({dow})")
    rhs.extend(controls)
    X = design_matrix(data, *terms, intercept=False)
    frame = data.join(X)
    formula = f"{outcome} ~ " + " + ".join(list(X.columns) + rhs)
    if offset is not None:
        fit_kwargs["offset"] = np.asarray(data[offset], dtype=float)

    funs = (cb.argvar.get("fun"), cb.arglag.get("fun"))
    if penalised is None:
        penalised = any(f in ("ps", "cr") for f in funs)
    if penalised:
        from .penalized import fit_pgam
        from .penalty import cbpen
        model = fit_pgam(formula, frame, {name: cbpen(cb)}, family=family, method=method, **fit_kwargs)
    else:
        model = fit_glm(formula, frame, family=family, **fit_kwargs)
    return DLNM(cb=cb, model=model, x=x, y=y, formula=formula, name=name, exposure=exposure,
                outcome=outcome, family=family, data=data)
