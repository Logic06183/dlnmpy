"""Attributable risk and minimum-risk exposure from a fitted DLNM.

These functions are not part of the R package ``dlnm``; they port the
reference R functions that Antonio Gasparrini publishes alongside his
papers and that most temperature-mortality studies copy:

- :func:`attrdl` ports ``attrdl.R`` (Gasparrini & Leone, *BMC Medical
  Research Methodology* 2014, "Attributable risk from distributed lag
  models"), computing attributable fractions and numbers from the backward
  or forward perspective, with empirical confidence intervals by parametric
  simulation of the coefficients.
- :func:`findmin` ports ``findmin.R`` (Tobías, Armstrong & Gasparrini,
  *Epidemiology* 2017, "Investigating uncertainty in the minimum mortality
  temperature"), locating the minimum of the overall cumulative
  exposure-response curve and its simulation-based uncertainty.

Built on top of them, :func:`mmt` and :func:`attr_table` give the summaries
reported in practice: the minimum mortality temperature with its interval
and percentile, and attributable numbers/fractions for total, cold and heat
(and, optionally, extreme and moderate ranges), all derived from one set of
simulated coefficients so that the components are mutually consistent.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._rcompat import quantile7
from .core import CrossBasis, OneBasis, onebasis
from .lag import lag_matrix, seqlag
from .predict import _resolve_coef, mkXpred, mkat

__all__ = ["attrdl", "findmin", "mmt", "attr_table", "simulate_coef", "MMTResult"]


# ----------------------------------------------------------------------------
def simulate_coef(coef, vcov, nsim: int = 5000, seed=None) -> np.ndarray:
    """Draw ``nsim`` coefficient vectors from N(coef, vcov).

    Uses the eigen-decomposition square root of ``vcov`` exactly as the R
    reference code does, so that a matrix of standard normals produced in R
    gives identical draws. Returns an array of shape ``(k, nsim)``.
    """
    coef = np.asarray(coef, dtype=float).ravel()
    vcov = np.asarray(vcov, dtype=float)
    rng = np.random.default_rng(seed)
    Z = rng.standard_normal((nsim, coef.size))
    return _coefsim_from_normals(coef, vcov, Z)


def _coefsim_from_normals(coef, vcov, Z) -> np.ndarray:
    """``coef + V^{1/2} Z'`` with the eigen square root (R: eigen$vectors %*%
    diag(sqrt(values)) %*% t(X)). ``Z`` is (nsim, k)."""
    w, U = np.linalg.eigh(vcov)
    # R's eigen() returns values in decreasing order; match that ordering so
    # that a given normal matrix maps to the same draws
    order = np.argsort(w)[::-1]
    w, U = w[order], U[:, order]
    w = np.maximum(w, 0)
    return coef[:, None] + U @ (np.sqrt(w)[:, None] * Z.T)


def _xpred_all(basis, at, cen):
    """Sum over lags of the (centred) prediction basis for the rows of ``at``."""
    predlag = seqlag(basis.lag) if isinstance(basis, CrossBasis) else np.array([0.0])
    n = at.shape[0]
    predvar = np.arange(n, dtype=float)
    X = mkXpred(basis, at, predvar, predlag, cen)
    Xall = np.zeros((n, X.shape[1]))
    for i in range(len(predlag)):
        Xall += X[i * n:(i + 1) * n]
    return Xall


# ----------------------------------------------------------------------------
def attrdl(x, basis, cases, model=None, coef=None, vcov=None, model_link=None,
           type: str = "af", dir: str = "back", tot: bool = True, cen=None,
           range=None, sim: bool = False, nsim: int = 5000, seed=None,
           coefsim=None, name: str = "cb"):
    """Attributable fraction or number from a DLNM (port of ``attrdl.R``).

    Parameters
    ----------
    x : array_like
        Exposure series, or (``dir="back"`` only) a matrix of lagged exposures.
    basis : CrossBasis
        The cross-basis built from ``x``.
    cases : array_like
        Daily cases (vector), or (``dir="forw"`` only) a matrix of future
        cases.
    model, coef, vcov, model_link
        As in :func:`dlnmpy.crosspred`. The link must be log or logit.
        ``coef``/``vcov`` may also be *reduced* estimates for the predictor
        space (from :func:`crossreduce`), in which case only ``dir="forw"``
        is allowed.
    type : {"af", "an"}
        Attributable fraction or attributable number.
    dir : {"back", "forw"}
        Backward perspective (risk at time t attributed to past exposures)
        or forward perspective (exposure at time t contributes to future
        cases).
    tot : bool
        Return the total over the series (True) or the daily contributions.
    cen : float
        Counterfactual reference exposure (e.g. the minimum mortality
        temperature). Falls back to the centring stored in ``basis``.
    range : (float, float), optional
        Only exposures within this range contribute; others are set to
        ``cen`` (null risk). Use it for cold/heat or extreme ranges.
    sim : bool
        With ``tot=True``, return ``nsim`` simulated totals (for empirical
        intervals) instead of the point estimate.
    seed, coefsim
        Random seed, or an explicit ``(k, nsim)`` matrix of simulated
        coefficients (e.g. from :func:`simulate_coef`) to share across calls.
    """
    if type not in ("an", "af"):
        raise ValueError("'type' must be 'an' or 'af'")
    out = _attr_core(x, basis, cases, model, coef, vcov, model_link, dir, tot, cen, range,
                     sim, nsim, seed, coefsim, name)
    if sim:
        return out["ansim"] if type == "an" else out["afsim"]
    return out["an"] if type == "an" else out["af"]


def _attr_core(x, basis, cases, model, coef, vcov, model_link, dir, tot, cen, range,
               sim, nsim, seed, coefsim, name):
    if dir not in ("back", "forw"):
        raise ValueError("'dir' must be 'back' or 'forw'")
    if not isinstance(basis, CrossBasis):
        raise TypeError("'basis' must be a CrossBasis")
    if cen is None:
        cen = basis.argvar.get("cen")
    if cen is None:
        raise ValueError("'cen' must be provided")
    cen = float(cen)

    x = np.asarray(x, dtype=float)
    if range is not None:
        x = x.copy()
        x[(x < range[0]) | (x > range[1])] = cen

    lag = basis.lag
    nlag = int(lag[1] - lag[0]) + 1
    if x.ndim == 1:
        at = lag_matrix(x, seqlag(lag)) if dir == "back" else np.tile(x[:, None], (1, nlag))
    else:
        if dir == "forw":
            raise ValueError("'x' must be a vector when dir='forw'")
        if x.shape[1] != nlag:
            raise ValueError("dimension of 'x' not compatible with 'basis'")
        at = x

    cases = np.asarray(cases, dtype=float)
    if cases.shape[0] != at.shape[0]:
        raise ValueError("'x' and 'cases' not consistent")
    if cases.ndim == 2 and cases.shape[1] > 1:
        if dir == "back":
            raise ValueError("'cases' must be a vector if dir='back'")
        if cases.shape[1] != nlag:
            raise ValueError("dimension of 'cases' not compatible")
        den = float(np.nansum(np.nanmean(cases, axis=1)))
        cases = np.mean(cases, axis=1)  # NaN if any lag missing, as in R rowMeans
    else:
        cases = cases.ravel()
        den = float(np.nansum(cases))
        if dir == "forw":
            fut = lag_matrix(cases, -seqlag(lag))
            cases = np.mean(fut, axis=1)

    # coefficients (full cross-basis or reduced predictor-space estimates)
    if model is not None:
        coef, vcov, model_link, _ = _resolve_coef(basis, model, None, None, model_link, name, "cb")
        if model_link not in ("log", "logit"):
            raise ValueError("'model' must have a log or logit link function")
    else:
        if coef is None or vcov is None:
            raise ValueError("At least 'model' or 'coef'-'vcov' must be provided")
        coef = np.asarray(coef, dtype=float).ravel()
        vcov = np.atleast_2d(np.asarray(vcov, dtype=float))
    typebasis = "one" if coef.size != basis.ncol else "cb"

    if typebasis == "cb":
        Xall = _xpred_all(basis, at, cen)
    else:
        if dir == "back":
            raise ValueError("only dir='forw' allowed for reduced estimates")
        av = dict(basis.argvar)
        fun = av.pop("fun")
        av.pop("cen", None)
        ob = onebasis(x, fun, **av)
        Xall = mkXpred(ob, x, x, np.array([0.0]), cen)
    if coef.size != Xall.shape[1]:
        raise ValueError("arguments 'basis' do not match 'model' or 'coef'-'vcov'")
    if vcov.shape != (coef.size, coef.size):
        raise ValueError("arguments 'coef' and 'vcov' do not match")

    af = 1 - np.exp(-(Xall @ coef))
    an = af * cases
    if tot:
        isna = np.isnan(an)
        af = float(np.sum(an[~isna]) / np.sum(cases[~isna]))
        an = af * den

    out = {"af": af, "an": an}
    if sim and not tot:
        raise ValueError("simulation samples only returned for tot=True")
    if sim:
        if coefsim is None:
            coefsim = simulate_coef(coef, vcov, nsim, seed)
        coefsim = np.asarray(coefsim, dtype=float)
        ani = (1 - np.exp(-(Xall @ coefsim))) * cases[:, None]
        ok = ~np.isnan(ani[:, 0])  # same missing pattern for every draw
        out["afsim"] = np.sum(ani[ok], axis=0) / np.sum(cases[ok])
        out["ansim"] = out["afsim"] * den
    return out


# ----------------------------------------------------------------------------
def findmin(basis, model=None, coef=None, vcov=None, at=None, from_=None, to=None,
            by=None, sim: bool = False, nsim: int = 5000, seed=None, coefsim=None,
            name: str = "cb"):
    """Minimum of the overall cumulative exposure-response curve (port of
    ``findmin.R``).

    Returns the predictor value at which the curve is lowest over the grid
    ``at`` (or ``from_``/``to``/``by``, default step 0.1 over the range). With
    ``sim=True`` returns ``nsim`` minima from simulated coefficients, whose
    quantiles give an empirical confidence interval.
    """
    if not isinstance(basis, (CrossBasis, OneBasis)):
        raise TypeError("'basis' must be a CrossBasis or OneBasis")
    one = isinstance(basis, OneBasis)
    lag = np.array([0, 0]) if one else basis.lag
    kind = "one" if one else "cb"
    coef, vcov, _, _ = _resolve_coef(basis, model, coef, vcov, None, name, kind)
    if by is None and at is None:
        by = 0.1
    at = mkat(at, from_, to, by, basis.range, lag, 1)
    if at.ndim == 2:
        raise ValueError("'at' must be a vector")
    predlag = seqlag(lag)
    X = mkXpred(basis, at, at, predlag, None)
    n = at.size
    Xall = np.zeros((n, X.shape[1]))
    for i in np.arange(len(predlag)):
        Xall += X[i * n:(i + 1) * n]
    pred = Xall @ coef
    if not sim:
        return float(at[int(np.argmin(pred))])
    if coefsim is None:
        coefsim = simulate_coef(coef, vcov, nsim, seed)
    predsim = Xall @ np.asarray(coefsim, dtype=float)
    return at[np.argmin(predsim, axis=0)]


# ----------------------------------------------------------------------------
@dataclass
class MMTResult:
    """Minimum mortality (risk) exposure with simulation-based uncertainty."""

    mmt: float
    low: float
    high: float
    percentile: float | None
    ci_level: float
    sims: np.ndarray

    def __repr__(self):
        pct = "" if self.percentile is None else f", percentile={self.percentile:.1f}"
        return (f"MMTResult(mmt={self.mmt:.2f}, {int(self.ci_level * 100)}% CI "
                f"({self.low:.2f}, {self.high:.2f}){pct})")


def mmt(basis, model=None, coef=None, vcov=None, x=None, from_=None, to=None, by=0.1,
        percentiles=(1, 99), ci_level: float = 0.95, nsim: int = 5000, seed=None,
        coefsim=None, name: str = "cb") -> MMTResult:
    """Minimum mortality temperature (or minimum-risk exposure) with an
    empirical confidence interval.

    The search is restricted to the ``percentiles`` of the exposure
    distribution ``x`` (default 1st to 99th, as is common practice, to avoid
    spurious minima in the sparse tails), unless ``from_``/``to`` are given.
    ``x`` is also used to report the percentile of the MMT.
    """
    if x is not None and (from_ is None or to is None):
        x = np.asarray(x, dtype=float)
        lo, hi = quantile7(x, np.asarray(percentiles) / 100)
        from_ = lo if from_ is None else from_
        to = hi if to is None else to
    kw = dict(model=model, coef=coef, vcov=vcov, from_=from_, to=to, by=by, name=name)
    point = findmin(basis, **kw)
    sims = findmin(basis, sim=True, nsim=nsim, seed=seed, coefsim=coefsim, **kw)
    a = (1 - ci_level) / 2
    low, high = np.quantile(sims, [a, 1 - a])
    pct = None
    if x is not None:
        xx = x[~np.isnan(x)]
        pct = float(100 * np.mean(xx <= point))
    return MMTResult(float(point), float(low), float(high), pct, ci_level, sims)


# ----------------------------------------------------------------------------
def attr_table(x, basis, cases, model=None, coef=None, vcov=None, model_link=None,
               cen=None, dir: str = "back", extreme_percentiles=(2.5, 97.5),
               ci_level: float = 0.95, nsim: int = 5000, seed=None, name: str = "cb"):
    """Attributable numbers and fractions for total, cold and heat (and
    extreme/moderate components), with empirical confidence intervals.

    ``cen`` is the reference (typically the MMT from :func:`mmt`). Cold is
    the range below ``cen``, heat above; "extreme" ranges lie beyond the
    given percentiles of ``x`` and "moderate" between them and ``cen``.
    All components use the same simulated coefficients. Note that cold and
    heat do not sum exactly to the total: the fraction is ``1 - exp(-sum)``
    over lags, which is not additive (the R reference behaves the same).
    Returns a pandas DataFrame with columns ``component``, ``range``,
    ``an``, ``an_low``, ``an_high``, ``af``, ``af_low``, ``af_high``.
    """
    import pandas as pd

    x = np.asarray(x, dtype=float)
    if cen is None:
        cen = basis.argvar.get("cen")
    if cen is None:
        raise ValueError("'cen' must be provided (e.g. the MMT)")
    cen = float(cen)
    xmin, xmax = float(np.nanmin(x)), float(np.nanmax(x))
    ranges = {"total": (xmin, xmax), "cold": (xmin, cen), "heat": (cen, xmax)}
    if extreme_percentiles is not None:
        plo, phi = quantile7(x, np.asarray(extreme_percentiles) / 100)
        ranges.update({
            "extreme cold": (xmin, float(plo)), "moderate cold": (float(plo), cen),
            "moderate heat": (cen, float(phi)), "extreme heat": (float(phi), xmax),
        })
    # one set of simulated coefficients shared by every component
    if model is not None:
        c, v, model_link, _ = _resolve_coef(basis, model, None, None, model_link, name, "cb")
    else:
        c, v = np.asarray(coef, float).ravel(), np.atleast_2d(np.asarray(vcov, float))
    coefsim = simulate_coef(c, v, nsim, seed)
    a = (1 - ci_level) / 2
    rows = []
    for comp, rng in ranges.items():
        o = _attr_core(x, basis, cases, None, c, v, model_link, dir, True, cen, rng, True, nsim, None, coefsim, name)
        rows.append({"component": comp, "range": rng, "an": o["an"],
                     "an_low": float(np.quantile(o["ansim"], a)), "an_high": float(np.quantile(o["ansim"], 1 - a)),
                     "af": o["af"], "af_low": float(np.quantile(o["afsim"], a)), "af_high": float(np.quantile(o["afsim"], 1 - a))})
    return pd.DataFrame(rows)
