"""Simulation-based (parametric bootstrap) uncertainty and model selection.

The delta-method intervals in :class:`~dlnmpy.predict.CrossPred` are exact
on the scale of the linear predictor and adequate for relative risks. For
derived quantities that are non-linear in the coefficients (a minimum, a
ratio of two curves, an attributable fraction, a difference between two
lags) the standard route is to draw coefficient vectors from
``N(coef, vcov)`` and push each draw through the prediction. These helpers
make that a one-liner for any function of the coefficients and give
empirical intervals, following the approach used for attributable risk
(Gasparrini & Leone 2014) and the MMT (Tobías et al. 2017).

Also here: QAIC for over-dispersed count models (Peng, Dominici & Louis
2006; Gasparrini et al. 2010), used to compare cross-basis specifications.
"""

from __future__ import annotations

import numpy as np
from scipy.special import gammaln

from .attribution import simulate_coef
from .predict import CrossPred, crosspred

__all__ = ["bootstrap", "simulate_pred", "empirical_ci", "qaic", "model_grid"]


def bootstrap(fn, coef, vcov, nsim: int = 1000, seed=None, coefsim=None):
    """Evaluate ``fn(coef_draw)`` for ``nsim`` draws from N(coef, vcov).

    Returns an array of shape ``(nsim, *fn(coef).shape)``. Pass ``coefsim``
    (``(k, nsim)``, e.g. from :func:`dlnmpy.attribution.simulate_coef`) to
    share one set of draws across several quantities.
    """
    coef = np.asarray(coef, dtype=float).ravel()
    if coefsim is None:
        coefsim = simulate_coef(coef, vcov, nsim, seed)
    coefsim = np.asarray(coefsim, dtype=float)
    out = [np.asarray(fn(coefsim[:, j]), dtype=float) for j in range(coefsim.shape[1])]
    return np.stack(out)


def empirical_ci(samples, level: float = 0.95, axis: int = 0):
    """Percentile interval of bootstrap samples: ``(low, high)``."""
    a = (1 - level) / 2
    return np.quantile(samples, a, axis=axis), np.quantile(samples, 1 - a, axis=axis)


def simulate_pred(basis, coef, vcov, nsim: int = 1000, seed=None, model_link=None, **kwargs):
    """Simulated overall cumulative and lag-specific predictions.

    Runs :func:`dlnmpy.crosspred` for ``nsim`` coefficient draws and returns
    a dict with ``allfit`` (``nsim x nvar``) and ``matfit`` (``nsim x nvar x
    nlag``) samples on the linear-predictor scale, the point prediction
    ``pred`` and ``coefsim``. Extra keyword arguments (``at``, ``cen``,
    ``by``, ``lag``, ``bylag`` ...) go to ``crosspred``.
    """
    pred = crosspred(basis, coef=coef, vcov=vcov, model_link=model_link, **kwargs)
    coefsim = simulate_coef(pred.coef, pred.vcov, nsim, seed)
    allfit, matfit = [], []
    for j in range(nsim):
        p = crosspred(basis, coef=coefsim[:, j], vcov=pred.vcov, model_link=model_link, **kwargs)
        allfit.append(p.allfit)
        matfit.append(p.matfit)
    return {"pred": pred, "coefsim": coefsim, "allfit": np.stack(allfit), "matfit": np.stack(matfit)}


def bootstrap_ci(pred: CrossPred, nsim: int = 1000, seed=None, level: float = 0.95, basis=None, **kwargs):
    """Empirical (percentile) intervals for the overall cumulative RR of a
    ``CrossPred``; useful as a check on the delta-method intervals when the
    curve is used on a scale where they are not exact."""
    if basis is None:
        raise ValueError("pass the basis the prediction was made from")
    sim = simulate_pred(basis, pred.coef, pred.vcov, nsim=nsim, seed=seed, model_link=pred.model_link,
                        at=pred.predvar if pred.at_matrix is None else pred.at_matrix, cen=pred.cen,
                        lag=pred.lag, bylag=pred.bylag, **kwargs)
    f = np.exp if pred.is_exp else (lambda a: a)
    low, high = empirical_ci(f(sim["allfit"]), level)
    return {"fit": f(pred.allfit), "low": low, "high": high, "samples": f(sim["allfit"])}


# ----------------------------------------------------------------------------
def qaic(results) -> float:
    """Quasi-AIC for an over-dispersed Poisson GLM fitted with statsmodels:
    ``-2 * loglik(Poisson) + 2 * phi * k`` with ``phi`` the estimated
    dispersion (Gasparrini, Armstrong & Kenward 2010; Peng et al. 2006).
    For a fit with ``scale = 1`` this reduces to the AIC.

    The Poisson log-likelihood is evaluated at the fitted values rather than
    read off ``results.llf``: for a quasi family statsmodels divides the
    log-likelihood by the estimated dispersion, so ``results.llf`` is
    ``loglik / phi``. Using it would scale ``-2 * loglik`` by ``1 / phi``,
    and because phi differs between candidate models the error does not
    cancel in a comparison -- it favours the more over-dispersed (typically
    the under-fitted) model. This is also why R's ``logLik()`` returns NA
    for a quasipoisson glm and the reference implementations compute
    ``sum(dpois(y, fitted, log = TRUE))`` by hand.
    """
    phi = float(getattr(results, "scale", 1.0))
    k = int(np.size(results.params))
    mu = np.asarray(results.fittedvalues, dtype=float)
    y = np.asarray(results.model.endog, dtype=float)
    llf = float(np.sum(y * np.log(mu) - mu - gammaln(y + 1.0)))
    return -2 * llf + 2 * phi * k


def model_grid(specs, fit_fn, criterion=qaic):
    """Fit a set of model specifications and rank them by a criterion.

    ``specs`` is a list of dicts (or any objects) passed one at a time to
    ``fit_fn(spec) -> results``; ``criterion(results)`` (default
    :func:`qaic`) scores each. Returns a DataFrame sorted by the criterion
    with the spec fields as columns, ``criterion`` and ``delta`` (difference
    from the best), and the fitted results in a ``results`` column.
    """
    import pandas as pd
    rows = []
    for spec in specs:
        res = fit_fn(spec)
        row = dict(spec) if isinstance(spec, dict) else {"spec": spec}
        row["criterion"] = float(criterion(res))
        row["results"] = res
        rows.append(row)
    df = pd.DataFrame(rows).sort_values("criterion").reset_index(drop=True)
    df["delta"] = df["criterion"] - df["criterion"].iloc[0]
    return df
