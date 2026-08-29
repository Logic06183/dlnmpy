"""Matplotlib plots for ``CrossPred`` and ``CrossReduce`` objects.

The R package offers "overall", "slices", "contour" and "3d" plot types.
The same are provided here through :func:`plot_crosspred`; each returns the
matplotlib ``Axes`` so the figure can be customised further.
"""

from __future__ import annotations

import numpy as np

__all__ = ["plot_crosspred", "plot_crossreduce"]


def _mpl():
    import matplotlib.pyplot as plt
    return plt


def _ci(ax, x, low, high, ci: str, **kw):
    if ci == "area":
        ax.fill_between(x, low, high, color=kw.pop("color", "0.85"), alpha=kw.pop("alpha", 1.0),
                        linewidth=0, **kw)
    elif ci == "lines":
        ax.plot(x, low, color=kw.pop("color", "0.5"), linestyle=kw.pop("linestyle", "--"), **kw)
        ax.plot(x, high, color="0.5", linestyle="--")
    elif ci == "bars":
        ax.vlines(x, low, high, color=kw.pop("color", "0.5"), **kw)


def plot_crosspred(pred, ptype=None, var=None, lag=None, ci: str = "area",
                   ci_level=None, cumul: bool = False, exp=None, ax=None,
                   xlab=None, ylab=None, zlab=None, title=None, **kwargs):
    """Plot a ``CrossPred`` object.

    Parameters
    ----------
    ptype : {"overall", "slices", "contour", "3d"}, optional
        Defaults to "slices" if ``var``/``lag`` given, "overall" for
        unlagged bases, otherwise "3d".
    var, lag : float or list, optional
        Values at which to cut slices (lag-response at ``var``,
        exposure-response at ``lag``).
    ci : {"area", "lines", "bars", "n"}
    exp : bool, optional
        Exponentiate (default: yes for log/logit links).
    xlab, ylab, zlab : str, optional
        Axis labels, with the same meaning as in R: for "overall" and
        "slices", ``ylab`` is the outcome; for "contour" and "3d", ``ylab``
        is the lag axis and ``zlab`` the outcome (colour bar or z axis).
    """
    plt = _mpl()
    if ptype is None:
        if var is not None or lag is not None:
            ptype = "slices"
        elif pred.lag[1] - pred.lag[0] == 0:
            ptype = "overall"
        else:
            ptype = "3d"
    ci_level = pred.ci_level if ci_level is None else ci_level
    from scipy.stats import norm
    z = norm.ppf(1 - (1 - ci_level) / 2)
    e = pred.is_exp if exp is None else bool(exp)
    f = np.exp if e else (lambda a: a)
    noeff = 1.0 if e else 0.0

    matfit, matse, lags = pred.matfit, pred.matse, pred.predlag
    if cumul:
        if pred.cumfit is None:
            raise ValueError("cumulative outcomes require cumul=True in crosspred()")
        matfit, matse, lags = pred.cumfit, pred.cumse, np.arange(pred.lag[0], pred.lag[1] + 1)

    if ptype == "overall":
        ax = ax or plt.gca()
        _ci(ax, pred.predvar, f(pred.allfit - z * pred.allse), f(pred.allfit + z * pred.allse), ci)
        ax.plot(pred.predvar, f(pred.allfit), **{"color": "C3", **kwargs})
        ax.axhline(noeff, color="0.3", linewidth=0.8)
        ax.set_xlabel(xlab or "Var"); ax.set_ylabel(ylab or "Outcome")
        if title: ax.set_title(title)
        return ax

    if ptype == "slices":
        vars_ = [] if var is None else list(np.atleast_1d(var))
        lags_ = [] if lag is None else list(np.atleast_1d(lag))
        if not vars_ and not lags_:
            raise ValueError("at least 'var' or 'lag' must be provided when ptype='slices'")
        n = len(vars_) + len(lags_)
        if ax is None:
            ncols = int(bool(vars_)) + int(bool(lags_))
            nrows = max(len(vars_), len(lags_))
            fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.2 * nrows), squeeze=False)
            axes_lag = [axes[i, 0] for i in range(len(lags_))] if lags_ else []
            axes_var = [axes[i, ncols - 1] for i in range(len(vars_))] if vars_ else []
            for row in axes:
                for a in row:
                    a.set_visible(False)
        else:
            if n > 1:
                raise ValueError("pass ax only for a single slice")
            axes_lag = [ax] if lags_ else []
            axes_var = [ax] if vars_ else []
        for a, l in zip(axes_lag, lags_):
            a.set_visible(True)
            j = pred._lag_index(l)
            _ci(a, pred.predvar, f(matfit[:, j] - z * matse[:, j]), f(matfit[:, j] + z * matse[:, j]), ci)
            a.plot(pred.predvar, f(matfit[:, j]), **{"color": "C3", **kwargs})
            a.axhline(noeff, color="0.3", linewidth=0.8)
            a.set_xlabel(xlab or "Var"); a.set_ylabel(ylab or "Outcome")
            a.set_title(title or f"Lag = {l:g}")
        for a, v in zip(axes_var, vars_):
            a.set_visible(True)
            i = pred._var_index(v)
            _ci(a, lags, f(matfit[i] - z * matse[i]), f(matfit[i] + z * matse[i]), ci)
            a.plot(lags, f(matfit[i]), **{"color": "C3", **kwargs})
            a.axhline(noeff, color="0.3", linewidth=0.8)
            a.set_xlabel(xlab or "Lag"); a.set_ylabel(ylab or "Outcome")
            a.set_title(title or f"Var = {v:g}")
        return (axes_lag + axes_var)[0] if n == 1 else axes_lag + axes_var

    if ptype == "contour":
        if pred.lag[1] == 0:
            raise ValueError("contour plot not conceivable for unlagged associations")
        ax = ax or plt.gca()
        zmat = f(matfit)
        levels = kwargs.pop("levels", 20)
        from matplotlib.colors import TwoSlopeNorm
        vmin, vmax = float(zmat.min()), float(zmat.max())
        norm_ = TwoSlopeNorm(vcenter=noeff, vmin=min(vmin, noeff - 1e-9), vmax=max(vmax, noeff + 1e-9))
        cf = ax.contourf(pred.predvar, lags, zmat.T, levels=levels, cmap=kwargs.pop("cmap", "RdBu_r"), norm=norm_)
        plt.colorbar(cf, ax=ax, label=zlab or "Outcome")
        ax.set_xlabel(xlab or "Var"); ax.set_ylabel(ylab or "Lag")
        if title: ax.set_title(title)
        return ax

    if ptype == "3d":
        if pred.lag[1] - pred.lag[0] == 0:
            raise ValueError("3D plot not conceivable for unlagged associations")
        if ax is None:
            fig = plt.figure(figsize=(7, 5.5))
            ax = fig.add_subplot(111, projection="3d")
        X, Y = np.meshgrid(pred.predvar, lags, indexing="ij")
        phi, theta = kwargs.pop("phi", 30), kwargs.pop("theta", 210)
        ax.plot_surface(X, Y, f(matfit), cmap=kwargs.pop("cmap", "viridis"),
                        linewidth=0.2, edgecolor="k", antialiased=True, **kwargs)
        ax.view_init(elev=phi, azim=theta)
        ax.set_xlabel(xlab or "Var"); ax.set_ylabel(ylab or "Lag"); ax.set_zlabel(zlab or "Outcome")
        if title: ax.set_title(title)
        return ax

    raise ValueError("ptype must be one of 'overall', 'slices', 'contour', '3d'")


def plot_crossreduce(red, ci: str = "area", ci_level=None, exp=None, ax=None,
                     xlab=None, ylab=None, title=None, **kwargs):
    """Plot a ``CrossReduce`` object."""
    plt = _mpl()
    from scipy.stats import norm
    ci_level = red.ci_level if ci_level is None else ci_level
    z = norm.ppf(1 - (1 - ci_level) / 2)
    e = red.is_exp if exp is None else bool(exp)
    f = np.exp if e else (lambda a: a)
    ax = ax or plt.gca()
    x = red.x
    _ci(ax, x, f(red.fit - z * red.se), f(red.fit + z * red.se), ci)
    ax.plot(x, f(red.fit), **{"color": "C3", **kwargs})
    ax.axhline(1.0 if e else 0.0, color="0.3", linewidth=0.8)
    ax.set_xlabel(xlab or ("Lag" if red.type == "var" else "Var"))
    ax.set_ylabel(ylab or "Outcome")
    if title: ax.set_title(title)
    return ax
