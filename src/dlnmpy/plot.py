"""Matplotlib plots for ``CrossPred`` and ``CrossReduce`` objects.

The R package offers "overall", "slices", "contour" and "3d" plot types; the
same are provided here through :func:`plot_crosspred` (or ``pred.plot``),
plus :func:`overlay_slices` (several curves on one axis, as with R's
``lines()``) and :func:`summary_figure` (the standard four-panel figure).
Each function returns the matplotlib ``Axes`` so figures can be customised.

Two themes are available through :func:`set_theme`: ``"journal"`` (default:
greyscale, print-safe, neutral sans-serif, thin marks, hairline gridlines,
no top/right spines) and ``"colour"`` (a validated colour-blind-safe palette:
one blue for estimates, a blue/red diverging map with a neutral grey midpoint
for surfaces, categorical hues in a fixed order for overlays). Every colour
can still be overridden through keyword arguments.
"""

from __future__ import annotations

import numpy as np

__all__ = ["plot_crosspred", "plot_crossreduce", "overlay_slices", "summary_figure",
           "set_style", "set_theme", "THEMES"]

# --- themes -------------------------------------------------------------------
THEMES = {
    "journal": {
        "line": "#111111", "band": "#dcdcdc", "ref": "#666666", "grid": "#e8e8e8",
        "series": ["#111111", "#5a5a5a", "#9a9a9a", "#c4c4c4"],
        "series_ls": ["-", "--", "-.", ":"],
        # greyscale cannot carry a sign, so surfaces use one monotone ramp
        # (light = low, dark = high) with the null-effect contour drawn on top
        "cmap_div": ["#f7f7f7", "#c6c6c6", "#8c8c8c", "#4d4d4d", "#1a1a1a"], "cmap_mode": "sequential",
        "cmap_seq": "Greys", "surface": "Greys",
    },
    "colour": {
        "line": "#2a78d6", "band": "#cde2fb", "ref": "#666666", "grid": "#ebebe8",
        "series": ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"],
        "series_ls": ["-"] * 8,
        "cmap_div": ["#0d366b", "#2a78d6", "#9ec5f4", "#f0efec", "#f3a9a5", "#e34948", "#7a1f1f"], "cmap_mode": "diverging",
        "cmap_seq": "Blues", "surface": "Blues",
    },
}
_theme = dict(THEMES["journal"])


def set_theme(name: str = "journal", **overrides):
    """Select the plotting theme (``"journal"`` or ``"colour"``) and apply the
    matching rcParams. Keyword overrides replace individual entries
    (``line``, ``band``, ``ref``, ``grid``, ``series``, ``cmap_div`` ...)."""
    global _theme
    if name not in THEMES:
        raise ValueError(f"theme must be one of {sorted(THEMES)}")
    _theme = {**THEMES[name], **overrides}
    set_style()
    return dict(_theme)


def set_style(font_size: float = 9.0):
    """Apply a journal-figure style (neutral sans-serif, thin axes, hairline
    solid gridlines, no top/right spines) to matplotlib's rcParams."""
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.size": font_size, "axes.titlesize": font_size + 1,
        "axes.titleweight": "normal", "axes.labelsize": font_size, "xtick.labelsize": font_size - 1,
        "ytick.labelsize": font_size - 1, "legend.fontsize": font_size - 1, "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 3, "ytick.major.size": 3,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
        "grid.color": _theme["grid"], "grid.linewidth": 0.6, "grid.linestyle": "-", "axes.axisbelow": True,
        "legend.frameon": False, "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
        "lines.linewidth": 1.3, "axes.edgecolor": "#444444", "xtick.color": "#444444", "ytick.color": "#444444",
        "axes.labelcolor": "#111111", "text.color": "#111111",
    })


def _mpl():
    import matplotlib.pyplot as plt
    return plt


def _finish(ax, xlab, ylab, title, noeff):
    ax.axhline(noeff, color=_theme["ref"], linewidth=0.6, linestyle=(0, (2, 2)))
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(True, color=_theme["grid"], linewidth=0.6)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    if title:
        ax.set_title(title, loc="left")


def _ci(ax, x, low, high, ci: str, **kw):
    if ci == "area":
        ax.fill_between(x, low, high, color=kw.pop("color", _theme["band"]), alpha=kw.pop("alpha", 1.0),
                        linewidth=0, **kw)
    elif ci == "lines":
        c = kw.pop("color", _theme["ref"])
        ax.plot(x, low, color=c, linestyle=kw.pop("linestyle", "--"), linewidth=0.8, **kw)
        ax.plot(x, high, color=c, linestyle="--", linewidth=0.8)
    elif ci == "bars":
        ax.vlines(x, low, high, color=kw.pop("color", _theme["ref"]), linewidth=0.8, **kw)


def _diverging_cmap(noeff, vmin, vmax, name="dlnm_div"):
    """Two-hue diverging colormap with a neutral midpoint pinned at ``noeff``."""
    from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
    cmap = LinearSegmentedColormap.from_list(name, _theme["cmap_div"])
    if _theme.get("cmap_mode") == "sequential":
        return cmap, Normalize(vmin=vmin, vmax=vmax)
    norm = TwoSlopeNorm(vcenter=noeff, vmin=min(vmin, noeff - 1e-9), vmax=max(vmax, noeff + 1e-9))
    return cmap, norm


def _scale(pred, exp):
    e = pred.is_exp if exp is None else bool(exp)
    return (np.exp if e else (lambda a: a)), (1.0 if e else 0.0), e


def _z(pred, ci_level):
    from scipy.stats import norm
    ci_level = pred.ci_level if ci_level is None else ci_level
    return norm.ppf(1 - (1 - ci_level) / 2)


# ------------------------------------------------------------------------------
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
        exposure-response at ``lag``); several values give a panel grid.
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
    z = _z(pred, ci_level)
    f, noeff, e = _scale(pred, exp)
    matfit, matse, lags = pred.matfit, pred.matse, pred.predlag
    if cumul:
        if pred.cumfit is None:
            raise ValueError("cumulative outcomes require cumul=True in crosspred()")
        matfit, matse, lags = pred.cumfit, pred.cumse, np.arange(pred.lag[0], pred.lag[1] + 1)
    outcome = ylab or ("RR" if e else "Effect")

    def _lag_col(value):
        """Column of ``matfit`` for lag ``value``. The cumulative outcomes have
        one column per integer lag, not one per ``pred.predlag``: resolving them
        against predlag silently returns the wrong lag whenever bylag != 1 (R
        sets bylag to 1 under cumul and indexes cumfit by name)."""
        if not cumul:
            return pred._lag_index(value)
        idx = np.nonzero(np.isclose(lags, value))[0]
        if idx.size == 0:
            raise ValueError(f"'lag'={value:g} must be one of the integer lags "
                             f"{lags[0]:g}..{lags[-1]:g} for cumulative outcomes")
        return int(idx[0])

    if ptype == "overall":
        ax = ax or plt.gca()
        _ci(ax, pred.predvar, f(pred.allfit - z * pred.allse), f(pred.allfit + z * pred.allse), ci)
        ax.plot(pred.predvar, f(pred.allfit), **{"color": _theme["line"], **kwargs})
        _finish(ax, xlab or "Var", outcome, title, noeff)
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
            fig, axes = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 2.6 * nrows), squeeze=False)
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
        for a, lg in zip(axes_lag, lags_):
            a.set_visible(True)
            j = _lag_col(lg)
            _ci(a, pred.predvar, f(matfit[:, j] - z * matse[:, j]), f(matfit[:, j] + z * matse[:, j]), ci)
            a.plot(pred.predvar, f(matfit[:, j]), **{"color": _theme["line"], **kwargs})
            _finish(a, xlab or "Var", outcome, title or f"Lag {lg:g}", noeff)
        for a, v in zip(axes_var, vars_):
            a.set_visible(True)
            i = pred._var_index(v)
            _ci(a, lags, f(matfit[i] - z * matse[i]), f(matfit[i] + z * matse[i]), ci)
            a.plot(lags, f(matfit[i]), **{"color": _theme["line"], **kwargs})
            _finish(a, xlab or "Lag", outcome, title or f"Var = {v:g}", noeff)
        if ax is None:
            fig.tight_layout()
        return (axes_lag + axes_var)[0] if n == 1 else axes_lag + axes_var

    if ptype == "contour":
        if pred.lag[1] == 0:
            raise ValueError("contour plot not conceivable for unlagged associations")
        ax = ax or plt.gca()
        zmat = f(matfit)
        levels = kwargs.pop("levels", 24)
        cmap, norm_ = _diverging_cmap(noeff, float(zmat.min()), float(zmat.max()))
        cf = ax.contourf(pred.predvar, lags, zmat.T, levels=levels, cmap=kwargs.pop("cmap", cmap), norm=norm_)
        ax.contour(pred.predvar, lags, zmat.T, levels=[noeff], colors=[_theme["ref"]], linewidths=0.6)
        from matplotlib.ticker import MaxNLocator
        cb = plt.colorbar(cf, ax=ax, pad=0.02, fraction=0.05, ticks=MaxNLocator(6))
        cb.set_label(zlab or outcome)
        cb.outline.set_linewidth(0.5)
        ax.grid(False)
        ax.set_xlabel(xlab or "Var")
        ax.set_ylabel(ylab or "Lag")
        if title:
            ax.set_title(title, loc="left")
        return ax

    if ptype == "3d":
        if pred.lag[1] - pred.lag[0] == 0:
            raise ValueError("3D plot not conceivable for unlagged associations")
        if ax is None:
            fig = plt.figure(figsize=(6, 4.8))
            ax = fig.add_subplot(111, projection="3d")
        X, Y = np.meshgrid(pred.predvar, lags, indexing="ij")
        phi, theta = kwargs.pop("phi", 30), kwargs.pop("theta", 210)
        zmat = f(matfit)
        cmap, norm_ = _diverging_cmap(noeff, float(zmat.min()), float(zmat.max()))
        ax.plot_surface(X, Y, zmat, cmap=kwargs.pop("cmap", cmap), norm=norm_,
                        linewidth=0.15, edgecolor="#555555", antialiased=True, **kwargs)
        ax.view_init(elev=phi, azim=theta)
        ax.set_xlabel(xlab or "Var")
        ax.set_ylabel(ylab or "Lag")
        ax.set_zlabel(zlab or outcome)
        ax.xaxis.pane.fill = ax.yaxis.pane.fill = ax.zaxis.pane.fill = False
        ax.grid(True)
        if title:
            ax.set_title(title, loc="left")
        return ax

    raise ValueError("ptype must be one of 'overall', 'slices', 'contour', '3d'")


# ------------------------------------------------------------------------------
def overlay_slices(pred, var=None, lag=None, ci: str = "n", ci_level=None, exp=None, ax=None,
                   xlab=None, ylab=None, title=None, labels=None, legend: bool = True):
    """Several lag-response (``var=[...]``) or exposure-response (``lag=[...]``)
    curves on one axis, with fixed-order colours (or line styles in the
    journal theme), direct labels at the line ends and a legend, as one would
    do in R with ``plot(...)`` followed by ``lines(...)``.
    """
    plt = _mpl()
    if (var is None) == (lag is None):
        raise ValueError("give either 'var' or 'lag' (a list of values)")
    z = _z(pred, ci_level)
    f, noeff, e = _scale(pred, exp)
    ax = ax or plt.gca()
    values = list(np.atleast_1d(var if var is not None else lag))
    if len(values) > len(_theme["series"]):
        raise ValueError(f"at most {len(_theme['series'])} curves can be overlaid; use small multiples instead")
    for k, v in enumerate(values):
        if var is not None:
            i = pred._var_index(v)
            x, fit, se = pred.predlag, pred.matfit[i], pred.matse[i]
        else:
            j = pred._lag_index(v)
            x, fit, se = pred.predvar, pred.matfit[:, j], pred.matse[:, j]
        col, ls = _theme["series"][k], _theme["series_ls"][k]
        lab = labels[k] if labels else f"{v:g}"
        if ci == "area":
            ax.fill_between(x, f(fit - z * se), f(fit + z * se), color=col, alpha=0.12, linewidth=0)
        elif ci == "lines":
            ax.plot(x, f(fit - z * se), color=col, linestyle=":", linewidth=0.7)
            ax.plot(x, f(fit + z * se), color=col, linestyle=":", linewidth=0.7)
        ax.plot(x, f(fit), color=col, linestyle=ls, label=lab)
        if not legend:  # direct labels at the line ends when no legend is drawn
            ax.annotate(lab, xy=(x[-1], f(fit)[-1]), xytext=(4, 0), textcoords="offset points",
                        va="center", fontsize=7, color="#333333")
    xl = xlab or ("Lag" if var is not None else "Var")
    _finish(ax, xl, ylab or ("RR" if e else "Effect"), title, noeff)
    if legend and len(values) > 1:
        ax.legend(title=("Var" if var is not None else "Lag"), loc="best")
    return ax


def summary_figure(pred, var=None, lag=None, xlab="Var", ylab=None, figsize=(7.2, 5.6), exp=None):
    """The standard four-panel figure: overall cumulative curve, contour of
    the exposure-lag-response surface, lag-response curves at ``var`` values
    and exposure-response curves at ``lag`` values (overlaid)."""
    plt = _mpl()
    f, noeff, e = _scale(pred, exp)
    outcome = ylab or ("RR" if e else "Effect")
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    plot_crosspred(pred, "overall", ax=axes[0, 0], xlab=xlab, ylab=outcome, title="a  Overall cumulative", exp=exp)
    plot_crosspred(pred, "contour", ax=axes[0, 1], xlab=xlab, ylab="Lag", zlab=outcome,
                   title="b  Exposure-lag-response surface", exp=exp)
    if var is None:
        q = np.quantile(pred.predvar, [0.01, 0.99])
        var = [float(pred.predvar[np.argmin(np.abs(pred.predvar - q[0]))]),
               float(pred.predvar[np.argmin(np.abs(pred.predvar - q[1]))])]
    if lag is None:
        lag = [float(pred.predlag[0]), float(pred.predlag[min(len(pred.predlag) - 1, len(pred.predlag) // 4)])]
    ax = overlay_slices(pred, var=var, ci="area", ax=axes[1, 0], xlab="Lag", ylab=outcome,
                        title="c  Lag-response at selected values", exp=exp)
    if ax.get_legend() is not None:  # only drawn for more than one value
        ax.get_legend().set_title(xlab)
    overlay_slices(pred, lag=lag, ci="area", ax=axes[1, 1], xlab=xlab, ylab=outcome,
                   title="d  Exposure-response at selected lags", exp=exp)
    fig.tight_layout()
    return fig, axes


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
    ax.plot(x, f(red.fit), **{"color": _theme["line"], **kwargs})
    _finish(ax, xlab or ("Lag" if red.type == "var" else "Var"), ylab or ("RR" if e else "Effect"), title, 1.0 if e else 0.0)
    return ax
