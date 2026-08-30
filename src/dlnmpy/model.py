"""Glue between fitted regression models and the prediction functions.

The prediction functions only need three things from a fitted model: the
coefficients of the basis columns, their covariance matrix and the link
function. This module extracts them from statsmodels results objects (GLM,
GEE, OLS, discrete models, mixed models) and from plain ``dict``/tuple
inputs so that any modelling library can be used.
"""

from __future__ import annotations

import re

import numpy as np

__all__ = ["extract_coef_vcov", "get_link", "fit_glm", "fit_clogit", "design_matrix"]


def _params_index(model):
    params = getattr(model, "params", None)
    if params is None:
        raise TypeError("'model' must expose .params (a statsmodels results object), "
                        "or pass 'coef' and 'vcov' directly")
    names = list(getattr(params, "index", []))
    if not names:
        raise ValueError("model parameters are not named; build the design matrix with "
                         "named columns (see CrossBasis.to_dataframe) or pass coef/vcov")
    return np.asarray(params, dtype=float), names


def _cov(model) -> np.ndarray:
    if hasattr(model, "cov_params"):
        cov = model.cov_params()
    elif hasattr(model, "vcov"):
        cov = model.vcov
    else:
        raise TypeError("cannot extract a covariance matrix from 'model'")
    return np.asarray(cov, dtype=float)


def _match(names: list[str], name: str, kind: str, ncol: int) -> list[int]:
    """Indices (in model order, then sorted by basis column) of the
    coefficients belonging to the basis called ``name``."""
    esc = re.escape(name)
    if ncol == 1:
        exact = [i for i, n in enumerate(names) if n == name]
        if exact:
            return exact
    # the separator between the basis name and the column label is optional
    # (R pastes them directly: "cbv1.l1"; to_dataframe uses "cb_v1_l1"), but
    # it must be a separator: without this anchor the name "cb" also matches
    # the columns of a second basis called "cb2".
    if kind == "cb":
        pat = re.compile(rf"^{esc}[._]?v(\d+)[._]l(\d+)$")
    else:
        pat = re.compile(rf"^{esc}[._]?b(\d+)$")
    found = []
    for i, n in enumerate(names):
        m = pat.match(n)
        if m:
            key = tuple(int(g) for g in m.groups())
            found.append((key, i))
    found.sort()
    return [i for _, i in found]


def extract_coef_vcov(model, name: str, kind: str, ncol: int):
    """Return ``(coef, vcov)`` for the ``ncol`` columns of basis ``name``."""
    params, names = _params_index(model)
    idx = _match(names, name, kind, ncol)
    if len(idx) != ncol:
        raise ValueError(
            f"found {len(idx)} coefficients matching basis '{name}' but the basis has "
            f"{ncol} columns. Name the design columns with to_dataframe(name=...) and "
            f"pass the same name, or supply coef/vcov directly")
    cov = _cov(model)
    return params[idx], cov[np.ix_(idx, idx)]


def get_link(model, model_link: str | None = None) -> str | None:
    """Best-effort extraction of the link function name from a results object."""
    if model_link is not None:
        return model_link
    if model is None:
        return None
    if isinstance(getattr(model, "link", None), str):  # e.g. PenalizedGLMResults
        return model.link
    m = getattr(model, "model", model)
    fam = getattr(m, "family", None)
    link = getattr(fam, "link", None)
    if link is not None:
        cname = type(link).__name__.lower()
        # compound names first: "log" is a substring of "cloglog" and
        # "loglog", so testing it earlier reports those as a log link and
        # makes crosspred exponentiate them into relative risks
        for cand in ("cloglog", "loglog", "logit", "log", "identity", "probit", "inverse"):
            if cand in cname:
                return cand
        return cname
    cname = type(m).__name__.lower()
    if cname in ("logit", "conditionallogit", "conditionalmnlogit"):
        return "logit"
    if cname in ("poisson", "negativebinomial", "phreg", "conditionalpoisson"):
        return "log"
    if cname in ("ols", "wls", "gls", "mixedlm"):
        return "identity"
    return None


def design_matrix(data, *terms, intercept: bool = True):
    """Column-bind a DataFrame with any number of basis objects/arrays.

    ``terms`` may be ``CrossBasis``/``OneBasis`` objects (use their
    ``to_dataframe`` names), tuples ``(name, obj)`` or DataFrames.
    """
    import pandas as pd
    parts = []
    if intercept:
        parts.append(pd.DataFrame({"Intercept": np.ones(len(data))}, index=data.index))
    for t in terms:
        if isinstance(t, tuple):
            name, obj = t
            df = obj.to_dataframe(name)
        elif hasattr(t, "to_dataframe"):
            df = t.to_dataframe()
        else:
            df = pd.DataFrame(t)
        df.index = data.index
        parts.append(df)
    return pd.concat(parts, axis=1)


def aliased_columns(X, tol: float = 1e-7) -> list[int]:
    """Indices of linearly dependent columns of ``X``, detected the way R's
    ``lm.fit``/``glm.fit`` do (LINPACK ``dqrdc2``): a column is aliased when
    its norm after orthogonalisation against the previous kept columns falls
    below ``tol`` times its original norm."""
    X = np.asarray(X, dtype=float)
    n, p = X.shape
    Q = np.zeros((n, 0))
    aliased = []
    for j in range(p):
        v = X[:, j]
        nrm0 = np.linalg.norm(v)
        if nrm0 == 0:
            aliased.append(j)
            continue
        r = v - Q @ (Q.T @ v)
        r = r - Q @ (Q.T @ r)  # re-orthogonalise for stability
        nrm = np.linalg.norm(r)
        if nrm < tol * nrm0:
            aliased.append(j)
        else:
            Q = np.column_stack((Q, r / nrm))
    return aliased


def fit_glm(formula: str, data, family: str = "quasipoisson", drop_aliased: bool = True,
            offset=None, exposure=None, freq_weights=None, var_weights=None,
            **kwargs):
    """Fit a GLM with statsmodels' formula API.

    ``family`` is one of ``"quasipoisson"`` (Poisson with Pearson dispersion,
    as in R), ``"poisson"``, ``"gaussian"``, ``"binomial"``,
    ``"quasibinomial"``. Rows with missing values are dropped, as R's ``glm``
    does by default. With ``drop_aliased=True`` linearly dependent design
    columns are removed before fitting (R reports them as ``NA``
    coefficients), which keeps the residual degrees of freedom, and hence
    the dispersion and standard errors, identical to R's.

    ``offset``, ``exposure``, ``freq_weights`` and ``var_weights`` are
    arguments of the model rather than of the fit, and are passed on as such;
    ``offset`` is R's ``offset()`` term. They are named explicitly because
    ``GLM.fit()`` accepts arbitrary keywords and would silently discard them.
    """
    import statsmodels.api as sm
    import statsmodels.formula.api as smf

    fam = family.lower()
    scale = None
    if fam == "quasipoisson":
        family_obj, scale = sm.families.Poisson(), "X2"
    elif fam == "poisson":
        family_obj = sm.families.Poisson()
    elif fam in ("gaussian", "normal"):
        family_obj = sm.families.Gaussian()
    elif fam in ("binomial", "logistic"):
        family_obj = sm.families.Binomial()
    elif fam == "quasibinomial":
        family_obj, scale = sm.families.Binomial(), "X2"
    else:
        raise ValueError(f"unknown family '{family}'")
    model = smf.glm(formula, data=data, family=family_obj, missing="drop",
                    offset=offset, exposure=exposure,
                    freq_weights=freq_weights, var_weights=var_weights)
    if drop_aliased:
        bad = aliased_columns(model.exog)
        if bad:
            import pandas as pd
            keep = [i for i in range(model.exog.shape[1]) if i not in bad]
            exog = pd.DataFrame(model.exog[:, keep], columns=[model.exog_names[i] for i in keep])
            endog = pd.Series(model.endog, name=model.endog_names)
            dropped = [model.exog_names[i] for i in bad]
            # these have already been subset by missing="drop" above
            model = sm.GLM(endog, exog, family=family_obj,
                           offset=getattr(model, "offset", None),
                           exposure=getattr(model, "exposure", None),
                           freq_weights=getattr(model, "freq_weights", None),
                           var_weights=getattr(model, "var_weights", None))
            model.aliased = dropped  # names of the columns R would report as NA
    # R's glm.fit solves the IRLS weighted least squares by QR; statsmodels
    # defaults to a pseudo-inverse, which differs when the design is close to
    # singular. Use QR unless the caller says otherwise.
    kwargs.setdefault("wls_method", "qr")
    kwargs.setdefault("tol", 1e-10)
    res = model.fit(scale=scale, **kwargs) if scale else model.fit(**kwargs)
    if scale == "X2":
        # R divides the Pearson chi-square by n - p with p the number of
        # non-aliased coefficients; statsmodels uses a numerical rank, which
        # can be smaller for ill-conditioned designs. Re-fit with the R
        # dispersion so that standard errors agree.
        p = res.params.size
        disp = float(res.pearson_chi2 / (res.nobs - p))
        if not np.isclose(disp, res.scale):
            res = model.fit(scale=disp, start_params=res.params, **kwargs)
    res.aliased = getattr(model, "aliased", [])
    return res


def fit_clogit(y, X, groups, **kwargs):
    """Conditional logistic regression (matched case-control designs), the
    Python counterpart of ``survival::clogit``.

    Uses statsmodels' ``ConditionalLogit`` fitted by Newton-Raphson with a
    tight tolerance. The default BFGS optimiser can stop early on flat
    likelihoods and give coefficients that differ from R's at the second
    decimal, so it is not used.
    """
    from statsmodels.discrete.conditional_models import ConditionalLogit

    kwargs.setdefault("method", "newton")
    kwargs.setdefault("tol", 1e-12)
    kwargs.setdefault("maxiter", 500)
    kwargs.setdefault("disp", 0)
    res = ConditionalLogit(np.asarray(y), X, groups=np.asarray(groups)).fit(**kwargs)
    # statsmodels approximates the Hessian numerically (about 1e-4 relative
    # error in the covariance). Replace it with Richardson-extrapolated central
    # differences of the analytic score, which agrees with R's information
    # matrix to about 1e-10.
    b = np.asarray(res.params, dtype=float)
    cov = np.linalg.inv(-_score_jacobian(res.model.score, b))
    inner = getattr(res, "_results", res)  # results may be wrapped
    inner.normalized_cov_params = cov
    inner.cov_params_default = cov
    return res


def _score_jacobian(score, b, h0: float = 1e-4) -> np.ndarray:
    p = b.size
    H = np.zeros((p, p))
    for j in range(p):
        h = h0 * max(1.0, abs(b[j]))
        e = np.zeros(p)
        e[j] = h
        d1 = (score(b + e) - score(b - e)) / (2 * h)
        e[j] = h / 2
        d2 = (score(b + e) - score(b - e)) / h
        H[:, j] = (4 * d2 - d1) / 3
    return (H + H.T) / 2
