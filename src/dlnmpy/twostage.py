"""Bridges between the two stages of a multi-location DLNM analysis.

Stage 1 reduces each location's DLNM to one-dimensional coefficients
(:func:`dlnmpy.crossreduce`); stage 2 pools them with a multivariate
meta-analysis, which lives in the separate package ``mixmetapy`` (as
``mixmeta`` is separate from ``dlnm`` in R). These helpers move coefficients
out of stage 1 and pooled or BLUP coefficients back into dlnmpy's prediction
machinery. They take and return plain arrays, so neither needs ``mixmetapy``;
any pooling method that gives a coefficient vector and a covariance matrix
works with :func:`predict_reduced`.
"""

from __future__ import annotations

import numpy as np

from .core import CrossBasis, OneBasis
from .predict import CrossPred, crosspred

__all__ = ["stack_reduced", "predict_reduced"]


def stack_reduced(reduced: list) -> tuple[np.ndarray, list]:
    """Stack the coefficients and covariances of a list of ``CrossReduce``
    objects (one per location) into the ``y`` and ``S`` inputs of
    ``mixmetapy.mixmeta``."""
    y = np.vstack([r.coef for r in reduced])
    S = [r.vcov for r in reduced]
    return y, S


def predict_reduced(basis, coef, vcov, at=None, from_=None, to=None, by=None, cen=None,
                    model_link: str = "log", ci_level: float = 0.95) -> CrossPred:
    """Exposure-response curve from reduced (one-dimensional) coefficients,
    e.g. pooled or BLUP coefficients from ``mixmetapy.mixmeta``.

    ``basis`` is a ``CrossBasis`` whose predictor-space specification
    (``argvar``) defines the one-dimensional basis, or a ``OneBasis``.
    """
    if isinstance(basis, CrossBasis):
        av = dict(basis.argvar)
        fun = av.pop("fun")
        av.pop("cen", None)
        ob = OneBasis(matrix=np.zeros((1, basis.df[0])), fun=fun, attrs=av, range=basis.range)
    elif isinstance(basis, OneBasis):
        ob = basis
    else:
        raise TypeError("'basis' must be a CrossBasis or OneBasis")
    return crosspred(ob, coef=coef, vcov=vcov, model_link=model_link, at=at, from_=from_,
                     to=to, by=by, cen=cen, ci_level=ci_level)
