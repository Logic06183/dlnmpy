"""dlnmpy: distributed lag non-linear models in Python.

A faithful port of the R package ``dlnm`` (Gasparrini) built from the
published methodology and the R source, and validated numerically against
R on the ``chicagoNMMAPS`` examples.

Typical workflow::

    import dlnmpy as dl
    df = dl.datasets.chicago_nmmaps()
    cb = dl.crossbasis(df["temp"], lag=21, argvar={"fun": "ns", "df": 4},
                       arglag={"fun": "ns", "df": 4})
    X = df.join(cb.to_dataframe("cb"))
    fit = dl.fit_glm("death ~ " + " + ".join(cb.to_dataframe("cb").columns)
                     + " + dow + ns_time", X, family="quasipoisson")
    pred = dl.crosspred(cb, fit, cen=21, by=1, name="cb")
    pred.plot("overall")
"""

from . import attribution, basis, datasets, meta, plot, uncertainty
from ._rcompat import pretty, quantile7
from .basis import bs, cr, integer, lin, ns, poly, ps, strata, thr
from .core import CrossBasis, OneBasis, crossbasis, onebasis
from .knots import equalknots, logknots
from .lag import exphist, lag_matrix, mklag, seqlag
from .model import design_matrix, extract_coef_vcov, fit_clogit, fit_glm, get_link
from .penalty import cbpen
from .predict import CrossPred, CrossReduce, crosspred, crossreduce
from .attribution import attr_table, attrdl, findmin, mmt
from .meta import MixMeta, mixmeta, predict_reduced, stack_reduced
from .penalized import PenalizedGLMResults, fit_pgam, fit_pglm
from .uncertainty import bootstrap, empirical_ci, model_grid, qaic, simulate_pred

__version__ = "0.4.1"

__all__ = [
    "onebasis", "crossbasis", "crosspred", "crossreduce", "exphist", "logknots",
    "equalknots", "cbpen", "OneBasis", "CrossBasis", "CrossPred", "CrossReduce",
    "lin", "poly", "strata", "thr", "integer", "ns", "bs", "ps", "cr",
    "fit_pgam", "fit_pglm", "PenalizedGLMResults", "bootstrap", "empirical_ci", "model_grid", "qaic", "simulate_pred",
    "lag_matrix", "mklag", "seqlag", "fit_glm", "fit_clogit", "design_matrix", "extract_coef_vcov",
    "get_link", "pretty", "quantile7", "basis", "datasets", "attribution", "plot", "uncertainty",
    "attrdl", "findmin", "mmt", "attr_table", "meta", "mixmeta", "MixMeta", "predict_reduced", "stack_reduced",
]
