"""Deprecated: multivariate meta-analysis moved to ``mixmetapy`` in 0.8.0.

The ``mixmeta`` port is now its own package, as ``mixmeta`` and ``dlnm`` are
separate in R. Install it with the ``twostage`` extra::

    pip install dlnmpy[twostage]

and import from the new places::

    from mixmetapy import mixmeta, MixMeta, vech, xpnd
    from dlnmpy import stack_reduced, predict_reduced

This module re-exports those names with a ``DeprecationWarning`` during 0.8.x
and 0.9.x and will be removed in 1.0.
"""

import warnings

from . import twostage as _twostage

try:
    import mixmetapy as _mixmetapy
except ImportError as e:
    raise ImportError(
        "dlnmpy.meta moved to the separate package mixmetapy in dlnmpy 0.8.0. "
        "Install it with:  pip install dlnmpy[twostage]   "
        "and then use 'from mixmetapy import mixmeta'."
    ) from e

_MOVED = {"mixmeta": _mixmetapy, "MixMeta": _mixmetapy, "vech": _mixmetapy, "xpnd": _mixmetapy,
          "stack_reduced": _twostage, "predict_reduced": _twostage}
_NEW_HOME = {_mixmetapy: "mixmetapy", _twostage: "dlnmpy"}

__all__ = list(_MOVED)


def __getattr__(name):
    if name in _MOVED:
        home = _MOVED[name]
        warnings.warn(f"dlnmpy.meta.{name} is deprecated and will be removed in dlnmpy 1.0; "
                      f"use {_NEW_HOME[home]}.{name}", DeprecationWarning, stacklevel=2)
        return getattr(home, name)
    raise AttributeError(f"module 'dlnmpy.meta' has no attribute {name!r}")
