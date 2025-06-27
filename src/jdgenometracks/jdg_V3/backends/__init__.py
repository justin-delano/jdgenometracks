"""
Backend implementations for jdgenometracks V3.

This module provides backend implementations for different plotting libraries,
all inheriting from a common BaseBackend interface.

Available backends:
- MatplotlibBackend: For matplotlib-based plotting
- PlotlyBackend: For plotly-based plotting
"""

from .backend_factory import (
    BackendFactory,
    create_backend,
    get_available_backends,
    register_backend,
)
from .base_backend import BaseBackend
from .matplotlib_backend import MatplotlibBackend
from .plotly_backend import PlotlyBackend

__all__ = [
    "BaseBackend",
    "MatplotlibBackend",
    "PlotlyBackend",
    "BackendFactory",
    "create_backend",
    "get_available_backends",
    "register_backend",
]
