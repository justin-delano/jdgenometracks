"""
Backend factory for jdgenometracks V3.

This module provides a factory pattern for creating and managing different
plotting backends (matplotlib, plotly, etc.).

Author: Assistant
Date: 2024
"""

from __future__ import annotations

from typing import Dict, Optional, Type

from .base_backend import BaseBackend
from .matplotlib_backend import MatplotlibBackend
from .plotly_backend import PlotlyBackend


class BackendFactory:
    """
    Factory class for creating plotting backends.

    This factory manages the registration and creation of different backend
    implementations, providing a unified interface for backend selection.
    """

    # Backend factory constants
    DEFAULT_BACKEND = "plotly"  # Default backend to use
    MATPLOTLIB_ALIASES = ["matplotlib", "mpl", "pyplot"]  # Aliases for matplotlib
    PLOTLY_ALIASES = ["plotly", "plot.ly", "dash"]  # Aliases for plotly

    def __init__(self):
        """Initialize the backend factory with default backends."""
        self._backends: Dict[str, Type[BaseBackend]] = {}
        self._register_default_backends()

    def _register_default_backends(self) -> None:
        """Register the default backend implementations."""
        # Register matplotlib backend with all its aliases
        for alias in self.MATPLOTLIB_ALIASES:
            self._backends[alias] = MatplotlibBackend

        # Register plotly backend with all its aliases
        for alias in self.PLOTLY_ALIASES:
            self._backends[alias] = PlotlyBackend

    def register_backend(self, name: str, backend_class: Type[BaseBackend]) -> None:
        """
        Register a new backend implementation.

        Args:
            name: Name/alias for the backend
            backend_class: Backend class that inherits from BaseBackend

        Raises:
            TypeError: If backend_class doesn't inherit from BaseBackend
        """
        if not issubclass(backend_class, BaseBackend):
            raise TypeError(
                f"Backend class must inherit from BaseBackend, got {backend_class}"
            )

        self._backends[name.lower()] = backend_class

    def create_backend(
        self, backend_name: Optional[str] = None, **kwargs
    ) -> BaseBackend:
        """
        Create a backend instance.

        Args:
            backend_name: Name of the backend to create (uses default if None)
            **kwargs: Additional arguments to pass to backend constructor

        Returns:
            Backend instance

        Raises:
            ValueError: If backend name is not recognized
        """
        if backend_name is None:
            backend_name = self.DEFAULT_BACKEND

        backend_key = backend_name.lower()
        if backend_key not in self._backends:
            available = list(self._backends.keys())
            raise ValueError(
                f"Unknown backend '{backend_name}'. Available backends: {available}"
            )

        backend_class = self._backends[backend_key]
        return backend_class(**kwargs)

    def get_available_backends(self) -> Dict[str, str]:
        """
        Get information about available backends.

        Returns:
            Dictionary mapping backend names to their class names
        """
        return {
            name: backend_class.__name__
            for name, backend_class in self._backends.items()
        }

    def is_backend_available(self, backend_name: str) -> bool:
        """
        Check if a backend is available.

        Args:
            backend_name: Name of the backend to check

        Returns:
            True if backend is available, False otherwise
        """
        return backend_name.lower() in self._backends

    def get_default_backend(self) -> str:
        """
        Get the name of the default backend.

        Returns:
            Name of the default backend
        """
        return self.DEFAULT_BACKEND

    def set_default_backend(self, backend_name: str) -> None:
        """
        Set the default backend.

        Args:
            backend_name: Name of the backend to set as default

        Raises:
            ValueError: If backend name is not recognized
        """
        if not self.is_backend_available(backend_name):
            available = list(self._backends.keys())
            raise ValueError(
                f"Cannot set unknown backend '{backend_name}' as default. "
                f"Available backends: {available}"
            )

        self.DEFAULT_BACKEND = backend_name.lower()


# Global backend factory instance
_backend_factory = BackendFactory()


def create_backend(backend_name: Optional[str] = None, **kwargs) -> BaseBackend:
    """
    Convenience function to create a backend using the global factory.

    Args:
        backend_name: Name of the backend to create
        **kwargs: Additional arguments to pass to backend constructor

    Returns:
        Backend instance
    """
    return _backend_factory.create_backend(backend_name, **kwargs)


def get_available_backends() -> Dict[str, str]:
    """
    Convenience function to get available backends from the global factory.

    Returns:
        Dictionary mapping backend names to their class names
    """
    return _backend_factory.get_available_backends()


def register_backend(name: str, backend_class: Type[BaseBackend]) -> None:
    """
    Convenience function to register a backend with the global factory.

    Args:
        name: Name/alias for the backend
        backend_class: Backend class that inherits from BaseBackend
    """
    _backend_factory.register_backend(name, backend_class)
