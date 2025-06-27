"""
Abstract base backend class for jdgenometracks V3.

This module provides the foundational BaseBackend class that all specific backend
implementations (matplotlib, plotly, etc.) inherit from. It defines the interface
for figure creation, track plotting, and output generation.

Author: Assistant
Date: 2024
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple


@dataclass
class BaseBackend(ABC):
    """
    Abstract base class for all plotting backends.

    This class provides the interface and common functionality that all backend types
    must implement. Each backend is responsible for creating figures, managing subplots,
    and rendering tracks using its specific plotting library.

    Attributes:
        figure: The backend-specific figure object
        backend_name: Name of the backend implementation
        figure_options: Dictionary containing backend-specific figure options
    """

    # Base backend constants
    DEFAULT_BACKEND_NAME = "base"  # Default backend name
    DEFAULT_FIGURE_HEIGHT = 8.0  # Default figure height in inches
    DEFAULT_FIGURE_WIDTH = 10.0  # Default figure width in inches
    MIN_TRACK_HEIGHT = 0.1  # Minimum track height
    MAX_TRACK_HEIGHT = 10.0  # Maximum track height

    figure: Optional[Any] = None
    backend_name: str = field(default=DEFAULT_BACKEND_NAME)
    figure_options: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Initialize backend after dataclass creation."""
        # Validate figure options
        self._validate_figure_options()

    def _validate_figure_options(self) -> None:
        """Validate backend-specific figure options."""
        # Basic validation - subclasses can override for more specific validation
        if not isinstance(self.figure_options, dict):
            raise TypeError("figure_options must be a dictionary")

    @abstractmethod
    def create_figure(
        self,
        num_rows: int,
        num_cols: int,
        height_props: List[float],
        width_props: List[float],
        row_titles: List[str],
        column_titles: List[str],
        **kwargs,
    ) -> Any:
        """
        Create a figure with subplots.

        Args:
            num_rows: Number of rows in the subplot grid
            num_cols: Number of columns in the subplot grid
            height_props: Relative heights of each row
            width_props: Relative widths of each column
            row_titles: Titles for each row
            column_titles: Titles for each column
            **kwargs: Additional backend-specific options

        Returns:
            Backend-specific figure object
        """
        pass

    @abstractmethod
    def get_subplot(self, row: int, col: int) -> Any:
        """
        Get a specific subplot for plotting.

        Args:
            row: Row index (0-based)
            col: Column index (0-based)

        Returns:
            Backend-specific subplot/axes object
        """
        pass

    @abstractmethod
    def set_subplot_title(self, row: int, col: int, title: str) -> None:
        """
        Set title for a specific subplot.

        Args:
            row: Row index (0-based)
            col: Column index (0-based)
            title: Title text
        """
        pass

    @abstractmethod
    def set_figure_title(self, title: str) -> None:
        """
        Set the overall figure title.

        Args:
            title: Figure title text
        """
        pass

    @abstractmethod
    def finalize_figure(self) -> None:
        """
        Finalize the figure (adjust layout, etc.).
        """
        pass

    @abstractmethod
    def save_figure(self, filename: str, **kwargs) -> None:
        """
        Save the figure to file.

        Args:
            filename: Output filename
            **kwargs: Backend-specific save options
        """
        pass

    @abstractmethod
    def show_figure(self) -> None:
        """
        Display the figure.
        """
        pass

    def validate_subplot_range(self, row: int, col: int) -> None:
        """
        Validate that subplot coordinates are within the figure bounds.

        Args:
            row: Row index to validate
            col: Column index to validate

        Raises:
            ValueError: If coordinates are out of bounds
        """
        if self.figure is None:
            raise ValueError("Figure has not been created yet")

        # Subclasses should override this with specific validation logic
        if row < 0 or col < 0:
            raise ValueError(
                f"Subplot coordinates must be non-negative, got ({row}, {col})"
            )

    def get_figure_dimensions(self) -> Tuple[float, float]:
        """
        Get the current figure dimensions.

        Returns:
            Tuple of (width, height) in inches
        """
        # Default implementation - subclasses should override
        return (self.DEFAULT_FIGURE_WIDTH, self.DEFAULT_FIGURE_HEIGHT)

    def set_figure_dimensions(self, width: float, height: float) -> None:
        """
        Set the figure dimensions.

        Args:
            width: Figure width in inches
            height: Figure height in inches
        """
        # Validate dimensions
        if width <= 0 or height <= 0:
            raise ValueError("Figure dimensions must be positive")

        # Store in figure options for subclasses to use
        self.figure_options.update({"width": width, "height": height})

    def get_backend_info(self) -> Dict[str, Any]:
        """
        Get information about this backend.

        Returns:
            Dictionary containing backend information
        """
        return {
            "name": self.backend_name,
            "figure_created": self.figure is not None,
            "figure_options": self.figure_options.copy(),
        }

    def __repr__(self) -> str:
        """String representation of the backend."""
        return f"{self.__class__.__name__}(backend_name='{self.backend_name}')"
