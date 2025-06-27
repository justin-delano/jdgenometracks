"""
Matplotlib backend implementation for jdgenometracks V3.

This module provides the MatplotlibBackend class that implements the backend
interface using matplotlib for figure creation and rendering.

Author: Assistant
Date: 2024
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np

from ..utils.constants import MatplotlibConstants
from .base_backend import BaseBackend


@dataclass
class MatplotlibBackend(BaseBackend):
    """
    Matplotlib implementation of the backend interface.

    This backend uses matplotlib for creating figures, managing subplots,
    and rendering genomic tracks.

    Attributes:
        figure: The matplotlib Figure object
        axes: 2D array of matplotlib Axes objects
        backend_name: Name of backend ("matplotlib")
    """

    # Matplotlib backend constants
    BACKEND_NAME = "matplotlib"  # Backend identifier
    DEFAULT_DPI = 100  # Default DPI for figure creation
    DEFAULT_FIGSIZE = (10.0, 8.0)  # Default figure size (width, height)
    DEFAULT_CONSTRAINED_LAYOUT = True  # Use constrained layout by default

    axes: Optional[np.ndarray] = None

    def __post_init__(self):
        """Initialize matplotlib backend."""
        self.backend_name = self.BACKEND_NAME
        super().__post_init__()

    def create_figure(
        self,
        num_rows: int,
        num_cols: int,
        height_props: List[float],
        width_props: List[float],
        row_titles: List[str],
        column_titles: List[str],
        **kwargs,
    ) -> matplotlib.figure.Figure:
        """
        Create a matplotlib figure with subplots.

        Args:
            num_rows: Number of rows in the subplot grid
            num_cols: Number of columns in the subplot grid
            height_props: Relative heights of each row
            width_props: Relative widths of each column
            row_titles: Titles for each row
            column_titles: Titles for each column
            **kwargs: Additional matplotlib options

        Returns:
            matplotlib Figure object
        """
        # Extract matplotlib-specific options
        figsize = kwargs.get("figsize", self.DEFAULT_FIGSIZE)
        dpi = kwargs.get("dpi", self.DEFAULT_DPI)
        constrained_layout = kwargs.get(
            "constrained_layout", self.DEFAULT_CONSTRAINED_LAYOUT
        )

        # Convert width/height to figsize if provided
        if "width" in kwargs and "height" in kwargs:
            # Convert pixels to inches using a reasonable DPI (72 points per inch is common)
            dpi_scale = 72.0
            width_inches = kwargs["width"] / dpi_scale
            height_inches = kwargs["height"] / dpi_scale
            figsize = (width_inches, height_inches)

        # Override with any options stored in figure_options
        figsize = self.figure_options.get("figsize", figsize)
        dpi = self.figure_options.get("dpi", dpi)
        constrained_layout = self.figure_options.get(
            "constrained_layout", constrained_layout
        )

        # Create figure and subplots
        self.figure, self.axes = plt.subplots(
            num_rows,
            num_cols,
            figsize=figsize,
            gridspec_kw={
                "height_ratios": height_props,
                "width_ratios": width_props,
                "hspace": kwargs.get("hspace", 0.3),
                "wspace": kwargs.get("wspace", 0.2),
            },
            dpi=dpi,
            constrained_layout=constrained_layout,
            sharex=kwargs.get("sharex", MatplotlibConstants.DEFAULT_SHAREX),
            sharey=kwargs.get("sharey", MatplotlibConstants.DEFAULT_SHAREY),
        )

        # Ensure axes is always 2D for consistent indexing
        if num_rows == 1 and num_cols == 1:
            self.axes = np.array([[self.axes]])
        elif num_rows == 1:
            self.axes = self.axes.reshape(1, -1)
        elif num_cols == 1:
            self.axes = self.axes.reshape(-1, 1)

        # Set subplot titles
        self._set_subplot_titles(column_titles, row_titles)

        return self.figure

    def get_subplot(self, row: int, col: int) -> plt.Axes:
        """
        Get a specific matplotlib Axes object.

        Args:
            row: Row index (0-based)
            col: Column index (0-based)

        Returns:
            matplotlib Axes object
        """
        self.validate_subplot_range(row, col)
        if self.axes is None:
            raise ValueError("Subplots have not been created yet")
        return self.axes[row, col]

    def set_subplot_title(self, row: int, col: int, title: str) -> None:
        """
        Set title for a specific subplot.

        Args:
            row: Row index (0-based)
            col: Column index (0-based)
            title: Title text
        """
        ax = self.get_subplot(row, col)
        ax.set_title(title)

    def set_figure_title(self, title: str) -> None:
        """
        Set the overall figure title.

        Args:
            title: Figure title text
        """
        if self.figure is None:
            raise ValueError("Figure has not been created yet")

        fontsize = self.figure_options.get("title_fontsize", 16)
        self.figure.suptitle(title, fontsize=fontsize)

    def finalize_figure(self) -> None:
        """
        Finalize the matplotlib figure (adjust layout, etc.).
        """
        if self.figure is None:
            raise ValueError("Figure has not been created yet")

        # Apply any final adjustments
        if not self.figure_options.get(
            "constrained_layout", self.DEFAULT_CONSTRAINED_LAYOUT
        ):
            self.figure.tight_layout()

    def save_figure(self, filename: str, **kwargs) -> None:
        """
        Save the matplotlib figure to file.

        Args:
            filename: Output filename
            **kwargs: matplotlib-specific save options (dpi, bbox_inches, etc.)
        """
        if self.figure is None:
            raise ValueError("Figure has not been created yet")

        # Set default save options
        save_kwargs = {
            "dpi": kwargs.get("dpi", self.DEFAULT_DPI),
            "bbox_inches": kwargs.get("bbox_inches", "tight"),
            "facecolor": kwargs.get("facecolor", "white"),
            "edgecolor": kwargs.get("edgecolor", "none"),
        }
        save_kwargs.update(kwargs)

        self.figure.savefig(filename, **save_kwargs)

    def show_figure(self) -> None:
        """
        Display the matplotlib figure.
        """
        if self.figure is None:
            raise ValueError("Figure has not been created yet")

        plt.show()

    def validate_subplot_range(self, row: int, col: int) -> None:
        """
        Validate that subplot coordinates are within the matplotlib figure bounds.

        Args:
            row: Row index to validate
            col: Column index to validate

        Raises:
            ValueError: If coordinates are out of bounds
        """
        super().validate_subplot_range(row, col)

        if self.axes is None:
            raise ValueError("Subplots have not been created yet")

        max_rows, max_cols = self.axes.shape
        if row >= max_rows or col >= max_cols:
            raise ValueError(
                f"Subplot coordinates ({row}, {col}) out of bounds. "
                f"Figure has {max_rows} rows and {max_cols} columns."
            )

    def get_figure_dimensions(self) -> Tuple[float, float]:
        """
        Get the current matplotlib figure dimensions.

        Returns:
            Tuple of (width, height) in inches
        """
        if self.figure is None:
            return self.DEFAULT_FIGSIZE

        return self.figure.get_size_inches()

    def set_figure_dimensions(self, width: float, height: float) -> None:
        """
        Set the matplotlib figure dimensions.

        Args:
            width: Figure width in inches
            height: Figure height in inches
        """
        super().set_figure_dimensions(width, height)

        if self.figure is not None:
            self.figure.set_size_inches(width, height)

    def _set_subplot_titles(
        self, column_titles: List[str], row_titles: List[str]
    ) -> None:
        """
        Set titles for matplotlib subplots.

        Args:
            column_titles: Titles for each column
            row_titles: Titles for each row
        """
        if self.axes is None:
            return

        # Set column titles (top row)
        for col_idx, title in enumerate(column_titles):
            if title and col_idx < self.axes.shape[1]:
                self.axes[0, col_idx].set_title(title)

        # Set row titles (leftmost column as y-labels)
        for row_idx, title in enumerate(row_titles):
            if title and row_idx < self.axes.shape[0]:
                self.axes[row_idx, 0].set_ylabel(
                    title, rotation=0, ha="right", va="center"
                )

    def close_figure(self) -> None:
        """Close the matplotlib figure and free memory."""
        if self.figure is not None:
            plt.close(self.figure)
            self.figure = None
            self.axes = None
