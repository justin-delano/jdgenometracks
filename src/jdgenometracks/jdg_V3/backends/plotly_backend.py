"""
Plotly backend implementation for jdgenometracks V3.

This module provides the PlotlyBackend class that implements the backend
interface using plotly for figure creation and rendering.

Author: Assistant
Date: 2024
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import plotly.graph_objects as go
import plotly.subplots as sp

from ..utils.constants import PlotlyConstants
from .base_backend import BaseBackend


@dataclass
class PlotlyBackend(BaseBackend):
    """
    Plotly implementation of the backend interface.

    This backend uses plotly for creating figures, managing subplots,
    and rendering genomic tracks.

    Attributes:
        figure: The plotly Figure object
        subplot_specs: Subplot specifications for plotly
        backend_name: Name of backend ("plotly")
    """

    # Plotly backend constants
    BACKEND_NAME = "plotly"  # Backend identifier
    DEFAULT_WIDTH = 1000  # Default figure width in pixels
    DEFAULT_HEIGHT = 600  # Default figure height in pixels
    DEFAULT_VERTICAL_SPACING = 0.02  # Default vertical spacing between subplots
    DEFAULT_HORIZONTAL_SPACING = 0.02  # Default horizontal spacing between subplots

    subplot_specs: Optional[List[List[dict]]] = None
    num_rows: int = 0
    num_cols: int = 0

    def __post_init__(self):
        """Initialize plotly backend."""
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
    ) -> go.Figure:
        """
        Create a plotly figure with subplots.

        Args:
            num_rows: Number of rows in the subplot grid
            num_cols: Number of columns in the subplot grid
            height_props: Relative heights of each row
            width_props: Relative widths of each column
            row_titles: Titles for each row
            column_titles: Titles for each column
            **kwargs: Additional plotly options

        Returns:
            plotly Figure object
        """
        # Store dimensions for validation
        self.num_rows = num_rows
        self.num_cols = num_cols

        # Extract plotly-specific options
        vertical_spacing = kwargs.get("vertical_spacing", self.DEFAULT_VERTICAL_SPACING)
        horizontal_spacing = kwargs.get(
            "horizontal_spacing", self.DEFAULT_HORIZONTAL_SPACING
        )

        # Override with any options stored in figure_options
        vertical_spacing = self.figure_options.get("vertical_spacing", vertical_spacing)
        horizontal_spacing = self.figure_options.get(
            "horizontal_spacing", horizontal_spacing
        )

        # Create subplot specifications
        self.subplot_specs = [[{} for _ in range(num_cols)] for _ in range(num_rows)]

        # Create figure with subplots
        self.figure = sp.make_subplots(
            rows=num_rows,
            cols=num_cols,
            row_heights=height_props,
            column_widths=width_props,
            subplot_titles=column_titles if len(column_titles) > 0 else None,
            vertical_spacing=vertical_spacing,
            horizontal_spacing=horizontal_spacing,
            shared_xaxes=kwargs.get(
                "shared_xaxes", PlotlyConstants.DEFAULT_SHARED_XAXES
            ),
            shared_yaxes=kwargs.get(
                "shared_yaxes", PlotlyConstants.DEFAULT_SHARED_YAXES
            ),
            specs=self.subplot_specs,
        )

        # Set figure dimensions
        width = kwargs.get("width", self.DEFAULT_WIDTH)
        height = kwargs.get("height", self.DEFAULT_HEIGHT)

        # Override with figure_options
        width = self.figure_options.get("width", width)
        height = self.figure_options.get("height", height)

        self.figure.update_layout(
            width=width,
            height=height,
            autosize=kwargs.get("autosize", PlotlyConstants.DEFAULT_AUTOSIZE),
            showlegend=kwargs.get(
                "showlegend", PlotlyConstants.DEFAULT_SHOWLEGEND_EMPTY
            ),
            # Set background colors
            plot_bgcolor=PlotlyConstants.DEFAULT_PLOT_BGCOLOR,
            paper_bgcolor=PlotlyConstants.DEFAULT_PAPER_BGCOLOR,
            # Set tight margins for genomic plots
            margin=dict(
                l=PlotlyConstants.DEFAULT_MARGIN_LEFT,
                r=PlotlyConstants.DEFAULT_MARGIN_RIGHT,
                t=PlotlyConstants.DEFAULT_MARGIN_TOP,
                b=PlotlyConstants.DEFAULT_MARGIN_BOTTOM,
                pad=PlotlyConstants.DEFAULT_MARGIN_PAD,
            ),
        )

        # Set minimal global defaults only (no axis hiding here)
        self._configure_minimal_defaults()
        

        # Set row titles (as y-axis titles for leftmost subplot in each row)
        self._set_row_titles(row_titles)

        return self.figure

    def get_subplot(self, row: int, col: int) -> Tuple[int, int]:
        """
        Get subplot coordinates for plotly (1-based indexing).

        Args:
            row: Row index (0-based)
            col: Column index (0-based)

        Returns:
            Tuple of (row, col) in 1-based indexing for plotly
        """
        self.validate_subplot_range(row, col)
        return (row + 1, col + 1)  # Plotly uses 1-based indexing

    def set_subplot_title(self, row: int, col: int, title: str) -> None:
        """
        Set title for a specific subplot.

        Args:
            row: Row index (0-based)
            col: Column index (0-based)
            title: Title text
        """
        if self.figure is None:
            raise ValueError("Figure has not been created yet")

        # Plotly subplot titles are set during figure creation
        # This method can be used to update them if needed
        subplot_idx = row * self.num_cols + col
        if subplot_idx < len(self.figure.layout.annotations):
            self.figure.layout.annotations[subplot_idx].text = title

    def set_figure_title(self, title: str) -> None:
        """
        Set the overall figure title.

        Args:
            title: Figure title text
        """
        if self.figure is None:
            raise ValueError("Figure has not been created yet")

        fontsize = self.figure_options.get("title_fontsize", 16)
        self.figure.update_layout(
            title=dict(
                text=title,
                font=dict(size=fontsize),
                x=0.5,  # Center the title
                xanchor="center",
            )
        )

    def finalize_figure(self) -> None:
        """
        Finalize the plotly figure (apply final layout adjustments).
        """
        if self.figure is None:
            raise ValueError("Figure has not been created yet")

        # Apply any final layout adjustments
        # Most layout is handled during creation for plotly
        pass

    def save_figure(self, filename: str, **kwargs) -> None:
        """
        Save the plotly figure to file.

        Args:
            filename: Output filename
            **kwargs: plotly-specific save options (format, scale, etc.)
        """
        if self.figure is None:
            raise ValueError("Figure has not been created yet")

        # Determine format from filename extension if not provided
        file_format = kwargs.get("format")
        if file_format is None:
            if filename.endswith(".html"):
                file_format = "html"
            elif filename.endswith(".png"):
                file_format = "png"
            elif filename.endswith(".jpg") or filename.endswith(".jpeg"):
                file_format = "jpg"
            elif filename.endswith(".pdf"):
                file_format = "pdf"
            elif filename.endswith(".svg"):
                file_format = "svg"
            else:
                file_format = "html"  # Default format

        if file_format == "html":
            self.figure.write_html(filename, **kwargs)
        else:
            # For static formats, use write_image
            save_kwargs = {
                "format": file_format,
                "scale": kwargs.get("scale", 1),
                "width": kwargs.get("width", self.figure.layout.width),
                "height": kwargs.get("height", self.figure.layout.height),
            }
            save_kwargs.update(kwargs)
            self.figure.write_image(filename, **save_kwargs)

    def show_figure(self) -> None:
        """
        Display the plotly figure.
        """
        if self.figure is None:
            raise ValueError("Figure has not been created yet")

        self.figure.show()

    def validate_subplot_range(self, row: int, col: int) -> None:
        """
        Validate that subplot coordinates are within the plotly figure bounds.

        Args:
            row: Row index to validate
            col: Column index to validate

        Raises:
            ValueError: If coordinates are out of bounds
        """
        super().validate_subplot_range(row, col)

        if row >= self.num_rows or col >= self.num_cols:
            raise ValueError(
                f"Subplot coordinates ({row}, {col}) out of bounds. "
                f"Figure has {self.num_rows} rows and {self.num_cols} columns."
            )

    def get_figure_dimensions(self) -> Tuple[float, float]:
        """
        Get the current plotly figure dimensions.

        Returns:
            Tuple of (width, height) in pixels
        """
        if self.figure is None:
            return (self.DEFAULT_WIDTH, self.DEFAULT_HEIGHT)

        layout = self.figure.layout
        width = layout.width or self.DEFAULT_WIDTH
        height = layout.height or self.DEFAULT_HEIGHT
        return (width, height)

    def set_figure_dimensions(self, width: float, height: float) -> None:
        """
        Set the plotly figure dimensions.

        Args:
            width: Figure width in pixels
            height: Figure height in pixels
        """
        super().set_figure_dimensions(width, height)

        if self.figure is not None:
            self.figure.update_layout(width=width, height=height)

    def _set_row_titles(self, row_titles: List[str]) -> None:
        """
        Set row titles for plotly subplots.

        Args:
            row_titles: Titles for each row
        """
        if self.figure is None:
            return

        # Set row titles as y-axis titles for the leftmost subplot in each row
        for row_idx, title in enumerate(row_titles):
            if title and row_idx < self.num_rows:
                plotly_row = row_idx + 1  # 1-based indexing
                self.figure.update_yaxes(
                    title_text=title, row=plotly_row, col=1  # Leftmost column
                )

    def _configure_minimal_defaults(self) -> None:
        """
        Configure minimal default settings for cleaner genomic plots.
        
        Only sets basic aesthetic defaults like gridlines and background.
        Individual tracks are responsible for their own axis configuration.
        """
        if self.figure is None:
            return

        # Apply minimal default settings - just turn off gridlines and zero lines
        self.figure.update_xaxes(
            showgrid=PlotlyConstants.DEFAULT_SHOWGRID_X,
            zeroline=PlotlyConstants.DEFAULT_ZEROLINE_X,
        )
        self.figure.update_yaxes(
            showgrid=PlotlyConstants.DEFAULT_SHOWGRID_Y,
            zeroline=PlotlyConstants.DEFAULT_ZEROLINE_Y,
        )

    def add_trace(self, trace: go.Scatter, row: int, col: int) -> None:
        """
        Add a trace to a specific subplot.

        Args:
            trace: Plotly trace object
            row: Row index (0-based)
            col: Column index (0-based)
        """
        if self.figure is None:
            raise ValueError("Figure has not been created yet")

        plotly_row, plotly_col = self.get_subplot(row, col)
        self.figure.add_trace(trace, row=plotly_row, col=plotly_col)
