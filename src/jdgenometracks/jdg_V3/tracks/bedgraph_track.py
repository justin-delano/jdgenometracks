"""
BedGraph track implementation for continuous genomic data visualization.

This module provides the BedGraphTrack class for visualizing continuous genomic
data from bedGraph format files. It supports lines, bars, and points plot types
with both matplotlib and Plotly backends.

Author: Assistant
Date: 2024
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import pandas as pd
import plotly.graph_objects as go
from matplotlib.axes import Axes

from jdgenometracks.jdg_V3.tracks.base_track import BaseTrack
from jdgenometracks.jdg_V3.utils.options import translate


@dataclass
class BedGraphTrack(BaseTrack):
    """
    A class for plotting continuous genomic data from BedGraph files.

    BedGraph tracks visualize continuous genomic data as lines, bars, or points.
    Supports standard bedGraph format with chromosome, start, end, and value columns.
    """

    # BedGraph-specific constants
    BEDGRAPH_COLUMNS = ["chrom", "chromStart", "chromEnd", "value", "name"]

    # Plot type constants
    SUPPORTED_PLOT_TYPES = ["lines", "bars", "points"]
    DEFAULT_PLOT_TYPE = "lines"

    # Data format constants
    MIN_BEDGRAPH_COLUMNS = 4  # Minimum required columns for BedGraph format
    MAX_BEDGRAPH_COLUMNS = 5  # Maximum expected columns for BedGraph format

    # Plotting constants
    DEFAULT_BAR_WIDTH_FALLBACK = 1000  # Default bar width when calculation fails
    DEFAULT_FILL_ALPHA = 0.3  # Default transparency for fill areas
    MIN_POINTS_FOR_BAR_WIDTH = 1  # Minimum points needed for bar width calculation
    MIDPOINT_DIVISOR = 2  # Divisor for midpoint calculations

    def __post_init__(self):
        """Initialize the BedGraph track."""
        super().__post_init__()

        # Validate plot type
        plot_type = self.track_options.get("plot.type", self.DEFAULT_PLOT_TYPE)
        if plot_type not in self.SUPPORTED_PLOT_TYPES:
            raise ValueError(
                f"Invalid plot type '{plot_type}'. Supported types: {self.SUPPORTED_PLOT_TYPES}"
            )

    def load_data(self, file_path: Optional[str] = None) -> pd.DataFrame:
        """
        Load BedGraph format data from file.

        Args:
            file_path: Path to BedGraph file (uses self.file_path if None)

        Returns:
            DataFrame with standardized BedGraph columns

        Raises:
            FileNotFoundError: If file doesn't exist
            ValueError: If file format is invalid
        """
        # Use provided file_path or fallback to instance file_path
        file_to_load = file_path or self.file_path
        if file_to_load is None:
            raise ValueError("No file path provided for data loading")

        # Use common file reading method
        data = self._read_tabbed_file(
            file_path=file_to_load,
            expected_columns=self.BEDGRAPH_COLUMNS,
            min_cols=self.MIN_BEDGRAPH_COLUMNS,
            max_cols=None,  # Allow any number of columns, ignore extras
        )

        # Apply data type conversions
        data = self._apply_bedgraph_dtypes(data)

        return data

    def validate_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Validate and clean BedGraph data.

        Args:
            data: Raw BedGraph data

        Returns:
            Validated and cleaned data

        Raises:
            ValueError: If data is invalid
        """
        # Use common validation methods
        self._validate_required_columns(
            data, ["chrom", "chromStart", "chromEnd", "value"]
        )
        self._validate_coordinate_ranges(data)
        self._validate_numeric_columns(data, ["value"])
        return data

    def _apply_bedgraph_dtypes(self, data: pd.DataFrame) -> pd.DataFrame:
        """Apply appropriate data types to BedGraph columns."""
        # Ensure we have the required columns
        if len(data.columns) < self.MIN_BEDGRAPH_COLUMNS:
            raise ValueError(
                f"BedGraph file must have at least {self.MIN_BEDGRAPH_COLUMNS} columns"
            )

        # Apply data types for required columns
        dtype_map = {"chrom": str, "chromStart": int, "chromEnd": int, "value": float}
        return self._apply_column_dtypes(data, dtype_map)

    def calculate_mid_points(self, data: pd.DataFrame) -> pd.Series:
        """
        Calculate the midpoints of genomic regions for the x-axis in plots.

        Args:
            data: The cleaned data

        Returns:
            The midpoints of the genomic regions
        """
        return (data["chromStart"] + data["chromEnd"]) / self.MIDPOINT_DIVISOR

    def plot_matplotlib(self, ax: Axes, start: int, end: int, **kwargs) -> None:
        """
        Plot the BedGraph data using matplotlib.

        Args:
            ax: The matplotlib axis to plot on
            start: Start position of the region to plot
            end: End position of the region to plot
            **kwargs: Additional plotting parameters
        """
        # Load data if not already loaded and filter for the region
        if self.data is None:
            self.data = self.load_data()

        cleaned_data = self._filter_by_coordinates(self.data, start, end)

        if cleaned_data.empty:
            return

        mid_points = self.calculate_mid_points(cleaned_data)
        y_values = cleaned_data["value"].astype(float)

        # Translate options for matplotlib
        mpl_opts = translate(self.track_options, target="mpl")

        # Extract specific option categories
        line_opts = mpl_opts.get("line", {})
        marker_opts = mpl_opts.get("marker", {})
        fill_opts = mpl_opts.get("fill", {})
        legend_opts = mpl_opts.get("legend", {})

        # Get plot type configuration
        plot_type = self.track_options.get("plot.type", "lines")

        # Create the plot based on type
        if plot_type == "lines":
            ax.plot(mid_points, y_values, label=self.track_name, **line_opts)
            # Add fill if requested
            if self.track_options.get("fill.enabled", False):
                # Use alpha from options if provided, otherwise use default
                alpha_value = fill_opts.get("alpha", self.DEFAULT_FILL_ALPHA)
                fill_opts_clean = {
                    k: v for k, v in fill_opts.items() if k not in ["alpha"]
                }
                # Use type: ignore to suppress false positive type error
                ax.fill_between(mid_points, y_values, alpha=alpha_value, **fill_opts_clean)  # type: ignore
        elif plot_type == "bars":
            bar_width = (
                mid_points.diff().median()
                if len(mid_points) > self.MIN_POINTS_FOR_BAR_WIDTH
                else self.DEFAULT_BAR_WIDTH_FALLBACK
            )
            ax.bar(
                mid_points,
                y_values,
                width=bar_width,
                label=self.track_name,
                **marker_opts,
            )
        elif plot_type == "points":
            ax.scatter(mid_points, y_values, label=self.track_name, **marker_opts)

        # Add legend if requested
        if self.show_legend:
            ax.legend(**legend_opts)

        # Add horizontal lines
        self._add_horizontal_lines_matplotlib(ax)

        # Configure axis styling (following original MPL plotter pattern)
        ax.spines["bottom"].set_visible(False)  # Hide bottom spine/border
        ax.xaxis.set_tick_params(bottom=False)  # Hide x-axis ticks

        # Set custom y-axis limits if specified
        ymin = self.track_options.get("ymin", None)
        ymax = self.track_options.get("ymax", None)
        if ymin is not None:
            ax.set_ylim(bottom=ymin)
        if ymax is not None:
            ax.set_ylim(top=ymax)

    def plot_plotly(
        self, fig: go.Figure, row: int, col: int, start: int, end: int, **kwargs
    ) -> None:
        """
        Plot the BedGraph data using Plotly.

        Args:
            fig: The Plotly figure to add the plot to
            row: The row number in the subplot grid
            col: The column number in the subplot grid
            start: Start position of the region to plot
            end: End position of the region to plot
            **kwargs: Additional plotting parameters
        """
        # Load data if not already loaded and filter for the region
        if self.data is None:
            self.data = self.load_data()

        cleaned_data = self._filter_by_coordinates(self.data, start, end)

        if cleaned_data.empty:
            return

        mid_points = self.calculate_mid_points(cleaned_data)
        y_values = cleaned_data["value"].astype(float)

        # Translate options for Plotly
        plotly_opts = translate(self.track_options, target="plotly")

        # Get plot type configuration
        plot_type = self.track_options.get("plot.type", "lines")

        # Extract line/marker options
        line_opts = plotly_opts.get("line", {})
        marker_opts = plotly_opts.get("marker", {})

        # Create the appropriate trace
        if plot_type == "lines":
            trace = go.Scatter(
                x=mid_points,
                y=y_values,
                mode="lines",
                name=self.track_name or "BedGraph Track",
                line=line_opts,
                showlegend=self.show_legend,
            )
            # Add fill if requested
            if self.track_options.get("fill.enabled", False):
                trace.fill = "tozeroy"

        elif plot_type == "bars":
            trace = go.Bar(
                x=mid_points,
                y=y_values,
                name=self.track_name or "BedGraph Track",
                marker=marker_opts,
                showlegend=self.show_legend,
            )
        elif plot_type == "points":
            trace = go.Scatter(
                x=mid_points,
                y=y_values,
                mode="markers",
                name=self.track_name or "BedGraph Track",
                marker=marker_opts,
                showlegend=self.show_legend,
            )

        fig.add_trace(trace, row=row, col=col)

        # Add horizontal lines
        self._add_horizontal_lines_plotly(fig, row, col)

        # Configure axes for this track - follow V2 pattern
        fig.update_xaxes(
            showline=False,        # Hide x-axis line (like V2)
            row=row,
            col=col,
        )
        fig.update_yaxes(
            linecolor="black",     # Show y-axis with black line (like V2)
            row=row,
            col=col,
        )
