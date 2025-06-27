"""
X-axis track implementation for genomic coordinate display.

This module provides the XAxisTrack class for displaying genomic coordinates
with intelligent tick formatting. It automatically adapts between bp, Kb, and Mb
units for optimal readability.

Author: Assistant
Date: 2024
"""

from dataclasses import dataclass
from typing import Optional

import pandas as pd
import plotly.graph_objects as go
from matplotlib.axes import Axes
from matplotlib.ticker import FixedLocator, Formatter

from ..utils.coordinates import GenomicTickCalculator, format_genomic_ticks
from ..utils.options import translate
from .base_track import BaseTrack


class BPFormatter(Formatter):
    """
    Unified intelligent formatter for genomic coordinates with adaptive tick placement.

    This formatter uses the GenomicTickCalculator to determine optimal tick
    locations and formatting based on the visible genomic range. It automatically
    adapts between bp, Kb, and Mb units for optimal readability.
    """

    # BPFormatter constants
    DEFAULT_SCALE_FACTOR = 1  # Default scale factor for formatting
    DEFAULT_UNIT_LABEL = "b"  # Default unit label
    TICK_TOLERANCE = 0.5  # Tolerance for matching tick locations
    DECIMAL_PLACES = 1  # Number of decimal places for non-integer values

    def __init__(self):
        self.tick_locations = []
        self.tick_labels = []
        self.unit_label = self.DEFAULT_UNIT_LABEL
        self.scale_factor = self.DEFAULT_SCALE_FACTOR

    def __call__(self, x, pos=None):
        """Return the formatted tick value."""
        # Find the closest tick location to format it properly
        if not self.tick_labels:
            return f"{int(x):,}"

        # Find the index of this tick
        tick_idx = None
        for i, loc in enumerate(self.tick_locations):
            if abs(loc - x) < self.TICK_TOLERANCE:  # Close enough
                tick_idx = i
                break

        if tick_idx is not None and tick_idx < len(self.tick_labels):
            return self.tick_labels[tick_idx]
        else:
            # Fallback formatting
            scaled = x / self.scale_factor
            if self.scale_factor == self.DEFAULT_SCALE_FACTOR:
                return f"{int(x):,}"
            elif scaled == int(scaled):
                return f"{int(scaled):,}"
            else:
                return f"{scaled:.{self.DECIMAL_PLACES}f}"

    def set_range(self, vmin, vmax):
        """Set the visible range and calculate optimal ticks."""
        # Calculate optimal ticks for this range
        self.tick_locations, self.unit_label, self.scale_factor = (
            GenomicTickCalculator.calculate_optimal_ticks(int(vmin), int(vmax))
        )

        # Format the tick labels
        self.tick_labels = format_genomic_ticks(
            self.tick_locations, self.unit_label, self.scale_factor
        )

        # Set the locations for matplotlib
        self.locs = self.tick_locations


@dataclass
class XAxisTrack(BaseTrack):
    """
    Track class for displaying genomic X-axis with intelligent tick formatting.

    Uses a unified BPFormatter that automatically adapts tick placement and
    formatting based on the genomic range being displayed, supporting both
    matplotlib and Plotly backends.

    Attributes:
        track_options: Dictionary containing styling and display options
    """

    # XAxisTrack constants
    TRACK_TYPE = "axis"  # Type identifier for axis tracks
    DEFAULT_TRACK_NAME = "x_axis"  # Default name for axis tracks
    DEFAULT_FONT_SIZE = 12  # Default font size for axis labels
    TICK_FONT_SCALE = 0.8  # Scale factor for tick label font size
    DEFAULT_GRID_WIDTH = 1  # Default grid line width
    DEFAULT_GRID_COLOR = "lightgray"  # Default grid line color
    DEFAULT_AXIS_TYPE = "verbose"  # Default axis type
    DEFAULT_UNIT_LABEL = "bp"  # Default unit label for axis

    def __post_init__(self):
        """Initialize for axis track (no data needed)."""
        # Set track type
        self.track_type = self.TRACK_TYPE

        # Axis tracks don't need data or file paths
        self.data = pd.DataFrame()
        self._is_loaded = True

        # Set default track name if not provided
        if self.track_name is None:
            self.track_name = self.DEFAULT_TRACK_NAME

        # Don't call super().__post_init__() to avoid validation requirements

    def load_data(self, file_path: Optional[str] = None) -> pd.DataFrame:
        """
        Load data for axis track (always empty).

        Args:
            file_path: Not used for axis tracks

        Returns:
            Empty DataFrame as axis tracks don't display data
        """
        return pd.DataFrame()

    def validate_data(self, data: pd.DataFrame) -> bool:
        """
        Validate data for axis track.

        Args:
            data: DataFrame to validate (should be empty)

        Returns:
            True if data is empty (as expected for axis tracks)
        """
        return len(data) == 0

    def subset_data_for_region(
        self, data: pd.DataFrame, start: int, end: int
    ) -> pd.DataFrame:
        """
        Subset data for axis track (always empty).

        Args:
            data: Input DataFrame (should be empty)
            start: Start position (not used)
            end: End position (not used)

        Returns:
            Same empty DataFrame
        """
        return data

    def add_verbose_axis_mpl(self, ax: Axes, chromosome: str) -> None:
        """
        Add a verbose X-axis using matplotlib with intelligent tick formatting.

        Args:
            ax: The matplotlib axis to add the X-axis to
            chromosome: The chromosome label to display
        """
        # Translate options for matplotlib
        mpl_opts = translate(self.track_options, target="mpl")
        text_opts = mpl_opts.get("text", {})

        # Get axis limits to calculate optimal ticks
        xlim = ax.get_xlim()

        # Create and configure the formatter
        formatter = BPFormatter()
        formatter.set_range(xlim[0], xlim[1])

        # Set up the ticks
        ax.xaxis.set_major_locator(FixedLocator(formatter.tick_locations))
        ax.xaxis.set_major_formatter(formatter)

        # Configure axis appearance
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.yaxis.set_visible(False)

        # Add chromosome label
        font_size = self.track_options.get("font_size", self.DEFAULT_FONT_SIZE)
        if self.track_options.get("verbose_label", True):
            label_text = f"Chromosome {chromosome} ({formatter.unit_label})"
            ax.set_xlabel(label_text, fontsize=font_size, **text_opts)

        # Style the tick labels
        ax.tick_params(
            axis="x",
            which="major",
            bottom=True,
            top=False,
            labelbottom=True,
            labelsize=font_size * self.TICK_FONT_SCALE,
        )

    def add_verbose_axis_plotly(
        self,
        fig: go.Figure,
        row: int,
        col: int,
        chromosome: str,
        x_range: Optional[tuple] = None,
    ) -> None:
        """
        Add a verbose X-axis using Plotly with intelligent tick formatting.

        Args:
            fig: The Plotly figure to add the X-axis to
            row: The row number in the subplot grid
            col: The column number in the subplot grid
            chromosome: The chromosome label to display
            x_range: Tuple of (min, max) x-axis range
        """
        # Translate options for Plotly (if needed for future use)

        # Calculate optimal ticks
        if x_range:
            tick_locations, unit_label, scale_factor = (
                GenomicTickCalculator.calculate_optimal_ticks(x_range[0], x_range[1])
            )
            tick_labels = format_genomic_ticks(tick_locations, unit_label, scale_factor)
        else:
            # Use default ticks if no range provided
            tick_locations = None
            tick_labels = None
            unit_label = self.DEFAULT_UNIT_LABEL

        if self.track_options.get("axis.show_chromosome", True):
            axis_title = f"Chromosome {chromosome} ({unit_label})" if unit_label else f"Chromosome {chromosome}"
        else:
            axis_title = None

        # Add empty trace to establish the subplot (required by Plotly)
        fig.add_trace(go.Scatter(x=[], y=[], showlegend=False), row=row, col=col)

        # Configure x-axis update parameters
        xaxis_params = {
            "title_text": axis_title,
            "tickmode": "array",  # Use array mode for custom ticks
            "tickvals": tick_locations,
            "ticktext": tick_labels,
            "ticks": "outside",     # Show tick marks outside
            "tickangle": 0,
            "ticklen": 5,           # Length of tick marks
            "tickwidth": 1,         # Width of tick marks
            "linecolor": "black",   # Color of axis line
            "tickcolor": "black",   # Color of tick marks
            "linewidth": 1,         # Width of axis line
            "row": row,
            "col": col,
        }

        # Set range if provided
        if x_range:
            xaxis_params["range"] = [x_range[0], x_range[1]]

        # Update the x-axis
        fig.update_xaxes(**xaxis_params)

        # Hide y-axis for this track (axis tracks don't need y-axis)
        fig.update_yaxes(
            visible=False,
            row=row,
            col=col,
        )

    def plot_matplotlib(
        self, ax: Axes, start: int, end: int, chromosome: str = "Unknown", **kwargs
    ) -> None:
        """
        Plot the X-axis using matplotlib.

        Args:
            ax: The matplotlib axis to add the X-axis to
            start: Start position (not used for axis display)
            end: End position (not used for axis display)
            chromosome: The chromosome name to display
            **kwargs: Additional plotting parameters
        """
        axis_type = self.track_options.get("axis.type", self.DEFAULT_AXIS_TYPE)

        if axis_type == "verbose":
            self.add_verbose_axis_mpl(ax, chromosome)
        else:
            # Simple axis - just show ticks
            ax.tick_params(axis="x", which="major", bottom=True)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.yaxis.set_visible(False)

    def plot_plotly(
        self,
        fig: go.Figure,
        row: int,
        col: int,
        start: int,
        end: int,
        chromosome: str = "Unknown",
        **kwargs,
    ) -> None:
        """
        Plot the X-axis using Plotly.

        Args:
            fig: The Plotly figure to add the plot to
            row: The row number in the subplot grid
            col: The column number in the subplot grid
            start: Start position
            end: End position
            chromosome: The chromosome name to display
            **kwargs: Additional plotting parameters
        """
        axis_type = self.track_options.get("axis.type", self.DEFAULT_AXIS_TYPE)
        x_range = (start, end)

        if axis_type == "verbose":
            self.add_verbose_axis_plotly(fig, row, col, chromosome, x_range)
        else:
            # Simple axis - add empty trace and configure x-axis properly
            fig.add_trace(go.Scatter(x=[], y=[], showlegend=False), row=row, col=col)
            
            # Configure x-axis with explicit visibility (no longer overriding backend defaults)
            fig.update_xaxes(
                showgrid=True,
                showline=True,       # Explicitly show x-axis line
                showticklabels=True, # Explicitly show tick labels
                ticks="outside",     # Show tick marks outside
                ticklen=5,           # Length of tick marks
                tickwidth=1,         # Width of tick marks
                tickcolor="black",   # Color of tick marks
                linecolor="black",   # Color of axis line
                linewidth=1,         # Width of axis line
                row=row,
                col=col,
            )

        # Hide y-axis for this track
        fig.update_yaxes(
            visible=False,
            row=row,
            col=col,
        )

        # Hide y-axis for this track
        fig.update_yaxes(
            visible=False,
            row=row,
            col=col,
        )
