"""
BED track implementation for genomic region visualization.

This module provides the BedTrack class for visualizing genomic regions
from BED format files. It supports both matplotlib and Plotly backends
with comprehensive styling options.

Author: Assistant
Date: 2024
"""

from dataclasses import dataclass
from typing import Dict, Optional

import pandas as pd
import plotly.graph_objects as go
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle

from jdgenometracks.tracks.base_track import BaseTrack
from jdgenometracks.utils.options import translate


@dataclass
class BedTrack(BaseTrack):
    """
    A class for plotting genomic regions from BED files.

    BED tracks visualize genomic intervals as rectangles, with optional
    labels and color coding. Supports standard BED format with up to
    12 columns and various styling options.
    """

    # BED-specific constants
    BED_COLUMNS = [
        "chrom",
        "chromStart",
        "chromEnd",
        "name",
        "score",
        "strand",
        "thickStart",
        "thickEnd",
        "itemRGB",
        "blockCount",
        "blockSizes",
        "blockStarts",
    ]

    # BED plotting constants
    DEFAULT_RECT_HEIGHT = 1.0
    DEFAULT_RECT_PADDING = 0.0
    DEFAULT_BED_Y_VALUE = 0.0

    # Color and styling constants
    RGB_COLOR_MAX = 255  # Maximum RGB color value
    RGB_COMPONENT_COUNT = 3  # Expected number of RGB components
    DEFAULT_FALLBACK_COLOR = "blue"  # Default color when parsing fails
    DEFAULT_FILL_ALPHA = 0.3  # Default transparency for fill areas

    # Layout and positioning constants
    RECT_PADDING_MULTIPLIER = 2  # Multiplier for rectangle padding calculations
    PLOTLY_COORDS_PER_SHAPE = 5  # Number of coordinate points per Plotly shape
    MIDPOINT_DIVISOR = 2  # Divisor for midpoint calculations

    def __post_init__(self):
        """Initialize the BED track."""
        super().__post_init__()
        # Data loading is now handled by load_data() method when needed

    def load_data(self, file_path: Optional[str] = None) -> pd.DataFrame:
        """
        Load BED format data from file.

        Args:
            file_path: Path to BED file (uses self.file_path if None)

        Returns:
            DataFrame with standardized BED columns

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
            expected_columns=self.BED_COLUMNS,
            max_cols=None,  # Allow any number of columns, ignore extras
        )

        # Apply data type conversions
        data = self._apply_bed_dtypes(data)

        return data

    def validate_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Validate and clean BED data.

        Args:
            data: Raw BED data

        Returns:
            Validated and cleaned data

        Raises:
            ValueError: If data is invalid
        """
        # Use common validation methods
        self._validate_required_columns(data, ["chrom", "chromStart", "chromEnd"])
        self._validate_coordinate_ranges(data)
        return data

    def _apply_bed_dtypes(self, data: pd.DataFrame) -> pd.DataFrame:
        """Apply appropriate data types to BED columns."""
        dtype_map = {"chromStart": int, "chromEnd": int, "score": float, "chrom": str}
        return self._apply_column_dtypes(data, dtype_map)

    def plot_matplotlib(self, ax: Axes, start: int, end: int, **kwargs) -> None:
        """
        Plot the BED data using matplotlib.

        Draws rectangles for genomic regions with optional labels and styling.

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

        # Translate options for matplotlib
        mpl_opts = translate(self.track_options, target="mpl")

        # Extract specific option categories
        rect_opts = mpl_opts.get("marker", {})
        text_opts = mpl_opts.get("text", {})
        legend_opts = mpl_opts.get("legend", {})

        # Get layout parameters
        rect_height = self.track_options.get("rect.height", self.DEFAULT_RECT_HEIGHT)
        rect_padding = self.track_options.get("rect.padding", self.DEFAULT_RECT_PADDING)

        # Calculate y-levels for overlapping regions
        track_y_levels = self._calculate_y_levels(
            cleaned_data, rect_height, rect_padding
        )

        # Plot each region
        max_y = 0
        for region_idx, (_, region) in enumerate(cleaned_data.iterrows()):
            # Use itemRGB column for color if available and requested
            current_rect_opts = rect_opts.copy()
            if (
                "itemRGB" in cleaned_data.columns
                and self.track_options.get("use_color_column", False)
                and pd.notna(region["itemRGB"])
            ):
                try:
                    rgb_values = [
                        int(val) / self.RGB_COLOR_MAX
                        for val in str(region["itemRGB"]).split(",")
                    ]
                    if len(rgb_values) == self.RGB_COMPONENT_COUNT:
                        current_rect_opts["color"] = rgb_values
                except (ValueError, AttributeError):
                    pass  # Use default color if parsing fails

            # Get y-level for this region
            y = track_y_levels.get(region_idx, self.DEFAULT_BED_Y_VALUE)
            max_y = max(max_y, y + rect_height)

            # Create and add rectangle
            region_rect = Rectangle(
                (region["chromStart"], y + rect_padding),
                region["chromEnd"] - region["chromStart"],
                rect_height - self.RECT_PADDING_MULTIPLIER * rect_padding,
                label=region.get("name", f"Region_{region_idx}"),
                **current_rect_opts,
            )
            ax.add_patch(region_rect)

            # Add labels if requested
            label_alignment = self.track_options.get("label.alignment", None)
            region_name = region.get("name", "")

            if label_alignment and region_name:
                if label_alignment == "above":
                    ax.text(
                        (region["chromEnd"] + region["chromStart"])
                        / self.MIDPOINT_DIVISOR,
                        y + rect_height + rect_padding,
                        region_name,
                        **text_opts,
                    )
                elif label_alignment == "left":
                    ax.text(
                        region["chromStart"],
                        y + rect_height / self.MIDPOINT_DIVISOR,
                        region_name,
                        **text_opts,
                    )
                elif label_alignment == "right":
                    ax.text(
                        region["chromEnd"],
                        y + rect_height / self.MIDPOINT_DIVISOR,
                        region_name,
                        **text_opts,
                    )

        # Style the axes
        ax.xaxis.set_tick_params(bottom=False)
        ax.yaxis.set_tick_params(left=False, labelleft=False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

        # Set y-limits
        if self.track_options.get("use_global_max", False):
            ax.set_ylim(0, max_y + rect_padding)
        else:
            ax.autoscale(enable=True, axis="y")

        # Add legend if requested
        if self.show_legend:
            ax.legend(**legend_opts)

        # Add horizontal lines
        self._add_horizontal_lines_matplotlib(ax)

    def _calculate_y_levels(
        self, data: pd.DataFrame, rect_height: float, rect_padding: float
    ) -> Dict[int, float]:
        """
        Calculate non-overlapping y-levels for BED regions.

        Args:
            data: BED data
            rect_height: Height of rectangles
            rect_padding: Padding between rectangles

        Returns:
            Dictionary mapping region index to y-level
        """
        if data.empty:
            return {}

        # Sort by start position for level assignment
        sorted_data = data.sort_values("chromStart").reset_index(drop=True)
        levels = {}
        level_ends = []  # Track the end position of each level

        for idx, (_, region) in enumerate(sorted_data.iterrows()):
            start, end = region["chromStart"], region["chromEnd"]

            # Find the first level where this region can fit
            assigned_level = 0
            for level, level_end in enumerate(level_ends):
                if start >= level_end:  # No overlap
                    level_ends[level] = end
                    assigned_level = level
                    break
            else:
                # Need a new level
                assigned_level = len(level_ends)
                level_ends.append(end)

            # Calculate y-position for this level
            y_pos = assigned_level * (rect_height + rect_padding)
            levels[idx] = y_pos

        return levels

    def plot_plotly(
        self, fig: go.Figure, row: int, col: int, start: int, end: int, **kwargs
    ) -> None:
        """
        Plot the BED data using Plotly.

        Creates rectangular shapes for genomic regions with optional
        labels and interactive features.

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

        # Translate options for Plotly
        plotly_opts = translate(self.track_options, target="plotly")

        # Get layout parameters
        rect_height = self.track_options.get("rect.height", self.DEFAULT_RECT_HEIGHT)
        rect_padding = self.track_options.get("rect.padding", self.DEFAULT_RECT_PADDING)

        # Calculate y-levels for overlapping regions
        track_y_levels = self._calculate_y_levels(
            cleaned_data, rect_height, rect_padding
        )

        # Add individual traces for each region (like V2) to ensure proper hover text
        for region_idx, (_, region) in enumerate(cleaned_data.iterrows()):
            # Get y-level for this region
            y = track_y_levels.get(region_idx, self.DEFAULT_BED_Y_VALUE)

            # Get region name
            region_name = region.get("name", f"Region_{region_idx}")

            # Handle colors
            if (
                "itemRGB" in cleaned_data.columns
                and self.track_options.get("use_color_column", False)
                and pd.notna(region["itemRGB"])
            ):
                try:
                    rgb_values = [int(val) for val in str(region["itemRGB"]).split(",")]
                    if len(rgb_values) == self.RGB_COMPONENT_COUNT:
                        fillcolor = f"rgb({rgb_values[0]}, {rgb_values[1]}, {rgb_values[2]})"
                        linecolor = fillcolor
                    else:
                        fillcolor = plotly_opts.get("marker", {}).get("color", self.DEFAULT_FALLBACK_COLOR)
                        linecolor = fillcolor
                except (ValueError, AttributeError):
                    fillcolor = plotly_opts.get("marker", {}).get("color", self.DEFAULT_FALLBACK_COLOR)
                    linecolor = fillcolor
            else:
                fillcolor = plotly_opts.get("marker", {}).get("color", self.DEFAULT_FALLBACK_COLOR)
                linecolor = fillcolor

            # Create individual trace for this region (like V2)
            fig.add_trace(
                go.Scatter(
                    x=[
                        region["chromStart"],
                        region["chromStart"],
                        region["chromEnd"],
                        region["chromEnd"],
                        region["chromStart"],
                    ],
                    y=[
                        y + rect_padding,
                        y + rect_height - rect_padding,
                        y + rect_height - rect_padding,
                        y + rect_padding,
                        y + rect_padding,
                    ],
                    mode="lines",
                    fill="toself",
                    fillcolor=fillcolor,
                    line=dict(color=linecolor),
                    name=region_name,  # Individual name for each region
                    showlegend=self.show_legend,
                ),
                row=row,
                col=col,
            )

        # Add horizontal lines
        self._add_horizontal_lines_plotly(fig, row, col)

        # Configure axes for this track - follow V2 pattern exactly
        fig.update_xaxes(
            row=row,
            col=col,
        )
        fig.update_yaxes(
            showticklabels=False,  # Hide y-axis tick labels (like V2)
            row=row,
            col=col,
        )
