"""
Layout optimization and grid management for jdgenometracks V3.

Provides efficient algorithms for:
- Assigning y-coordinates to genomic regions to prevent overlaps
- Managing subplot grids and positioning
- Optimizing track layout and spacing

Consolidates layout_optimizer.py functionality with improved algorithms.
"""

import heapq
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from .constants import PlotConstants


def assign_y_levels_sweep_line(
    regions_df: pd.DataFrame, rect_height: float = 1.0, rect_padding: float = 0.0
) -> List[float]:
    """
    Assign y-levels to BED regions using a sweep line algorithm to prevent overlaps.

    This algorithm is more efficient than the sparse matrix approach:
    - Time complexity: O(n log n) where n is the number of regions
    - Space complexity: O(n) instead of O(genomic_range_width)
    - Optimal y-level assignment: Uses minimal number of y-levels needed

    Args:
        regions_df: DataFrame with 'chromStart' and 'chromEnd' columns
        rect_height: Height of each rectangle
        rect_padding: Padding between rectangles

    Returns:
        List of y-values for each region in the same order as input DataFrame
    """
    if regions_df.empty:
        return []

    # Create events for start and end of each region
    events = []
    for idx, region in regions_df.iterrows():
        events.append((region["chromStart"], "start", idx))
        events.append((region["chromEnd"], "end", idx))

    # Sort events by position, with 'end' events before 'start' events at same position
    # This ensures that when a region ends at position X and another starts at X,
    # the y-level is freed before being potentially reused
    events.sort(key=lambda x: (x[0], x[1] == "start"))

    # Track active regions and their y-levels
    active_regions = {}  # region_idx -> y_level
    available_y_levels = []  # min-heap of available y-levels
    next_y_level = 0.0

    # Result array to store y-level for each region
    y_levels = [0.0] * len(regions_df)

    for position, event_type, region_idx in events:
        if event_type == "start":
            # Assign y-level to this region
            if available_y_levels:
                y_level = heapq.heappop(available_y_levels)
            else:
                y_level = next_y_level
                next_y_level += rect_height + rect_padding

            active_regions[region_idx] = y_level
            y_levels[region_idx] = y_level

        else:  # event_type == "end"
            # Free up the y-level used by this region
            if region_idx in active_regions:
                y_level = active_regions.pop(region_idx)
                heapq.heappush(available_y_levels, y_level)

    return y_levels


def calculate_grid_dimensions(tracks: List[Any]) -> Tuple[int, int]:
    """
    Calculate the grid dimensions needed for a list of tracks.

    Args:
        tracks: List of track objects with subplot_x and subplot_y attributes

    Returns:
        Tuple of (max_rows, max_cols)
    """
    if not tracks:
        return 1, 1

    max_rows = max(track.subplot_y for track in tracks if track is not None) + 1
    max_cols = max(track.subplot_x for track in tracks if track is not None) + 1

    return max_rows, max_cols


def normalize_layout_proportions(proportions: List[float]) -> List[float]:
    """
    Normalize a list of proportions to sum to 1.0.

    Args:
        proportions: List of proportion values

    Returns:
        Normalized proportions that sum to 1.0
    """
    if not proportions:
        return []

    total = sum(proportions)
    if total <= 0:
        # If all proportions are 0 or negative, return equal proportions
        return [1.0 / len(proportions)] * len(proportions)

    return [prop / total for prop in proportions]


def create_subplot_grid_spec(
    max_rows: int,
    max_cols: int,
    height_props: List[float],
    width_props: List[float],
    vertical_spacing: float = PlotConstants.VERTICAL_SPACING,
    horizontal_spacing: float = PlotConstants.HORIZONTAL_SPACING,
) -> Dict[str, Any]:
    """
    Create subplot grid specification for both matplotlib and plotly.

    Args:
        max_rows: Number of rows in grid
        max_cols: Number of columns in grid
        height_props: Relative heights for each row
        width_props: Relative widths for each column
        vertical_spacing: Spacing between rows
        horizontal_spacing: Spacing between columns

    Returns:
        Dictionary with grid specification parameters
    """
    # Normalize proportions
    norm_height_props = normalize_layout_proportions(height_props)
    norm_width_props = normalize_layout_proportions(width_props)

    return {
        "rows": max_rows,
        "cols": max_cols,
        "row_heights": norm_height_props,
        "column_widths": norm_width_props,
        "vertical_spacing": vertical_spacing,
        "horizontal_spacing": horizontal_spacing,
    }


def optimize_track_spacing(
    tracks: List[Any], column_regions: List[Optional[str]] = None
) -> Dict[int, float]:
    """
    Optimize vertical spacing between tracks in the same column.

    Args:
        tracks: List of track objects
        column_regions: Optional list of genomic regions per column

    Returns:
        Dictionary mapping column index to optimal spacing
    """
    max_rows, max_cols = calculate_grid_dimensions(tracks)
    column_spacing = {}

    for col_idx in range(max_cols):
        # Count tracks in this column
        col_tracks = [
            track
            for track in tracks
            if track is not None and track.subplot_x == col_idx
        ]

        # Calculate optimal spacing based on track count and types
        if len(col_tracks) <= 1:
            spacing = PlotConstants.VERTICAL_SPACING
        else:
            # More tracks = tighter spacing
            spacing = max(
                PlotConstants.VERTICAL_SPACING * 0.5,
                PlotConstants.VERTICAL_SPACING / len(col_tracks),
            )

        column_spacing[col_idx] = spacing

    return column_spacing


class LayoutManager:
    """Manages layout calculations and optimizations for track plots."""

    def __init__(
        self, tracks: List[Any], figure_options: Optional[Dict[str, Any]] = None
    ):
        self.tracks = tracks
        self.figure_options = figure_options or {}
        self.max_rows, self.max_cols = calculate_grid_dimensions(tracks)

    def get_layout_parameters(self) -> Dict[str, Any]:
        """
        Get all layout parameters with defaults applied.

        Returns:
            Dictionary with complete layout specification
        """
        # Get parameters from figure options with defaults
        height_props = self.figure_options.get(
            "height_props", [PlotConstants.DEFAULT_HEIGHT_PROP] * self.max_rows
        )
        width_props = self.figure_options.get(
            "width_props", [PlotConstants.DEFAULT_WIDTH_PROP] * self.max_cols
        )
        row_titles = self.figure_options.get("row_titles", [""] * self.max_rows)
        column_titles = self.figure_options.get("column_titles", [""] * self.max_cols)

        # Validate layout
        from .validation import validate_grid_layout

        validate_grid_layout(
            self.tracks, height_props, width_props, row_titles, column_titles
        )

        # Create grid specification
        grid_spec = create_subplot_grid_spec(
            self.max_rows, self.max_cols, height_props, width_props
        )

        return {
            "max_rows": self.max_rows,
            "max_cols": self.max_cols,
            "height_props": height_props,
            "width_props": width_props,
            "row_titles": row_titles,
            "column_titles": column_titles,
            "grid_spec": grid_spec,
            "total_height": self.figure_options.get("total_height"),
            "total_width": self.figure_options.get("total_width"),
            "plot_title": self.figure_options.get("plot_title"),
        }

    def optimize_bed_track_layout(
        self, column_regions: List[Optional[str]]
    ) -> Dict[int, Dict[int, float]]:
        """
        Optimize y-level assignments for BED tracks to prevent overlaps.

        Args:
            column_regions: List of genomic regions for each column

        Returns:
            Dictionary mapping column_idx -> {track_id: {region_idx: y_level}}
        """
        bed_y_levels = {}

        # Process each column separately
        for col_idx in range(self.max_cols):
            column_y_levels = self._optimize_column_bed_layout(col_idx, column_regions)
            bed_y_levels[col_idx] = column_y_levels

        return bed_y_levels

    def _optimize_column_bed_layout(
        self, col_idx: int, column_regions: List[Optional[str]]
    ) -> Dict[int, Dict[int, float]]:
        """Optimize BED track layout for a single column."""
        # Collect all BED regions for this column across all tracks
        all_regions = []
        track_region_mapping = (
            {}
        )  # Maps global region index to (track, local_region_index)

        region_counter = 0
        for track in self.tracks:
            if track is None or track.subplot_x != col_idx:
                continue

            # Check if this is a BED track
            if (
                hasattr(track, "__class__")
                and "bed" in track.__class__.__name__.lower()
            ):
                # Get formatted data for this track
                if hasattr(track, "format_data"):
                    region_str = (
                        column_regions[col_idx]
                        if col_idx < len(column_regions)
                        else None
                    )
                    cleaned_data = track.format_data(
                        subset_region=region_str, axis_shift=None
                    )

                    # Add regions to global list
                    for local_idx, (_, region) in enumerate(cleaned_data.iterrows()):
                        all_regions.append(
                            {
                                "chromStart": region["chromStart"],
                                "chromEnd": region["chromEnd"],
                            }
                        )
                        track_region_mapping[region_counter] = (track, local_idx)
                        region_counter += 1

        if not all_regions:
            return {}

        # Convert to DataFrame for the optimizer
        regions_df = pd.DataFrame(all_regions)

        # Assign y-levels globally using sweep line algorithm
        global_y_levels = assign_y_levels_sweep_line(
            regions_df,
            rect_height=PlotConstants.DEFAULT_RECT_HEIGHT,
            rect_padding=PlotConstants.DEFAULT_RECT_PADDING,
        )

        # Map y-levels back to individual tracks
        column_y_levels = {}
        for global_idx, y_level in enumerate(global_y_levels):
            track, local_idx = track_region_mapping[global_idx]
            track_id = id(track)  # Use track object ID as key
            if track_id not in column_y_levels:
                column_y_levels[track_id] = {}
            column_y_levels[track_id][local_idx] = y_level

        return column_y_levels
