from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from jdgenometracks.tracks import XAxisTrack

from .config import PlotDefaults, get_default_column_regions
from .tracks.BedTrack import BedTrack
from .tracks.GenomeTrack import GenomeTrack
from .utils import TrackUtils


@dataclass
class BasePlotter(ABC):
    """
    Base class for generating multi-track genomic plots.

    This class contains shared functionality between matplotlib and plotly plotters,
    including input validation, track management, and common plotting logic.

    Attributes:
        tracks: List of tracks to plot, can contain instances of GenomeTrack or its subclasses.
    """

    tracks: Union[List[GenomeTrack], Sequence[Any]]

    def __post_init__(self) -> None:
        """
        Validates tracks and performs common initialization.
        """
        if not isinstance(self.tracks, (list, tuple)):
            raise TypeError(f"Tracks must be a list or tuple, got {type(self.tracks)}.")
        if len(self.tracks) == 0:
            raise ValueError("Tracks cannot be empty.")

        # Validate tracks have required methods
        self._validate_track_methods()

    @abstractmethod
    def _validate_track_methods(self) -> None:
        """
        Validate that tracks have required plotting methods for the specific backend.
        Must be implemented by subclasses.
        """
        pass

    def _validate_inputs(
        self,
        max_rows: int,
        max_cols: int,
        height_props: Union[List[float], Sequence[float]],
        row_titles: Union[List[str], Sequence[str]],
        width_props: Union[List[float], Sequence[float]],
        column_titles: Union[List[str], Sequence[str]],
    ) -> None:
        """
        Validates the inputs for the plot, ensuring the lengths match the expected dimensions.

        Args:
            max_rows: Maximum number of rows in the plot grid.
            max_cols: Maximum number of columns in the plot grid.
            height_props: Proportions of the figure height for each row.
            row_titles: Titles for each row in the figure.
            width_props: Proportions of the figure width for each column.
            column_titles: Titles for each column in the figure.

        Raises:
            AssertionError: If the lengths of inputs don't match the expected number of rows or columns.
        """
        assert (
            len(height_props) == max_rows
        ), f"Number of height_props should equal {max_rows}"
        assert (
            len(row_titles) == max_rows
        ), f"Number of row_titles should equal {max_rows}"
        assert (
            len(width_props) == max_cols
        ), f"Number of width_props should equal {max_cols}"
        assert (
            len(column_titles) == max_cols
        ), f"Number of column_titles should equal {max_cols}"

    def _prepare_plot_data(
        self,
        fig_options: Optional[Dict[str, Any]] = None,
    ) -> Tuple[int, int, Dict[str, Any]]:
        """
        Prepares common data needed for plotting.

        Args:
            fig_options: Optional figure configuration options.

        Returns:
            Tuple containing:
            - max_rows: Maximum number of rows
            - max_cols: Maximum number of columns
            - plot_data: Dictionary with all prepared plotting data
        """
        if fig_options is None:
            fig_options = {}

        # Calculate grid dimensions
        max_rows = max(track.subplot_y for track in self.tracks) + 1
        max_cols = max(track.subplot_x for track in self.tracks) + 1

        # Get configuration options with defaults
        column_regions = fig_options.get(
            "column_regions", get_default_column_regions(max_cols)
        )

        # Validate column regions
        if (
            not isinstance(column_regions, (list, tuple))
            or len(column_regions) != max_cols
        ):
            raise ValueError(
                "column_regions must be a list/tuple with length equal to number of columns in tracks."
            )

        # Get layout options with defaults
        height_props = fig_options.get(
            "height_props", [PlotDefaults.DEFAULT_HEIGHT_PROP for _ in range(max_rows)]
        )
        row_titles = fig_options.get("row_titles", ["" for _ in range(max_rows)])
        width_props = fig_options.get(
            "width_props", [PlotDefaults.DEFAULT_WIDTH_PROP for _ in range(max_cols)]
        )
        column_titles = fig_options.get("column_titles", ["" for _ in range(max_cols)])
        total_height = fig_options.get("total_height", None)
        total_width = fig_options.get("total_width", None)
        plot_title = fig_options.get("plot_title", None)
        relative_x_axis = fig_options.get("relative_x_axis", False)

        # Validate inputs
        self._validate_inputs(
            max_rows, max_cols, height_props, row_titles, width_props, column_titles
        )

        # Calculate column limits
        column_limits = self._calculate_column_limits(
            max_cols, column_regions, relative_x_axis
        )

        # Preprocess BED tracks to assign y-levels globally
        bed_y_levels = self._preprocess_bed_tracks(column_limits, column_regions)

        plot_data = {
            "column_regions": column_regions,
            "height_props": height_props,
            "row_titles": row_titles,
            "width_props": width_props,
            "column_titles": column_titles,
            "total_height": total_height,
            "total_width": total_width,
            "plot_title": plot_title,
            "relative_x_axis": relative_x_axis,
            "column_limits": column_limits,
            "bed_y_levels": bed_y_levels,
        }

        return max_rows, max_cols, plot_data

    def _preprocess_bed_tracks(
        self,
        column_limits: list[dict[str, Any]],
        column_regions: list[str | None] | tuple[str | None, ...],
    ) -> dict[str, dict[int, float]]:
        """
        Preprocess all BED tracks to assign y-levels globally, ensuring no overlaps
        within or between tracks in the same column.

        Args:
            column_limits: List of column limit dictionaries
            column_regions: List of region strings for each column

        Returns:
            Dictionary mapping column index to track y-level assignments
        """
        from .layout_optimizer import assign_y_levels_sweep_line

        bed_y_levels = {}

        # Process each column separately
        for col_idx in range(len(column_limits)):
            # Collect all BED regions for this column across all tracks
            all_regions = []
            track_region_mapping = (
                {}
            )  # Maps global region index to (track, local_region_index)

            region_counter = 0
            for track in self.tracks:
                if track is None:
                    continue

                # Check if this is a BED track using class name to avoid import issues
                is_bed_track = (
                    isinstance(track, BedTrack)
                    or track.__class__.__name__ == "BedTrack"
                )

                if not is_bed_track or track.subplot_x != col_idx:
                    continue

                # Get formatted data for this track
                cleaned_data = track.format_data(
                    subset_region=column_regions[col_idx], axis_shift=None
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
                bed_y_levels[col_idx] = {}
                continue

            # Convert to DataFrame for the optimizer
            import pandas as pd

            regions_df = pd.DataFrame(all_regions)

            # Assign y-levels globally using sweep line algorithm
            global_y_levels = assign_y_levels_sweep_line(
                regions_df,
                rect_height=PlotDefaults.DEFAULT_RECT_HEIGHT,
                rect_padding=PlotDefaults.DEFAULT_RECT_PADDING,
            )

            # Map y-levels back to individual tracks
            column_y_levels = {}
            for global_idx, y_level in enumerate(global_y_levels):
                track, local_idx = track_region_mapping[global_idx]
                track_id = id(track)  # Use track object ID as key
                if track_id not in column_y_levels:
                    column_y_levels[track_id] = {}
                column_y_levels[track_id][local_idx] = y_level

            bed_y_levels[col_idx] = column_y_levels

        return bed_y_levels

    def _calculate_column_limits(
        self,
        max_cols: int,
        column_regions: Union[List[Optional[str]], Sequence[Optional[str]]],
        relative_x_axis: bool,
    ) -> List[Dict[str, Any]]:
        """
        Calculate limits for each column of tracks.

        Args:
            max_cols: Number of columns in the plot grid.
            column_regions: List of genomic regions for each column.
            relative_x_axis: Whether to use relative x-axis coordinates.

        Returns:
            List of dictionaries containing column limit information.
        """
        column_limits = []
        for column in range(max_cols):
            col_tracks = [
                track
                for track in self.tracks
                if track is not None and track.subplot_x == column
            ]
            chromosome, xmin, xmax, max_bed_regions, axis_shift = (
                TrackUtils.get_column_limits(
                    col_tracks, column_regions[column], relative_x_axis
                )
            )
            column_limits.append(
                {
                    "chromosome": chromosome,
                    "xmin": xmin,
                    "xmax": xmax,
                    "max_bed_regions": max_bed_regions,
                    "axis_shift": axis_shift,
                }
            )
        return column_limits

    def _plot_track(
        self,
        track: GenomeTrack,
        plot_data: Dict[str, Any],
        backend_specific_args: Dict[str, Any],
    ) -> None:
        """
        Plot a single track using backend-specific logic.

        Args:
            track: The track to plot.
            plot_data: Common plotting data from _prepare_plot_data.
            backend_specific_args: Backend-specific arguments (figure, axes, etc.).
        """
        column_limits = plot_data["column_limits"]
        column_regions = plot_data["column_regions"]
        bed_y_levels = plot_data["bed_y_levels"]

        try:
            # Use class name check to handle different import paths
            if isinstance(track, BedTrack) or track.__class__.__name__ == "BedTrack":
                # Get precomputed y-levels for this track
                track_y_levels = bed_y_levels.get(track.subplot_x, {}).get(
                    id(track), {}
                )
                self._plot_bed_track(
                    track,  # type: ignore
                    column_limits,
                    column_regions,
                    track_y_levels,
                    backend_specific_args,
                )
            elif isinstance(track, XAxisTrack):
                self._plot_axis_track(track, column_limits, backend_specific_args)
            else:
                self._plot_generic_track(track, column_regions, backend_specific_args)

            # Add horizontal lines
            self._add_track_hlines(track, backend_specific_args)

            # Update axis properties
            self._update_track_axes(track, column_limits, backend_specific_args)

        except Exception as e:
            raise RuntimeError(
                f"Error plotting track at grid ({track.subplot_y},{track.subplot_x}): {e}"
            )

    @abstractmethod
    def _plot_bed_track(
        self,
        track: BedTrack,
        column_limits: List[Dict[str, Any]],
        column_regions: Union[List[Optional[str]], Sequence[Optional[str]]],
        track_y_levels: Dict[int, float],
        backend_args: Dict[str, Any],
    ) -> None:
        """Plot a BED track using backend-specific methods."""
        pass

    @abstractmethod
    def _plot_axis_track(
        self,
        track: XAxisTrack,
        column_limits: List[Dict[str, Any]],
        backend_args: Dict[str, Any],
    ) -> None:
        """Plot an axis track using backend-specific methods."""
        pass

    @abstractmethod
    def _plot_generic_track(
        self,
        track: GenomeTrack,
        column_regions: Union[List[Optional[str]], Sequence[Optional[str]]],
        backend_args: Dict[str, Any],
    ) -> None:
        """Plot a generic track using backend-specific methods."""
        pass

    @abstractmethod
    def _add_track_hlines(
        self,
        track: GenomeTrack,
        backend_args: Dict[str, Any],
    ) -> None:
        """Add horizontal lines to a track using backend-specific methods."""
        pass

    @abstractmethod
    def _update_track_axes(
        self,
        track: GenomeTrack,
        column_limits: List[Dict[str, Any]],
        backend_args: Dict[str, Any],
    ) -> None:
        """Update track axes properties using backend-specific methods."""
        pass

    @abstractmethod
    def plot_all_tracks(
        self,
        fig_options: Optional[Dict[str, Any]] = None,
    ) -> Any:
        """
        Plot all tracks. Must be implemented by subclasses.

        Args:
            fig_options: Optional figure configuration options.

        Returns:
            The figure object (type depends on backend).
        """
        pass
