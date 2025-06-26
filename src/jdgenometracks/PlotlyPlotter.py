from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import plotly.graph_objects as go
import plotly.subplots as ps

from jdgenometracks.tracks import XAxisTrack

from .base_plotter import BasePlotter
from .config import PlotDefaults, PlotlyDefaults
from .tracks.BedTrack import BedTrack
from .tracks.GenomeTrack import GenomeTrack


@dataclass
class PlotlyPlotter(BasePlotter):
    """
    A class for generating multi-track genomic plots using Plotly.

    Attributes:
        tracks: List of tracks to plot, can contain instances of GenomeTrack or its subclasses.
    """

    def _validate_track_methods(self) -> None:
        """
        Validate that tracks have required Plotly plotting methods.
        """
        for track in self.tracks:
            if track is not None and not hasattr(track, "plot_plotly"):
                raise TypeError(
                    f"Each track must have a 'plot_plotly' method (got {type(track)})."
                )

    def _initialize_subplots(
        self,
        max_rows: int,
        max_cols: int,
        height_props: list[float],
        width_props: list[float],
        row_titles: list[str],
        column_titles: list[str],
        fig_options: dict | None = None,
    ) -> go.Figure:
        """
        Initializes the Plotly subplots figure with the given layout properties.

        Args:
            max_rows: Number of rows in the subplot grid.
            max_cols: Number of columns in the subplot grid.
            height_props: Heights of each row as a proportion of the total figure height.
            width_props: Widths of each column as a proportion of the total figure width.
            row_titles: Titles for each row.
            column_titles: Titles for each column.
            fig_options: Optional figure configuration options.

        Returns:
            go.Figure: A Plotly figure with subplots initialized.
        """
        if fig_options is None:
            fig_options = {}
        subplots = ps.make_subplots(
            rows=max_rows,
            cols=max_cols,
            shared_xaxes=fig_options.get(
                "shared_xaxes", PlotlyDefaults.DEFAULT_SHARED_XAXES
            ),
            shared_yaxes=fig_options.get(
                "shared_yaxes", PlotlyDefaults.DEFAULT_SHARED_YAXES
            ),
            row_heights=height_props,
            row_titles=row_titles,
            column_widths=width_props,
            column_titles=column_titles,
            vertical_spacing=fig_options.get(
                "vertical_spacing", PlotDefaults.VERTICAL_SPACING
            ),
            horizontal_spacing=fig_options.get(
                "horizontal_spacing", PlotDefaults.HORIZONTAL_SPACING
            ),
        )
        return subplots

    def plot_all_tracks(
        self,
        fig_options: dict | None = None,
    ) -> go.Figure:
        """
        Plots all tracks into a single Plotly figure with optional customization.

        Args:
            fig_options: Optional figure configuration options.

        Returns:
            go.Figure: The complete Plotly figure with all tracks.
        """
        # Use base class to prepare common data
        max_rows, max_cols, plot_data = self._prepare_plot_data(fig_options)

        # Initialize the subplots figure
        subplots = self._initialize_subplots(
            max_rows,
            max_cols,
            plot_data["height_props"],
            plot_data["width_props"],
            plot_data["row_titles"],
            plot_data["column_titles"],
            fig_options=fig_options,
        )

        # Set layout properties for the entire figure
        subplots.update_layout(
            autosize=PlotlyDefaults.DEFAULT_AUTOSIZE,
            height=plot_data["total_height"],
            width=plot_data["total_width"],
            plot_bgcolor=(
                fig_options.get("plot_bgcolor", PlotDefaults.DEFAULT_PLOT_BGCOLOR)
                if fig_options
                else PlotDefaults.DEFAULT_PLOT_BGCOLOR
            ),
            margin=(
                fig_options.get("margin", PlotDefaults.DEFAULT_MARGIN)
                if fig_options
                else PlotDefaults.DEFAULT_MARGIN
            ),
            title=plot_data["plot_title"],
        )

        # Plot each track using the base class method
        for track in self.tracks:
            if track is None:
                continue

            backend_args = {
                "figure": subplots,
                "row": track.subplot_y + 1,
                "col": track.subplot_x + 1,
            }

            self._plot_track(track, plot_data, backend_args)

        return subplots

    def _plot_bed_track(
        self,
        track: BedTrack,
        column_limits: list[dict[str, Any]],
        column_regions: list[str | None],
        track_y_levels: dict[int, float],
        backend_args: dict[str, Any],
    ) -> None:
        """Plot a BED track using Plotly."""
        track.plot_plotly(
            backend_args["figure"],
            backend_args["row"],
            backend_args["col"],
            track_y_levels,
            subset_region=column_regions[track.subplot_x],
            xmin=column_limits[track.subplot_x]["xmin"],
        )

    def _plot_axis_track(
        self,
        track: XAxisTrack,
        column_limits: list[dict[str, Any]],
        backend_args: dict[str, Any],
    ) -> None:
        """Plot an axis track using Plotly."""
        column_limit = column_limits[track.subplot_x]
        track.plot_plotly(
            backend_args["figure"],
            backend_args["row"],
            backend_args["col"],
            chromosome=column_limit["chromosome"],
            xmin=column_limit["xmin"],
            xmax=column_limit["xmax"],
        )

    def _plot_generic_track(
        self,
        track: GenomeTrack,
        column_regions: list[str | None],
        backend_args: dict[str, Any],
    ) -> None:
        """Plot a generic track using Plotly."""
        track.plot_plotly(
            backend_args["figure"],
            backend_args["row"],
            backend_args["col"],
            subset_region=column_regions[track.subplot_x],
        )

    def _add_track_hlines(
        self,
        track: GenomeTrack,
        backend_args: dict[str, Any],
    ) -> None:
        """Add horizontal lines to a track using Plotly."""
        track.add_hlines_plotly(
            backend_args["figure"],
            backend_args["row"],
            backend_args["col"],
        )

    def _update_track_axes(
        self,
        track: GenomeTrack,
        column_limits: list[dict[str, Any]],
        backend_args: dict[str, Any],
    ) -> None:
        """Update track axes properties using Plotly."""
        xmin = column_limits[track.subplot_x]["xmin"]
        xmax = column_limits[track.subplot_x]["xmax"]

        backend_args["figure"].update_xaxes(
            range=[xmin, xmax],
            row=backend_args["row"],
            col=backend_args["col"],
        )
        backend_args["figure"].update_yaxes(
            row=backend_args["row"],
            col=backend_args["col"],
            **track.track_options.get("yaxis", {}),
        )
