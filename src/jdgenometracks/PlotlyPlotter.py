from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import plotly.graph_objects as go
import plotly.subplots as ps
import scipy.sparse as sp

from jdgenometracks.tracks import XAxisTrack

from .tracks.BedTrack import BedTrack
from .tracks.GenomeTrack import GenomeTrack
from .utils import TrackUtils
from .utils_units import convert_to_pixels, parse_size_with_units


@dataclass
class PlotlyPlotter:
    """
    A class for generating multi-track genomic plots using Plotly.

    Attributes:
        tracks (np.ndarray): Array of tracks to plot, can contain instances of GenomeTrack or its subclasses.
        total_height (float): The total height of the figure in inches.
    """

    tracks: list[GenomeTrack]

    def _validate_inputs(
        self,
        max_rows: int,
        max_cols: int,
        height_props: list[float],
        row_titles: list[str],
        width_props: list[float],
        column_titles: list[str],
    ) -> None:
        """
        Validates the inputs for the plot, ensuring the lengths match the expected dimensions.

        Args:
            height_props (list[float]): Proportions of the figure height for each row.
            row_titles (list[str]): Titles for each row in the figure.
            width_props (list[float]): Proportions of the figure width for each column.
            column_titles (list[str]): Titles for each column in the figure.

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

    def _initialize_subplots(
        self,
        max_rows,
        max_cols,
        height_props: list[float],
        width_props: list[float],
        row_titles: list[str],
        column_titles: list[str],
        fig_options: dict | None = None,
    ) -> go.Figure:
        """
        Initializes the Plotly subplots figure with the given layout properties.

        Args:
            height_props (list[float]): Heights of each row as a proportion of the total figure height.
            width_props (list[float]): Widths of each column as a proportion of the total figure width.
            row_titles (list[str]): Titles for each row.
            column_titles (list[str]): Titles for each column.
            num_distinct_rows (int): The number of distinct rows in the plot.
            vertical_spacing (float): Vertical spacing between rows.
            horizontal_spacing (float): Horizontal spacing between columns.
            shared_xaxes (str): Determines if x-axes are shared across columns (can be 'columns', 'all', or None).
            shared_yaxes (str): Determines if y-axes are shared across rows (can be 'rows', 'all', or None).

        Returns:
            go.Figure: A Plotly figure with subplots initialized.
        """
        if fig_options is None:
            fig_options = {}
        subplots = ps.make_subplots(
            rows=max_rows,
            cols=max_cols,
            shared_xaxes=fig_options.get("shared_xaxes", "columns"),
            shared_yaxes=fig_options.get("shared_yaxes", "rows"),
            row_heights=height_props,
            row_titles=row_titles,
            column_widths=width_props,
            column_titles=column_titles,
            vertical_spacing=fig_options.get("vertical_spacing", 0.02),
            horizontal_spacing=fig_options.get("horizontal_spacing", 0.05),
        )
        return subplots

    def plot_all_tracks(
        self,
        fig_options: dict | None = None,
    ) -> go.Figure:
        """
        Plots all tracks into a single Plotly figure with optional customization.

        Args:


        Returns:
            go.Figure: The complete Plotly figure with all tracks.
        """
        """
        Plots all tracks in a single matplotlib figure.
        """
        if fig_options is None:
            fig_options = {}

        max_rows = max(track.subplot_y for track in self.tracks) + 1
        max_cols = max(track.subplot_x for track in self.tracks) + 1
        column_regions = fig_options.get(
            "column_regions", [None for _ in range(max_cols)]
        )

        if (
            not isinstance(column_regions, (list, tuple))
            or len(column_regions) != max_cols
        ):
            raise ValueError(
                "column_regions must be a list/tuple with length equal to number of columns in tracks."
            )

        height_props = fig_options.get("height_props", [1 for _ in range(max_rows)])
        row_titles = fig_options.get("row_titles", ["" for _ in range(max_rows)])
        width_props = fig_options.get("width_props", [1 for _ in range(max_cols)])
        column_titles = fig_options.get("column_titles", ["" for _ in range(max_cols)])
        total_height = fig_options.get("total_height", None)
        total_width = fig_options.get("total_width", None)
        plot_title = fig_options.get("plot_title", None)
        relative_x_axis = fig_options.get("relative_x_axis", False)

        self._validate_inputs(
            max_rows, max_cols, height_props, row_titles, width_props, column_titles
        )
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
        bed_region_coverages = [
            sp.csr_matrix((1, column["xmax"] - column["xmin"]))
            for column in column_limits
        ]
        # Initialize the subplots figure with custom spacings and shared axes settings
        subplots = self._initialize_subplots(
            max_rows,
            max_cols,
            height_props,
            width_props,
            row_titles,
            column_titles,
            fig_options=fig_options,
        )

        # Set layout properties for the entire figure
        subplots.update_layout(
            autosize=True,
            height=total_height,
            width=total_width,
            plot_bgcolor=fig_options.get("plot_bgcolor", "white"),
            margin=fig_options.get("margin", dict(l=0.1, r=0.1, t=50, b=20, pad=4)),
            title=plot_title,
        )

        # Plot each track
        for track in self.tracks:
            if track is None:
                continue
            if not hasattr(track, "plot_plotly"):
                raise TypeError(f"Track {track} does not have a 'plot_plotly' method.")
            try:
                if isinstance(track, BedTrack):
                    bed_region_coverages[track.subplot_x] = track.plot_plotly(
                        subplots,
                        track.subplot_y + 1,
                        track.subplot_x + 1,
                        bed_region_coverages[track.subplot_x],
                        subset_region=column_regions[track.subplot_x],
                        xmin=column_limits[track.subplot_x]["xmin"],
                    )
                elif isinstance(track, XAxisTrack):
                    track.plot_plotly(
                        subplots,
                        track.subplot_y + 1,
                        track.subplot_x + 1,
                        chromosome=column_limits[track.subplot_x]["chromosome"],
                    )
                else:
                    track.plot_plotly(
                        subplots,
                        track.subplot_y + 1,
                        track.subplot_x + 1,
                        subset_region=column_regions[track.subplot_x],
                    )
                track.add_hlines_plotly(
                    subplots, track.subplot_y + 1, track.subplot_x + 1
                )
                subplots.update_xaxes(
                    range=[xmin, xmax],
                    row=track.subplot_y + 1,
                    col=track.subplot_x + 1,
                )
                subplots.update_yaxes(
                    row=track.subplot_y + 1,
                    col=track.subplot_x + 1,
                    **track.track_options.get("yaxis", {}),
                )

            except Exception as e:
                raise RuntimeError(
                    f"Error plotting track at grid ({track.subplot_y},{track.subplot_x}): {e}"
                )

        return subplots
