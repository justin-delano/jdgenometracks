from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Sequence

import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp

from jdgenometracks.tracks import XAxisTrack

from .tracks.BedTrack import BedTrack
from .utils import TrackUtils
from .utils_units import convert_to_inches, parse_size_with_units


@dataclass
class MPLPlotter:
    """
    A class for generating multi-track genomic plots using Matplotlib.

    Attributes:
        tracks (Sequence[Sequence[Sequence[Any]]]): 2D grid of subplots, each containing a sequence of tracks.
        total_height (float): Total height of the plot figure in inches.
        total_width (float): Total width of the plot figure in inches.
    """

    tracks: Sequence[Any]

    def __post_init__(self) -> None:
        """
        Ensures tracks are stored as a grid and validates input types.
        """
        if not isinstance(self.tracks, (list, tuple)):
            raise TypeError(f"Tracks must be a list or tuple, got {type(self.tracks)}.")
        if len(self.tracks) == 0:
            raise ValueError("Tracks cannot be empty.")
        for track in self.tracks:
            if track is not None and not hasattr(track, "plot_mpl"):
                raise TypeError(
                    f"Each track must have a 'plot_mpl' method (got {type(track)})."
                )

    def _validate_inputs(
        self,
        max_rows: int,
        max_cols: int,
        height_props: Sequence[float],
        row_titles: Sequence[str],
        width_props: Sequence[float],
        column_titles: Sequence[str],
    ) -> None:
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

    def _set_subplot_titles(
        self,
        subplots: Any,  # np.ndarray, but numpy not imported at top level
        column_titles: Sequence[str],
        row_titles: Sequence[str],
    ) -> None:

        for col_idx, title in enumerate(column_titles):
            subplots[0, col_idx].set_title(title)
        for row_idx, title in enumerate(row_titles):
            subplots[row_idx, 0].set_ylabel(title)

    def plot_all_tracks(
        self,
        fig_options: Optional[Dict[str, Any]] = None,
    ) -> tuple[
        matplotlib.figure.Figure, Any
    ]:  # np.ndarray, but numpy not imported at top level
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

        try:
            fig, axes = plt.subplots(
                max_rows,
                max_cols,
                figsize=(total_width, total_height),
                sharex="col",
                gridspec_kw={
                    "height_ratios": height_props,
                    "width_ratios": width_props,
                },
                layout="constrained",
            )
        except Exception as e:
            raise RuntimeError(f"Error creating subplots: {e}")

        axes = np.atleast_2d(axes).T

        self._set_subplot_titles(axes, column_titles, row_titles)
        if plot_title:
            fig.suptitle(
                str(plot_title), fontsize=fig_options.get("suptitle.fontsize", 16)
            )
        for track in self.tracks:
            if track is None:
                continue
            if not hasattr(track, "plot_mpl"):
                raise TypeError(f"Track {track} does not have a 'plot_mpl' method.")
            try:
                if isinstance(track, BedTrack):
                    bed_region_coverages[track.subplot_x] = track.plot_mpl(
                        axes[track.subplot_y, track.subplot_x],
                        bed_region_coverages[track.subplot_x],
                        subset_region=column_regions[track.subplot_x],
                        xmin=column_limits[track.subplot_x]["xmin"],
                    )
                elif isinstance(track, XAxisTrack):
                    track.plot_mpl(
                        axes[track.subplot_y, track.subplot_x],
                        chromosome=column_limits[track.subplot_x]["chromosome"],
                    )
                else:
                    track.plot_mpl(
                        axes[track.subplot_y, track.subplot_x],
                        subset_region=column_regions[track.subplot_x],
                    )
                axes[track.subplot_y, track.subplot_x].set_xlim(
                    column_limits[track.subplot_x]["xmin"],
                    column_limits[track.subplot_x]["xmax"],
                )
                track.add_hlines_mpl(axes[track.subplot_y, track.subplot_x])
                axes[track.subplot_y, track.subplot_x].spines["top"].set_visible(False)
                axes[track.subplot_y, track.subplot_x].spines["right"].set_visible(
                    False
                )
            except Exception as e:
                raise RuntimeError(
                    f"Error plotting track at grid ({track.subplot_y},{track.subplot_x}): {e}"
                )
        return fig, axes
