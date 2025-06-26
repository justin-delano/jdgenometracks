from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional, Sequence

import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np

from jdgenometracks.tracks import XAxisTrack

from .base_plotter import BasePlotter
from .config import MPLDefaults, PlotDefaults
from .tracks.BedTrack import BedTrack
from .tracks.GenomeTrack import GenomeTrack


@dataclass
class MPLPlotter(BasePlotter):
    """
    A class for generating multi-track genomic plots using Matplotlib.

    Attributes:
        tracks: Sequence of tracks to plot, can contain instances of GenomeTrack or its subclasses.
    """

    def _validate_track_methods(self) -> None:
        """
        Validate that tracks have required Matplotlib plotting methods.
        """
        for track in self.tracks:
            if track is not None and not hasattr(track, "plot_mpl"):
                raise TypeError(
                    f"Each track must have a 'plot_mpl' method (got {type(track)})."
                )

    def _set_subplot_titles(
        self,
        subplots: Any,  # np.ndarray, but numpy not imported at top level
        column_titles: Sequence[str],
        row_titles: Sequence[str],
    ) -> None:
        """Set titles for matplotlib subplots."""
        for col_idx, title in enumerate(column_titles):
            subplots[0, col_idx].set_title(title)
        for row_idx, title in enumerate(row_titles):
            subplots[row_idx, 0].set_ylabel(title)

    def plot_all_tracks(
        self,
        fig_options: Optional[Dict[str, Any]] = None,
    ) -> tuple[matplotlib.figure.Figure, Any]:
        """
        Plots all tracks in a single matplotlib figure.

        Args:
            fig_options: Optional figure configuration options.

        Returns:
            Tuple of (matplotlib Figure, axes array).
        """
        # Use base class to prepare common data
        max_rows, max_cols, plot_data = self._prepare_plot_data(fig_options)

        # Create matplotlib figure
        try:
            fig, axes = plt.subplots(
                max_rows,
                max_cols,
                figsize=(plot_data["total_width"], plot_data["total_height"]),
                sharex=MPLDefaults.DEFAULT_SHAREX,
                gridspec_kw={
                    "height_ratios": plot_data["height_props"],
                    "width_ratios": plot_data["width_props"],
                },
                layout=MPLDefaults.DEFAULT_LAYOUT,
            )
        except Exception as e:
            raise RuntimeError(f"Error creating subplots: {e}")

        axes = np.atleast_2d(axes).T

        # Set subplot titles
        self._set_subplot_titles(
            axes, plot_data["column_titles"], plot_data["row_titles"]
        )

        # Set figure title if provided
        if plot_data["plot_title"]:
            fontsize = (
                fig_options.get(
                    "suptitle.fontsize", PlotDefaults.DEFAULT_SUPTITLE_FONTSIZE
                )
                if fig_options
                else PlotDefaults.DEFAULT_SUPTITLE_FONTSIZE
            )
            fig.suptitle(str(plot_data["plot_title"]), fontsize=fontsize)

        # Plot each track using the base class method
        for track in self.tracks:
            if track is None:
                continue

            backend_args = {
                "axes": axes,
                "ax": axes[track.subplot_y, track.subplot_x],
            }

            self._plot_track(track, plot_data, backend_args)

            # Hide matplotlib spines
            axes[track.subplot_y, track.subplot_x].spines["top"].set_visible(
                MPLDefaults.HIDE_TOP_SPINE
            )
            axes[track.subplot_y, track.subplot_x].spines["right"].set_visible(
                MPLDefaults.HIDE_RIGHT_SPINE
            )

        return fig, axes

    def _plot_bed_track(
        self,
        track: BedTrack,
        column_limits: list[dict[str, Any]],
        column_regions: list[str | None],
        track_y_levels: dict[int, float],
        backend_args: dict[str, Any],
    ) -> None:
        """Plot a BED track using Matplotlib."""
        track.plot_mpl(
            backend_args["ax"],
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
        """Plot an axis track using Matplotlib."""
        column_limit = column_limits[track.subplot_x]
        track.plot_mpl(
            backend_args["ax"],
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
        """Plot a generic track using Matplotlib."""
        track.plot_mpl(
            backend_args["ax"],
            subset_region=column_regions[track.subplot_x],
        )

    def _add_track_hlines(
        self,
        track: GenomeTrack,
        backend_args: dict[str, Any],
    ) -> None:
        """Add horizontal lines to a track using Matplotlib."""
        track.add_hlines_mpl(backend_args["ax"])

    def _update_track_axes(
        self,
        track: GenomeTrack,
        column_limits: list[dict[str, Any]],
        backend_args: dict[str, Any],
    ) -> None:
        """Update track axes properties using Matplotlib."""
        xmin = column_limits[track.subplot_x]["xmin"]
        xmax = column_limits[track.subplot_x]["xmax"]

        backend_args["ax"].set_xlim(xmin, xmax)
