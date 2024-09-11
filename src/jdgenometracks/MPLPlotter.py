from __future__ import annotations

from dataclasses import dataclass

import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
from matplotlib.axes import Axes

from .tracks.BedTrack import BedTrack
from .tracks.GenomeTrack import GenomeTrack
from .utils import TrackUtils


@dataclass
class MPLPlotter:
    """
    A class for generating multi-track genomic plots using Matplotlib.

    Attributes:
        tracks (np.ndarray): Array of GenomeTrack or its subclasses to be plotted.
        total_height (float): Total height of the plot figure in inches.
    """

    tracks: np.ndarray
    total_height: float
    total_width: float

    def __post_init__(self):
        """
        Ensures tracks are stored as a numpy array for consistency.
        """
        if not isinstance(self.tracks, np.ndarray):
            self.tracks = np.array(self.tracks)

    def plot_single_track(
        self, subplot: Axes, track: GenomeTrack, **kwargs
    ) -> np.ndarray | None:
        """
        Plots a single track on a subplot.

        Args:
            subplot (Axes): The matplotlib subplot to plot on.
            track (GenomeTrack): The genomic track to plot.
            **kwargs: Additional keyword arguments to pass to the plot method.

        Returns:
            sp.spmatrix | None: Bed region coverage if applicable, otherwise None.
        """
        if isinstance(track, BedTrack):
            bed_region_coverage = track.plot_mpl(subplot, **kwargs)
        else:
            track.plot_mpl(subplot, **kwargs)
            bed_region_coverage = None

        track.add_hlines_mpl(subplot)
        subplot.spines["top"].set_visible(False)
        subplot.spines["right"].set_visible(False)
        return bed_region_coverage

    def _validate_inputs(
        self,
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
        num_distinct_rows = sum(
            [not track.share_with_previous for track in self.tracks[:, 0]]
        )

        assert (
            len(height_props) == num_distinct_rows
        ), f"Number of height_props should equal {num_distinct_rows}"
        assert (
            len(row_titles) == num_distinct_rows
        ), f"Number of row_titles should equal {num_distinct_rows}"
        assert (
            len(width_props) == self.tracks.shape[1]
        ), f"Number of width_props should equal {self.tracks.shape[1]}"
        assert (
            len(column_titles) == self.tracks.shape[1]
        ), f"Number of column_titles should equal {self.tracks.shape[1]}"

    def _set_subplot_titles(
        self, subplots: np.ndarray, column_titles: list[str], row_titles: list[str]
    ) -> None:
        """
        Sets the titles for the rows and columns of the subplots.

        Args:
            subplots (np.ndarray): The array of subplots.
            column_titles (list[str]): The titles for the columns.
            row_titles (list[str]): The titles for the rows.
        """
        for col_idx, title in enumerate(column_titles or []):
            subplots[0, col_idx].set_title(title)

        for row_idx, title in enumerate(row_titles or []):
            subplots[row_idx, 0].set_ylabel(title)

    def plot_all_tracks(
        self,
        plot_title: str | None = None,
        column_regions: list[str | None] | None = None,
        total_height: float | None = None,
        total_width: float | None = None,
        height_props: list[float] | None = None,
        row_titles: list[str] | None = None,
        width_props: list[float] | None = None,
        column_titles: list[str] | None = None,
        relative_x_axis: bool = False,
        show_fig: bool = False,
    ) -> tuple[matplotlib.figure.Figure, np.ndarray]:
        """
        Plots all tracks in a single matplotlib figure.

        Args:
            plot_title (str | None): Title for the entire plot.
            column_regions (list[str | None] | None): Genomic regions for subsetting in each column.
            total_height (float | None): Total figure height in inches.
            height_props (list[float] | None): Proportions of figure height used by each track.
            row_titles (list[str] | None): Titles for each track row.
            width_props (list[float] | None): Proportions of figure width used by each column.
            column_titles (list[str] | None): Titles for each track column.
            relative_x_axis (bool): If True, x-axes start at 0bp.
            show_fig (bool): If True, display the figure.

        Returns:
            tuple[matplotlib.figure.Figure, np.ndarray]: The matplotlib figure and the subplots array.
        """
        # Reshape tracks if they are a 1D array
        if self.tracks.ndim == 1:
            self.tracks = self.tracks.reshape(-1, 1)

        # Default column regions to None
        if column_regions is None:
            column_regions = [None for _ in range(self.tracks.shape[1])]

        # Calculate number of distinct rows (tracks not sharing axes with the previous row)
        num_distinct_rows = sum(
            [not track.share_with_previous for track in self.tracks[:, 0]]
        )

        # Default layout properties
        height_props = height_props or TrackUtils.get_height_props(self.tracks)
        row_titles = row_titles or [""] * num_distinct_rows
        width_props = width_props or [1 / self.tracks.shape[1]] * self.tracks.shape[1]
        column_titles = column_titles or [""] * self.tracks.shape[1]
        total_height = total_height or self.total_height
        total_width = total_width or self.total_width

        # Validate inputs
        self._validate_inputs(height_props, row_titles, width_props, column_titles)

        # Create figure and subplots
        fig, axes = plt.subplots(
            num_distinct_rows,
            self.tracks.shape[1],
            figsize=(total_width, total_height),
            sharex="col",
            gridspec_kw={
                "height_ratios": height_props,
                "width_ratios": width_props,
            },
            layout="constrained",
        )

        # Ensure subplots is always a 2D array
        axes = np.atleast_2d(axes).T

        # Set row and column titles
        self._set_subplot_titles(axes, column_titles, row_titles)

        # Set plot title if provided
        if plot_title:
            fig.suptitle(plot_title, fontsize=16)

        # Plot each track in the subplots
        for data_col in range(self.tracks.shape[1]):
            plot_row = 0

            xmin, xmax, extra_options = TrackUtils.get_col_limits(
                self.tracks[:, data_col], column_regions[data_col], relative_x_axis
            )

            # Initialize bed region coverage for BedTracks
            bed_region_coverage = sp.csr_matrix((1, xmax - xmin))

            for data_row in range(self.tracks.shape[0]):
                track = self.tracks[data_row, data_col]

                if track is None:
                    continue

                # If the track shares its x-axis with the previous track, adjust row index
                if track.share_with_previous:
                    plot_row -= 1
                else:
                    bed_region_coverage = sp.csr_matrix((1, xmax - xmin))

                # Plot BedTrack with bed region coverage
                if isinstance(track, BedTrack):
                    extra_options["bed_region_coverage"] = bed_region_coverage
                    bed_region_coverage = self.plot_single_track(
                        axes[plot_row, data_col],
                        track,
                        region=column_regions[data_col],
                        **extra_options,
                    )
                else:
                    self.plot_single_track(
                        axes[plot_row, data_col],
                        track,
                        region=column_regions[data_col],
                        **extra_options,
                    )
                axes[plot_row, data_col].set_xlim(xmin, xmax)
                plot_row += 1

        # Show the figure if required
        if show_fig:
            fig.show()

        return fig, axes
