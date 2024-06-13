from __future__ import annotations

from dataclasses import dataclass

import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from matplotlib.axes import Axes

from .tracks.BedTrack import BedTrack
from .tracks.GenomeTrack import GenomeTrack
from .utils import _get_col_limits, _get_height_props


@dataclass
class MPLPlotter:
    tracks: npt.NDArray[GenomeTrack]  # type: ignore
    total_height: float

    def plot_single_track(self, subplot: Axes, track: GenomeTrack, **kwargs) -> None:
        track.plot_mpl(subplot, **kwargs)
        track.add_hlines_mpl(subplot, **kwargs)
        subplot.spines["top"].set_visible(False)
        subplot.spines["right"].set_visible(False)

    def plot_all_tracks(
        self,
        plot_title: str | None = None,
        column_regions: list[str | None] | None = None,
        total_height: float | None = None,
        height_props: list[float] | None = None,
        row_titles: list[str] | None = None,
        width_props: list[float] | None = None,
        column_titles: list[str] | None = None,
        relative_x_axis: bool = False,
        show_fig: bool = False,
    ) -> None | tuple[matplotlib.figure.Figure, np.ndarray]:
        """Plot all saved tracks in a single figure

        Args:
            plot_title (str | None, optional): Adds a title to the plot. Defaults to None.
            column_regions (list[str  |  None] | None, optional): Allows subsetting of each column to a genomic region in format "chr:start-end". Defaults to None.
            total_height (float | None, optional): Total figure height. Defaults to None.
            height_props (list[float] | None, optional): Proportion of figure height used by each track. Should have length equal to number of track rows being plotted. Defaults to None.
            row_titles (list[str] | None, optional): Titles for each track row. Should have length equal to number of track rows being plotted. Defaults to None.
            width_props (list[float] | None, optional): Proportion of figure width used by each column. Should have length equal to number of track columns being plotted. Defaults to None.
            column_titles (list[str] | None, optional): Titles for each track column. Should have length equal to number of track columns being plotted. Defaults to None.
            relative_x_axis (bool, optional): If true, x-axes start at 0bp. Defaults to False.
            show_fig (bool, optional): If true, figure is shown as well as returned. Defaults to False.

        Returns:
            tuple[matplotlib.figure.Figure, np.ndarray]: Matplotlib (fig, ax) tuple.
        """

        if column_regions is None:
            column_regions = [None for _ in range(self.tracks.shape[1])]

        num_distinct_rows = sum(
            [not track.share_with_previous for track in self.tracks[:, 0]]
        )

        if self.tracks.ndim == 1:
            self.tracks = self.tracks.reshape(-1, 1)

        height_props = height_props or _get_height_props(self.tracks)
        row_titles = row_titles or [""] * num_distinct_rows
        width_props = width_props or [1 / self.tracks.shape[1]] * self.tracks.shape[1]
        column_titles = column_titles or [""] * self.tracks.shape[1]
        assert (
            len(height_props) == num_distinct_rows
        ), f"Number of height_props should be equal to {num_distinct_rows}"
        assert (
            len(row_titles) == num_distinct_rows
        ), f"Number of row_titles should be equal to {num_distinct_rows}"
        assert (
            len(width_props) == self.tracks.shape[1]
        ), f"Number of width_props should be equal to {self.tracks.shape[1]}"
        assert (
            len(column_titles) == self.tracks.shape[1]
        ), f"Number of column_titles should be equal to {self.tracks.shape[1]}"
        if total_height is None:
            total_height = self.total_height

        fig, subplots = plt.subplots(
            num_distinct_rows,
            self.tracks.shape[1],
            figsize=(20, total_height),
            sharex="col",
            gridspec_kw={
                "height_ratios": height_props,
                "width_ratios": width_props,
            },
            layout="constrained",
        )

        if subplots.ndim == 1:
            subplots = subplots.reshape(-1, 1)
        if column_titles or row_titles:
            for data_col, title in enumerate(column_titles or []):
                subplots[0, data_col].set_title(title)
            for data_row, title in enumerate(row_titles or []):
                subplots[data_row, 0].set_ylabel(title)

        if plot_title is not None:
            fig.suptitle(plot_title, fontsize=16)

        for data_col in range(self.tracks.shape[1]):
            plot_row = 0
            plot_col = data_col

            xmin, xmax, extra_options = _get_col_limits(
                self.tracks[:, data_col], column_regions[data_col], relative_x_axis
            )

            for data_row in range(self.tracks.shape[0]):
                track = self.tracks[data_row, data_col]

                if track is None:
                    continue

                if track.share_with_previous:
                    plot_row -= 1

                if isinstance(track, BedTrack):
                    extra_options["offset"] = len(
                        subplots[plot_row, plot_col].patches  # type: ignore
                    )

                self.plot_single_track(
                    subplots[plot_row, plot_col],  # type: ignore
                    track,
                    region=column_regions[data_col],
                    **extra_options,
                )
                subplots[plot_row, plot_col].set_xlim(xmin, xmax)
                plot_row += 1

        if show_fig:
            fig.show()
        return fig, subplots
