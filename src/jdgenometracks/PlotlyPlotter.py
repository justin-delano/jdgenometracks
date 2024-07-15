from __future__ import annotations
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import plotly.graph_objects as go
import plotly.subplots as ps
import scipy as sp
from .tracks.BedTrack import BedTrack
from .tracks.GenomeTrack import GenomeTrack
from .utils import _get_col_limits, _get_height_props


@dataclass
class PlotlyPlotter:
    tracks: npt.ArrayLike[GenomeTrack]
    total_height: float

    def __post_init__(self):
        self.tracks = np.array(self.tracks) if not isinstance(self.tracks, np.ndarray) else self.tracks

    def plot_single_track(
        self, subplots: go.Figure, track: GenomeTrack, row: int, col: int, **kwargs
    ) -> None:
        if isinstance(track, BedTrack):
            bed_region_coverage = track.plot_plotly(subplots, row, col, **kwargs)
        else:
            track.plot_plotly(subplots, row, col, **kwargs)
            bed_region_coverage = None

        track.add_hlines_plotly(subplots, row, col, **kwargs)
        return bed_region_coverage

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
    ) -> go.Figure:
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
            go.Figure: Complete plotly figure
        """
        if self.tracks.ndim == 1:
            self.tracks = self.tracks.reshape(-1, 1)

        if column_regions is None:
            column_regions = [None for _ in range(self.tracks.shape[1])]

        num_distinct_rows = sum(
            [not track.share_with_previous for track in self.tracks[:, 0]]
        )

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

        subplots = ps.make_subplots(
            rows=num_distinct_rows,
            cols=self.tracks.shape[1],
            shared_xaxes="columns",  # type: ignore
            row_heights=height_props,
            row_titles=row_titles,
            column_widths=width_props,
            column_titles=column_titles,
            vertical_spacing=0.02,
            horizontal_spacing=0.05,
        )
        if plot_title is not None:
            subplots.update_layout(
                title_text=plot_title, title_font_size=16, title_x=0.5
            )
        subplots.update_layout(
            autosize=True,
            height=48 * total_height,
            plot_bgcolor="white",
            margin=dict(l=0.1, r=0.1, t=50, b=20, pad=4),
        )

        for data_col in range(self.tracks.shape[1]):
            plot_row = 1
            plot_col = data_col + 1

            xmin, xmax, extra_options = _get_col_limits(
                self.tracks[:, data_col], column_regions[data_col], relative_x_axis
            )

            bed_region_coverage = sp.sparse.csr_matrix((1, xmax - xmin))
            for data_row in range(self.tracks.shape[0]):
                track = self.tracks[data_row, data_col]

                if track is None:
                    continue

                if track.share_with_previous:
                    plot_row -= 1
                else:
                    bed_region_coverage = sp.sparse.csr_matrix((1, xmax - xmin))

                if isinstance(track, BedTrack):
                    extra_options["bed_region_coverage"] = bed_region_coverage
                    bed_region_coverage = self.plot_single_track(
                    subplots,
                    track,
                    row=plot_row,
                    col=plot_col,
                    region=column_regions[data_col],
                    **extra_options,
                    )
                else:
                    self.plot_single_track(
                    subplots,
                    track,
                    row=plot_row,
                    col=plot_col,
                    region=column_regions[data_col],
                    **extra_options,
                )

                subplots.update_xaxes(
                    range=[xmin, xmax],
                    tickformatstops=[
                        dict(dtickrange=[None, 1e3], value=","),
                        dict(dtickrange=[1e3, 7e5], value=".4s"),
                        dict(dtickrange=[7e5, None], value=".4s"),
                    ],
                    ticksuffix="b",
                    row=plot_row,
                    col=plot_col,
                )
                plot_row += 1

        if show_fig:
            subplots.show()
        return subplots
