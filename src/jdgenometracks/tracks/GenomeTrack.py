from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
import plotly.graph_objects as go
from matplotlib.axes import Axes


@dataclass
class GenomeTrack:
    # Required
    file_path: str
    track_type: str
    track_name: str
    data: pd.DataFrame
    # Style
    main_color: str = "blue"
    secondary_color: str = "red"
    # Text
    axis_font_size: int = 12
    text_font_size: int = 12
    # Horizontal Lines
    hlines: list | None = None
    hline_color: str = "black"
    hline_width: int = 1
    # Structure
    height_prop: float | None = None
    share_with_previous: bool = False
    x_axis_type: str | None = None
    x_axis_interval: int | None = None
    axis_ticks: bool = False

    def add_hlines_mpl(self, ax: Axes, **kwargs):
        if self.hlines is None:
            pass
        else:
            for hline in self.hlines:
                ax.axhline(hline, color=self.hline_color, lw=self.hline_width)

    def add_hlines_plotly(self, fig: go.Figure, row: int, col: int, **kwargs):
        if self.hlines is None:
            pass
        else:
            for hline in self.hlines:
                fig.add_hline(
                    y=hline,
                    line=dict(color=self.hline_color, width=self.hline_width),
                    row=row,  # type: ignore
                    col=col,  # type: ignore
                )

    def format_data(self, subset_region: str | None = None, axis_shift: int = 0):
        if subset_region is not None:
            chrom, start, end = (
                subset_region.split(":")[0],
                int(subset_region.split(":")[1].split("-")[0]),
                int(subset_region.split(":")[1].split("-")[1]),
            )
            formatted_data = self.data[
                (self.data["chrom"] == chrom)
                & (
                    (self.data["chromStart"].between(start, end))
                    | (self.data["chromEnd"].between(start, end))
                )
            ].reset_index(drop=True)
        else:
            formatted_data = self.data

        if axis_shift:
            formatted_data["chromStart"] -= axis_shift
            formatted_data["chromEnd"] -= axis_shift

        return formatted_data

    def plot_mpl(self, ax: Axes, **kwargs):
        pass

    def plot_plotly(self, fig: go.Figure, row: int, col: int, **kwargs):
        pass
