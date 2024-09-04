from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import plotly.graph_objects as go
from matplotlib.axes import Axes

from .GenomeTrack import GenomeTrack


@dataclass
class BedGraphTrack(GenomeTrack):
    ymin: float | None = None
    ymax: float | None = None
    plot_type: str = "lines"
    secondary_color: str = "red"
    plot_line_width: int = 2
    plot_point_sizes: Iterable[float] | None = None
    scatter_marker_mpl: str = "."
    scatter_marker_plotly: str = "0"
    fill_color: str = "rgba(0,0,0,0)"
    alpha: float = 1

    def plot_mpl(self, ax: Axes, **kwargs):

        subset_region = kwargs.get("subset_region", None)
        axis_shift = kwargs.get("axis_shift", 0)
        if subset_region is not None:
            xdata = self.format_data(subset_region=subset_region, axis_shift=axis_shift)
        else:
            xdata = self.format_data(axis_shift=axis_shift)

        mid_points = (xdata["chromStart"] + xdata["chromEnd"]) / 2
        y_values = xdata["value"]

        if self.plot_type == "lines":
            ax.plot(
                mid_points,
                y_values,
                color=self.main_color,
                lw=self.plot_line_width,
            )
            ax.fill_between(
                mid_points,
                y_values,
                0,
                color=self.main_color,
                alpha=0.5,
            )
        elif self.plot_type == "bars":
            ax.bar(mid_points, y_values, color=self.main_color)
        elif self.plot_type == "points":
            ax.scatter(
                mid_points,
                y_values,
                color=self.main_color,
                s=self.plot_point_sizes,  # type: ignore
                marker=self.scatter_marker_mpl,
            )

        if self.ymin is not None:
            ax.set_ylim(bottom=self.ymin)
        if self.ymax is not None:
            ax.set_ylim(top=self.ymax)

        ax.spines["bottom"].set_visible(False)
        ax.xaxis.set_tick_params(bottom=False)

    def plot_plotly(
        self,
        fig: go.Figure,
        row: int,
        col: int,
        **kwargs,
    ):
        subset_region = kwargs.get("subset_region", None)
        axis_shift = kwargs.get("axis_shift", 0)
        if subset_region is not None:
            xdata = self.format_data(subset_region=subset_region, axis_shift=axis_shift)
        else:
            xdata = self.format_data(axis_shift=axis_shift)

        mid_points = (xdata["chromStart"] + xdata["chromEnd"]) / 2
        y_values = xdata["value"]

        if self.plot_type == "lines":
            fig.add_trace(
                go.Scatter(
                    x=mid_points,
                    y=y_values,
                    mode="lines",
                    line=dict(color=self.main_color, width=self.plot_line_width),
                    fill="tozeroy",
                    fillcolor=self.fill_color,
                    opacity=self.alpha,
                    showlegend=False,
                    name=self.track_name,
                ),
                row=row,
                col=col,
            )
        elif self.plot_type == "bars":
            fig.add_trace(
                go.Bar(
                    x=mid_points,
                    y=y_values,
                    marker=dict(color=self.main_color, line_color=self.main_color),
                    showlegend=False,
                    name=self.track_name,
                ),
                row=row,
                col=col,
            )
            fig.update_layout(bargap=0)

        elif self.plot_type == "points":
            fig.add_trace(
                go.Scatter(
                    x=mid_points,
                    y=y_values,
                    mode="markers",
                    marker=dict(
                        color=self.main_color,
                        line_color=self.main_color,
                        size=self.plot_point_sizes,
                        symbol=self.scatter_marker_plotly,
                        opacity=self.alpha,
                        line_width=self.plot_line_width,
                    ),
                    showlegend=False,
                    name=self.track_name,
                ),
                row=row,
                col=col,
            )

        fig.update_xaxes(showline=False, row=row, col=col)
        fig.update_yaxes(linecolor="black", row=row, col=col)
        if self.ymin is not None or self.ymax is not None:
            fig.update_yaxes(range=[self.ymin, self.ymax], row=row, col=col)
