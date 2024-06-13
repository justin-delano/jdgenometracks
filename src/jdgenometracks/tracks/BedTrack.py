from dataclasses import dataclass

import plotly.graph_objects as go
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle

from .GenomeTrack import GenomeTrack


@dataclass
class BedTrack(GenomeTrack):
    rect_border_width: int = 2
    rect_height: float = 1
    rect_padding: float = 0.2
    use_global_max: bool = False
    use_offset: bool = True

    def plot_mpl(self, ax: Axes, **kwargs):

        offset = 0 if not self.use_offset else kwargs["offset"]
        y = offset

        subset_region = kwargs.get("subset_region", None)
        axis_shift = kwargs.get("axis_shift", 0)
        if subset_region is not None:
            xdata = self.format_data(subset_region=subset_region, axis_shift=axis_shift)
        else:
            xdata = self.format_data(axis_shift=axis_shift)

        for idx, region in xdata.iterrows():
            if "itemRGB" in region:
                color = [int(val) / 255 for val in region["itemRGB"].split(",")]
            else:
                color = self.main_color

            region_rect = Rectangle(
                (region["chromStart"], y + self.rect_padding),
                region["chromEnd"] - region["chromStart"],
                self.rect_height - 2 * self.rect_padding,
                color=color,
                lw=self.rect_border_width,
            )
            ax.add_patch(region_rect)

            ax.text(
                (region["chromEnd"] + region["chromStart"]) / 2,
                y + self.rect_height - self.rect_padding,
                region["name"],
                fontsize=self.text_font_size,
                color="black",
                va="bottom",
                ha="center",
            )

            if (
                idx < xdata.shape[0] - 1  # type: ignore
                and region["chromEnd"] >= xdata.iloc[idx + 1, :]["chromStart"]  # type: ignore
            ):
                y += 1
            else:
                y = max(offset, y - 1)

        ax.xaxis.set_tick_params(bottom=False)
        ax.yaxis.set_tick_params(left=False, labelleft=False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        if self.use_global_max:
            ax.set_ylim(0, max(kwargs.get("max_regions", 1), 1) * self.rect_height)

    def plot_plotly(
        self,
        fig: go.Figure,
        row: int,
        col: int,
        **kwargs,
    ):
        offset = 0 if not self.use_offset else kwargs["offset"]
        y = offset

        subset_region = kwargs.get("subset_region", None)
        axis_shift = kwargs.get("axis_shift", 0)
        if subset_region is not None:
            xdata = self.format_data(subset_region=subset_region, axis_shift=axis_shift)
        else:
            xdata = self.format_data(axis_shift=axis_shift)

        for idx, region in xdata.iterrows():
            if "itemRGB" in region:
                color = f"rgb({region['itemRGB']})"
            else:
                color = self.main_color

            fig.add_trace(
                go.Scatter(
                    x=[
                        region["chromStart"],
                        region["chromStart"],
                        region["chromEnd"],
                        region["chromEnd"],
                        region["chromStart"],
                    ],
                    y=[
                        y + self.rect_padding,
                        y + self.rect_height - 2 * self.rect_padding,
                        y + self.rect_height - 2 * self.rect_padding,
                        y + self.rect_padding,
                        y + self.rect_padding,
                    ],
                    mode="lines",
                    fill="toself",
                    fillcolor=color,
                    line=dict(color=color),
                    showlegend=False,
                    name=region["name"],
                ),
                row=row,
                col=col,
            )

            if (
                idx < xdata.shape[0] - 1  # type: ignore
                and region["chromEnd"] >= xdata.iloc[idx + 1, :]["chromStart"]  # type: ignore
            ):
                y += 1
            else:
                y = max(offset, y - 1)

        fig.update_xaxes(
            showline=False,
            linewidth=self.rect_border_width,
            linecolor="black",
            row=row,
            col=col,
        )

        fig.update_yaxes(showticklabels=False, row=row, col=col)
        if self.use_global_max:
            fig.update_yaxes(
                range=[0, kwargs.get("max_regions", 1) * 1.1 * self.rect_height],
                row=row,
                col=col,
            )
