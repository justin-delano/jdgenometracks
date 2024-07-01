from dataclasses import dataclass

import plotly.graph_objects as go
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle
import scipy as sp
import numpy as np
from .GenomeTrack import GenomeTrack


@dataclass
class BedTrack(GenomeTrack):
    rect_border_width: int = 2
    rect_height: float = 1
    rect_padding: float = 0.2
    use_global_max: bool = False
    use_offset: bool = True

    def plot_mpl(self, ax: Axes, **kwargs):

        bed_region_coverage = kwargs.get("bed_region_coverage", 0)

        subset_region = kwargs.get("subset_region", None)
        subset_region_start  = int(subset_region.split(":")[1].split("-")[0])

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

            try:
                y = bed_region_coverage[0, max(0, region["chromStart"]-subset_region_start):min(region["chromEnd"]-subset_region_start, bed_region_coverage.shape[1])].max()+1
            except ValueError as e:
                y=0

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

            new_coverage = np.zeros(bed_region_coverage.shape)
            new_coverage[0, max(0, region["chromStart"]-subset_region_start):min(region["chromEnd"]-subset_region_start, bed_region_coverage.shape[1])] = y
            bed_region_coverage = bed_region_coverage.maximum(sp.sparse.csr_matrix(new_coverage))


        ax.xaxis.set_tick_params(bottom=False)
        ax.yaxis.set_tick_params(left=False, labelleft=False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        if self.use_global_max:
            ax.set_ylim(0, max(kwargs.get("max_regions", 1), 1) * self.rect_height)

        return bed_region_coverage

    def plot_plotly(
        self,
        fig: go.Figure,
        row: int,
        col: int,
        **kwargs,
    ):
        bed_region_coverage = kwargs.get("bed_region_coverage", 0)
        subset_region = kwargs.get("subset_region", None)
        subset_region_start  = int(subset_region.split(":")[1].split("-")[0])
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

            try:
                y = bed_region_coverage[0, max(0, region["chromStart"]-subset_region_start):min(region["chromEnd"]-subset_region_start, bed_region_coverage.shape[1])].max()+1
            except ValueError as e:
                y=0

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
            new_coverage = np.zeros(bed_region_coverage.shape)
            new_coverage[0, max(0, region["chromStart"]-subset_region_start):min(region["chromEnd"]-subset_region_start, bed_region_coverage.shape[1])] = y
            bed_region_coverage = bed_region_coverage.maximum(sp.sparse.csr_matrix(new_coverage))


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

        return bed_region_coverage