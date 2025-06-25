from dataclasses import dataclass, field
from typing import Literal, Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import scipy.sparse as sp
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle

from jdgenometracks.option_mapping import translate

from .GenomeTrack import GenomeTrack


@dataclass
class BedTrack(GenomeTrack):
    """
    A class for plotting genomic regions from BED files using matplotlib or Plotly.
    """

    track_options: dict = field(default_factory=dict)

    def __post_init__(self):
        super().__post_init__()
        if self.data is None:
            self.read_data()

    def read_data(self):
        """
        Reads BED data from the file and assigns it to the `data` attribute. If the file is empty,
        initializes an empty DataFrame with predefined columns.
        """
        if self.file_path is None:
            raise ValueError("File path must be provided.")

        try:
            data = pd.read_csv(self.file_path, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            data = pd.DataFrame(
                columns=["chrom", "chromStart", "chromEnd", "name", "score"]
            )

        # Define possible BED columns and adjust based on data size
        possible_bed_columns = [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRGB",
        ]
        data.columns = possible_bed_columns[: len(data.columns)]
        self.data = self.set_df_col_dtype(data)

    def get_subset_region_start(self, subset_region: Optional[str], **kwargs) -> int:
        """
        Extracts the start position from the subset region, or returns the minimum x value.

        Args:
            subset_region (Optional[str]): Region in "chrom:start-end" format.

        Returns:
            int: The starting position of the region.
        """

    def update_coverage(
        self, region, subset_region_start, bed_region_coverage, y
    ) -> np.ndarray:
        """
        Updates the coverage matrix for a region by extending coverage in that range.

        Args:
            region: The region of interest from the BED file.
            subset_region_start: The start of the subset region.
            bed_region_coverage: The existing coverage matrix.
            y: The y-value of the region for the coverage.

        Returns:
            np.ndarray: The updated coverage matrix.
        """
        new_coverage = np.zeros(bed_region_coverage.shape)
        new_coverage[
            0,
            max(0, region["chromStart"] - subset_region_start)
            - 1 : min(
                region["chromEnd"] - subset_region_start, bed_region_coverage.shape[1]
            )
            + 1,
        ] = y + self.track_options.get("rect.height", 1)
        return bed_region_coverage.maximum(sp.csr_matrix(new_coverage))

    def plot_mpl(
        self,
        ax: Axes,
        bed_region_coverage: np.ndarray,
        subset_region: Optional[str] = None,
        axis_shift: Optional[int] = None,
        xmin: Optional[int] = None,
    ):
        """
        Plots the BED data using matplotlib. Draws rectangles for genomic regions.

        Args:
            ax (Axes): The matplotlib axis to plot on.
            bed_region_coverage (np.ndarray): The coverage data matrix.
            subset_region (Optional[str]): The region to subset for plotting.
            axis_shift (Optional[int]): The genomic axis shift for visualization.

        Returns:
            np.ndarray: The updated coverage matrix after plotting.
        """
        if subset_region is not None:
            subset_region_start = int(subset_region.split(":")[1].split("-")[0])
        else:
            subset_region_start = xmin or 0

        cleaned_data = self.format_data(
            subset_region=subset_region, axis_shift=axis_shift
        )

        # Use unified options
        mpl_opts = translate(self.track_options, target="mpl")
        rect_opts = mpl_opts.get("marker", {})
        text_opts = mpl_opts.get("text", {})
        legend_opts = mpl_opts.get("legend", {})
        for idx, region in cleaned_data.iterrows():
            if "itemRGB" in region and self.track_options.get(
                "use_color_column", False
            ):
                rect_opts["color"] = [
                    int(val) / 255 for val in region["itemRGB"].split(",")
                ]
            try:
                y = bed_region_coverage[
                    0,
                    max(0, region["chromStart"] - subset_region_start)
                    + 1 : min(
                        region["chromEnd"] - subset_region_start,
                        bed_region_coverage.shape[1],
                    )
                    - 1,
                ].max()
            except ValueError:
                y = 0
            rect_height = self.track_options.get("rect.height", 1)
            rect_padding = self.track_options.get("rect.padding", 0)
            region_rect = Rectangle(
                (region["chromStart"], y + rect_padding),
                region["chromEnd"] - region["chromStart"],
                rect_height - 2 * rect_padding,
                label=region["name"],
                **rect_opts,
            )
            ax.add_patch(region_rect)
            label_alignment = self.track_options.get("label.alignment", False)
            if label_alignment == "above":
                ax.text(
                    (region["chromEnd"] + region["chromStart"]) / 2,
                    y + rect_height + rect_padding,
                    region["name"],
                    **text_opts,
                )
            elif label_alignment == "left":
                ax.text(
                    region["chromStart"],
                    y + rect_height / 2,
                    region["name"],
                    **text_opts,
                )
            elif label_alignment == "right":
                ax.text(
                    region["chromEnd"],
                    y + rect_height / 2,
                    region["name"],
                    **text_opts,
                )
            bed_region_coverage = self.update_coverage(
                region, subset_region_start, bed_region_coverage, y
            )
        ax.xaxis.set_tick_params(bottom=False)
        ax.yaxis.set_tick_params(left=False, labelleft=False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        if self.track_options.get("use_global_max", False):
            ax.set_ylim(0, bed_region_coverage.max() * 1.1)
        else:
            ax.autoscale(enable=True, axis="y")
        if self.show_legend:
            ax.legend(**legend_opts)
        return bed_region_coverage

    def plot_plotly(
        self,
        fig: go.Figure,
        row: int,
        col: int,
        bed_region_coverage: np.ndarray,
        subset_region: Optional[str] = None,
        axis_shift: Optional[int] = None,
        xmin: Optional[int] = None,
    ):
        """
        Plots the BED data using Plotly. Draws filled polygons for genomic regions.

        Args:
            fig (go.Figure): The Plotly figure to add traces to.
            row (int): The row index in the subplot grid.
            col (int): The column index in the subplot grid.
            bed_region_coverage (np.ndarray): The coverage data matrix.
            subset_region (Optional[str]): The region to subset for plotting.
            axis_shift (Optional[int]): The genomic axis shift for visualization.

        Returns:
            np.ndarray: The updated coverage matrix after plotting.
        """
        if subset_region is not None:
            subset_region_start = int(subset_region.split(":")[1].split("-")[0])
        else:
            subset_region_start = xmin or 0

        cleaned_data = self.format_data(
            subset_region=subset_region, axis_shift=axis_shift
        )

        plotly_opts = translate(self.track_options, target="plotly")
        marker_opts = plotly_opts.get("marker", {})
        fill_opts = plotly_opts.get("fill", {})
        text_opts = plotly_opts.get("text", {})
        legend_opts = plotly_opts.get("legend", {})
        line_opts = plotly_opts.get("line", {})
        xaxis_opts = plotly_opts.get("xaxis", {})
        yaxis_opts = plotly_opts.get("yaxis", {})
        for idx, region in cleaned_data.iterrows():
            if "itemRGB" in region and self.track_options.get(
                "use_color_column", False
            ):
                fill_opts["fillcolor"] = f"rgb({region['itemRGB']})"
            try:
                y = bed_region_coverage[
                    0,
                    max(0, region["chromStart"] - subset_region_start)
                    - 1 : min(
                        region["chromEnd"] - subset_region_start,
                        bed_region_coverage.shape[1],
                    )
                    + 1,
                ].max()
            except ValueError:
                y = 0
            rect_height = self.track_options.get("rect.height", 1)
            rect_padding = self.track_options.get("rect.padding", 0)
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
                        y + rect_padding,
                        y + rect_height - 2 * rect_padding,
                        y + rect_height - 2 * rect_padding,
                        y + rect_padding,
                        y + rect_padding,
                    ],
                    mode="lines",
                    fill="toself",
                    showlegend=self.show_legend,
                    name=region["name"],
                    **marker_opts,
                    **legend_opts,
                    **fill_opts,
                    line=line_opts,
                ),
                row=row,
                col=col,
            )
            label_alignment = self.track_options.get("label_alignment", False)
            if label_alignment == "above":
                fig.add_annotation(
                    x=(region["chromEnd"] + region["chromStart"]) / 2,
                    y=y + rect_height - rect_padding,
                    text=region["name"],
                    xanchor="center",
                    yanchor="bottom",
                    showarrow=False,
                    row=row,
                    col=col,
                    **text_opts,
                )
            elif label_alignment == "left":
                fig.add_annotation(
                    x=region["chromStart"],
                    y=y + rect_height / 2,
                    text=region["name"],
                    xanchor="right",
                    yanchor="middle",
                    showarrow=False,
                    row=row,
                    col=col,
                    **text_opts,
                )
            elif label_alignment == "right":
                fig.add_annotation(
                    x=region["chromEnd"],
                    y=y + rect_height / 2,
                    text=region["name"],
                    xanchor="left",
                    yanchor="middle",
                    showarrow=False,
                    row=row,
                    col=col,
                    **text_opts,
                )
            bed_region_coverage = self.update_coverage(
                region, subset_region_start, bed_region_coverage, y
            )

        fig.update_xaxes(row=row, col=col, **xaxis_opts)
        if self.track_options.get("use_global_max", False):
            fig.update_yaxes(
                range=[0, bed_region_coverage.max() * 1.1], row=row, col=col
            )
        fig.update_yaxes(
            row=row,
            col=col,
            showticklabels=False,
            **yaxis_opts,
        )
        return bed_region_coverage
