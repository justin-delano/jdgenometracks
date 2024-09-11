from dataclasses import dataclass, field
from typing import Literal

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import scipy.sparse as sp
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle

from .GenomeTrack import GenomeTrack


@dataclass
class BedTrack(GenomeTrack):
    """
    A class for plotting genomic regions from BED files using matplotlib or Plotly.

    Attributes:
        rect_height (float): Height of each rectangle representing a region.
        rect_padding (float): Padding around the rectangles.
        use_global_max (bool): Whether to use global maximum coverage for the y-axis.
        use_color_column (bool): Whether to color rectangles based on a column in the BED file.
        mpl_rect_options (dict): Matplotlib options for rectangles.
        mpl_text_options (dict): Matplotlib options for text labels.
        plotly_options (dict): Plotly options for the plotting traces.
    """

    rect_height: float = 1
    rect_padding: float = 0
    use_global_max: bool = False
    use_color_column: bool = False

    mpl_rect_options: dict = field(default_factory=dict)
    mpl_text_options: dict = field(default_factory=dict)
    mpl_text_alignment: Literal["left", "above", "right"] = "above"
    plotly_options: dict = field(default_factory=dict)

    def __post_init__(self):
        """Ensures data is loaded during initialization."""
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

    def get_subset_region_start(self, subset_region: str, **kwargs) -> int:
        """
        Extracts the start position from the subset region, or returns the minimum x value.

        Args:
            subset_region (str): Region in "chrom:start-end" format.

        Returns:
            int: The starting position of the region.
        """
        if subset_region is not None:
            return int(subset_region.split(":")[1].split("-")[0])
        return kwargs.get("xmin", 0)

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
        ] = (
            y + self.rect_height
        )
        return bed_region_coverage.maximum(sp.csr_matrix(new_coverage))

    def plot_mpl(self, ax: Axes, **kwargs):
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
        bed_region_coverage = kwargs.get("bed_region_coverage", np.zeros((1, 100)))
        subset_region = kwargs.get("subset_region", None)
        subset_region_start = self.get_subset_region_start(subset_region, **kwargs)
        axis_shift = kwargs.get("axis_shift", 0)

        cleaned_data = self.format_data(
            subset_region=subset_region, axis_shift=axis_shift
        )

        for idx, region in cleaned_data.iterrows():
            # Set color if itemRGB is present and use_color_column is enabled
            if "itemRGB" in region and self.use_color_column:
                self.mpl_rect_options["color"] = [
                    int(val) / 255 for val in region["itemRGB"].split(",")
                ]

            # Calculate y-position from bed_region_coverage
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

            # Add rectangle patch
            region_rect = Rectangle(
                (region["chromStart"], y + self.rect_padding),
                region["chromEnd"] - region["chromStart"],
                self.rect_height - 2 * self.rect_padding,
                label=region["name"],
                **self.mpl_rect_options,
            )
            ax.add_patch(region_rect)

            # Add text label
            if self.mpl_text_alignment == "above":
                ax.text(
                    (region["chromEnd"] + region["chromStart"]) / 2,
                    y + self.rect_height + self.rect_padding,
                    region["name"],
                    **self.mpl_text_options,
                )
            elif self.mpl_text_alignment == "left":
                ax.text(
                    region["chromStart"],
                    y + self.rect_height / 2 - self.rect_padding,
                    region["name"],
                    **self.mpl_text_options,
                )
            elif self.mpl_text_alignment == "right":
                ax.text(
                    region["chromEnd"],
                    y + self.rect_height / 2 - self.rect_padding,
                    region["name"],
                    **self.mpl_text_options,
                )

            # Update coverage
            bed_region_coverage = self.update_coverage(
                region, subset_region_start, bed_region_coverage, y
            )

        ax.xaxis.set_tick_params(bottom=False)
        ax.yaxis.set_tick_params(left=False, labelleft=False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

        # if self.use_global_max:
        if self.use_global_max:
            ax.set_ylim(
                0, bed_region_coverage.max() * 1.1
            )  # Force the y-limit to the calculated maximum
        else:
            ax.autoscale(
                enable=True, axis="y"
            )  # Let Matplotlib auto-scale if use_global_max is False
        if self.show_legend:
            ax.legend()

        return bed_region_coverage

    def plot_plotly(self, fig: go.Figure, row: int, col: int, **kwargs):
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
        bed_region_coverage = kwargs.get("bed_region_coverage", np.zeros((1, 100)))
        subset_region = kwargs.get("subset_region", None)
        subset_region_start = self.get_subset_region_start(subset_region, **kwargs)
        axis_shift = kwargs.get("axis_shift", 0)

        cleaned_data = self.format_data(
            subset_region=subset_region, axis_shift=axis_shift
        )

        for idx, region in cleaned_data.iterrows():
            # Set color if itemRGB is present and use_color_column is enabled
            if "itemRGB" in region and self.use_color_column:
                self.plotly_options["fillcolor"] = f"rgb({region['itemRGB']})"

            # Calculate y-position from bed_region_coverage
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

            # Add polygon trace
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
                    fill="toself",
                    marker=dict(color="rgba(0,0,0,0)"),
                    showlegend=self.show_legend,
                    name=region["name"],
                    **self.plotly_options,
                ),
                row=row,
                col=col,
            )

            # Update coverage
            bed_region_coverage = self.update_coverage(
                region, subset_region_start, bed_region_coverage, y
            )

        if "showticklabels" not in self.plotly_yaxis_options:
            fig.update_yaxes(
                showticklabels=False, row=row, col=col, **self.plotly_yaxis_options
            )
        else:
            fig.update_yaxes(row=row, col=col, **self.plotly_yaxis_options)

        if self.use_global_max:
            fig.update_yaxes(
                range=[0, bed_region_coverage.max() * 1.1], row=row, col=col
            )

        return bed_region_coverage
