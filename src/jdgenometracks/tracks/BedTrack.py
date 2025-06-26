from dataclasses import dataclass, field
from typing import Optional

import pandas as pd
import plotly.graph_objects as go
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle

from jdgenometracks.config import FileDefaults, PlotDefaults
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
                columns=FileDefaults.BED_COLUMNS[:5]  # Use first 5 columns as minimum
            )

        # Define possible BED columns and adjust based on data size
        data.columns = FileDefaults.BED_COLUMNS[: len(data.columns)]
        self.data = self.set_df_col_dtype(data)

    def plot_mpl(
        self,
        ax: Axes,
        track_y_levels: dict[int, float],
        subset_region: Optional[str] = None,
        axis_shift: Optional[int] = None,
        xmin: Optional[int] = None,
    ) -> None:
        """
        Plots the BED data using matplotlib. Draws rectangles for genomic regions.

        Args:
            ax (Axes): The matplotlib axis to plot on.
            track_y_levels: Dictionary mapping region index to y-level.
            subset_region (Optional[str]): The region to subset for plotting.
            axis_shift (Optional[int]): The genomic axis shift for visualization.
            xmin (Optional[int]): Minimum x-coordinate.
        """
        cleaned_data = self.format_data(
            subset_region=subset_region, axis_shift=axis_shift
        )

        # Use unified options
        mpl_opts = translate(self.track_options, target="mpl")
        rect_opts = mpl_opts.get("marker", {})
        text_opts = mpl_opts.get("text", {})
        legend_opts = mpl_opts.get("legend", {})

        rect_height = self.track_options.get(
            "rect.height", PlotDefaults.DEFAULT_RECT_HEIGHT
        )
        rect_padding = self.track_options.get(
            "rect.padding", PlotDefaults.DEFAULT_RECT_PADDING
        )

        max_y = 0
        for region_idx, (_, region) in enumerate(cleaned_data.iterrows()):
            if "itemRGB" in cleaned_data.columns and self.track_options.get(
                "use_color_column", False
            ):
                rect_opts["color"] = [
                    int(val) / 255 for val in region["itemRGB"].split(",")
                ]

            # Use precomputed y-level
            y = track_y_levels.get(region_idx, PlotDefaults.DEFAULT_BED_Y_VALUE)
            max_y = max(max_y, y + rect_height)

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

        ax.xaxis.set_tick_params(bottom=False)
        ax.yaxis.set_tick_params(left=False, labelleft=False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

        if self.track_options.get("use_global_max", False):
            ax.set_ylim(0, max_y + rect_padding)
        else:
            ax.autoscale(enable=True, axis="y")

        if self.show_legend:
            ax.legend(**legend_opts)

    def plot_plotly(
        self,
        fig: go.Figure,
        row: int,
        col: int,
        track_y_levels: dict[int, float],
        subset_region: Optional[str] = None,
        axis_shift: Optional[int] = None,
        xmin: Optional[int] = None,
    ) -> None:
        """
        Plots the BED data using Plotly. Draws filled polygons for genomic regions.

        Args:
            fig (go.Figure): The Plotly figure to add traces to.
            row (int): The row index in the subplot grid.
            col (int): The column index in the subplot grid.
            track_y_levels: Dictionary mapping region index to y-level.
            subset_region (Optional[str]): The region to subset for plotting.
            axis_shift (Optional[int]): The genomic axis shift for visualization.
            xmin (Optional[int]): Minimum x-coordinate.
        """
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

        rect_height = self.track_options.get("rect.height", 1)
        rect_padding = self.track_options.get("rect.padding", 0)

        max_y = 0
        for region_idx, (_, region) in enumerate(cleaned_data.iterrows()):
            if "itemRGB" in cleaned_data.columns and self.track_options.get(
                "use_color_column", False
            ):
                fill_opts["fillcolor"] = f"rgb({region['itemRGB']})"

            # Use precomputed y-level
            y = track_y_levels.get(region_idx, 0)
            max_y = max(max_y, y + rect_height)

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

        fig.update_xaxes(row=row, col=col, **xaxis_opts)
        if self.track_options.get("use_global_max", False):
            fig.update_yaxes(range=[0, max_y + rect_padding], row=row, col=col)
        fig.update_yaxes(
            row=row,
            col=col,
            showticklabels=False,
            **yaxis_opts,
        )
