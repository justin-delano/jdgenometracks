from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional, Union

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from matplotlib.axes import Axes

from jdgenometracks.config import ErrorMessages, FileDefaults, TrackDefaults

from .GenomeTrack import GenomeTrack

# Define common types for scalar values and array-like structures
Scalar_Type = Union[int, float]
ArrayLike_Type = Union[Iterable[Scalar_Type], np.ndarray]

# Global dictionary for plot types and their specific functions for MPL and Plotly
PLOT_TYPES = {
    "lines": {
        "mpl": {"plot": Axes.plot, "fill": "tozeroy"},
        "plotly": {"plot": go.Scatter, "mode": "lines", "fill": "tozeroy"},
    },
    "bars": {
        "mpl": {"plot": Axes.bar},
        "plotly": {"plot": go.Bar},
    },
    "points": {
        "mpl": {"plot": Axes.scatter},
        "plotly": {"plot": go.Scatter, "mode": "markers"},
    },
}


@dataclass
class BedGraphTrack(GenomeTrack):
    """
    A class for plotting genomic tracks from bedGraph data using matplotlib or Plotly.

    Attributes:
        plotly_options (dict): Customization options for Plotly plots.
    """

    track_options: dict = field(default_factory=dict)

    def __post_init__(self):
        super().__post_init__()
        if (
            self.track_options.get("plot.type")
            not in TrackDefaults.SUPPORTED_PLOT_TYPES
        ):
            raise ValueError(
                ErrorMessages.INVALID_PLOT_TYPE.format(
                    plot_type=self.track_options.get("plot.type"),
                    supported_types=TrackDefaults.SUPPORTED_PLOT_TYPES,
                )
            )
        if self.data is None:
            self.read_data()

    def read_data(self):
        """
        Reads bedGraph-like data from a file, and sets the data attribute.
        """
        if self.file_path is None:
            raise ValueError("File path must be provided.")

        try:
            data = pd.read_csv(self.file_path, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            data = pd.DataFrame(columns=FileDefaults.BEDGRAPH_COLUMNS)

        # Ensure the dataframe contains valid columns
        data = data.iloc[:, : len(FileDefaults.BEDGRAPH_COLUMNS)]
        data.columns = FileDefaults.BEDGRAPH_COLUMNS[: len(data.columns)]
        self.data = self.set_df_col_dtype(data)

    def get_cleaned_data(
        self, subset_region: Optional[str] = None, axis_shift: int = 0
    ) -> pd.DataFrame:
        """
        Retrieves and formats the data, applying a subset region and/or axis shift if provided.

        Args:
            subset_region (Optional[tuple]): Tuple of (start, end) coordinates to subset the data.
            axis_shift (int): Shift the genomic coordinates for visualization.

        Returns:
            pd.DataFrame: The cleaned subset of data ready for plotting.
        """
        return self.format_data(subset_region=subset_region, axis_shift=axis_shift)

    def calculate_mid_points(self, data: pd.DataFrame) -> pd.Series:
        """
        Calculates the midpoints of genomic regions for the x-axis in plots.

        Args:
            data (pd.DataFrame): The cleaned data.

        Returns:
            np.ndarray: The midpoints of the genomic regions.
        """
        return (data["chromStart"] + data["chromEnd"]) / 2

    def plot_mpl(
        self,
        ax: Axes,
        subset_region: Optional[str] = None,
        axis_shift: int = 0,
    ):
        """
        Plots the genomic track using matplotlib.

        Args:
            ax (Axes): The matplotlib axis to plot on.
            subset_region (Optional[tuple]): Subset the data to a specific region.
            axis_shift (Optional[int]): Shift the genomic coordinates for plotting.
            ymin (Optional[float]): Minimum y-axis value.
            ymax (Optional[float]): Maximum y-axis value.

        Raises:
            ValueError: If no data is available for the selected region.
        """
        from jdgenometracks.option_mapping import translate

        cleaned_data = self.get_cleaned_data(
            subset_region=subset_region, axis_shift=axis_shift
        )

        if cleaned_data.empty:
            raise ValueError("No data available for the selected region.")

        mid_points = self.calculate_mid_points(cleaned_data)
        y_values = cleaned_data["value"].astype(float)

        # Use unified options
        mpl_opts = translate(self.track_options, target="mpl")
        rect_opts = mpl_opts.get("marker", {})
        fill_opts = mpl_opts.get("fill", {})
        legend_opts = mpl_opts.get("legend", {})
        plot_fn = PLOT_TYPES[self.track_options.get("plot.type", "lines")]["mpl"][
            "plot"
        ]
        fill_between = PLOT_TYPES[self.track_options.get("plot.type", "lines")][
            "mpl"
        ].get("fill", False)

        plot_fn(ax, mid_points, y_values, label=self.track_name, **rect_opts)
        if fill_between:
            ax.fill_between(mid_points, y_values, 0, **fill_opts)  # type: ignore
        ymin = self.track_options.get("ymin", None)
        ymax = self.track_options.get("ymax", None)
        if ymin is not None:
            ax.set_ylim(bottom=ymin)
        if ymax is not None:
            ax.set_ylim(top=ymax)
        ax.spines["bottom"].set_visible(False)
        ax.xaxis.set_tick_params(bottom=False)
        if self.show_legend:
            ax.legend(**legend_opts)

    def plot_plotly(
        self,
        fig: go.Figure,
        row: int,
        col: int,
        subset_region: Optional[str] = None,
        axis_shift: int = 0,
    ):
        """
        Plots the genomic track using Plotly.

        Args:
            fig (go.Figure): The Plotly figure to add the trace to.
            row (int): The row number in the Plotly subplot.
            col (int): The column number in the Plotly subplot.
            subset_region (Optional[tuple]): Subset the data to a specific region.
            axis_shift (Optional[int]): Shift the genomic coordinates for plotting.

        Raises:
            ValueError: If no data is available for the selected region.
        """
        from jdgenometracks.option_mapping import translate

        cleaned_data = self.get_cleaned_data(
            subset_region=subset_region, axis_shift=axis_shift
        )

        if cleaned_data.empty:
            raise ValueError("No data available for the selected region.")

        mid_points = self.calculate_mid_points(cleaned_data)
        y_values = cleaned_data["value"]

        # Use unified options
        plotly_opts = translate(self.track_options, target="plotly")
        marker_opts = plotly_opts.get("marker", {})
        fill_opts = plotly_opts.get("fill", {})
        legend_opts = plotly_opts.get("legend", {})
        xaxis_opts = plotly_opts.get("xaxis", {})
        yaxis_opts = plotly_opts.get("yaxis", {})
        plot_fn = PLOT_TYPES[self.track_options.get("plot.type", "lines")]["plotly"][
            "plot"
        ]
        plot_mode = PLOT_TYPES[self.track_options.get("plot.type", "lines")][
            "plotly"
        ].get("mode")
        fill_area = PLOT_TYPES[self.track_options.get("plot.type", "lines")][
            "plotly"
        ].get("fill")
        plot_type = self.track_options.get("plot.type", "lines")
        if plot_type == "bars":
            trace = plot_fn(
                x=mid_points,
                y=y_values,
                name=self.track_name,
                showlegend=self.show_legend,
                marker=marker_opts,
                **legend_opts,
            )
            fig.update_layout(bargap=0)
        else:
            trace = plot_fn(
                x=mid_points,
                y=y_values,
                mode=plot_mode,
                fill=fill_area,
                name=self.track_name,
                showlegend=self.show_legend,
                marker=marker_opts,
                **fill_opts,
                **legend_opts,
            )
        fig.add_trace(trace, row=row, col=col)
        fig.update_xaxes(showline=False, row=row, col=col, **xaxis_opts)
        fig.update_yaxes(linecolor="black", row=row, col=col, **yaxis_opts)
        ymin = self.track_options.get("ymin", None)
        ymax = self.track_options.get("ymax", None)
        if ymin is not None or ymax is not None:
            fig.update_yaxes(range=[ymin, ymax], row=row, col=col)
