from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Literal, Optional, Union

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from matplotlib.axes import Axes

from .GenomeTrack import GenomeTrack

# Define common types for scalar values and array-like structures
Scalar_Type = Union[int, float]
ArrayLike_Type = Union[Iterable[Scalar_Type], np.ndarray]

# Global dictionary for plot types and their specific functions for MPL and Plotly
PLOT_TYPES = {
    "lines": {
        "mpl": {"plot": Axes.plot, "fill": True},
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
        ymin (Optional[float]): Minimum y-axis limit for the plot.
        ymax (Optional[float]): Maximum y-axis limit for the plot.
        plot_type (Literal): Type of plot ('lines', 'bars', 'points').
        mpl_options (dict): Customization options for matplotlib plots.
        plotly_options (dict): Customization options for Plotly plots.
    """

    ymin: Optional[float] = None
    ymax: Optional[float] = None
    plot_type: Literal["lines", "bars", "points"] = "lines"
    mpl_plot_options: dict = field(default_factory=dict)
    mpl_fill_options: dict = field(default_factory=dict)
    plotly_options: dict = field(default_factory=dict)

    def __post_init__(self):
        """Initializes the class and ensures the plot_type is valid."""
        super().__post_init__()
        if self.plot_type not in PLOT_TYPES:
            raise ValueError(
                f"Invalid plot_type: {self.plot_type}. Must be one of {list(PLOT_TYPES.keys())}."
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
            data = pd.DataFrame(
                columns=["chrom", "chromStart", "chromEnd", "value", "name"]
            )

        # Ensure the dataframe contains valid columns
        possible_bedgraph_columns = ["chrom", "chromStart", "chromEnd", "value", "name"]
        data = data.iloc[:, : len(possible_bedgraph_columns)]
        data.columns = possible_bedgraph_columns[: len(data.columns)]
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

    def plot_mpl(self, ax: Axes, **kwargs):
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
        subset_region = kwargs.get("subset_region", None)
        axis_shift = kwargs.get("axis_shift", 0)
        cleaned_data = self.get_cleaned_data(
            subset_region=subset_region, axis_shift=axis_shift
        )

        if cleaned_data.empty:
            raise ValueError("No data available for the selected region.")

        mid_points = self.calculate_mid_points(cleaned_data)
        y_values = cleaned_data["value"].astype(float)

        plot_fn = PLOT_TYPES[self.plot_type]["mpl"]["plot"]
        fill_between = PLOT_TYPES[self.plot_type]["mpl"].get("fill", False)

        plot_fn(
            ax, mid_points, y_values, label=self.track_name, **self.mpl_plot_options
        )
        if fill_between:
            ax.fill_between(mid_points, y_values, 0, **self.mpl_fill_options)  # type: ignore

        # Set y-axis limits
        ymin = kwargs.get("ymin", self.ymin)
        ymax = kwargs.get("ymax", self.ymax)
        if ymin is not None:
            ax.set_ylim(bottom=ymin)
        if ymax is not None:
            ax.set_ylim(top=ymax)

        # Customize appearance
        ax.spines["bottom"].set_visible(False)
        ax.xaxis.set_tick_params(bottom=False)
        if self.show_legend:
            ax.legend()

    def plot_plotly(self, fig: go.Figure, row: int, col: int, **kwargs):
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
        subset_region = kwargs.get("subset_region", None)
        axis_shift = kwargs.get("axis_shift", 0)
        cleaned_data = self.get_cleaned_data(
            subset_region=subset_region, axis_shift=axis_shift
        )

        if cleaned_data.empty:
            raise ValueError("No data available for the selected region.")

        mid_points = self.calculate_mid_points(cleaned_data)
        y_values = cleaned_data["value"]

        plot_fn = PLOT_TYPES[self.plot_type]["plotly"]["plot"]
        plot_mode = PLOT_TYPES[self.plot_type]["plotly"].get("mode")
        fill_area = PLOT_TYPES[self.plot_type]["plotly"].get("fill")

        if self.plot_type == "bars":
            trace = plot_fn(
                x=mid_points,
                y=y_values,
                name=self.track_name,
                showlegend=self.show_legend,
                **self.plotly_options,
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
                **self.plotly_options,
            )

        fig.add_trace(trace, row=row, col=col)
        fig.update_xaxes(showline=False, row=row, col=col)
        fig.update_yaxes(linecolor="black", row=row, col=col)

        # Set y-axis limits
        ymin = kwargs.get("ymin", self.ymin)
        ymax = kwargs.get("ymax", self.ymax)
        if ymin is not None or ymax is not None:
            fig.update_yaxes(range=[ymin, ymax], row=row, col=col)
