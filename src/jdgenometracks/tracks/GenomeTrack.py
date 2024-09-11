from __future__ import annotations

import os
from dataclasses import dataclass, field

import pandas as pd
import plotly.graph_objects as go
from matplotlib.axes import Axes


@dataclass
class GenomeTrack:
    """
    A base class for genomic track visualization, providing support for handling
    BED and BedGraph files, including horizontal line annotations, and managing
    matplotlib and Plotly-based plotting.

    Attributes:
        file_path (str | None): Path to the input file for the track.
        track_type (str | None): Type of track ('bed', 'bedgraph').
        track_name (str | None): Name of the track.
        data (pd.DataFrame | None): Data for the track, if preloaded.
        axis_font_size (int): Font size for the axis labels.
        text_font_size (int): Font size for the text annotations.
        hlines (list | None): Horizontal lines to be drawn on the plot.
        mpl_hline_options (dict): Matplotlib options for horizontal lines.
        plotly_hline_options (dict): Plotly options for horizontal lines.
        height_prop (float | None): Proportional height of the track.
        share_with_previous (bool): Whether the track shares the x-axis with the previous track.
        x_axis_type (str | None): Type of the x-axis (e.g., 'linear').
        x_axis_interval (int | None): Interval for x-axis ticks.
        axis_ticks (bool): Whether to show axis ticks.
        show_legend (bool): Whether to display the legend in the plot.
    """

    file_path: str | None = None
    track_type: str | None = None
    track_name: str | None = None
    data: pd.DataFrame | None = None

    # Horizontal Lines
    hlines: list[float] | None = None
    mpl_hline_options: dict = field(default_factory=dict)
    plotly_hline_options: dict = field(default_factory=dict)

    # X-axis properties
    show_xaxis_ticks: bool = False
    mpl_xaxis_options: dict = field(default_factory=dict)
    plotly_xaxis_options: dict = field(default_factory=dict)

    # Y-axis properties
    show_yaxis_ticks: bool = True
    mpl_yaxis_options: dict = field(default_factory=dict)
    plotly_yaxis_options: dict = field(default_factory=dict)

    # Structure
    height_prop: float | None = None
    share_with_previous: bool = False
    mpl_legend_options: dict = field(default_factory=dict)
    show_legend: bool = False

    def __post_init__(self):
        """
        Validates and initializes the instance after dataclass instantiation.
        Loads the track name and type based on the file path if not explicitly provided.
        """
        # Ensure data or file path is provided
        if self.data is None:
            assert (
                self.file_path is not None
            ), "Either data or file_path must be provided."

        # Infer track type from file extension if not provided
        if self.track_type is None and self.file_path is not None:
            self.track_type = self._infer_track_type(self.file_path)
            assert self.track_type in [
                "bed",
                "bedgraph",
            ], f"Unsupported track type: {self.track_type}"

        # Infer track name from the file name if not provided
        if self.track_name is None and self.file_path is not None:
            self.track_name = self._infer_track_name(self.file_path)

    def _infer_track_type(self, file_path: str) -> str:
        """
        Infers the track type based on the file extension.

        Args:
            file_path (str): Path to the file.

        Returns:
            str: The inferred track type ('bed', 'bedgraph').
        """
        extension = os.path.splitext(file_path)[1][1:].lower()
        if extension not in ["bed", "bedgraph"]:
            raise ValueError(f"Unsupported file extension: {extension}")
        return extension

    def _infer_track_name(self, file_path: str) -> str:
        """
        Infers the track name from the file name.

        Args:
            file_path (str): Path to the file.

        Returns:
            str: The inferred track name.
        """
        return os.path.basename(file_path).split(".")[0]

    @staticmethod
    def set_df_col_dtype(data: pd.DataFrame) -> pd.DataFrame:
        """
        Sets the appropriate data types for common columns in genomic data files.

        Args:
            data (pd.DataFrame): The DataFrame to set column data types for.

        Returns:
            pd.DataFrame: The DataFrame with updated column data types.
        """
        dtypes = {
            "chrom": str,
            "chromStart": int,
            "chromEnd": int,
            "name": str,
            "score": float,
            "strand": str,
            "thickStart": int,
            "thickEnd": int,
            "itemRGB": str,
            "value": float,
        }
        for column in data.columns:
            if column in dtypes:
                data[column] = data[column].astype(dtypes[column])
        return data

    def add_hlines_mpl(self, ax: Axes):
        """
        Adds horizontal lines to a matplotlib plot.

        Args:
            ax (Axes): The matplotlib axis to add the horizontal lines to.
        """
        if self.hlines:
            for hline in self.hlines:
                ax.axhline(hline, **self.mpl_hline_options)

    def add_hlines_plotly(self, fig: go.Figure, row: int, col: int):
        """
        Adds horizontal lines to a Plotly plot.

        Args:
            fig (go.Figure): The Plotly figure to add the horizontal lines to.
            row (int): The row index in the subplot.
            col (int): The column index in the subplot.
        """
        if self.hlines:
            for hline in self.hlines:
                fig.add_hline(
                    y=hline,
                    row=row,  # type: ignore
                    col=col,  # type: ignore
                    **self.plotly_hline_options,
                )

    def format_data(
        self, subset_region: str | None = None, axis_shift: int = 0
    ) -> pd.DataFrame:
        """
        Formats the data by subsetting and applying an axis shift, if specified.

        Args:
            subset_region (str | None): Genomic region in the format 'chrom:start-end'.
            axis_shift (int): Number of base pairs to shift the genomic coordinates.

        Returns:
            pd.DataFrame: The formatted data.
        """
        assert self.data is not None, "Data must be provided to format"
        formatted_data = self.data

        # Subset the data to a specified region
        if subset_region:
            chrom, start, end = self._parse_region(subset_region)
            formatted_data = self.data[
                (self.data["chrom"] == chrom)
                & (
                    (self.data["chromStart"].between(start, end))
                    | (self.data["chromEnd"].between(start, end))
                )
            ].reset_index(drop=True)

        # Apply axis shift if specified
        if axis_shift:
            formatted_data["chromStart"] -= axis_shift
            formatted_data["chromEnd"] -= axis_shift

        return formatted_data

    def _parse_region(self, region: str):
        """
        Parses a genomic region in the format 'chrom:start-end'.

        Args:
            region (str): Genomic region in 'chrom:start-end' format.

        Returns:
            tuple: A tuple containing chromosome (str), start (int), and end (int).
        """
        chrom, positions = region.split(":")
        start, end = map(int, positions.split("-"))
        return chrom, start, end

    def plot_mpl(self, ax: Axes, **kwargs):
        """
        Placeholder for plotting with matplotlib. Should be implemented by subclasses.

        Args:
            ax (Axes): The matplotlib axis to plot on.
        """
        raise NotImplementedError("Subclasses should implement this method.")

    def plot_plotly(self, fig: go.Figure, row: int, col: int, **kwargs):
        """
        Placeholder for plotting with Plotly. Should be implemented by subclasses.

        Args:
            fig (go.Figure): The Plotly figure to plot on.
            row (int): The row index in the subplot.
            col (int): The column index in the subplot.
        """
        raise NotImplementedError("Subclasses should implement this method.")
