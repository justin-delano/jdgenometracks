from __future__ import annotations

from typing import Any, Optional, Sequence, Tuple

import pandas as pd
from matplotlib import axis
from matplotlib.pyplot import subplot

from .config import ErrorMessages, FileDefaults, PlotDefaults, TrackDefaults
from .tracks.BedGraphTrack import BedGraphTrack
from .tracks.BedTrack import BedTrack
from .tracks.GenomeTrack import GenomeTrack
from .tracks.SpacerTrack import SpacerTrack
from .tracks.XAxisTrack import XAxisTrack


class TrackFactory:
    @staticmethod
    def create_track(
        file_path: str | None = None,
        track_type: str | None = None,
        track_options: dict[str, Any] = {},
        **kwargs: Any,
    ) -> GenomeTrack:
        """
        Creates a GenomeTrack object based on the file path and type.
        Allows either inferring the track type from the file extension or using the provided `track_type`.

        Args:
            file_path (str): Path to the file containing track data.
            track_type (str, optional): Explicit track type input. If not provided, it is inferred from the file name.
            **kwargs: Additional keyword arguments for track creation.

        Returns:
            GenomeTrack: The created track object.

        Raises:
            ValueError: If required arguments are missing or invalid.
        """
        # Validate input
        if track_type is None and file_path is None:
            raise ValueError("Either 'track_type' or 'file_path' must be provided.")
        if track_type is not None and not isinstance(track_type, str):
            raise TypeError("track_type must be a string if provided.")
        if file_path is not None and not isinstance(file_path, str):
            raise TypeError("file_path must be a string if provided.")

        # Handle special track types (axis, spacer) without file_path
        if track_type == "axis":
            track_name = kwargs.pop(
                "track_name", FileDefaults.DEFAULT_AXIS_TRACK_NAME
            )  # Use pop to remove it from kwargs
            return XAxisTrack(
                data=pd.DataFrame(),
                file_path="Axis",
                track_name=track_name,
                track_type="Axis",
                track_options=track_options,
                **kwargs,
            )
        elif track_type == "spacer":
            track_name = kwargs.pop(
                "track_name", FileDefaults.DEFAULT_SPACER_TRACK_NAME
            )  # Use pop to remove it from kwargs
            return SpacerTrack(
                data=pd.DataFrame(),
                file_path="Spacer",
                track_name=track_name,
                track_type="Spacer",
                track_options=track_options,
                **kwargs,
            )

        # For other track types, file_path is required
        if file_path is None:
            raise ValueError(ErrorMessages.MISSING_FILE_PATH)

        # Infer track type and track name if not provided explicitly
        track_type = track_type or TrackFactory._infer_track_type(file_path)
        track_name = kwargs.pop("track_name", TrackFactory._infer_track_name(file_path))

        # Load the data and create the track
        data = TrackFactory._load_data(file_path, track_type)
        return TrackFactory._create_track_with_data(
            file_path,
            track_type,
            track_name,
            data,
            track_options=track_options,
            **kwargs,
        )

    @staticmethod
    def _infer_track_type(file_path: str) -> str:
        """
        Infers the track type based on the file extension.
        Defaults to 'bed' or 'bedgraph' but can be extended for more types.

        Args:
            file_path (str): The path to the file.

        Returns:
            str: The inferred track type.

        Raises:
            NotImplementedError: If the file extension is not supported.
        """
        extension = file_path.split(".")[-1].lower()
        if extension == "bed":
            return "bed"
        elif extension == "bedgraph":
            return "bedgraph"
        else:
            raise NotImplementedError(f"Unsupported file extension: {extension}")

    @staticmethod
    def _infer_track_name(file_path: str) -> str:
        """
        Infers the track name based on the file name.

        Args:
            file_path (str): The path to the file.

        Returns:
            str: The inferred track name.
        """
        return file_path.split("/")[-1].split(".")[0]

    @staticmethod
    def _load_data(file_path: str, track_type: str) -> pd.DataFrame:
        """
        Loads data from the specified file.

        Args:
            file_path (str): The path to the file.
            track_type (str): The type of track ('bed', 'bedgraph').

        Returns:
            pd.DataFrame: The loaded data.

        Raises:
            pd.errors.EmptyDataError: If the file is empty.
        """
        try:
            return pd.read_csv(file_path, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            return TrackFactory._create_empty_dataframe(track_type)

    @staticmethod
    def _create_empty_dataframe(track_type: str) -> pd.DataFrame:
        """
        Creates an empty DataFrame based on the track type.

        Args:
            track_type (str): The type of track ('bed', 'bedgraph').

        Returns:
            pd.DataFrame: An empty DataFrame with appropriate columns.
        """
        if track_type == "bed":
            return pd.DataFrame(
                columns=["chrom", "chromStart", "chromEnd", "name", "score"]
            )
        elif track_type == "bedgraph":
            return pd.DataFrame(
                columns=["chrom", "chromStart", "chromEnd", "value", "name"]
            )
        else:
            raise NotImplementedError(f"Track type '{track_type}' is not supported")

    @staticmethod
    def prepare_data(track_type: str, data: pd.DataFrame) -> pd.DataFrame:
        """
        Prepares data for a given track type ('bed', 'bedgraph').

        Args:
            track_type (str): The type of track ('bed', 'bedgraph').
            data (pd.DataFrame): The input data.

        Returns:
            pd.DataFrame: The prepared data.
        """
        if track_type == "bed":
            return TrackFactory._prepare_bed_data(data)
        elif track_type == "bedgraph":
            return TrackFactory._prepare_bedgraph_data(data)
        else:
            raise NotImplementedError(f"Track type '{track_type}' is not supported")

    @staticmethod
    def _create_track_with_data(
        file_path: str,
        track_type: str,
        track_name: str,
        data: pd.DataFrame,
        subplot_x: int = 0,
        subplot_y: int = 0,
        track_options: dict[str, Any] = {},
        show_legend: bool = False,
    ) -> GenomeTrack:
        """
        Creates a track with the loaded data.

        Args:
            file_path (str): Path to the file.
            track_type (str): The type of track ('bed', 'bedgraph').
            track_name (str): The name of the track.
            data (pd.DataFrame): The loaded data.
            **kwargs: Additional arguments to pass to the track.

        Returns:
            GenomeTrack: The constructed track object.
        """

        prepared_data = TrackFactory.prepare_data(track_type, data)
        if track_type == "bed":
            return BedTrack(
                file_path=file_path,
                track_type=track_type,
                track_name=track_name,
                data=prepared_data,
                track_options=track_options,
                show_legend=show_legend,
                subplot_x=subplot_x,
                subplot_y=subplot_y,
            )
        elif track_type == "bedgraph":
            return BedGraphTrack(
                file_path=file_path,
                track_type=track_type,
                track_name=track_name,
                data=prepared_data,
                track_options=track_options,
                show_legend=show_legend,
                subplot_x=subplot_x,
                subplot_y=subplot_y,
            )
        else:
            raise NotImplementedError(f"Track type '{track_type}' is not supported")

    @staticmethod
    def _prepare_bed_data(data: pd.DataFrame) -> pd.DataFrame:
        """
        Prepares BED data by assigning column names and setting the correct data types.

        Args:
            data (pd.DataFrame): The input data.

        Returns:
            pd.DataFrame: Prepared BED data.
        """
        bed_columns = [
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
        }
        data.columns = bed_columns[: len(data.columns)]
        for column in data.columns:
            data[column] = data[column].astype(dtypes.get(column, object))
        return data

    @staticmethod
    def _prepare_bedgraph_data(data: pd.DataFrame) -> pd.DataFrame:
        """
        Prepares BedGraph data by assigning column names and setting the correct data types.

        Args:
            data (pd.DataFrame): The input data.

        Returns:
            pd.DataFrame: Prepared BedGraph data.
        """
        bedgraph_columns = ["chrom", "chromStart", "chromEnd", "value", "name"]
        if len(data.columns) > len(bedgraph_columns):
            data = data.iloc[:, : len(bedgraph_columns)]

        dtypes = {
            "chrom": str,
            "chromStart": int,
            "chromEnd": int,
            "value": float,
            "name": str,
        }
        data.columns = bedgraph_columns[: len(data.columns)]
        for column in data.columns:
            data[column] = data[column].astype(dtypes.get(column, object))
        return data


class TrackUtils:
    """
    A utility class containing helper methods for track-related operations.
    """

    @staticmethod
    def get_height_props(tracks: list[list[Any]]) -> list[float]:
        """
        Computes the height proportions for the given tracks (list of lists).

        Args:
            tracks (list[list[Any]]): 2D list of tracks.

        Returns:
            list[float]: List of height proportions for each track row.
        """
        return [
            (
                row[0].height_prop
                if row[0] is not None
                and hasattr(row[0], "height_prop")
                and row[0].height_prop is not None
                else PlotDefaults.DEFAULT_HEIGHT_PROP
            )
            for row in tracks
            if not (
                row[0] is not None
                and hasattr(row[0], "share_with_previous")
                and row[0].share_with_previous
            )
        ]

    @staticmethod
    def get_xlim_bedlim(
        column_tracks: Sequence[Any], column_region: str | None
    ) -> tuple[int, int, int]:
        """
        Determines the x-axis limits and the maximum number of BED regions for a list of tracks (column).

        Args:
            tracks (list[Any]): List of tracks (column).
            column_region (str | None): The genomic region.

        Returns:
            tuple[int, int, int]: The minimum x-axis value, maximum x-axis value, and maximum number of BED regions.
        """
        max_bed_regions = PlotDefaults.DEFAULT_MAX_BED_REGIONS
        xmin = xmax = None
        for track in column_tracks:
            if not hasattr(track, "data") or track.data is None or track.data.empty:
                continue
            formatted_data = track.format_data(subset_region=column_region)
            if formatted_data.empty:
                continue
            xmin = min(
                xmin or formatted_data["chromStart"].min(),
                formatted_data["chromStart"].min(),
            )
            xmax = max(
                xmax or formatted_data["chromEnd"].max(),
                formatted_data["chromEnd"].max(),
            )
            if isinstance(track, BedTrack):
                max_bed_regions = max(max_bed_regions, formatted_data.shape[0])
        assert xmin is not None and xmax is not None, ErrorMessages.NO_DATA_TO_PLOT
        return xmin, xmax, max_bed_regions

    @staticmethod
    def get_column_limits(
        column_tracks: Sequence[Any],
        region: Optional[str] = None,
        relative_x_axis: bool = False,
    ):
        """
        Determines the x-axis limits and other options for a column of tracks (list).

        Args:
            column_tracks: A column of tracks.
            region (Optional[str]): Genomic region in 'chrom:start-end' format.
            relative_x_axis (bool): Whether the x-axis should start at 0.

        Returns:
            tuple[int, int, dict]: x-axis min, x-axis max, and a dictionary of additional options.
        """
        chromosome = (
            region.split(":")[0]
            if region
            else (
                column_tracks[0].data.iloc[0, 0]
                if hasattr(column_tracks[0], "data") and not column_tracks[0].data.empty
                else None
            )
        )
        xmin, xmax, max_bed_regions = TrackUtils.get_xlim_bedlim(column_tracks, region)
        axis_shift = PlotDefaults.DEFAULT_AXIS_SHIFT
        if relative_x_axis:
            axis_shift = xmin - PlotDefaults.RELATIVE_X_AXIS_OFFSET
            xmax, xmin = (
                xmax - xmin + PlotDefaults.RELATIVE_X_AXIS_OFFSET,
                PlotDefaults.RELATIVE_X_AXIS_START,
            )
        return chromosome, xmin, xmax, max_bed_regions, axis_shift
