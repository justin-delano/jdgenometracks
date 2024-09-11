from __future__ import annotations

import numpy as np
import pandas as pd

from .tracks.BedGraphTrack import BedGraphTrack
from .tracks.BedTrack import BedTrack
from .tracks.GenomeTrack import GenomeTrack
from .tracks.SpacerTrack import SpacerTrack
from .tracks.XAxisTrack import XAxisTrack


class TrackFactory:
    @staticmethod
    def create_track(
        file_path: str | None = None, track_type: str | None = None, **kwargs
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
        """
        # Handle special track types (axis, spacer) without file_path
        if track_type == "axis":
            track_name = kwargs.pop(
                "track_name", "Axis"
            )  # Use pop to remove it from kwargs
            return XAxisTrack("Axis", "Axis", track_name, pd.DataFrame(), **kwargs)
        elif track_type == "spacer":
            track_name = kwargs.pop(
                "track_name", "Spacer"
            )  # Use pop to remove it from kwargs
            return SpacerTrack("Spacer", "Spacer", track_name, pd.DataFrame(), **kwargs)

        # For other track types, file_path is required
        if file_path is None:
            raise ValueError(
                "file_path is required for track types other than 'axis' or 'spacer'"
            )

        # Infer track type and track name if not provided explicitly
        track_type = track_type or TrackFactory._infer_track_type(file_path)
        track_name = kwargs.pop("track_name", TrackFactory._infer_track_name(file_path))

        # Load the data and create the track
        data = TrackFactory._load_data(file_path, track_type)
        return TrackFactory._create_track_with_data(
            file_path, track_type, track_name, data, **kwargs
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
    def _create_track_with_data(
        file_path: str, track_type: str, track_name: str, data: pd.DataFrame, **kwargs
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
        if track_type == "bed":
            return BedTrack(
                file_path,
                track_type,
                track_name,
                TrackFactory._prepare_bed_data(data),
                **kwargs,
            )
        elif track_type == "bedgraph":
            return BedGraphTrack(
                file_path,
                track_type,
                track_name,
                TrackFactory._prepare_bedgraph_data(data),
                **kwargs,
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
    def get_height_props(tracks: np.ndarray) -> list[float]:
        """
        Computes the height proportions for the given tracks.

        Args:
            tracks (np.ndarray): Array of tracks.

        Returns:
            list[float]: List of height proportions for each track.
        """
        return [
            track.height_prop if track.height_prop is not None else 1
            for track in tracks[:, 0]
            if not track.share_with_previous
        ]

    @staticmethod
    def get_xlim_bedlim(
        tracks: np.ndarray, column_region: str | None
    ) -> tuple[int, int, int]:
        """
        Determines the x-axis limits and the maximum number of BED regions.

        Args:
            tracks (np.ndarray): Array of tracks.
            column_region (str | None): The genomic region.

        Returns:
            tuple[int, int, int]: The minimum x-axis value, maximum x-axis value, and maximum number of BED regions.
        """
        max_bed_regions = 0
        xmin = xmax = None

        for track in tracks:
            if track.data is None or track.data.empty:
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

        assert xmin is not None and xmax is not None, "No data available to plot."
        return xmin, xmax, max_bed_regions

    @staticmethod
    def get_col_limits(
        tracks_col: np.ndarray, column_region: str | None, relative_x_axis: bool
    ) -> tuple[int, int, dict]:
        """
        Determines the x-axis limits and other options for a column of tracks.

        Args:
            tracks_col (np.ndarray): A column of tracks.
            column_region (str | None): Genomic region in 'chrom:start-end' format.
            relative_x_axis (bool): Whether the x-axis should start at 0.

        Returns:
            tuple[int, int, dict]: x-axis min, x-axis max, and a dictionary of additional options.
        """
        extra_options = {}
        chromosome = (
            column_region.split(":")[0]
            if column_region
            else tracks_col[0].data.iloc[0, 0]
        )

        extra_options.update({"chromosome": chromosome, "max_regions": 0})

        xmin, xmax, max_bed_regions = TrackUtils.get_xlim_bedlim(
            tracks_col, column_region
        )
        extra_options["max_regions"] = max_bed_regions

        if relative_x_axis:
            extra_options["axis_shift"] = xmin - 1
            xmax, xmin = xmax - xmin + 1, 0

        extra_options["xmin"], extra_options["xmax"] = xmin, xmax

        return xmin, xmax, extra_options
