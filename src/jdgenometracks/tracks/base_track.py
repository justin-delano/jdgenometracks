"""
Base track class for jdgenometracks V3.

This module provides the foundational BaseTrack class that all specific track
types inherit from. It defines the interface for data loading, validation,
and common track functionality.

Author: Assistant
Date: 2024
"""

from __future__ import annotations

import os
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Union

import pandas as pd
import plotly.graph_objects as go
from matplotlib.axes import Axes

from ..utils.constants import TrackConstants
from ..utils.validation import validate_genomic_region


@dataclass
class BaseTrack(ABC):
    """
    Abstract base class for all genomic tracks.

    This class provides the interface and common functionality that all track types
    must implement. Each track type is responsible for loading its own data format
    and implementing its own plotting methods.

    Attributes:
        file_path: Path to the data file (optional if data is provided directly)
        data: Pre-loaded DataFrame (optional if file_path is provided)
        track_name: Human-readable name for the track
        track_type: Type of track (e.g., 'bed', 'bedgraph')
        track_options: Dictionary containing all track styling and behavior options
        show_legend: Whether to show legend for this track
        hlines: List of y-values for horizontal reference lines
        subplot_x: Column index in subplot grid (0 = first column)
        subplot_y: Row index in subplot grid (0 = first row)
        _is_loaded: Internal flag indicating if data has been loaded
    """

    # Base track constants
    DEFAULT_TRACK_HEIGHT = 1.0  # Default track height in track units
    DEFAULT_MIN_COLUMNS = 3  # Default minimum columns for empty DataFrame

    # File extension parsing constants
    EXTENSION_START_INDEX = 1  # Index to start extracting extension (skip the '.')
    FILENAME_COMPONENT_INDEX = 0  # Index for filename component when splitting by '.'

    # Horizontal line styling constants
    DEFAULT_HLINE_COLOR = "gray"  # Default horizontal line color
    DEFAULT_HLINE_STYLE_MPL = "--"  # Default matplotlib line style
    DEFAULT_HLINE_STYLE_PLOTLY = "dash"  # Default Plotly line style
    DEFAULT_HLINE_Y_POSITION = 0  # Default y-position for horizontal lines

    file_path: Optional[str] = None
    data: Optional[Union[pd.DataFrame, List[Dict[str, Any]]]] = None
    track_name: Optional[str] = None
    track_type: Optional[str] = None
    track_options: Dict[str, Any] = field(default_factory=dict)
    show_legend: bool = False
    hlines: List[float] = field(default_factory=list)
    subplot_x: int = 0  # Column index in subplot grid (default: single column layout)
    subplot_y: int = 0  # Row index in subplot grid (default: first row)
    _is_loaded: bool = field(default=False, init=False)

    def __post_init__(self):
        """
        Validates and initializes the instance after dataclass instantiation.
        """
        # Ensure either data or file path is provided
        if self.data is None and self.file_path is None:
            raise ValueError("Either data or file_path must be provided.")

        # Convert data to DataFrame if provided as list of dictionaries
        if self.data is not None:
            if isinstance(self.data, list):
                # Convert list of dictionaries to DataFrame
                self.data = pd.DataFrame(self.data)
            elif not isinstance(self.data, pd.DataFrame):
                raise ValueError("Data must be a pandas DataFrame or list of dictionaries")
            
            self._is_loaded = True

        # Infer track type from file extension if not provided
        if self.track_type is None and self.file_path is not None:
            self.track_type = self._infer_track_type(self.file_path)

        # Validate track type if provided
        if self.track_type is not None:
            if self.track_type not in TrackConstants.SUPPORTED_TRACK_TYPES:
                raise ValueError(f"Unsupported track type: {self.track_type}")

        # Infer track name from the file name if not provided
        if self.track_name is None and self.file_path is not None:
            self.track_name = self._infer_track_name(self.file_path)
        elif self.track_name is None:
            self.track_name = f"{self.track_type or 'unknown'}_track"

    def _infer_track_type(self, file_path: str) -> str:
        """
        Infer the track type based on the file extension.

        Args:
            file_path: Path to the file

        Returns:
            The inferred track type ('bed', 'bedgraph', etc.)

        Raises:
            ValueError: If the file extension is not supported
        """
        extension = os.path.splitext(file_path)[self.EXTENSION_START_INDEX][
            self.EXTENSION_START_INDEX :
        ].lower()

        # Map file extensions to track types
        extension_map = {
            "bed": "bed",
            "bedgraph": "bedgraph",
            "bg": "bedgraph",
            "bdg": "bedgraph",
            "wig": "bedgraph",
            "wiggle": "bedgraph",
        }

        if extension in extension_map:
            return extension_map[extension]
        else:
            raise ValueError(f"Unsupported file extension: {extension}")

    def _infer_track_name(self, file_path: str) -> str:
        """
        Infer the track name based on the file name.

        Args:
            file_path: Path to the file

        Returns:
            The inferred track name
        """
        return os.path.splitext(os.path.basename(file_path))[
            self.FILENAME_COMPONENT_INDEX
        ]

    @abstractmethod
    def load_data(self, file_path: Optional[str] = None) -> pd.DataFrame:
        """
        Load data from the file path.

        Each track type must implement this method to handle its specific
        data format (BED, BedGraph, etc.).

        Args:
            file_path: Path to data file (uses self.file_path if None)

        Returns:
            DataFrame containing the loaded and processed data

        Raises:
            FileNotFoundError: If file doesn't exist
            ValueError: If file format is invalid
        """
        pass

    @abstractmethod
    def validate_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Validate and clean the loaded data.

        Each track type should implement validation specific to its format.

        Args:
            data: Raw loaded data

        Returns:
            Validated and cleaned data

        Raises:
            ValueError: If data is invalid
        """
        pass

    def get_data(self) -> pd.DataFrame:
        """
        Get the track data, loading it if necessary.

        Returns:
            DataFrame containing the track data
        """
        if not self._is_loaded:
            if self.file_path is None:
                raise ValueError("No file path provided for data loading")

            self.data = self.load_data()
            self.data = self.validate_data(self.data)
            self._is_loaded = True

        # Ensure data is not None and is a DataFrame at this point
        if self.data is None:
            raise ValueError("Data loading failed - no data available")
        
        # At this point, data should always be a DataFrame due to __post_init__ conversion
        # Type assertion to help the type checker understand this
        assert isinstance(self.data, pd.DataFrame), "Data should be DataFrame after initialization"
        return self.data

    def subset_by_region(self, region: str) -> pd.DataFrame:
        """
        Subset the track data by a genomic region.

        Args:
            region: Genomic region string (e.g., "chr1:1000-2000")

        Returns:
            DataFrame containing only data within the specified region
        """
        # Validate and parse the region
        chr_name, start, end = validate_genomic_region(region)

        # Get the data
        data = self.get_data()

        if data.empty:
            return data

        # Default implementation assumes standard genomic coordinate columns
        # Subclasses can override this for track-specific logic
        return self._filter_by_coordinates(data, start, end)

    # Common helper methods for subclasses

    def _read_tabbed_file(
        self,
        file_path: str,
        expected_columns: List[str],
        min_cols: Optional[int] = None,
        max_cols: Optional[int] = None,
    ) -> pd.DataFrame:
        """
        Read a tab-separated file with common error handling.

        Args:
            file_path: Path to the file
            expected_columns: List of expected column names
            min_cols: Minimum number of required columns
            max_cols: Maximum number of allowed columns

        Returns:
            DataFrame with assigned column names

        Raises:
            FileNotFoundError: If file doesn't exist
            ValueError: If file format is invalid
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

        try:
            # Read tab-separated file without headers
            data = pd.read_csv(file_path, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            # Return empty DataFrame with proper columns
            return pd.DataFrame(
                columns=expected_columns[: min_cols or self.DEFAULT_MIN_COLUMNS]
            )
        except Exception as e:
            raise ValueError(f"Error reading file {file_path}: {e}")

        # Validate number of columns
        num_cols = len(data.columns)
        if min_cols and num_cols < min_cols:
            raise ValueError(
                f"File must have at least {min_cols} columns, got {num_cols}"
            )
        if max_cols and num_cols > max_cols:
            raise ValueError(f"File has too many columns: {num_cols} (max {max_cols})")

        # Assign column names - generate names for extra columns if needed
        if num_cols <= len(expected_columns):
            data.columns = expected_columns[:num_cols]
        else:
            # Generate additional column names for extra columns
            column_names = list(expected_columns)
            for i in range(len(expected_columns), num_cols):
                column_names.append(f"extra_col_{i}")
            data.columns = column_names

        return data

    def _validate_required_columns(
        self, data: pd.DataFrame, required_cols: List[str]
    ) -> None:
        """
        Validate that required columns are present in the data.

        Args:
            data: DataFrame to validate
            required_cols: List of required column names

        Raises:
            ValueError: If required columns are missing
        """
        if data.empty:
            return  # Empty data is valid

        missing_cols = [col for col in required_cols if col not in data.columns]
        if missing_cols:
            raise ValueError(f"Data missing required columns: {missing_cols}")

    def _validate_coordinate_ranges(
        self,
        data: pd.DataFrame,
        start_col: str = "chromStart",
        end_col: str = "chromEnd",
    ) -> None:
        """
        Validate that genomic coordinates are valid (start < end).

        Args:
            data: DataFrame to validate
            start_col: Name of start coordinate column
            end_col: Name of end coordinate column

        Raises:
            ValueError: If invalid coordinates are found
        """
        if data.empty:
            return

        if start_col in data.columns and end_col in data.columns:
            invalid_coords = data[start_col] >= data[end_col]
            if invalid_coords.any():
                raise ValueError("Data contains invalid coordinates (start >= end)")

    def _validate_numeric_columns(
        self, data: pd.DataFrame, numeric_cols: List[str]
    ) -> None:
        """
        Validate that specified columns contain numeric values.

        Args:
            data: DataFrame to validate
            numeric_cols: List of columns that should be numeric

        Raises:
            ValueError: If non-numeric values are found
        """
        if data.empty:
            return

        for col in numeric_cols:
            if col in data.columns:
                non_numeric = pd.isna(pd.to_numeric(data[col], errors="coerce"))
                if non_numeric.any():
                    raise ValueError(f"Column '{col}' contains non-numeric values")

    def _apply_column_dtypes(
        self, data: pd.DataFrame, dtype_map: Dict[str, Any]
    ) -> pd.DataFrame:
        """
        Apply data types to DataFrame columns.

        Args:
            data: DataFrame to modify
            dtype_map: Dictionary mapping column names to data types

        Returns:
            DataFrame with applied data types
        """
        for col, dtype in dtype_map.items():
            if col in data.columns:
                if dtype in [int, float]:
                    data[col] = pd.to_numeric(data[col], errors="coerce")
                else:
                    data[col] = data[col].astype(dtype)
        return data

    def _filter_by_coordinates(
        self,
        data: pd.DataFrame,
        start: int,
        end: int,
        start_col: str = "chromStart",
        end_col: str = "chromEnd",
    ) -> pd.DataFrame:
        """
        Filter data by genomic coordinates.

        Args:
            data: DataFrame to filter
            start: Region start coordinate
            end: Region end coordinate
            start_col: Name of start coordinate column
            end_col: Name of end coordinate column

        Returns:
            Filtered DataFrame
        """
        if data.empty:
            return data

        # Filter by coordinates - keep regions that overlap with the query region
        filtered_data = data[(data[end_col] > start) & (data[start_col] < end)].copy()

        return filtered_data

    def _add_horizontal_lines_matplotlib(self, ax: Axes) -> None:
        """Add horizontal lines to matplotlib plot if specified in options."""
        hlines = self.track_options.get("hlines", [])
        if not hlines:
            return

        for hline in hlines:
            # Use get() instead of pop() to avoid modifying the original dict
            y_pos = hline.get("y", self.DEFAULT_HLINE_Y_POSITION)
            color = hline.get("color", self.DEFAULT_HLINE_COLOR)
            linestyle = hline.get("linestyle", self.DEFAULT_HLINE_STYLE_MPL)

            # Create a copy of hline without the special keys for **kwargs
            other_opts = {
                k: v for k, v in hline.items() if k not in ["y", "color", "linestyle"]
            }
            ax.axhline(y=y_pos, color=color, linestyle=linestyle, **other_opts)

    def _add_horizontal_lines_plotly(self, fig: go.Figure, row: int, col: int) -> None:
        """Add horizontal lines to Plotly plot if specified in options."""
        hlines = self.track_options.get("hlines", [])
        if not hlines:
            return

        for hline in hlines:
            # Use get() instead of pop() to avoid modifying the original dict
            y_pos = hline.get("y", self.DEFAULT_HLINE_Y_POSITION)
            color = hline.get("color", self.DEFAULT_HLINE_COLOR)
            linestyle = hline.get("linestyle", self.DEFAULT_HLINE_STYLE_PLOTLY)

            # Create a copy of hline without the special keys for **kwargs
            other_opts = {
                k: v for k, v in hline.items() if k not in ["y", "color", "linestyle"]
            }

            fig.add_hline(
                y=y_pos,
                row=row,  # type: ignore
                col=col,  # type: ignore
                line_color=color,
                line_dash=linestyle,
                **other_opts,
            )

    @abstractmethod
    def plot_matplotlib(self, ax: Axes, start: int, end: int, **kwargs) -> None:
        """
        Plot the track using matplotlib.

        Args:
            ax: Matplotlib axes object
            start: Start position of the region to plot
            end: End position of the region to plot
            **kwargs: Additional plotting parameters
        """
        pass

    @abstractmethod
    def plot_plotly(
        self, fig: go.Figure, row: int, col: int, start: int, end: int, **kwargs
    ) -> None:
        """
        Plot the track using Plotly.

        Args:
            fig: Plotly figure object
            row: Row in subplot grid
            col: Column in subplot grid
            start: Start position of the region to plot
            end: End position of the region to plot
            **kwargs: Additional plotting parameters
        """
        pass

    def get_height(self) -> float:
        """
        Get the track height from options.

        Returns:
            Track height in track units
        """
        return self.track_options.get("height", self.DEFAULT_TRACK_HEIGHT)

    def get_genomic_bounds(self) -> Optional[Dict[str, Union[str, int]]]:
        """
        Get the genomic bounds (chromosome, start, end) from the track data.

        Returns:
            Dictionary with 'chrom', 'start', 'end' keys, or None if no genomic data
        """
        # Load data if not already loaded
        if not self._is_loaded:
            try:
                self.data = self.load_data()
            except Exception:
                # If data can't be loaded, return None
                return None

        if self.data is None:
            return None

        # Ensure data is a DataFrame
        if isinstance(self.data, list):
            if not self.data:  # Empty list
                return None
            df = pd.DataFrame(self.data)
        else:
            df = self.data
            if df.empty:
                return None

        # Check if the required genomic columns exist
        if not all(col in df.columns for col in ['chrom', 'chromStart', 'chromEnd']):
            return None

        # Get the first chromosome (assuming all data is from the same chromosome)
        chrom = df['chrom'].iloc[0] if len(df['chrom']) > 0 else None
        if chrom is None:
            return None

        # Get min start and max end positions
        min_start = df['chromStart'].min()
        max_end = df['chromEnd'].max()

        return {
            'chrom': str(chrom),
            'start': int(min_start),
            'end': int(max_end)
        }

    def __repr__(self) -> str:
        """String representation of the track."""
        return f"{self.__class__.__name__}(name='{self.track_name}', type='{self.track_type}')"
