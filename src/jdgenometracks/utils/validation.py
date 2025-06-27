"""
Input validation utilities for jdgenometracks V3.

Provides validation functions for:
- Configuration dictionaries and options
- Track specifications
- File paths and data
- Plot parameters
"""

import os
from typing import Any, Dict, List, Optional, Union

import pandas as pd

from .constants import ConfigSchema, ErrorMessages, FileConstants, TrackConstants


class ConfigValidator:
    """Validates configuration dictionaries and options."""

    @staticmethod
    def validate_track_config(track_config: Dict[str, Any]) -> None:
        """
        Validate a single track configuration.

        Args:
            track_config: Track configuration dictionary

        Raises:
            ValueError: If configuration is invalid
        """
        if not isinstance(track_config, dict):
            raise ValueError("Track configuration must be a dictionary")

        # Check required fields
        if "track_type" not in track_config:
            raise ValueError("Track configuration must include 'track_type'")

        track_type = track_config["track_type"]
        if not validate_track_type(track_type):
            raise ValueError(f"Invalid track_type: {track_type}")

        # Check file_path requirement for file-based tracks
        if track_type in ["bed", "bedgraph"] and "file_path" not in track_config:
            raise ValueError(f"Track type '{track_type}' requires 'file_path'")

        # Validate file_path if provided
        if "file_path" in track_config:
            file_path = track_config["file_path"]
            if not isinstance(file_path, str):
                raise ValueError("file_path must be a string")
            if not os.path.exists(file_path):
                raise ValueError(f"File not found: {file_path}")

        # Validate subplot coordinates
        subplot_x = track_config.get("subplot_x", 0)
        subplot_y = track_config.get("subplot_y", 0)
        if not isinstance(subplot_x, int) or subplot_x < 0:
            raise ValueError("subplot_x must be a non-negative integer")
        if not isinstance(subplot_y, int) or subplot_y < 0:
            raise ValueError("subplot_y must be a non-negative integer")

        # Validate track options if provided
        if "options" in track_config:
            ConfigValidator.validate_track_options(track_config["options"])

    @staticmethod
    def validate_track_options(options: Dict[str, Any]) -> None:
        """
        Validate track-specific options.

        Args:
            options: Track options dictionary

        Raises:
            ValueError: If options are invalid
        """
        if not isinstance(options, dict):
            raise ValueError("Track options must be a dictionary")

        # Validate plot type if specified
        if "plot.type" in options:
            plot_type = options["plot.type"]
            if not validate_plot_type(plot_type):
                raise ValueError(f"Invalid plot.type: {plot_type}")

        # Validate label alignment if specified
        if "label.alignment" in options:
            alignment = options["label.alignment"]
            if not validate_label_alignment(alignment):
                raise ValueError(f"Invalid label.alignment: {alignment}")

        # Validate numeric options
        numeric_options = [
            "rect.height",
            "rect.padding",
            "ymin",
            "ymax",
            "line.width",
            "marker.size",
            "fill.opacity",
        ]
        for option in numeric_options:
            if option in options:
                value = options[option]
                if not isinstance(value, (int, float)):
                    raise ValueError(f"Option '{option}' must be numeric")
                if value < 0:
                    raise ValueError(f"Option '{option}' must be non-negative")

    @staticmethod
    def validate_figure_options(figure_options: Dict[str, Any]) -> None:
        """
        Validate figure-level options.

        Args:
            figure_options: Figure options dictionary

        Raises:
            ValueError: If options are invalid
        """
        if not isinstance(figure_options, dict):
            raise ValueError("Figure options must be a dictionary")

        # Validate backend
        backend = figure_options.get("backend", "plotly")
        if backend not in ["matplotlib", "plotly"]:
            raise ValueError(f"Invalid backend: {backend}")

        # Validate dimensions if provided
        for dim in ["total_height", "total_width"]:
            if dim in figure_options:
                value = figure_options[dim]
                if value is not None and (
                    not isinstance(value, (int, float)) or value <= 0
                ):
                    raise ValueError(f"{dim} must be a positive number or None")

        # Validate proportions if provided
        for prop in ["height_props", "width_props"]:
            if prop in figure_options:
                props = figure_options[prop]
                if props is not None:
                    if not isinstance(props, list) or not all(
                        isinstance(p, (int, float)) for p in props
                    ):
                        raise ValueError(f"{prop} must be a list of numbers or None")
                    if any(p <= 0 for p in props):
                        raise ValueError(f"All values in {prop} must be positive")

        # Validate titles if provided
        for title in ["row_titles", "column_titles"]:
            if title in figure_options:
                titles = figure_options[title]
                if titles is not None:
                    if not isinstance(titles, list) or not all(
                        isinstance(t, str) for t in titles
                    ):
                        raise ValueError(f"{title} must be a list of strings or None")

    @staticmethod
    def validate_full_config(config: Dict[str, Any]) -> None:
        """
        Validate a complete configuration dictionary.

        Args:
            config: Full configuration dictionary

        Raises:
            ValueError: If configuration is invalid
        """
        if not isinstance(config, dict):
            raise ValueError("Configuration must be a dictionary")

        # Check required tracks
        if "tracks" not in config:
            raise ValueError("Configuration must include 'tracks'")

        tracks = config["tracks"]
        if not isinstance(tracks, list) or len(tracks) == 0:
            raise ValueError("'tracks' must be a non-empty list")

        # Validate each track
        for i, track_config in enumerate(tracks):
            try:
                ConfigValidator.validate_track_config(track_config)
            except ValueError as e:
                raise ValueError(f"Error in track {i+1}: {e}")

        # Validate figure options if provided
        if "figure_options" in config:
            ConfigValidator.validate_figure_options(config["figure_options"])


class DataValidator:
    """Validates data files and content."""

    @staticmethod
    def validate_file_exists(file_path: str) -> None:
        """
        Validate that a file exists and is readable.

        Args:
            file_path: Path to the file

        Raises:
            FileNotFoundError: If file doesn't exist
            PermissionError: If file is not readable
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

        if not os.path.isfile(file_path):
            raise ValueError(f"Path is not a file: {file_path}")

        if not os.access(file_path, os.R_OK):
            raise PermissionError(f"File is not readable: {file_path}")

    @staticmethod
    def validate_file_extension(file_path: str, expected_type: str) -> None:
        """
        Validate that a file has the expected extension.

        Args:
            file_path: Path to the file
            expected_type: Expected file type ("bed", "bedgraph")

        Raises:
            ValueError: If extension doesn't match expected type
        """
        ext = os.path.splitext(file_path)[1].lower()

        if expected_type == "bed" and ext not in [".bed"]:
            raise ValueError(f"Expected .bed file, got: {ext}")
        elif expected_type == "bedgraph" and ext not in [".bedgraph", ".bg"]:
            raise ValueError(f"Expected .bedgraph or .bg file, got: {ext}")

    @staticmethod
    def validate_file_size(file_path: str, max_size_mb: float = 100.0) -> None:
        """
        Validate that a file is not too large.

        Args:
            file_path: Path to the file
            max_size_mb: Maximum allowed size in MB

        Raises:
            ValueError: If file is too large
        """
        if not os.path.exists(file_path):
            return  # Will be caught by validate_file_exists

        size_bytes = os.path.getsize(file_path)
        size_mb = size_bytes / (1024 * 1024)

        if size_mb > max_size_mb:
            raise ValueError(f"File too large: {size_mb:.1f}MB (max: {max_size_mb}MB)")


class ParameterValidator:
    """Validates plotting parameters and options."""

    @staticmethod
    def validate_color_value(color: Any) -> None:
        """
        Validate a color value.

        Args:
            color: Color value to validate

        Raises:
            ValueError: If color is invalid
        """
        if isinstance(color, str):
            # Basic validation for string colors
            if len(color) == 0:
                raise ValueError("Color string cannot be empty")
            # Could add more sophisticated color validation here
        elif isinstance(color, (list, tuple)):
            # RGB/RGBA tuple validation
            if len(color) not in [3, 4]:
                raise ValueError("Color tuple must have 3 (RGB) or 4 (RGBA) values")
            if not all(isinstance(c, (int, float)) for c in color):
                raise ValueError("Color tuple values must be numeric")
            if not all(0 <= c <= 1 for c in color):
                raise ValueError("Color tuple values must be between 0 and 1")
        else:
            raise ValueError("Color must be a string or RGB/RGBA tuple")


# Basic type validations moved from constants.py
def validate_plot_type(plot_type: str) -> bool:
    """Validate that a plot type is supported."""
    return plot_type in TrackConstants.SUPPORTED_PLOT_TYPES


def validate_label_alignment(alignment: str) -> bool:
    """Validate that a label alignment is supported."""
    return alignment in TrackConstants.LABEL_ALIGNMENTS


def validate_track_type(track_type: str) -> bool:
    """Validate that a track type is supported."""
    return track_type in TrackConstants.SUPPORTED_TRACK_TYPES


def get_file_type_from_extension(filepath: str) -> Optional[str]:
    """Get the file type from a file extension."""
    ext = os.path.splitext(filepath)[1].lower()
    return FileConstants.SUPPORTED_EXTENSIONS.get(ext)


# Data format validations moved from data_loading.py
def validate_bed_data(data: pd.DataFrame) -> None:
    """
    Validate BED format data.

    Args:
        data: BED DataFrame to validate

    Raises:
        ValueError: If data is invalid
    """
    if data.empty:
        return  # Empty data is valid

    # Check required columns
    required_cols = ["chrom", "chromStart", "chromEnd"]
    missing_cols = [col for col in required_cols if col not in data.columns]
    if missing_cols:
        raise ValueError(f"BED data missing required columns: {missing_cols}")

    # Check coordinate validity
    if "chromStart" in data.columns and "chromEnd" in data.columns:
        invalid_coords = data["chromStart"] >= data["chromEnd"]
        if invalid_coords.any():
            raise ValueError("BED data contains invalid coordinates (start >= end)")


def validate_bedgraph_data(data: pd.DataFrame) -> None:
    """
    Validate BedGraph format data.

    Args:
        data: BedGraph DataFrame to validate

    Raises:
        ValueError: If data is invalid
    """
    if data.empty:
        return  # Empty data is valid

    # Check required columns
    required_cols = ["chrom", "chromStart", "chromEnd", "value"]
    missing_cols = [col for col in required_cols if col not in data.columns]
    if missing_cols:
        raise ValueError(f"BedGraph data missing required columns: {missing_cols}")

    # Check coordinate validity
    invalid_coords = data["chromStart"] >= data["chromEnd"]
    if invalid_coords.any():
        raise ValueError("BedGraph data contains invalid coordinates (start >= end)")

    # Check that values are numeric
    if "value" in data.columns:
        try:
            pd.to_numeric(data["value"])
        except (ValueError, TypeError):
            raise ValueError("BedGraph data contains non-numeric values")


# Options validation moved from options.py
def validate_unified_options(options: Dict[str, Any]) -> None:
    """
    Validate a unified options dictionary.

    Args:
        options: Options dictionary to validate

    Raises:
        ValueError: If options are invalid
    """
    if not isinstance(options, dict):
        raise ValueError("Options must be a dictionary")

    # Check for valid option categories
    valid_categories = {
        "line",
        "marker",
        "fill",
        "text",
        "legend",
        "hline",
        "plot",
        "rect",
        "label",
        "axis",
        "figure",
    }

    for key in options.keys():
        if "." in key:
            category = key.split(".")[0]
            if category not in valid_categories:
                print(f"Warning: Unknown option category '{category}' in key '{key}'")


# Grid validation (consolidated from layout.py and existing validate_grid_dimensions)
def validate_grid_layout(
    tracks: List[Any],
    height_props: Optional[List[float]] = None,
    width_props: Optional[List[float]] = None,
    row_titles: Optional[List[str]] = None,
    column_titles: Optional[List[str]] = None,
) -> None:
    """
    Validate that layout parameters match the grid dimensions.

    Consolidates validation from both layout.py and existing validate_grid_dimensions.

    Args:
        tracks: List of track objects
        height_props: Height proportions for each row
        width_props: Width proportions for each column
        row_titles: Titles for each row
        column_titles: Titles for each column

    Raises:
        ValueError: If dimensions don't match
    """
    if not tracks:
        raise ValueError("Track list cannot be empty")

    # Calculate actual grid dimensions
    max_row = max(t.subplot_y for t in tracks if t is not None)
    max_col = max(t.subplot_x for t in tracks if t is not None)
    num_rows = max_row + 1
    num_cols = max_col + 1

    # Validate height proportions
    if height_props is not None:
        if len(height_props) != num_rows:
            raise ValueError(
                f"height_props length ({len(height_props)}) must match number of rows ({num_rows})"
            )
        if any(h <= 0 for h in height_props):
            raise ValueError("All height proportions must be positive")

    # Validate width proportions
    if width_props is not None:
        if len(width_props) != num_cols:
            raise ValueError(
                f"width_props length ({len(width_props)}) must match number of columns ({num_cols})"
            )
        if any(w <= 0 for w in width_props):
            raise ValueError("All width proportions must be positive")

    # Validate row titles
    if row_titles is not None:
        if len(row_titles) != num_rows:
            raise ValueError(
                f"row_titles length ({len(row_titles)}) must match number of rows ({num_rows})"
            )

    # Validate column titles
    if column_titles is not None:
        if len(column_titles) != num_cols:
            raise ValueError(
                f"column_titles length ({len(column_titles)}) must match number of columns ({num_cols})"
            )


def validate_plotting_config(config: Dict[str, Any]) -> None:
    """
    Comprehensive validation of a plotting configuration.

    Args:
        config: Complete configuration dictionary

    Raises:
        ValueError: If any part of the configuration is invalid
    """
    # Validate overall structure
    ConfigValidator.validate_full_config(config)

    # Validate all referenced files exist
    for track_config in config["tracks"]:
        if "file_path" in track_config:
            DataValidator.validate_file_exists(track_config["file_path"])

    # Additional cross-validation could go here
    print("Configuration validation passed")


def sanitize_options(options: Dict[str, Any]) -> Dict[str, Any]:
    """
    Sanitize and normalize options dictionary.

    Args:
        options: Raw options dictionary

    Returns:
        Sanitized options dictionary
    """
    sanitized = {}

    for key, value in options.items():
        # Convert string values to appropriate types where possible
        if isinstance(value, str):
            # Try to convert numeric strings
            try:
                if "." in value:
                    sanitized[key] = float(value)
                else:
                    sanitized[key] = int(value)
            except ValueError:
                # Keep as string if conversion fails
                sanitized[key] = value
        else:
            sanitized[key] = value

    return sanitized


def validate_genomic_region(region: str) -> tuple[str, int, int]:
    """
    Parse and validate a genomic region in the format 'chrom:start-end'.

    Args:
        region: Genomic region in 'chrom:start-end' format

    Returns:
        A tuple containing chromosome (str), start (int), and end (int)

    Raises:
        ValueError: If the region format is invalid
    """
    if not isinstance(region, str):
        raise ValueError("Region must be a string")

    if ":" not in region:
        raise ValueError("Region must contain ':'")

    try:
        chrom, positions = region.split(":", 1)
        if "-" not in positions:
            raise ValueError("Region must contain '-' between coordinates")

        start_str, end_str = positions.split("-", 1)
        start, end = int(start_str), int(end_str)

        if start < 0 or end < 0:
            raise ValueError("Coordinates must be non-negative")

        if start >= end:
            raise ValueError("Start coordinate must be less than end coordinate")

        return chrom.strip(), start, end

    except ValueError as e:
        if "invalid literal for int()" in str(e):
            raise ValueError("Coordinates must be integers") from e
        raise
