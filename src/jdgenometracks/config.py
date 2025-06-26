"""
Configuration constants and default values for jdgenometracks.

This module contains all magic numbers, default values, and configuration
constants used throughout the jdgenometracks package.
"""

from __future__ import annotations

from typing import List

# =============================================================================
# VERSION INFORMATION
# =============================================================================

VERSION = "0.2.1"
PYTHON_MIN_VERSION = "3.8"

# =============================================================================
# PLOTTING DEFAULTS
# =============================================================================


class PlotDefaults:
    """Default values for plotting operations."""

    # Figure dimensions and spacing
    VERTICAL_SPACING = 0.02
    HORIZONTAL_SPACING = 0.05
    DEFAULT_HEIGHT_PROP = 1.0
    DEFAULT_WIDTH_PROP = 1.0

    # Margins and padding
    DEFAULT_MARGIN = {"l": 0.1, "r": 0.1, "t": 50, "b": 20, "pad": 4}

    # Font sizes
    DEFAULT_SUPTITLE_FONTSIZE = 16

    # Background colors
    DEFAULT_PLOT_BGCOLOR = "white"

    # Axis and coordinate defaults
    DEFAULT_SUBPLOT_X = 0
    DEFAULT_SUBPLOT_Y = 0
    DEFAULT_AXIS_SHIFT = 0
    RELATIVE_X_AXIS_START = 0
    RELATIVE_X_AXIS_OFFSET = 1

    # BED track specific
    BED_Y_LIMIT_MULTIPLIER = 1.1
    DEFAULT_BED_Y_VALUE = 0
    DEFAULT_MAX_BED_REGIONS = 0

    # Rectangle defaults for BED tracks
    DEFAULT_RECT_HEIGHT = 1
    DEFAULT_RECT_PADDING = 0


# =============================================================================
# UNIT CONVERSION CONSTANTS
# =============================================================================


class UnitConversion:
    """Constants for unit conversions."""

    # Pixels per inch (standard web/screen DPI)
    PIXELS_PER_INCH = 96.0

    # Metric conversions
    CM_PER_INCH = 2.54
    MM_PER_INCH = 25.4

    # Points (typography)
    POINTS_PER_INCH = 72.0


# =============================================================================
# GENOMIC COORDINATE FORMATTING
# =============================================================================


class GenomicFormatting:
    """Constants for genomic coordinate formatting."""

    # Base pair unit thresholds
    BASE_PAIR_THRESHOLD = 1e3
    KILOBASE_THRESHOLD = 7e5

    # Unit exponents
    BASE_EXPONENT = 0
    KILOBASE_EXPONENT = 3
    MEGABASE_EXPONENT = 6

    # Unit labels
    BASE_UNIT = "b"
    KILOBASE_UNIT = "Kb"
    MEGABASE_UNIT = "Mb"

    # Formatting precision
    DEFAULT_DECIMAL_PRECISION = 10.0
    PRECISION_TOLERANCE = 1e-8
    MAX_SIGNIFICANT_FIGURES = 3
    THRESHOLD_MULTIPLIER = 1e-3


# =============================================================================
# FILE FORMAT DEFAULTS
# =============================================================================


class FileDefaults:
    """Default values for file format handling."""

    # Supported file extensions and their types
    SUPPORTED_EXTENSIONS = {".bed": "bed", ".bedgraph": "bedgraph", ".bg": "bedgraph"}

    # Default track names
    DEFAULT_AXIS_TRACK_NAME = "Axis"
    DEFAULT_SPACER_TRACK_NAME = "Spacer"

    # Column names for different file formats
    BED_COLUMNS = [
        "chrom",
        "chromStart",
        "chromEnd",
        "name",
        "score",
        "strand",
        "thickStart",
        "thickEnd",
        "itemRGB",
        "blockCount",
        "blockSizes",
        "blockStarts",
    ]

    BEDGRAPH_COLUMNS = ["chrom", "chromStart", "chromEnd", "value", "name"]

    # Data types for columns
    BED_COLUMN_DTYPES = {
        "chrom": str,
        "chromStart": int,
        "chromEnd": int,
        "name": str,
        "score": float,
        "strand": str,
        "thickStart": int,
        "thickEnd": int,
        "itemRGB": str,
        "blockCount": int,
        "blockSizes": str,
        "blockStarts": str,
    }

    BEDGRAPH_COLUMN_DTYPES = {
        "chrom": str,
        "chromStart": int,
        "chromEnd": int,
        "value": float,
        "name": str,
    }


# =============================================================================
# TRACK CONFIGURATION
# =============================================================================


class TrackDefaults:
    """Default values for track configuration."""

    # Track options
    DEFAULT_SHOW_LEGEND = False
    DEFAULT_USE_COLOR_COLUMN = False

    # Plot types for bedgraph tracks
    SUPPORTED_PLOT_TYPES = ["lines", "bars", "points"]
    DEFAULT_PLOT_TYPE = "lines"

    # Label alignment options
    LABEL_ALIGNMENTS = ["above", "center", "below"]
    DEFAULT_LABEL_ALIGNMENT = "center"


# =============================================================================
# MATPLOTLIB SPECIFIC DEFAULTS
# =============================================================================


class MPLDefaults:
    """Matplotlib-specific default values."""

    # Axis sharing
    DEFAULT_SHAREX = "col"
    DEFAULT_SHAREY = False

    # Layout
    DEFAULT_LAYOUT = "constrained"

    # Spine visibility
    HIDE_TOP_SPINE = False
    HIDE_RIGHT_SPINE = False


# =============================================================================
# PLOTLY SPECIFIC DEFAULTS
# =============================================================================


class PlotlyDefaults:
    """Plotly-specific default values."""

    # Axis sharing
    DEFAULT_SHARED_XAXES = "columns"
    DEFAULT_SHARED_YAXES = "rows"

    # Autosize
    DEFAULT_AUTOSIZE = True

    # Show legend for empty traces
    DEFAULT_SHOWLEGEND_EMPTY = False


# =============================================================================
# ERROR MESSAGES
# =============================================================================


class ErrorMessages:
    """Standard error messages used throughout the package."""

    MISSING_DATA_OR_PATH = "Either data or file_path must be provided."
    UNSUPPORTED_TRACK_TYPE = "Unsupported track type: {track_type}"
    INVALID_FILE_PATH = "file_path must be a string if provided."
    INVALID_TRACK_TYPE = "track_type must be a string if provided."
    MISSING_FILE_PATH = (
        "file_path is required for track types other than 'axis' or 'spacer'"
    )
    INVALID_COLUMN_REGIONS = "column_regions must be a list/tuple with length equal to number of columns in tracks."
    NO_DATA_TO_PLOT = "No data available to plot."
    MISSING_PLOT_METHOD = "Track {track} does not have a 'plot_{backend}' method."
    INVALID_PLOT_TYPE = (
        "Invalid plot.type: {plot_type}. Must be one of {supported_types}."
    )
    INVALID_SIZE_VALUE = "Invalid size value: {value}"
    INVALID_SIZE_STRING = "Invalid size string: {string}"
    UNSUPPORTED_UNIT_INCHES = "Unsupported unit for inches conversion: {unit}"
    UNSUPPORTED_UNIT_PIXELS = "Unsupported unit for pixel conversion: {unit}"


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================


def get_default_height_props(num_rows: int) -> List[float]:
    """Get default height proportions for a given number of rows."""
    return [PlotDefaults.DEFAULT_HEIGHT_PROP] * num_rows


def get_default_width_props(num_cols: int) -> List[float]:
    """Get default width proportions for a given number of columns."""
    return [PlotDefaults.DEFAULT_WIDTH_PROP] * num_cols


def get_default_row_titles(num_rows: int) -> List[str]:
    """Get default (empty) row titles for a given number of rows."""
    return [""] * num_rows


def get_default_column_titles(num_cols: int) -> List[str]:
    """Get default (empty) column titles for a given number of columns."""
    return [""] * num_cols


def get_default_column_regions(num_cols: int) -> List[None]:
    """Get default (None) column regions for a given number of columns."""
    return [None] * num_cols


# =============================================================================
# VALIDATION HELPERS
# =============================================================================


def validate_plot_type(plot_type: str) -> bool:
    """Validate that a plot type is supported."""
    return plot_type in TrackDefaults.SUPPORTED_PLOT_TYPES


def validate_label_alignment(alignment: str) -> bool:
    """Validate that a label alignment is supported."""
    return alignment in TrackDefaults.LABEL_ALIGNMENTS


def get_file_type_from_extension(filepath: str) -> str | None:
    """Get the file type from a file extension."""
    import os

    ext = os.path.splitext(filepath)[1].lower()
    return FileDefaults.SUPPORTED_EXTENSIONS.get(ext)
