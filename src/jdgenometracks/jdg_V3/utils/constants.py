"""
Core constants and configuration for jdgenometracks V3.

Consolidated configuration module containing all constants, defaults, and helper functions
organized by functional area.
"""

from __future__ import annotations

from typing import List, Optional

# =============================================================================
# VERSION INFORMATION
# =============================================================================

VERSION = "3.0.0"
PYTHON_MIN_VERSION = "3.8"

# =============================================================================
# PLOTTING CONSTANTS
# =============================================================================


class PlotConstants:
    """Core plotting constants and defaults."""

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

    # Grid positioning
    DEFAULT_SUBPLOT_X = 0
    DEFAULT_SUBPLOT_Y = 0
    DEFAULT_AXIS_SHIFT = 0
    RELATIVE_X_AXIS_START = 0
    RELATIVE_X_AXIS_OFFSET = 1

    # BED track specific
    BED_Y_LIMIT_MULTIPLIER = 1.1
    DEFAULT_BED_Y_VALUE = 0
    DEFAULT_MAX_BED_REGIONS = 0
    DEFAULT_RECT_HEIGHT = 1
    DEFAULT_RECT_PADDING = 0


# =============================================================================
# UNIT CONVERSION CONSTANTS
# =============================================================================


class UnitConstants:
    """Constants for unit conversions between different measurement systems."""

    # Screen/Web DPI standard
    PIXELS_PER_INCH = 96.0

    # Metric conversions
    CM_PER_INCH = 2.54
    MM_PER_INCH = 25.4

    # Typography
    POINTS_PER_INCH = 72.0


# =============================================================================
# GENOMIC COORDINATE CONSTANTS
# =============================================================================


class GenomicConstants:
    """Constants for genomic coordinate handling and formatting."""

    # Scale thresholds
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

    # Formatting
    DEFAULT_DECIMAL_PRECISION = 10.0
    PRECISION_TOLERANCE = 1e-8
    MAX_SIGNIFICANT_FIGURES = 3
    THRESHOLD_MULTIPLIER = 1e-3


# =============================================================================
# FILE FORMAT CONSTANTS
# =============================================================================


class FileConstants:
    """File format specifications and defaults."""

    # Supported formats
    SUPPORTED_EXTENSIONS = {".bed": "bed", ".bedgraph": "bedgraph", ".bg": "bedgraph"}

    # Default names
    DEFAULT_AXIS_TRACK_NAME = "Axis"
    DEFAULT_SPACER_TRACK_NAME = "Spacer"

    # BED format columns (12-column specification)
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

    # BedGraph format columns
    BEDGRAPH_COLUMNS = ["chrom", "chromStart", "chromEnd", "value", "name"]

    # Data type mappings
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
# TRACK CONFIGURATION CONSTANTS
# =============================================================================


class TrackConstants:
    """Track-specific configuration constants."""

    # General options
    DEFAULT_SHOW_LEGEND = False
    DEFAULT_USE_COLOR_COLUMN = False

    # Plot types
    SUPPORTED_PLOT_TYPES = ["lines", "bars", "points"]
    DEFAULT_PLOT_TYPE = "lines"

    # Label alignments
    LABEL_ALIGNMENTS = ["above", "center", "below", "left", "right"]
    DEFAULT_LABEL_ALIGNMENT = "center"

    # Track types
    SUPPORTED_TRACK_TYPES = ["bed", "bedgraph", "axis", "spacer"]


# =============================================================================
# BACKEND-SPECIFIC CONSTANTS
# =============================================================================


class MatplotlibConstants:
    """Matplotlib-specific constants."""

    # Axis sharing
    DEFAULT_SHAREX = "col"
    DEFAULT_SHAREY = False

    # Layout
    DEFAULT_LAYOUT = "constrained"

    # Spine visibility - hide top and right spines for cleaner genomic plots
    HIDE_TOP_SPINE = True
    HIDE_RIGHT_SPINE = True


class PlotlyConstants:
    """Plotly-specific constants."""

    # Axis sharing
    DEFAULT_SHARED_XAXES = "columns"
    DEFAULT_SHARED_YAXES = "rows"

    # Layout
    DEFAULT_AUTOSIZE = True
    DEFAULT_SHOWLEGEND_EMPTY = False
    
    # Background and appearance
    DEFAULT_PLOT_BGCOLOR = "white"
    DEFAULT_PAPER_BGCOLOR = "white"
    
    # Margins (in pixels) - tighter margins for genomic plots
    DEFAULT_MARGIN_LEFT = 60
    DEFAULT_MARGIN_RIGHT = 20
    DEFAULT_MARGIN_TOP = 80
    DEFAULT_MARGIN_BOTTOM = 60
    DEFAULT_MARGIN_PAD = 4
    
    # Grid lines - turn off by default for cleaner genomic plots
    DEFAULT_SHOWGRID_X = False
    DEFAULT_SHOWGRID_Y = False
    DEFAULT_ZEROLINE_X = False
    DEFAULT_ZEROLINE_Y = False


# =============================================================================
# ERROR MESSAGES
# =============================================================================


class ErrorMessages:
    """Centralized error message templates."""

    # Data and file errors
    MISSING_DATA_OR_PATH = "Either data or file_path must be provided."
    MISSING_FILE_PATH = (
        "file_path is required for track types other than 'axis' or 'spacer'"
    )
    INVALID_FILE_PATH = "file_path must be a string if provided."

    # Track errors
    UNSUPPORTED_TRACK_TYPE = "Unsupported track type: {track_type}"
    INVALID_TRACK_TYPE = "track_type must be a string if provided."
    MISSING_PLOT_METHOD = "Track {track} does not have a 'plot_{backend}' method."

    # Configuration errors
    INVALID_COLUMN_REGIONS = (
        "column_regions must be a list/tuple with length equal to number of columns."
    )
    INVALID_PLOT_TYPE = (
        "Invalid plot.type: {plot_type}. Must be one of {supported_types}."
    )

    # Data errors
    NO_DATA_TO_PLOT = "No data available to plot."

    # Unit conversion errors
    INVALID_SIZE_VALUE = "Invalid size value: {value}"
    INVALID_SIZE_STRING = "Invalid size string: {string}"
    UNSUPPORTED_UNIT_INCHES = "Unsupported unit for inches conversion: {unit}"
    UNSUPPORTED_UNIT_PIXELS = "Unsupported unit for pixel conversion: {unit}"


# =============================================================================
# DEFAULT VALUE GENERATORS
# =============================================================================


def get_default_height_props(num_rows: int) -> List[float]:
    """Generate default height proportions for a given number of rows."""
    return [PlotConstants.DEFAULT_HEIGHT_PROP] * num_rows


def get_default_width_props(num_cols: int) -> List[float]:
    """Generate default width proportions for a given number of columns."""
    return [PlotConstants.DEFAULT_WIDTH_PROP] * num_cols


def get_default_row_titles(num_rows: int) -> List[str]:
    """Generate default (empty) row titles for a given number of rows."""
    return [""] * num_rows


def get_default_column_titles(num_cols: int) -> List[str]:
    """Generate default (empty) column titles for a given number of columns."""
    return [""] * num_cols


def get_default_column_regions(num_cols: int) -> List[None]:
    """Generate default (None) column regions for a given number of columns."""
    return [None] * num_cols


# =============================================================================
# CONFIGURATION SCHEMA
# =============================================================================


class ConfigSchema:
    """Schema definitions for configuration validation."""

    # Required fields for different track types
    TRACK_REQUIRED_FIELDS = {
        "bed": ["track_type"],
        "bedgraph": ["track_type"],
        "axis": ["track_type"],
        "spacer": ["track_type"],
    }

    # Optional fields with defaults
    TRACK_OPTIONAL_FIELDS = {
        "file_path": None,
        "track_name": None,
        "subplot_x": PlotConstants.DEFAULT_SUBPLOT_X,
        "subplot_y": PlotConstants.DEFAULT_SUBPLOT_Y,
        "show_legend": TrackConstants.DEFAULT_SHOW_LEGEND,
        "options": {},
    }

    # Figure options schema
    FIGURE_OPTIONAL_FIELDS = {
        "backend": "plotly",
        "total_height": None,
        "total_width": None,
        "column_regions": None,
        "height_props": None,
        "width_props": None,
        "row_titles": None,
        "column_titles": None,
        "plot_title": None,
        "relative_x_axis": False,
    }
