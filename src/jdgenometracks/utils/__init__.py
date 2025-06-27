"""
Utility functions for jdgenometracks V3.

Consolidated utility modules providing:
- Constants and configuration values
- Unit conversion and size handling
- Genomic coordinate utilities and tick calculation
- Layout optimization and grid management
- Data loading and preprocessing
- Input validation
- Option mapping and translation
"""

# Import constants
from .constants import (
    PYTHON_MIN_VERSION,
    VERSION,
    ConfigSchema,
    ErrorMessages,
    FileConstants,
    GenomicConstants,
    MatplotlibConstants,
    PlotConstants,
    PlotlyConstants,
    TrackConstants,
    UnitConstants,
    get_default_column_regions,
    get_default_column_titles,
    get_default_height_props,
    get_default_row_titles,
    get_default_width_props,
)

# Import coordinate utilities
from .coordinates import (
    GenomicRegion,
    GenomicTickCalculator,
    calculate_axis_shift,
    format_genomic_coordinate,
    format_genomic_ticks,
    parse_region_string,
)

# Import data loading utilities
from .data_loading import DataLoader, DataPreprocessor, load_and_preprocess_file

# Import layout utilities
from .layout import (
    LayoutManager,
    assign_y_levels_sweep_line,
    calculate_grid_dimensions,
    create_subplot_grid_spec,
    normalize_layout_proportions,
    optimize_track_spacing,
)

# Import option mapping utilities
from .options import (
    INTERNAL_KEYS,
    PRESET_STYLES,
    UNIFIED_STYLE_MAP,
    create_fill_options,
    create_line_options,
    create_marker_options,
    create_text_options,
    flatten,
    get_internal_options,
    get_preset_style,
    merge_options,
    normalize_options,
    translate,
)

# Import unit conversion utilities
from .units import (
    SIZE_PRESETS,
    convert_to_inches,
    convert_to_pixels,
    get_preset_size,
    normalize_size_for_backend,
    parse_size_with_units,
)

# Import validation utilities (consolidated)
from .validation import (
    ConfigValidator,
    DataValidator,
    ParameterValidator,
    get_file_type_from_extension,
    sanitize_options,
    validate_bed_data,
    validate_bedgraph_data,
    validate_genomic_region,
    validate_grid_layout,
    validate_label_alignment,
    validate_plot_type,
    validate_plotting_config,
    validate_track_type,
    validate_unified_options,
)

__all__ = [
    # Constants and configuration
    "VERSION",
    "PYTHON_MIN_VERSION",
    "PlotConstants",
    "UnitConstants",
    "GenomicConstants",
    "FileConstants",
    "TrackConstants",
    "MatplotlibConstants",
    "PlotlyConstants",
    "ErrorMessages",
    "ConfigSchema",
    # Validation helpers
    "validate_plot_type",
    "validate_label_alignment",
    "validate_track_type",
    "get_file_type_from_extension",
    # Default generators
    "get_default_height_props",
    "get_default_width_props",
    "get_default_row_titles",
    "get_default_column_titles",
    "get_default_column_regions",
    # Unit conversion
    "parse_size_with_units",
    "convert_to_inches",
    "convert_to_pixels",
    "normalize_size_for_backend",
    "get_preset_size",
    "SIZE_PRESETS",
    # Coordinate utilities
    "GenomicRegion",
    "GenomicTickCalculator",
    "format_genomic_coordinate",
    "format_genomic_ticks",
    "calculate_axis_shift",
    "parse_region_string",
    # Layout utilities
    "assign_y_levels_sweep_line",
    "calculate_grid_dimensions",
    "validate_grid_layout",
    "normalize_layout_proportions",
    "create_subplot_grid_spec",
    "optimize_track_spacing",
    "LayoutManager",
    # Data loading
    "DataLoader",
    "DataPreprocessor",
    "load_and_preprocess_file",
    # Validation
    "ConfigValidator",
    "DataValidator",
    "ParameterValidator",
    "validate_plotting_config",
    "sanitize_options",
    # Option mapping
    "flatten",
    "translate",
    "get_internal_options",
    "normalize_options",
    "merge_options",
    "validate_unified_options",
    "create_line_options",
    "create_marker_options",
    "create_fill_options",
    "create_text_options",
    "get_preset_style",
    "UNIFIED_STYLE_MAP",
    "PRESET_STYLES",
    "INTERNAL_KEYS",
]
