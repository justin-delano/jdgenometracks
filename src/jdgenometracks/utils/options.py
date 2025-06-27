"""
Option mapping and translation utilities for jdgenometracks V3.

Consolidates option_mapping.py functionality with improvements:
- Unified style vocabulary for backend-agnostic plotting
- Translation between matplotlib and plotly options
- Option normalization and validation
"""

import math
from typing import Any, Callable, Dict, Mapping, Optional

# =============================================================================
# CONSTANTS AND HELPERS
# =============================================================================

SEP = "."  # Dotted-key separator used throughout


def flatten(
    data: Mapping[str, Any], *, sep: str = SEP, parent_key: str = ""
) -> Dict[str, Any]:
    """
    Return data with nested dictionaries flattened using dotted keys.

    Args:
        data: Dictionary to flatten
        sep: Separator for nested keys
        parent_key: Parent key prefix

    Returns:
        Flattened dictionary

    Examples:
        >>> flatten({'line': {'width': 2}})
        {'line.width': 2}
    """
    items: Dict[str, Any] = {}
    for key, value in data.items():
        new_key = f"{parent_key}{sep}{key}" if parent_key else key
        if isinstance(value, Mapping):
            items.update(flatten(value, sep=sep, parent_key=new_key))
        else:
            items[new_key] = value
    return items


def _apply_transform(value: Any, transform: Optional[Callable[[Any], Any]]) -> Any:
    """Apply a transformation function if provided."""
    return transform(value) if transform else value


# =============================================================================
# UNIFIED STYLE MAPPING
# =============================================================================

UNIFIED_STYLE_MAP: Dict[str, Mapping[str, Any]] = {
    # Marker styling (scatter plots, points)
    "marker.size": {
        "mpl": "s",  # marker area in points^2
        "plotly": "size",  # marker radius in px
        "transform": {
            "mpl→plotly": lambda area: math.sqrt(area),
            "plotly→mpl": lambda radius: radius**2,
        },
    },
    "marker.color": {"mpl": "c", "plotly": "color"},
    "marker.symbol": {"mpl": "marker", "plotly": "symbol"},
    "marker.edge_color": {"mpl": "edgecolors", "plotly": "line.color"},
    "marker.edge_width": {"mpl": "linewidths", "plotly": "line.width"},
    "marker.opacity": {"mpl": "alpha", "plotly": "opacity"},
    # Line styling
    "line.width": {"mpl": "linewidth", "plotly": "width"},
    "line.color": {"mpl": "color", "plotly": "color"},
    "line.dash": {
        "mpl": "linestyle",
        "plotly": "dash",
        "values": {
            "solid": {"mpl": "-", "plotly": "solid"},
            "dashed": {"mpl": "--", "plotly": "dash"},
            "dotted": {"mpl": ":", "plotly": "dot"},
            "dashdot": {"mpl": "-.", "plotly": "dashdot"},
        },
    },
    # Fill styling (areas, rectangles, polygons)
    "fill.color": {"mpl": "facecolor", "plotly": "fillcolor"},
    "fill.opacity": {"mpl": "alpha", "plotly": "opacity"},
    "fill.edge_color": {"mpl": "edgecolor", "plotly": "line.color"},
    "fill.edge_width": {"mpl": "linewidth", "plotly": "line.width"},
    # Text styling
    "text.font_family": {"mpl": "fontfamily", "plotly": "textfont.family"},
    "text.font_size": {"mpl": "fontsize", "plotly": "textfont.size"},
    "text.font_color": {"mpl": "color", "plotly": "textfont.color"},
    "text.angle": {"mpl": "rotation", "plotly": "textangle"},
    "text.vertical_alignment": {
        "mpl": "va",
        "plotly": "textposition",
        "values": {
            "top": {"mpl": "top", "plotly": "top center"},
            "center": {"mpl": "center", "plotly": "middle center"},
            "bottom": {"mpl": "bottom", "plotly": "bottom center"},
        },
    },
    "text.horizontal_alignment": {
        "mpl": "ha",
        "plotly": "textposition",
        "values": {
            "left": {"mpl": "left", "plotly": "middle left"},
            "center": {"mpl": "center", "plotly": "middle center"},
            "right": {"mpl": "right", "plotly": "middle right"},
        },
    },
    # Axis styling
    "xaxis.title": {"mpl": "xlabel", "plotly": "layout.xaxis.title.text"},
    "yaxis.title": {"mpl": "ylabel", "plotly": "layout.yaxis.title.text"},
    "xaxis.range": {"mpl": "xlim", "plotly": "layout.xaxis.range"},
    "yaxis.range": {"mpl": "ylim", "plotly": "layout.yaxis.range"},
    "xaxis.showline": {"mpl": "bottom", "plotly": "xaxis.showline"},
    "yaxis.showline": {"mpl": "left", "plotly": "yaxis.showline"},
    # Legend styling
    "legend.position": {
        "mpl": "loc",
        "plotly": "layout.legend.xanchor",
        "values": {
            "upper_right": {"mpl": "upper right", "plotly": "right"},
            "upper_left": {"mpl": "upper left", "plotly": "left"},
            "lower_right": {"mpl": "lower right", "plotly": "right"},
            "lower_left": {"mpl": "lower left", "plotly": "left"},
        },
    },
    "legend.show": {"mpl": "legend", "plotly": "showlegend"},
    # Title styling
    "title.text": {"mpl": "title", "plotly": "layout.title.text"},
    "title.font_size": {"mpl": "titlefontsize", "plotly": "layout.title.font.size"},
    # Rectangle/patch specific (for BED tracks)
    "rect.height": {"mpl": "height", "plotly": "height"},
    "rect.padding": {"mpl": "padding", "plotly": "padding"},
    # Plot-specific options
    "plot.type": {"mpl": "plot_type", "plotly": "plot_type"},
    "axis.type": {"mpl": "axis_type", "plotly": "axis_type"},
    "axis.show_chromosome": {"mpl": "show_chromosome", "plotly": "show_chromosome"},
    # Y-axis limits for tracks
    "ymin": {"mpl": "ymin", "plotly": "ymin"},
    "ymax": {"mpl": "ymax", "plotly": "ymax"},
    # Color column usage
    "use_color_column": {"mpl": "use_color_column", "plotly": "use_color_column"},
    "use_global_max": {"mpl": "use_global_max", "plotly": "use_global_max"},
    # Label alignment
    "label.alignment": {"mpl": "label_alignment", "plotly": "label_alignment"},
}


# Keys to ignore during translation (internal use only)
INTERNAL_KEYS = {
    "plot.type",
    "axis.type",
    "axis.show_chromosome",
    "rect.height",
    "rect.padding",
    "ymin",
    "ymax",
    "use_color_column",
    "use_global_max",
    "label.alignment",
    "fill.enabled",  # Internal flag to enable/disable fill
}


# =============================================================================
# TRANSLATION FUNCTIONS
# =============================================================================


def translate(
    unified_options: Dict[str, Any], target: str
) -> Dict[str, Dict[str, Any]]:
    """
    Translate unified options to backend-specific format.

    Args:
        unified_options: Dictionary with unified option keys
        target: Target backend ("mpl" or "plotly")

    Returns:
        Dictionary grouped by element type (marker, line, fill, etc.)

    Examples:
        >>> translate({"line.width": 2, "marker.color": "red"}, "mpl")
        {"line": {"linewidth": 2}, "marker": {"c": "red"}}
    """
    if target not in ["mpl", "plotly"]:
        raise ValueError(f"Target must be 'mpl' or 'plotly', got: {target}")

    result: Dict[str, Dict[str, Any]] = {}

    # Flatten input options to handle nested dictionaries
    flat_options = flatten(unified_options)

    for unified_key, value in flat_options.items():
        # Skip internal keys that aren't translated
        if unified_key in INTERNAL_KEYS:
            continue

        if unified_key not in UNIFIED_STYLE_MAP:
            # Unknown key - pass through as-is
            element_type = unified_key.split(".")[0] if "." in unified_key else "other"
            if element_type not in result:
                result[element_type] = {}
            result[element_type][unified_key] = value
            continue

        mapping = UNIFIED_STYLE_MAP[unified_key]

        # Get target key
        if target not in mapping:
            continue  # Skip if target not supported

        target_key = mapping[target]

        # Apply value transformation if needed
        if "values" in mapping and value in mapping["values"]:
            target_value = mapping["values"][value][target]
        elif "transform" in mapping:
            transform_key = f"{target}→{target}"  # No cross-conversion for now
            if transform_key in mapping["transform"]:
                target_value = mapping["transform"][transform_key](value)
            else:
                target_value = value
        else:
            target_value = value

        # Group by element type
        element_type = unified_key.split(".")[0]
        if element_type not in result:
            result[element_type] = {}

        # Handle nested target keys (like layout.xaxis.title.text)
        if "." in target_key:
            _set_nested_dict(result[element_type], target_key, target_value)
        else:
            result[element_type][target_key] = target_value

    return result


def _set_nested_dict(target_dict: Dict[str, Any], key_path: str, value: Any) -> None:
    """Set a value in a nested dictionary using dot notation."""
    keys = key_path.split(".")
    current = target_dict

    for key in keys[:-1]:
        if key not in current:
            current[key] = {}
        current = current[key]

    current[keys[-1]] = value


def get_internal_options(unified_options: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract internal options that are not translated to backends.

    Args:
        unified_options: Dictionary with unified option keys

    Returns:
        Dictionary with internal options only
    """
    flat_options = flatten(unified_options)
    return {k: v for k, v in flat_options.items() if k in INTERNAL_KEYS}


def normalize_options(options: Dict[str, Any]) -> Dict[str, Any]:
    """
    Normalize option keys and values.

    Args:
        options: Raw options dictionary

    Returns:
        Normalized options
    """
    normalized = {}

    for key, value in options.items():
        # Normalize key format (lowercase, dots for separators)
        norm_key = key.lower().replace("_", ".")

        # Normalize common values
        if isinstance(value, str):
            norm_value = value.lower().strip()

            # Normalize boolean strings
            if norm_value in ["true", "yes", "1"]:
                norm_value = True
            elif norm_value in ["false", "no", "0"]:
                norm_value = False
            else:
                norm_value = value  # Keep original case for non-boolean strings
        else:
            norm_value = value

        normalized[norm_key] = norm_value

    return normalized


def merge_options(*option_dicts: Dict[str, Any]) -> Dict[str, Any]:
    """
    Merge multiple option dictionaries with later ones taking precedence.

    Args:
        *option_dicts: Variable number of option dictionaries

    Returns:
        Merged options dictionary
    """
    merged = {}

    for options in option_dicts:
        if options:
            merged.update(options)

    return merged


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================


def create_line_options(
    color: str = "blue", width: float = 1.0, style: str = "solid"
) -> Dict[str, Any]:
    """Create standard line options."""
    return {"line.color": color, "line.width": width, "line.dash": style}


def create_marker_options(
    color: str = "blue", size: float = 6.0, symbol: str = "circle"
) -> Dict[str, Any]:
    """Create standard marker options."""
    return {"marker.color": color, "marker.size": size, "marker.symbol": symbol}


def create_fill_options(color: str = "blue", opacity: float = 0.7) -> Dict[str, Any]:
    """Create standard fill options."""
    return {"fill.color": color, "fill.opacity": opacity}


def create_text_options(
    font_size: float = 12.0, font_color: str = "black", font_family: str = "Arial"
) -> Dict[str, Any]:
    """Create standard text options."""
    return {
        "text.font_size": font_size,
        "text.font_color": font_color,
        "text.font_family": font_family,
    }


# =============================================================================
# PRESET OPTION COLLECTIONS
# =============================================================================

PRESET_STYLES = {
    "default": {
        "line.color": "blue",
        "line.width": 1.0,
        "marker.color": "blue",
        "marker.size": 6.0,
        "fill.color": "lightblue",
        "fill.opacity": 0.7,
    },
    "publication": {
        "line.color": "black",
        "line.width": 1.5,
        "marker.color": "black",
        "marker.size": 4.0,
        "fill.color": "gray",
        "fill.opacity": 0.5,
        "text.font_family": "Arial",
        "text.font_size": 10.0,
    },
    "vibrant": {
        "line.color": "red",
        "line.width": 2.0,
        "marker.color": "red",
        "marker.size": 8.0,
        "fill.color": "orange",
        "fill.opacity": 0.8,
    },
}


def get_preset_style(style_name: str) -> Dict[str, Any]:
    """
    Get a preset style configuration.

    Args:
        style_name: Name of the preset style

    Returns:
        Style options dictionary

    Raises:
        ValueError: If style name is not recognized
    """
    if style_name not in PRESET_STYLES:
        available = list(PRESET_STYLES.keys())
        raise ValueError(f"Unknown style '{style_name}'. Available: {available}")

    return PRESET_STYLES[style_name].copy()
