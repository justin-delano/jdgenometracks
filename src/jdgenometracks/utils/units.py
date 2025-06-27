"""
Unit conversion utilities for jdgenometracks V3.

Handles conversion between different measurement units (pixels, inches, cm, mm, points)
for both matplotlib and plotly backends.
"""

import re
from typing import Tuple

from .constants import ErrorMessages, UnitConstants


def parse_size_with_units(val, default_unit: str = "px") -> Tuple[float, str]:
    """
    Parse a size string with optional units (e.g., '800px', '10in').

    Args:
        val: Size value as string, int, or float
        default_unit: Unit to use if none specified

    Returns:
        Tuple of (value, unit)

    Raises:
        ValueError: If value cannot be parsed

    Examples:
        >>> parse_size_with_units("800px")
        (800.0, 'px')
        >>> parse_size_with_units(10.5)
        (10.5, 'px')
        >>> parse_size_with_units("5cm")
        (5.0, 'cm')
    """
    if isinstance(val, (int, float)):
        return float(val), default_unit

    if not isinstance(val, str):
        raise ValueError(ErrorMessages.INVALID_SIZE_VALUE.format(value=val))

    # Parse numeric value and optional unit
    match = re.match(r"([0-9.]+)\s*([a-zA-Z]*)", val)
    if not match:
        raise ValueError(ErrorMessages.INVALID_SIZE_STRING.format(string=val))

    value = float(match.group(1))
    unit = match.group(2) or default_unit
    return value, unit


def convert_to_inches(val: float, unit: str) -> float:
    """
    Convert a value with unit to inches (for matplotlib).

    Args:
        val: Numeric value
        unit: Unit string ('in', 'px', 'cm', 'mm', 'pt')

    Returns:
        Value converted to inches

    Raises:
        ValueError: If unit is not supported
    """
    if unit == "in":
        return val
    elif unit == "px":
        return val / UnitConstants.PIXELS_PER_INCH
    elif unit == "cm":
        return val / UnitConstants.CM_PER_INCH
    elif unit == "mm":
        return val / UnitConstants.MM_PER_INCH
    elif unit == "pt":
        return val / UnitConstants.POINTS_PER_INCH
    else:
        raise ValueError(ErrorMessages.UNSUPPORTED_UNIT_INCHES.format(unit=unit))


def convert_to_pixels(val: float, unit: str) -> float:
    """
    Convert a value with unit to pixels (for plotly).

    Args:
        val: Numeric value
        unit: Unit string ('px', 'in', 'cm', 'mm', 'pt')

    Returns:
        Value converted to pixels

    Raises:
        ValueError: If unit is not supported
    """
    if unit == "px":
        return val
    elif unit == "in":
        return val * UnitConstants.PIXELS_PER_INCH
    elif unit == "cm":
        return val * UnitConstants.PIXELS_PER_INCH / UnitConstants.CM_PER_INCH
    elif unit == "mm":
        return val * UnitConstants.PIXELS_PER_INCH / UnitConstants.MM_PER_INCH
    elif unit == "pt":
        return val * UnitConstants.PIXELS_PER_INCH / UnitConstants.POINTS_PER_INCH
    else:
        raise ValueError(ErrorMessages.UNSUPPORTED_UNIT_PIXELS.format(unit=unit))


def normalize_size_for_backend(size_str, backend: str) -> float:
    """
    Normalize a size string to the appropriate units for a backend.

    Args:
        size_str: Size string with units (e.g., "800px", "10in")
        backend: Target backend ("matplotlib" or "plotly")

    Returns:
        Normalized size value

    Examples:
        >>> normalize_size_for_backend("800px", "plotly")
        800.0
        >>> normalize_size_for_backend("8in", "matplotlib")
        8.0
    """
    value, unit = parse_size_with_units(size_str)

    if backend.lower() == "matplotlib":
        return convert_to_inches(value, unit)
    elif backend.lower() == "plotly":
        return convert_to_pixels(value, unit)
    else:
        raise ValueError(f"Unsupported backend: {backend}")


# Common size presets for convenience
SIZE_PRESETS = {
    "small": {"matplotlib": 6, "plotly": 400},
    "medium": {"matplotlib": 8, "plotly": 600},
    "large": {"matplotlib": 12, "plotly": 800},
    "xlarge": {"matplotlib": 16, "plotly": 1200},
}


def get_preset_size(preset: str, backend: str) -> float:
    """
    Get a preset size for a given backend.

    Args:
        preset: Size preset name ("small", "medium", "large", "xlarge")
        backend: Target backend ("matplotlib" or "plotly")

    Returns:
        Size value appropriate for the backend

    Raises:
        ValueError: If preset or backend is not recognized
    """
    if preset not in SIZE_PRESETS:
        raise ValueError(f"Unknown size preset: {preset}")

    backend_key = backend.lower()
    if backend_key not in SIZE_PRESETS[preset]:
        raise ValueError(f"Unsupported backend: {backend}")

    return SIZE_PRESETS[preset][backend_key]
