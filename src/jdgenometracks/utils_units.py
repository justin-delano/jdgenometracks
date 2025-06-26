import re

from .config import ErrorMessages, UnitConversion


def parse_size_with_units(val, default_unit="px"):
    """
    Parse a size string with optional units (e.g., '800px', '10in').
    Returns (value, unit). If no unit, uses default_unit.
    """
    if isinstance(val, (int, float)):
        return float(val), default_unit
    if not isinstance(val, str):
        raise ValueError(ErrorMessages.INVALID_SIZE_VALUE.format(value=val))
    m = re.match(r"([0-9.]+)\s*([a-zA-Z]*)", val)
    if not m:
        raise ValueError(ErrorMessages.INVALID_SIZE_STRING.format(string=val))
    value = float(m.group(1))
    unit = m.group(2) or default_unit
    return value, unit


def convert_to_inches(val, unit):
    """Convert a value with unit to inches (for matplotlib)."""
    if unit == "in":
        return val
    elif unit == "px":
        return (
            val / UnitConversion.PIXELS_PER_INCH
        )  # 96 px per inch is a common default
    elif unit == "cm":
        return val / UnitConversion.CM_PER_INCH
    elif unit == "mm":
        return val / UnitConversion.MM_PER_INCH
    else:
        raise ValueError(ErrorMessages.UNSUPPORTED_UNIT_INCHES.format(unit=unit))


def convert_to_pixels(val, unit):
    """Convert a value with unit to pixels (for plotly)."""
    if unit == "px":
        return val
    elif unit == "in":
        return val * UnitConversion.PIXELS_PER_INCH
    elif unit == "cm":
        return val * UnitConversion.PIXELS_PER_INCH / UnitConversion.CM_PER_INCH
    elif unit == "mm":
        return val * UnitConversion.PIXELS_PER_INCH / UnitConversion.MM_PER_INCH
    else:
        raise ValueError(ErrorMessages.UNSUPPORTED_UNIT_PIXELS.format(unit=unit))
