import re


def parse_size_with_units(val, default_unit="px"):
    """
    Parse a size string with optional units (e.g., '800px', '10in').
    Returns (value, unit). If no unit, uses default_unit.
    """
    if isinstance(val, (int, float)):
        return float(val), default_unit
    if not isinstance(val, str):
        raise ValueError(f"Invalid size value: {val}")
    m = re.match(r"([0-9.]+)\s*([a-zA-Z]*)", val)
    if not m:
        raise ValueError(f"Invalid size string: {val}")
    value = float(m.group(1))
    unit = m.group(2) or default_unit
    return value, unit


def convert_to_inches(val, unit):
    """Convert a value with unit to inches (for matplotlib)."""
    if unit == "in":
        return val
    elif unit == "px":
        return val / 96.0  # 96 px per inch is a common default
    elif unit == "cm":
        return val / 2.54
    elif unit == "mm":
        return val / 25.4
    else:
        raise ValueError(f"Unsupported unit for inches conversion: {unit}")


def convert_to_pixels(val, unit):
    """Convert a value with unit to pixels (for plotly)."""
    if unit == "px":
        return val
    elif unit == "in":
        return val * 96.0
    elif unit == "cm":
        return val * 96.0 / 2.54
    elif unit == "mm":
        return val * 96.0 / 25.4
    else:
        raise ValueError(f"Unsupported unit for pixel conversion: {unit}")
