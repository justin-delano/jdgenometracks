"""unified_plot_mapping.py
-----------------------
A unified style mapping and translator between a backend‑agnostic plotting
vocabulary and the two most common Python back‑ends: **Matplotlib** and
**Plotly**.

Revision 2025‑06‑23 (numpy‑free)
================================
* **Dependency change** – removed the optional NumPy requirement. All unit
  conversions now use the standard library’s ``math`` module.
* Behaviour is otherwise identical: ``translate()`` still returns nested
  dictionaries grouped by element type (``'marker'``, ``'line'``, etc.).

Public symbols
--------------
* ``UNIFIED_STYLE_MAP`` – master mapping dict for every unified key.
* ``flatten``            – flatten nested unified specs into dotted‑key form.
* ``translate``          – convert a unified spec to Matplotlib **or** Plotly
                           kwargs, grouped by element type.
"""

from __future__ import annotations

import math  # NumPy removed – we only need basic math
from typing import Any, Callable, Dict, Mapping, Optional

# ---------------------------------------------------------------------------
# Constants & helpers
# ---------------------------------------------------------------------------

SEP: str = "."  # dotted‑key separator used throughout


def flatten(
    data: Mapping[str, Any], *, sep: str = SEP, parent_key: str = ""
) -> Dict[str, Any]:
    """Return *data* with nested dictionaries flattened using dotted keys.

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
    return transform(value) if transform else value


# ---------------------------------------------------------------------------
# Master unified‑style → backend mapping dictionary
# ---------------------------------------------------------------------------
UNIFIED_STYLE_MAP: Dict[str, Mapping[str, Any]] = {
    # ---------------------------------------------------------------------
    # Marker (scatter, bubble)
    # ---------------------------------------------------------------------
    "marker.size": {
        "mpl": "s",  # marker *area* in points^2
        "plotly": "marker.size",  # marker *radius* in px
        "transform": {
            "mpl→plotly": lambda area: math.sqrt(area),
            "plotly→mpl": lambda radius: radius**2,
        },
    },
    "marker.color": {"mpl": "c", "plotly": "marker.color"},
    "marker.symbol": {"mpl": "marker", "plotly": "marker.symbol"},
    "marker.edge_color": {"mpl": "edgecolors", "plotly": "marker.line.color"},
    "marker.edge_width": {"mpl": "linewidths", "plotly": "marker.line.width"},
    "marker.opacity": {"mpl": "alpha", "plotly": "marker.opacity"},
    # ---------------------------------------------------------------------
    # Line (line plots or outlines)
    # ---------------------------------------------------------------------
    "line.width": {"mpl": "linewidth", "plotly": "line.width"},
    "line.color": {"mpl": "color", "plotly": "line.color"},
    "line.dash": {
        "mpl": "linestyle",
        "plotly": "line.dash",
        "values": {
            "solid": {"mpl": "-", "plotly": "solid"},
            "dashed": {"mpl": "--", "plotly": "dash"},
            "dotted": {"mpl": ":", "plotly": "dot"},
            "dashdot": {"mpl": "-.", "plotly": "dashdot"},
        },
    },
    # ---------------------------------------------------------------------
    # Fill (area, polygons, bars)
    # ---------------------------------------------------------------------
    "fill.color": {"mpl": "facecolor", "plotly": "fillcolor"},
    "fill.opacity": {"mpl": "alpha", "plotly": "opacity"},
    "fill.edge_color": {"mpl": "edgecolor", "plotly": "line.color"},
    "fill.edge_width": {"mpl": "linewidth", "plotly": "line.width"},
    # ---------------------------------------------------------------------
    # Text (labels & annotations)
    # ---------------------------------------------------------------------
    "text.font_family": {"mpl": "fontfamily", "plotly": "textfont.family"},
    "text.font_size": {"mpl": "fontsize", "plotly": "textfont.size"},
    "text.font_color": {"mpl": "color", "plotly": "textfont.color"},
    "text.angle": {"mpl": "rotation", "plotly": "textangle"},
    "text.vertical_alignment": {
        "mpl": "va",
        "plotly": "textposition",
    },
    # ---------------------------------------------------------------------
    # Axes (common subset only – extend if needed)
    # ---------------------------------------------------------------------
    "xaxis.title": {"mpl": "xlabel", "plotly": "layout.xaxis.title.text"},
    "yaxis.title": {"mpl": "ylabel", "plotly": "layout.yaxis.title.text"},
    "xaxis.range": {"mpl": "xlim", "plotly": "layout.xaxis.range"},
    "yaxis.range": {"mpl": "ylim", "plotly": "layout.yaxis.range"},
    # ---------------------------------------------------------------------
    # Legend & Title (layout‑level properties)
    # ---------------------------------------------------------------------
    "legend.position": {
        "mpl": "loc",
        "plotly": "layout.legend.xanchor",  # simplified – refine as needed
    },
    "title.text": {"mpl": "title", "plotly": "layout.title.text"},
    "title.font_size": {"mpl": "titlefontsize", "plotly": "layout.title.font.size"},
}

# keys to ignore which are used for internal purposes only
UNIFIED_STYLE_MAP_IGNORE_KEYS = {
    "plot.type",  # used to determine plot type (lines, bars, etc.)
    "rect.height",  # used for bar plots in Plotly
    "rect.padding",  # used for bar plots in Plotly
    "label.alignment",  # used for bed track labels in matplotlib
    "axis.type",  # used to determine axis type (verbose, simple, etc.)
    "axis.show_chromosome",  # used to determine if chromosome label should be shown
    "ymin",  # used to set y-axis minimum in matplotlib
    "ymax",  # used to set y-axis maximum in matplotlib
    "use_global_max",  # how to fit bars for bed tracks
    "hlines",  # used for horizontal lines in matplotlib
    "plot_bgcolor",
    "margin",
    "show_gridlines",  # used to toggle gridlines in Plotly
    
}

# ---------------------------------------------------------------------------
# Public translator – always returns *nested* dictionaries
# ---------------------------------------------------------------------------


def translate(
    spec: Mapping[str, Any],
    *,
    target: str = "mpl",
    mapping: Mapping[str, Mapping[str, Any]] = UNIFIED_STYLE_MAP,
    ignore_keys: set = UNIFIED_STYLE_MAP_IGNORE_KEYS,
    sep: str = SEP,
) -> Dict[str, Dict[str, Any]]:
    """Translate a unified *spec* into backend kwargs grouped by element.

    Parameters
    ----------
    spec    : nested mapping using the unified vocabulary.
    target  : ``'mpl'`` or ``'plotly'``.
    mapping : (optional) override mapping dictionary.
    sep     : dotted‑key separator (usually '.').

    Returns
    -------
    Dict[str, Dict[str, Any]]
        Dictionary whose keys are the *first* segment of each unified key
        (for example ``'marker'`` or ``'line'``) and whose values are a dict of
        backend‑specific kwargs suitable for that element.
    """

    if target not in {"mpl", "plotly"}:
        raise ValueError("target must be 'mpl' or 'plotly'")

    flat_spec = flatten(spec, sep=sep)
    output: Dict[str, Dict[str, Any]] = {}
    arrow = "mpl→plotly" if target == "plotly" else "plotly→mpl"

    for ukey, value in flat_spec.items():
        # If the key should be ignored, pass it through unchanged
        if ukey in ignore_keys:
            output[ukey] = value
            continue

        entry = mapping.get(ukey)
        if not entry:
            raise KeyError(f"No mapping defined for '{ukey}'")

        target_key = entry.get(target)
        if target_key is None:
            # unified option unsupported in this backend – silently skip
            continue

        # unit conversions (area<->radius etc.)
        tform_map = entry.get("transform", {})
        value = _apply_transform(value, tform_map.get(arrow))

        # discrete enum remapping (dash styles etc.)
        if "values" in entry and isinstance(value, str):
            value_map = entry["values"].get(value)
            if value_map and (target in value_map):
                value = value_map[target]

        group, *_ = ukey.split(sep)
        group_dict = output.setdefault(group, {})
        group_key = target_key.split(sep)[-1] if sep in target_key else target_key
        group_dict[group_key] = value

    return output
