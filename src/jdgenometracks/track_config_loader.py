import json
from typing import Any, Dict

from jdgenometracks import TrackFactory
from jdgenometracks.utils_units import (
    convert_to_inches,
    convert_to_pixels,
    parse_size_with_units,
)


def load_plot_config(json_path: str) -> Dict[str, Any]:
    """Load the high-level plot config from a JSON file."""
    try:
        with open(json_path, "r") as f:
            config = json.load(f)
    except Exception as e:
        raise RuntimeError(f"Failed to load plot config from {json_path}: {e}")
    if not isinstance(config, dict):
        raise ValueError("Plot config must be a JSON object (dict at top level).")
    if "tracks" not in config or not isinstance(config["tracks"], list):
        raise ValueError("Plot config must contain a 'tracks' list.")
    return config


def create_tracks_from_config(config: Dict[str, Any]):
    tracks = []
    for i, cfg in enumerate(config["tracks"]):
        track_kwargs = dict(cfg)

        if "subplot_x" not in cfg:
            track_kwargs["subplot_x"] = 0  # default subplot_x to 0 if not specified
        if "subplot_y" not in cfg:
            track_kwargs["subplot_y"] = i
        try:
            tracks.append(TrackFactory.create_track(**track_kwargs))
        except Exception as e:
            raise RuntimeError(f"Error creating track #{i+1}: {e}")
    if not tracks:
        raise ValueError("No tracks were created from the config.")

    return config.get("figure_options", {}), tracks


def plot_from_config(json_path: str, backend: str = "plotly") -> Any:
    """Load a plot config and render the plot using the specified backend and options."""
    config = load_plot_config(json_path)
    fig_options, tracks = create_tracks_from_config(config)

    # Parse units for total_height and total_width (default: px)
    def _parse_fig_size(fig_options, backend):
        for key, conv in [
            (
                "total_height",
                convert_to_pixels if backend == "plotly" else convert_to_inches,
            ),
            (
                "total_width",
                convert_to_pixels if backend == "plotly" else convert_to_inches,
            ),
        ]:
            if fig_options.get(key):
                val, unit = parse_size_with_units(fig_options[key], default_unit="px")
                fig_options[key] = conv(val, unit)

    if backend == "plotly":
        from jdgenometracks import PlotlyPlotter

        _parse_fig_size(fig_options, "plotly")
        plotter = PlotlyPlotter(tracks)
        try:
            fig = plotter.plot_all_tracks(
                fig_options=fig_options,
            )
        except Exception as e:
            raise RuntimeError(f"Error plotting with Plotly: {e}")
        return fig
    elif backend == "matplotlib":
        from jdgenometracks import MPLPlotter

        _parse_fig_size(fig_options, "matplotlib")
        plotter = MPLPlotter(tracks)
        try:
            fig, axes = plotter.plot_all_tracks(
                fig_options=fig_options,
            )
        except Exception as e:
            raise RuntimeError(f"Error plotting with Matplotlib: {e}")
        return fig, axes
    else:
        raise ValueError(f"Unknown backend: {backend}")
