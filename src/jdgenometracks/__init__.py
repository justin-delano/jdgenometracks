"""
jdgenometracks

A Python package for generating genome browser tracks using Matplotlib and Plotly.

Usage Example:
--------------
from jdgenometracks import MPLPlotter, BedTrack

# Create tracks
tracks = [BedTrack(...), ...]

# Plot with Matplotlib
plotter = MPLPlotter(tracks=tracks, total_height=6, total_width=10)
fig = plotter.plot()  # or plotter.plot_single_track(...)

# Plot with Plotly
from jdgenometracks import PlotlyPlotter
plotter = PlotlyPlotter(tracks=tracks, total_height=6, total_width=10)
fig = plotter.plot()  # or plotter.plot_single_track(...)

See the README or demo notebooks for more details.
"""

from .config import VERSION

__version__ = VERSION
__author__ = "Justin Delano"

# High-level convenience function for quick plotting
from typing import Any, Sequence

from .MPLPlotter import MPLPlotter  # noqa: F401
from .PlotlyPlotter import PlotlyPlotter  # noqa: F401
from .tracks import BedGraphTrack, BedTrack, SpacerTrack, XAxisTrack  # noqa: F401
from .utils import TrackFactory, TrackUtils  # noqa: F401


def plot_tracks(tracks: Sequence[Any], backend: str = "matplotlib", **kwargs):
    """
    Quickly plot tracks using the specified backend ("matplotlib" or "plotly").

    Args:
        tracks (Sequence): List or array of track objects.
        backend (str): Which backend to use ("matplotlib" or "plotly").
        **kwargs: Additional arguments passed to the plotter.

    Returns:
        The figure object (matplotlib Figure or plotly Figure).
    """
    import numpy as np

    tracks_array = np.array(tracks)
    if backend == "matplotlib":
        from .MPLPlotter import MPLPlotter

        plotter = MPLPlotter(
            tracks=tracks_array,
            total_height=kwargs.get("total_height", 6),
            total_width=kwargs.get("total_width", 10),
        )
        return plotter.plot_all_tracks(**kwargs)
    elif backend == "plotly":
        from .PlotlyPlotter import PlotlyPlotter

        plotter = PlotlyPlotter(
            tracks=tracks_array,
            total_height=kwargs.get("total_height", 600),
            total_width=kwargs.get("total_width", 1000),
        )
        return plotter.plot_all_tracks(**kwargs)
    else:
        raise ValueError(f"Unknown backend: {backend}. Use 'matplotlib' or 'plotly'.")
