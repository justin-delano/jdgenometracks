"""
Spacer track implementation for adding empty space between tracks.

This module provides the SpacerTrack class for adding empty space between
other tracks in genomic visualizations. It creates invisible tracks that
can be used for layout spacing.

Author: Assistant
Date: 2024
"""

from dataclasses import dataclass
from typing import Optional

import pandas as pd
import plotly.graph_objects as go
from matplotlib.axes import Axes

from ..utils.options import translate
from .base_track import BaseTrack


@dataclass
class SpacerTrack(BaseTrack):
    """
    Track class for adding empty space between plots.

    This track creates invisible space that can be used to separate
    other tracks visually. It doesn't display any data but takes up
    space in the layout.

    Attributes:
        track_options: Dictionary containing styling and display options
    """

    # Spacer track constants
    TRACK_TYPE = "spacer"  # Type identifier for spacer tracks
    DEFAULT_TRACK_NAME = "spacer"  # Default name for spacer tracks

    def __post_init__(self):
        """Initialize an empty DataFrame for this track."""
        # Set track type
        self.track_type = self.TRACK_TYPE

        # Initialize with empty data
        self.data = pd.DataFrame()
        self._is_loaded = True

        # Call parent post_init (but skip validation since we don't need file_path)
        if self.track_name is None:
            self.track_name = self.DEFAULT_TRACK_NAME

    def load_data(self, file_path: Optional[str] = None) -> pd.DataFrame:
        """
        Load data for spacer track (always empty).

        Returns:
            Empty DataFrame
        """
        return pd.DataFrame()

    def validate_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Validate spacer data (always valid).

        Args:
            data: Data to validate

        Returns:
            The same data (always valid)
        """
        return data

    def plot_matplotlib(self, ax: Axes, start: int, end: int, **kwargs) -> None:
        """
        Plot an empty space using matplotlib by hiding the axis.

        Args:
            ax: The matplotlib axis to hide
            **kwargs: Additional plotting parameters (unused)
        """
        # Translate options (if any) for matplotlib
        translate(self.track_options, target="mpl")

        # For a spacer, just hide the axis completely
        ax.set_visible(False)

    def plot_plotly(
        self, fig: go.Figure, row: int, col: int, start: int, end: int, **kwargs
    ) -> None:
        """
        Create plotly traces for the spacer track (empty).

        Args:
            **kwargs: Additional plotting parameters (unused)

        Returns:
            List containing empty scatter trace
        """
        # Spacer tracks don't add any plotly traces
        # Just ensure the subplot exists but is empty
        pass
