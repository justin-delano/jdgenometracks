from dataclasses import dataclass, field

import pandas as pd
import plotly.graph_objects as go
from matplotlib.axes import Axes

from jdgenometracks.option_mapping import translate

from .GenomeTrack import GenomeTrack


@dataclass
class SpacerTrack(GenomeTrack):
    """
    Track class for adding empty space between plots.
    """

    track_options: dict = field(default_factory=dict)

    def __post_init__(self):
        """Initialize an empty DataFrame for this track."""
        self.data = pd.DataFrame()

    def plot_mpl(self, ax: Axes):
        """
        Plot an empty space using matplotlib by hiding the axis.

        Args:
            ax (Axes): The matplotlib axis to hide.
        """
        # Use unified options (if any)
        translate(self.track_options, target="mpl")
        # For a spacer, just hide the axis (no options needed)
        ax.set_visible(False)

    def plot_plotly(self, fig: go.Figure, row: int, col: int):
        """
        Plot an empty space using Plotly by hiding the axes.

        Args:
            fig (go.Figure): The Plotly figure to add the empty space to.
            row (int): The row number in the subplot grid.
            col (int): The column number in the subplot grid.
        """
        translate(self.track_options, target="plotly")
        fig.add_trace(go.Scatter(x=[], y=[], showlegend=False), row=row, col=col)
        fig.update_xaxes(visible=False, row=row, col=col)
        fig.update_yaxes(visible=False, row=row, col=col)
