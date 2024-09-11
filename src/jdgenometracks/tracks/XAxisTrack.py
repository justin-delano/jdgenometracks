import math
from dataclasses import dataclass

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from matplotlib.axes import Axes
from matplotlib.ticker import Formatter

from .GenomeTrack import GenomeTrack


class BPFormatter(Formatter):
    """
    Custom formatter for formatting genomic axis tick values in base pairs (b), kilobases (Kb), or megabases (Mb).
    The choice of unit depends on the distance between the visible ticks.
    """

    def __init__(self):
        self.format = ""
        self.exponent = 0
        self.unit = ""

    def __call__(self, x, pos=None):
        """
        Return the formatted tick value *x* at position *pos*.

        Args:
            x (float): The tick value.
            pos (int, optional): Tick position (unused).

        Returns:
            str: The formatted tick value with units.
        """
        xp = x / (10.0**self.exponent)
        if abs(xp) < 1e-8:
            xp = 0
        return (
            self.format.format(xp) + " " + self.unit
            if x == self.locs[-2]  # type: ignore
            else self.format.format(xp)
        )

    def set_locs(self, locs):
        """
        Set the locations of the ticks and update the formatting based on the visible range.

        Args:
            locs (list): List of tick locations.
        """
        self.locs = locs
        if len(self.locs) > 0:  # type: ignore
            self._set_unit()
            self._set_format()

    def _set_unit(self):
        """Set the appropriate unit (b, Kb, Mb) based on the tick locations."""
        if self.axis is None:
            return
        vmin, vmax = sorted(self.axis.get_view_interval())
        locs = np.asarray(self.locs)
        locs = locs[(vmin <= locs) & (locs <= vmax)]
        locs = np.abs(locs)

        range_diff = np.abs(locs[-1] - locs[0])
        if range_diff <= 1e3:
            self.exponent = 0
            self.unit = "b"
        elif range_diff <= 7e5:
            self.exponent = 3
            self.unit = "Kb"
        else:
            self.exponent = 6
            self.unit = "Mb"

    def _set_format(self):
        """Set the format string to format all tick labels."""
        if len(self.locs) < 2:  # type: ignore
            _locs = (
                [*self.locs, *self.axis.get_view_interval()]  # type: ignore
                if self.axis is not None
                else self.locs
            )
        else:
            _locs = self.locs

        locs = np.asarray(_locs) / 10.0**self.exponent
        loc_range = np.ptp(locs) or np.max(np.abs(locs)) or 1

        loc_range_oom = int(math.floor(math.log10(loc_range)))
        sigfigs = max(0, 3 - loc_range_oom)
        thresh = 1e-3 * 10**loc_range_oom

        while sigfigs >= 0:
            if np.abs(locs - np.round(locs, decimals=sigfigs)).max() < thresh:
                sigfigs -= 1
            else:
                break
        sigfigs += 1

        self.format = "{:,." + str(sigfigs) + "f}"


@dataclass
class XAxisTrack(GenomeTrack):
    """
    Track class for displaying the genomic X-axis.

    Attributes:
        axis_type (str): Specifies the type of axis (default is 'verbose').
        verbose_label (bool): Whether to show a verbose label (default is True).
        font_size (int): Font size for the axis label (default is 12).
        height_prop (float): Proportion of the height of the axis.
    """

    axis_type: str = "verbose"
    verbose_label: bool = True
    font_size: int = 12
    height_prop: float = 0.01

    def __post_init__(self):
        """Initialize an empty DataFrame for this track."""
        self.data = pd.DataFrame()

    def add_verbose_axis_mpl(self, ax: Axes, chromosome: str):
        """
        Add a verbose X-axis using matplotlib.

        Args:
            ax (Axes): The matplotlib axis to add the X-axis to.
            chromosome (str): The chromosome label to display.
        """
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(True)

        ax.yaxis.set_tick_params(left=False, right=False, labelleft=False)
        ax.xaxis.set_visible(True)
        ax.xaxis.set_major_formatter(BPFormatter())
        ax.xaxis.set_tick_params(
            bottom=True, labelsize=self.font_size, labelbottom=True
        )

        ax.set_xlabel(chromosome if self.verbose_label else "", fontsize=self.font_size)

    def plot_mpl(self, ax: Axes, chromosome: str, **kwargs):
        """
        Plot the X-axis using matplotlib.

        Args:
            ax (Axes): The matplotlib axis to plot on.
            chromosome (str): The chromosome label.
        """
        if self.axis_type == "verbose":
            self.add_verbose_axis_mpl(ax, chromosome)
        else:
            raise NotImplementedError(
                f"Axis type '{self.axis_type}' is not implemented."
            )

    def add_verbose_axis_plotly(
        self, fig: go.Figure, row: int, col: int, chromosome: str
    ):
        """
        Add a verbose X-axis using Plotly.

        Args:
            fig (go.Figure): The Plotly figure to add the X-axis to.
            row (int): The row number in the subplot grid.
            col (int): The column number in the subplot grid.
            chromosome (str): The chromosome label to display.
        """
        fig.add_trace(go.Scatter(x=[], y=[], showlegend=False), row=row, col=col)

        fig.update_xaxes(
            title_text=chromosome if self.verbose_label else "",
            linecolor="black",
            showticklabels=True,
            tickangle=0,
            row=row,
            col=col,
        )
        fig.update_yaxes(showticklabels=False, row=row, col=col)

    def plot_plotly(
        self, fig: go.Figure, row: int, col: int, chromosome: str, **kwargs
    ):
        """
        Plot the X-axis using Plotly.

        Args:
            fig (go.Figure): The Plotly figure to plot on.
            row (int): The row number in the subplot grid.
            col (int): The column number in the subplot grid.
            chromosome (str): The chromosome label.
        """
        if self.axis_type == "verbose":
            self.add_verbose_axis_plotly(fig, row, col, chromosome)
        else:
            raise NotImplementedError(
                f"Axis type '{self.axis_type}' is not implemented."
            )
