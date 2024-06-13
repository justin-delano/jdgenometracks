import math
from dataclasses import dataclass

import numpy as np
import plotly.graph_objects as go
from matplotlib.axes import Axes
from matplotlib.ticker import Formatter

from .GenomeTrack import GenomeTrack


class BPFormatter(Formatter):
    """
    Format tick values as pretty numbers and add as offset the unit

    The units are "b", "Kb", "Mb"
    The choice is made based on distance between extreme visible locs

    """

    def __init__(self):
        self.format = ""
        self.exponent = 0
        self.unit = ""

    def __call__(self, x, pos=None):
        """
        Return the format for tick value *x* at position *pos*.
        """
        xp = (x) / (10.0**self.exponent)
        if abs(xp) < 1e-8:
            xp = 0
        if len(self.locs) < 2 or x == self.locs[-2]:
            return self.format.format(xp) + " " + self.unit
        else:
            return self.format.format(xp)

    def set_locs(self, locs):
        # docstring inherited
        self.locs = locs
        if len(self.locs) > 0:
            self._set_unit()
            self._set_format()

    def _set_unit(self):
        # restrict to visible ticks
        if self.axis is None:
            return
        vmin, vmax = sorted(self.axis.get_view_interval())
        locs = np.asarray(self.locs)
        locs = locs[(vmin <= locs) & (locs <= vmax)]
        locs = np.abs(locs)
        if np.abs(locs[-1] - locs[0]) <= 1e3:
            self.exponent = 0
            self.unit = "b"
        elif np.abs(locs[-1] - locs[0]) <= 7e5:
            self.exponent = 3
            self.unit = "Kb"
        else:
            self.exponent = 6
            self.unit = "Mb"

    # This is adapted from ScalarFormatter
    def _set_format(self):
        # set the format string to format all the ticklabels
        if len(self.locs) < 2:
            # Temporarily augment the locations with the axis end points.
            if self.axis is not None:
                _locs = [*self.locs, *self.axis.get_view_interval()]
            else:
                _locs = self.locs
        else:
            _locs = self.locs
        locs = np.asarray(_locs) / 10.0**self.exponent
        loc_range = np.ptp(locs)
        # Curvilinear coordinates can yield two identical points.
        if loc_range == 0:
            loc_range = np.max(np.abs(locs))
        # Both points might be zero.
        if loc_range == 0:
            loc_range = 1
        if len(self.locs) < 2:
            # We needed the end points only for the loc_range calculation.
            locs = locs[:-2]
        loc_range_oom = int(math.floor(math.log10(loc_range)))
        # first estimate:
        sigfigs = max(0, 3 - loc_range_oom)
        # refined estimate:
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
    axis_type: str = "verbose"
    verbose_label: bool = True
    font_size: int = 12
    height_prop: float = 0.01

    def add_verbose_axis_mpl(self, ax: Axes, **kwargs):

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(True)
        ax.yaxis.set_tick_params(left=False, right=False, labelleft=False)

        ax.xaxis.set_visible(True)
        ax.xaxis.set_major_formatter(BPFormatter())
        ax.xaxis.set_tick_params(
            bottom=True, labelsize=self.axis_font_size, labelbottom=True
        )

        if self.verbose_label:
            ax.set_xlabel(kwargs["chromosome"], fontsize=self.font_size)
        else:
            ax.set_xlabel("")

    def plot_mpl(self, ax, **kwargs):
        if self.axis_type == "verbose":
            return self.add_verbose_axis_mpl(ax, **kwargs)
        else:
            raise NotImplementedError

    def add_verbose_axis_plotly(self, fig: go.Figure, row: int, col: int, **kwargs):
        fig.add_trace(go.Scatter(x=[], y=[], showlegend=False), row=row, col=col)

        if self.verbose_label:
            fig.update_xaxes(
                title_text=f"{kwargs['chromosome']}",
                linecolor="black",
                showticklabels=True,
                tickangle=0,
                row=row,
                col=col,
            )
        else:
            fig.update_xaxes(
                title_text="",
                linecolor="black",
                showticklabels=True,
                tickangle=0,
                row=row,
                col=col,
            )

        fig.update_yaxes(showticklabels=False, row=row, col=col)

    def plot_plotly(self, fig, row, col, **kwargs):
        if self.axis_type == "verbose":
            self.add_verbose_axis_plotly(fig, row, col, **kwargs)
        else:
            raise NotImplementedError


@dataclass
class SpacerTrack(GenomeTrack):
    height_prop: float = 0.1

    def plot_mpl(self, ax: Axes, **kwargs):
        ax.set_visible(False)

    def plot_plotly(self, fig: go.Figure, row: int, col: int, **kwargs):
        fig.add_trace(go.Scatter(x=[], y=[], showlegend=False), row=row, col=col)
        fig.update_xaxes(visible=False, row=row, col=col)
        fig.update_yaxes(visible=False, row=row, col=col)
