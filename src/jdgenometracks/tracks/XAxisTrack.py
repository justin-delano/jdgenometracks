from dataclasses import dataclass, field
from typing import Optional

import pandas as pd
import plotly.graph_objects as go
from matplotlib.axes import Axes
from matplotlib.ticker import FixedLocator, Formatter

from jdgenometracks.genomic_ticks import GenomicTickCalculator, format_genomic_ticks
from jdgenometracks.option_mapping import translate

from .GenomeTrack import GenomeTrack


class BPFormatter(Formatter):
    """
    Unified intelligent formatter for genomic coordinates with adaptive tick placement.

    This formatter uses the GenomicTickCalculator to determine optimal tick
    locations and formatting based on the visible genomic range. It automatically
    adapts between bp, Kb, and Mb units for optimal readability.
    """

    def __init__(self):
        self.tick_locations = []
        self.tick_labels = []
        self.unit_label = "b"
        self.scale_factor = 1

    def __call__(self, x, pos=None):
        """Return the formatted tick value."""
        # Find the closest tick location to format it properly
        if not self.tick_labels:
            return f"{int(x):,}"

        # Find the index of this tick
        tick_idx = None
        for i, loc in enumerate(self.tick_locations):
            if abs(loc - x) < 0.5:  # Close enough
                tick_idx = i
                break

        if tick_idx is not None and tick_idx < len(self.tick_labels):
            return self.tick_labels[tick_idx]
        else:
            # Fallback formatting
            scaled = x / self.scale_factor
            if self.scale_factor == 1:
                return f"{int(x):,}"
            elif scaled == int(scaled):
                return f"{int(scaled):,}"
            else:
                return f"{scaled:.1f}"

    def set_range(self, vmin, vmax):
        """Set the visible range and calculate optimal ticks."""
        # Calculate optimal ticks for this range
        self.tick_locations, self.unit_label, self.scale_factor = (
            GenomicTickCalculator.calculate_optimal_ticks(int(vmin), int(vmax))
        )

        # Format the tick labels
        self.tick_labels = format_genomic_ticks(
            self.tick_locations, self.unit_label, self.scale_factor
        )

        # Set the locations for matplotlib
        self.locs = self.tick_locations


@dataclass
class XAxisTrack(GenomeTrack):
    """
    Track class for displaying genomic X-axis with intelligent tick formatting.

    Uses a unified BPFormatter that automatically adapts tick placement and
    formatting based on the genomic range being displayed, supporting both
    matplotlib and Plotly backends.

    Attributes:
        axis_type (str): Specifies the type of axis (default is 'verbose').
        verbose_label (bool): Whether to show a verbose label (default is True).
        font_size (int): Font size for the axis label (default is 12).
        height_prop (float): Proportion of the height of the axis.
    """

    track_options: dict = field(default_factory=dict)

    def __post_init__(self):
        """Initialize an empty DataFrame for this track."""
        self.data = pd.DataFrame()

    def add_verbose_axis_mpl(self, ax: Axes, chromosome: str):
        """
        Add a verbose X-axis using matplotlib with intelligent tick formatting.

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

        # Use the unified intelligent BPFormatter
        formatter = BPFormatter()
        # Set the range to calculate optimal ticks
        vmin, vmax = ax.get_xlim()
        formatter.set_range(vmin, vmax)

        # Use fixed locator with calculated tick positions
        ax.xaxis.set_major_locator(FixedLocator(formatter.tick_locations))
        ax.xaxis.set_major_formatter(formatter)

        mpl_opts = translate(self.track_options, target="mpl")
        text_opts = mpl_opts.get("text", {})
        ax.xaxis.set_tick_params(bottom=True, labelbottom=True, **text_opts)
        ax.set_xlabel(
            chromosome if self.track_options.get("axis.show_chromosome", True) else "",
            **text_opts,
        )

    def plot_mpl(
        self,
        ax: Axes,
        chromosome: str,
        xmin: Optional[int] = None,
        xmax: Optional[int] = None,
        **kwargs,
    ):
        """
        Plot the X-axis using matplotlib with unified intelligent tick formatting.

        Args:
            ax (Axes): The matplotlib axis to plot on.
            chromosome (str): The chromosome label.
            xmin (Optional[int]): Minimum x-coordinate for intelligent tick placement.
            xmax (Optional[int]): Maximum x-coordinate for intelligent tick placement.
        """
        if self.track_options.get("axis.type") == "verbose":
            self.add_verbose_axis_mpl(ax, chromosome)

            # Apply unified intelligent tick placement if range is provided
            if xmin is not None and xmax is not None:
                formatter = BPFormatter()
                formatter.set_range(xmin, xmax)

                ax.xaxis.set_major_locator(FixedLocator(formatter.tick_locations))
                ax.xaxis.set_major_formatter(formatter)
        else:
            raise NotImplementedError(
                f"Axis type '{self.track_options.get('axis.type')}' is not implemented."
            )

    def add_verbose_axis_plotly(
        self,
        fig: go.Figure,
        row: int,
        col: int,
        chromosome: str,
        xmin: Optional[int] = None,
        xmax: Optional[int] = None,
    ):
        """
        Add a verbose X-axis using Plotly with intelligent tick formatting.

        Args:
            fig (go.Figure): The Plotly figure to add the X-axis to.
            row (int): The row number in the subplot grid.
            col (int): The column number in the subplot grid.
            chromosome (str): The chromosome label to display.
            xmin (int): Minimum x-coordinate for the axis range.
            xmax (int): Maximum x-coordinate for the axis range.
        """
        plotly_opts = translate(self.track_options, target="plotly")
        line_opts = plotly_opts.get("line", {})
        xaxis_opts = plotly_opts.get("xaxis", {})
        yaxis_opts = plotly_opts.get("yaxis", {})

        # Add empty trace
        fig.add_trace(go.Scatter(x=[], y=[], showlegend=False), row=row, col=col)

        # Configure base axis options
        axis_config = {
            "title_text": (
                chromosome
                if self.track_options.get("axis.show_chromosome", True)
                else ""
            ),
            "showticklabels": True,
            "tickangle": 0,
            **xaxis_opts,
            **line_opts,
        }

        # Apply intelligent tick placement using GenomicTickCalculator
        if xmin is not None and xmax is not None:
            tick_locations, unit_label, scale_factor = (
                GenomicTickCalculator.calculate_optimal_ticks(xmin, xmax)
            )

            # Format tick labels with unit suffix only on the last tick
            from jdgenometracks.genomic_ticks import format_genomic_ticks

            tick_labels = format_genomic_ticks(
                tick_locations, unit_label, scale_factor, show_unit_on_last=True
            )

            # Use explicit tick positioning with custom labels
            axis_config.update(
                {
                    "tickmode": "array",
                    "tickvals": tick_locations,
                    "ticktext": tick_labels,
                    "ticks": "outside",  # Show tick marks
                    "ticklen": 5,  # Length of tick marks
                    "tickwidth": 1,  # Width of tick marks
                    "tickcolor": "black",  # Color of tick marks
                }
            )

        fig.update_xaxes(row=row, col=col, **axis_config)
        fig.update_yaxes(showticklabels=False, row=row, col=col, **yaxis_opts)

    def plot_plotly(
        self,
        fig: go.Figure,
        row: int,
        col: int,
        chromosome: str,
        xmin: Optional[int] = None,
        xmax: Optional[int] = None,
        **kwargs,
    ):
        """
        Plot the X-axis using Plotly with unified intelligent tick formatting.

        Args:
            fig (go.Figure): The Plotly figure to plot on.
            row (int): The row number in the subplot grid.
            col (int): The column number in the subplot grid.
            chromosome (str): The chromosome label.
            xmin (Optional[int]): Minimum x-coordinate for intelligent tick placement.
            xmax (Optional[int]): Maximum x-coordinate for intelligent tick placement.
        """
        if self.track_options.get("axis.type") == "verbose":
            self.add_verbose_axis_plotly(fig, row, col, chromosome, xmin, xmax)
        else:
            raise NotImplementedError(
                f"Axis type '{self.track_options.get('axis.type')}' is not implemented."
            )
