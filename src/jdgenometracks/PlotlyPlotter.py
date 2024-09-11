from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import plotly.graph_objects as go
import plotly.subplots as ps
import scipy.sparse as sp

from .tracks.BedTrack import BedTrack
from .tracks.GenomeTrack import GenomeTrack
from .utils import TrackUtils


@dataclass
class PlotlyPlotter:
    """
    A class for generating multi-track genomic plots using Plotly.

    Attributes:
        tracks (np.ndarray): Array of tracks to plot, can contain instances of GenomeTrack or its subclasses.
        total_height (float): The total height of the figure in inches.
    """

    tracks: np.ndarray
    total_height: float
    total_width: float

    def __post_init__(self):
        """
        Converts tracks into a numpy array if not already one and performs initialization.
        """
        if not isinstance(self.tracks, np.ndarray):
            self.tracks = np.array(self.tracks)

    def plot_single_track(
        self, subplots: go.Figure, track: GenomeTrack, row: int, col: int, **kwargs
    ) -> np.ndarray | None:
        """
        Plots a single track on the given subplot.

        Parameters:
            subplots (go.Figure): The Plotly figure to plot on.
            track (GenomeTrack): The track to plot.
            row (int): The row index of the subplot to plot on.
            col (int): The column index of the subplot to plot on.
            **kwargs: Additional arguments to pass to the plot method.

        Returns:
            np.ndarray | None: The bed region coverage matrix if applicable, otherwise None.
        """

        if isinstance(track, BedTrack):
            bed_region_coverage = track.plot_plotly(subplots, row, col, **kwargs)
        else:
            track.plot_plotly(subplots, row, col, **kwargs)
            bed_region_coverage = None

        track.add_hlines_plotly(subplots, row, col)
        if track.show_yaxis_ticks:
            subplots.update_yaxes(
                row=row,
                col=col,
                **track.plotly_yaxis_options,
            )

        return bed_region_coverage

    def _validate_inputs(
        self,
        height_props: list[float],
        row_titles: list[str],
        width_props: list[float],
        column_titles: list[str],
    ) -> None:
        """
        Validates the inputs for the plot, ensuring the lengths match the expected dimensions.

        Args:
            height_props (list[float]): Proportions of the figure height for each row.
            row_titles (list[str]): Titles for each row in the figure.
            width_props (list[float]): Proportions of the figure width for each column.
            column_titles (list[str]): Titles for each column in the figure.

        Raises:
            AssertionError: If the lengths of inputs don't match the expected number of rows or columns.
        """
        num_distinct_rows = sum(
            [not track.share_with_previous for track in self.tracks[:, 0]]
        )

        assert (
            len(height_props) == num_distinct_rows
        ), f"Number of height_props should equal {num_distinct_rows}"
        assert (
            len(row_titles) == num_distinct_rows
        ), f"Number of row_titles should equal {num_distinct_rows}"
        assert (
            len(width_props) == self.tracks.shape[1]
        ), f"Number of width_props should equal {self.tracks.shape[1]}"
        assert (
            len(column_titles) == self.tracks.shape[1]
        ), f"Number of column_titles should equal {self.tracks.shape[1]}"

    def _initialize_subplots(
        self,
        height_props: list[float],
        width_props: list[float],
        row_titles: list[str],
        column_titles: list[str],
        num_distinct_rows: int,
        vertical_spacing: float = 0.02,
        horizontal_spacing: float = 0.05,
        shared_xaxes: str = "columns",
        shared_yaxes: str = "rows",
    ) -> go.Figure:
        """
        Initializes the Plotly subplots figure with the given layout properties.

        Args:
            height_props (list[float]): Heights of each row as a proportion of the total figure height.
            width_props (list[float]): Widths of each column as a proportion of the total figure width.
            row_titles (list[str]): Titles for each row.
            column_titles (list[str]): Titles for each column.
            num_distinct_rows (int): The number of distinct rows in the plot.
            vertical_spacing (float): Vertical spacing between rows.
            horizontal_spacing (float): Horizontal spacing between columns.
            shared_xaxes (str): Determines if x-axes are shared across columns (can be 'columns', 'all', or None).
            shared_yaxes (str): Determines if y-axes are shared across rows (can be 'rows', 'all', or None).

        Returns:
            go.Figure: A Plotly figure with subplots initialized.
        """
        subplots = ps.make_subplots(
            rows=num_distinct_rows,
            cols=self.tracks.shape[1],
            shared_xaxes=shared_xaxes,  # Customizable x-axis sharing # type: ignore
            shared_yaxes=shared_yaxes,  # Customizable y-axis sharing # type: ignore
            row_heights=height_props,
            row_titles=row_titles,
            column_widths=width_props,
            column_titles=column_titles,
            vertical_spacing=vertical_spacing,  # Customizable spacing
            horizontal_spacing=horizontal_spacing,  # Customizable spacing
        )
        return subplots

    def plot_all_tracks(
        self,
        column_regions: list[str | None] | None = None,
        total_height: float | None = None,
        total_width: float | None = None,
        height_props: list[float] | None = None,
        row_titles: list[str] | None = None,
        width_props: list[float] | None = None,
        column_titles: list[str] | None = None,
        shared_xaxes: str = "columns",
        shared_yaxes: str = "rows",
        margin: dict | None = None,
        vertical_spacing: float = 0.02,
        horizontal_spacing: float = 0.05,
        title_options: dict | None = None,
        font_options: dict | None = None,
        layout_options: dict | None = None,
        show_gridlines: bool = True,
        plot_bgcolor: str = "white",
        relative_x_axis: bool = False,
        show_fig: bool = False,
    ) -> go.Figure:
        """
        Plots all tracks into a single Plotly figure with optional customization.

        Args:
            column_regions (list[str | None] | None): Subset each column to a genomic region in format "chr:start-end".
            total_height (float | None): The total height of the figure.
            height_props (list[float] | None): Heights of each row as proportions of the total height.
            row_titles (list[str] | None): Titles for each row in the figure.
            width_props (list[float] | None): Widths of each column as proportions of the total width.
            column_titles (list[str] | None): Titles for each column in the figure.
            relative_x_axis (bool): If true, x-axes start at 0bp.
            show_fig (bool): If true, the figure is displayed.
            layout_options (dict | None): Additional options to customize Plotly layout (e.g., font, title, colors).
            margin (dict | None): Custom margins for the plot.
            vertical_spacing (float): Vertical spacing between subplots.
            horizontal_spacing (float): Horizontal spacing between subplots.
            title_options (dict | None): Options for customizing the title of the plot.
            font_options (dict | None): Font settings for the plot elements.
            x_axis_title (str | None): Custom title for the x-axis.
            y_axis_title (str | None): Custom title for the y-axis.
            show_gridlines (bool): Toggle gridlines on/off.
            plot_bgcolor (str): Background color for the entire plot.

        Returns:
            go.Figure: The complete Plotly figure with all tracks.
        """
        # Ensure the tracks array has 2 dimensions (rows and columns)
        if self.tracks.ndim == 1:
            self.tracks = self.tracks.reshape(-1, 1)

        # Default to no region subsetting for columns
        if column_regions is None:
            column_regions = [None] * self.tracks.shape[1]  # type: ignore
        assert column_regions is not None

        # Determine number of distinct rows (tracks that don't share axes with the previous row)
        num_distinct_rows = sum(
            not track.share_with_previous for track in self.tracks[:, 0]
        )

        # Set default properties for layout if not provided
        height_props = height_props or TrackUtils.get_height_props(self.tracks)
        row_titles = row_titles or [""] * num_distinct_rows
        width_props = width_props or [
            1 / self.tracks.shape[1] for _ in range(self.tracks.shape[1])
        ]
        column_titles = column_titles or ["" for _ in range(self.tracks.shape[1])]
        total_height = total_height or self.total_height
        total_width = total_width or self.total_width

        # Validations
        self._validate_inputs(height_props, row_titles, width_props, column_titles)

        # Initialize the subplots figure with custom spacings and shared axes settings
        subplots = self._initialize_subplots(
            height_props,
            width_props,
            row_titles,
            column_titles,
            num_distinct_rows,
            vertical_spacing=vertical_spacing,
            horizontal_spacing=horizontal_spacing,
            shared_xaxes=shared_xaxes,
            shared_yaxes=shared_yaxes,
        )

        # Set the plot title if provided
        if title_options:
            subplots.update_layout(**title_options)

        # Set layout properties for the entire figure
        subplots.update_layout(
            autosize=True,
            height=total_height,
            width=total_width,
            plot_bgcolor=plot_bgcolor,
            margin=margin or dict(l=0.1, r=0.1, t=50, b=20, pad=4),
            **(layout_options or {}),  # Allow additional layout options
        )

        # Plot each track
        for col_idx in range(self.tracks.shape[1]):
            xmin, xmax, extra_options = TrackUtils.get_col_limits(
                self.tracks[:, col_idx], column_regions[col_idx], relative_x_axis
            )
            bed_region_coverage = sp.csr_matrix((1, xmax - xmin))

            plot_row = 1  # Row counter
            for row_idx in range(self.tracks.shape[0]):
                track = self.tracks[row_idx, col_idx]

                if track is None:
                    continue

                if track.share_with_previous:
                    plot_row -= 1
                else:
                    bed_region_coverage = sp.csr_matrix((1, xmax - xmin))

                # Handle BedTrack separately for bed region coverage
                if isinstance(track, BedTrack):
                    extra_options["bed_region_coverage"] = bed_region_coverage
                    bed_region_coverage = self.plot_single_track(
                        subplots,
                        track,
                        plot_row,
                        col_idx + 1,
                        region=column_regions[col_idx],
                        **extra_options,
                    )
                else:
                    self.plot_single_track(
                        subplots,
                        track,
                        plot_row,
                        col_idx + 1,
                        region=column_regions[col_idx],
                        **extra_options,
                    )

                # Update the x-axes
                subplots.update_xaxes(
                    range=[xmin, xmax],
                    row=plot_row,
                    col=col_idx + 1,
                    showgrid=show_gridlines,  # Toggle gridlines
                )
                plot_row += 1

        # Final layout updates for fonts and axis titles
        if font_options:
            subplots.update_layout(**font_options)
            # subplots.update_annotations(**plotly_font_options)
            # subplots.update_xaxes(**plotly_font_options)

        # Optionally show the figure
        if show_fig:
            subplots.show()

        return subplots
