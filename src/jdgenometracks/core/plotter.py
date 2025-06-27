"""
Main plotting orchestrator for jdgenometracks V3.

This module provides high-level functions for creating genomic visualizations
using the backend architecture and track system. It orchestrates the creation
of figures, track layout, and rendering.

Author: Assistant
Date: 2024
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union

from ..backends import BaseBackend, create_backend
from ..tracks.base_track import BaseTrack
from ..utils.constants import MatplotlibConstants, PlotConstants
from ..utils.validation import validate_genomic_region


@dataclass
class PlotConfiguration:
    """
    Configuration for genomic plots.

    This dataclass encapsulates all the configuration options for creating
    genomic visualizations, including layout, styling, and track options.
    """

    # Plot configuration constants
    DEFAULT_FIGURE_HEIGHT = 8.0  # Default figure height in inches
    DEFAULT_FIGURE_WIDTH = 12.0  # Default figure width in inches
    MIN_TRACK_HEIGHT = 0.1  # Minimum relative track height
    MAX_TRACK_HEIGHT = 10.0  # Maximum relative track height

    # Required configuration
    tracks: List[BaseTrack] = field(default_factory=list)
    genomic_region: Optional[str] = None

    # Figure layout
    figure_width: float = DEFAULT_FIGURE_WIDTH
    figure_height: float = DEFAULT_FIGURE_HEIGHT
    track_heights: Optional[List[float]] = None
    row_heights: Optional[List[float]] = None  # Height proportions for subplot rows

    # Titles and labels
    figure_title: Optional[str] = None
    track_titles: Optional[List[str]] = None
    column_titles: Optional[List[str]] = None

    # Layout options
    vertical_spacing: float = PlotConstants.VERTICAL_SPACING
    horizontal_spacing: float = PlotConstants.HORIZONTAL_SPACING

    # Export options
    output_file: Optional[str] = None
    save_options: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate configuration after initialization."""
        self._validate_configuration()

    def _validate_configuration(self) -> None:
        """Validate the plot configuration."""
        if not self.tracks:
            raise ValueError("At least one track must be provided")

        if not all(isinstance(track, BaseTrack) for track in self.tracks):
            raise TypeError("All tracks must inherit from BaseTrack")

        if self.track_heights is not None:
            if len(self.track_heights) != len(self.tracks):
                raise ValueError(
                    f"Number of track heights ({len(self.track_heights)}) must match "
                    f"number of tracks ({len(self.tracks)})"
                )

            for height in self.track_heights:
                if not (self.MIN_TRACK_HEIGHT <= height <= self.MAX_TRACK_HEIGHT):
                    raise ValueError(
                        f"Track height {height} is outside valid range "
                        f"[{self.MIN_TRACK_HEIGHT}, {self.MAX_TRACK_HEIGHT}]"
                    )

        if self.track_titles is not None:
            if len(self.track_titles) != len(self.tracks):
                raise ValueError(
                    f"Number of track titles ({len(self.track_titles)}) must match "
                    f"number of tracks ({len(self.tracks)})"
                )

        if self.figure_width <= 0 or self.figure_height <= 0:
            raise ValueError("Figure dimensions must be positive")

    def get_track_heights(self) -> List[float]:
        """
        Get track heights, using defaults if not specified.

        Returns:
            List of relative track heights
        """
        if self.track_heights is not None:
            return self.track_heights.copy()

        # Use track-specific heights or default
        heights = []
        for track in self.tracks:
            track_height = track.get_height()
            heights.append(track_height)

        return heights

    def get_track_titles(self) -> List[str]:
        """
        Get track titles, using track names if not specified.

        Returns:
            List of track titles
        """
        if self.track_titles is not None:
            return self.track_titles.copy()

        # Use track names as titles
        return [
            track.track_name or f"Track {i+1}" for i, track in enumerate(self.tracks)
        ]


class GenomicPlotter:
    """
    Main plotter class for creating genomic visualizations.

    This class orchestrates the creation of genomic plots using the backend
    architecture and track system. It handles figure creation, track layout,
    and rendering coordination.
    """

    # Plotter constants
    DEFAULT_SUBPLOT_COLS = 1  # Default number of subplot columns
    DEFAULT_COLUMN_WIDTH = 1.0  # Default relative column width

    def __init__(self, config: PlotConfiguration, backend_name: str = "plotly", backend_options: Optional[Dict[str, Any]] = None):
        """
        Initialize the genomic plotter.

        Args:
            config: Plot configuration object
            backend_name: Name of the backend to use ("matplotlib" or "plotly")
            backend_options: Additional options for the backend
        """
        self.config = config
        self.backend_name = backend_name
        self.backend_options = backend_options or {}
        self.backend: Optional[BaseBackend] = None
        self.region_start: Optional[int] = None
        self.region_end: Optional[int] = None
        self.chromosome: Optional[str] = None

    def create_plot(self) -> BaseBackend:
        """
        Create the genomic plot.

        Returns:
            Backend with the created figure

        Raises:
            ValueError: If configuration is invalid or plotting fails
        """
        # Parse or infer genomic region
        self._infer_genomic_region()

        # Create backend
        self.backend = create_backend(
            self.backend_name, figure_options=self.backend_options
        )

        # Set up figure layout based on track subplot positions
        max_subplot_x = max(track.subplot_x for track in self.config.tracks) if self.config.tracks else 0
        max_subplot_y = max(track.subplot_y for track in self.config.tracks) if self.config.tracks else 0
        
        num_cols = max_subplot_x + 1
        num_rows = max_subplot_y + 1
        
        # Use row_heights if provided, otherwise use track_heights if they match num_rows, otherwise use defaults
        if self.config.row_heights is not None:
            raw_height_props = self.config.row_heights
        elif self.config.track_heights is not None and len(self.config.track_heights) == num_rows:
            raw_height_props = self.config.track_heights
        else:
            raw_height_props = [1.0] * num_rows
            
        # Convert proportions (should sum to 1.0) to ratios for matplotlib/plotly
        # If the values don't sum to 1.0, normalize them to proportions first
        total_height = sum(raw_height_props)
        if abs(total_height - 1.0) > 0.001:  # Not already normalized
            # These are relative values, convert to ratios by normalizing
            proportions = [h / total_height for h in raw_height_props]
        else:
            # Already proportions, use as-is
            proportions = raw_height_props
            
        # Convert proportions to height ratios for backend
        # Multiply by a scale factor to make differences more pronounced
        scale_factor = num_rows  # This makes the ratios more visible
        height_props = [p * scale_factor for p in proportions]
        width_props = [self.DEFAULT_COLUMN_WIDTH] * num_cols
        row_titles = self.config.get_track_titles() if self.config.track_titles else []
        
        # Use column titles from config if provided, otherwise empty titles
        if self.config.column_titles:
            column_titles = self.config.column_titles
        else:
            column_titles = [""] * num_cols

        # Create figure
        self.backend.create_figure(
            num_rows=num_rows,
            num_cols=num_cols,
            height_props=height_props,
            width_props=width_props,
            row_titles=row_titles,
            column_titles=column_titles,
            vertical_spacing=self.config.vertical_spacing,
            horizontal_spacing=self.config.horizontal_spacing,
            width=self.config.figure_width,
            height=self.config.figure_height,
        )

        # Set figure title
        if self.config.figure_title:
            self.backend.set_figure_title(self.config.figure_title)

        # Plot tracks
        self._plot_tracks()

        # Finalize figure
        self.backend.finalize_figure()

        # Save figure if requested
        if self.config.output_file:
            self.backend.save_figure(
                self.config.output_file, **self.config.save_options
            )

        return self.backend

    def _parse_genomic_region(self) -> None:
        """Parse the genomic region string."""
        if not self.config.genomic_region:
            return

        self.chromosome, self.region_start, self.region_end = validate_genomic_region(
            self.config.genomic_region
        )

    def _infer_genomic_region(self) -> None:
        """
        Infer genomic region from track data if no region is explicitly provided.
        
        This method examines all tracks with genomic data and finds the overall
        bounds (minimum start and maximum end) to set as the default region.
        """
        if self.config.genomic_region:
            # Region already provided, parse it
            self._parse_genomic_region()
            return

        # Collect genomic bounds from all tracks
        all_bounds = []
        chromosomes = set()

        for track in self.config.tracks:
            bounds = track.get_genomic_bounds()
            if bounds:
                all_bounds.append(bounds)
                chromosomes.add(bounds['chrom'])

        if not all_bounds:
            # No genomic data found, use defaults
            self.chromosome = "chr1"
            self.region_start = 0
            self.region_end = 1000000
            return

        if len(chromosomes) > 1:
            # Multiple chromosomes found, use the first one and warn
            print(f"Warning: Multiple chromosomes found in tracks: {chromosomes}. Using {list(chromosomes)[0]}")

        # Use the first chromosome
        self.chromosome = list(chromosomes)[0]

        # Find overall bounds
        self.region_start = min(bounds['start'] for bounds in all_bounds)
        self.region_end = max(bounds['end'] for bounds in all_bounds)

        # Set the genomic region string for consistency
        self.config.genomic_region = f"{self.chromosome}:{self.region_start}-{self.region_end}"

    def _plot_tracks(self) -> None:
        """Plot all tracks on their respective subplots."""
        if self.backend is None:
            raise ValueError("Backend must be created before plotting tracks")

        for track in self.config.tracks:
            # Use the track's subplot position
            row = track.subplot_y
            col = track.subplot_x

            # Get subplot for this track
            if self.backend_name.lower() in ["matplotlib", "mpl", "pyplot"]:
                subplot = self.backend.get_subplot(row, col)
                self._plot_track_matplotlib(track, subplot, row, col)
            else:
                subplot_coords = self.backend.get_subplot(row, col)
                self._plot_track_plotly(track, subplot_coords, row, col)

    def _plot_track_matplotlib(
        self, track: BaseTrack, subplot: Any, row: int, col: int
    ) -> None:
        """
        Plot a track using matplotlib backend.

        Args:
            track: Track to plot
            subplot: Matplotlib axes object
            row: Row index
            col: Column index
        """
        # Call track's matplotlib plotting method
        if hasattr(track, "plot_matplotlib"):
            if self.region_start is not None and self.region_end is not None:
                # Plot with specific region
                track.plot_matplotlib(
                    subplot,
                    start=self.region_start,
                    end=self.region_end,
                    chromosome=self.chromosome or "Unknown",
                )
            else:
                # Plot without region (track will handle defaults)
                track.plot_matplotlib(
                    subplot,
                    start=0,
                    end=1000000,  # Default large range
                    chromosome=self.chromosome or "Unknown",
                )
        else:
            raise ValueError(
                f"Track {track.track_name} does not support matplotlib plotting"
            )

        # Hide matplotlib spines for cleaner genomic plots
        if MatplotlibConstants.HIDE_TOP_SPINE:
            subplot.spines["top"].set_visible(False)
        if MatplotlibConstants.HIDE_RIGHT_SPINE:
            subplot.spines["right"].set_visible(False)

    def _plot_track_plotly(
        self, track: BaseTrack, subplot_coords: Tuple[int, int], row: int, col: int
    ) -> None:
        """
        Plot a track using plotly backend.

        Args:
            track: Track to plot
            subplot_coords: Plotly subplot coordinates (1-based)
            row: Row index (0-based)
            col: Column index (0-based)
        """
        if self.backend is None or self.backend.figure is None:
            raise ValueError("Backend and figure must be created before plotting")

        plotly_row, plotly_col = subplot_coords

        # Call track's plotly plotting method
        if hasattr(track, "plot_plotly"):
            if self.region_start is not None and self.region_end is not None:
                # Plot with specific region
                track.plot_plotly(
                    self.backend.figure,
                    plotly_row,
                    plotly_col,
                    start=self.region_start,
                    end=self.region_end,
                    chromosome=self.chromosome or "Unknown",
                )
            else:
                # Plot without region (track will handle defaults)
                track.plot_plotly(
                    self.backend.figure,
                    plotly_row,
                    plotly_col,
                    start=0,
                    end=1000000,  # Default large range
                    chromosome=self.chromosome or "Unknown",
                )
        else:
            raise ValueError(
                f"Track {track.track_name} does not support plotly plotting"
            )

    def show(self) -> None:
        """Display the plot."""
        if self.backend is None:
            raise ValueError("Plot must be created before showing")
        self.backend.show_figure()

    def save(self, filename: str, **kwargs) -> None:
        """
        Save the plot to file.

        Args:
            filename: Output filename
            **kwargs: Backend-specific save options
        """
        if self.backend is None:
            raise ValueError("Plot must be created before saving")
        self.backend.save_figure(filename, **kwargs)


# Convenience functions for easy plotting


def plot_genomic_tracks(
    tracks: List[BaseTrack],
    genomic_region: Optional[str] = None,
    backend: str = "plotly",
    figure_title: Optional[str] = None,
    output_file: Optional[str] = None,
    show_plot: bool = True,
    **kwargs,
) -> BaseBackend:
    """
    High-level function to plot genomic tracks.

    Args:
        tracks: List of tracks to plot
        genomic_region: Genomic region string (e.g., "chr1:1000-2000")
        backend: Backend to use ("matplotlib" or "plotly")
        figure_title: Title for the figure
        output_file: Optional output filename
        show_plot: Whether to display the plot
        **kwargs: Additional configuration options

    Returns:
        Backend with the created figure
    """
    # Create configuration
    config = PlotConfiguration(
        tracks=tracks,
        genomic_region=genomic_region,
        figure_title=figure_title,
        output_file=output_file,
        **kwargs,
    )

    # Create plotter and generate plot
    plotter = GenomicPlotter(config, backend_name=backend)
    backend_instance = plotter.create_plot()

    # Show plot if requested
    if show_plot:
        plotter.show()

    return backend_instance


def quick_plot(
    tracks: Union[BaseTrack, List[BaseTrack]],
    region: Optional[str] = None,
    backend: str = "plotly",
    title: Optional[str] = None,
    save_as: Optional[str] = None,
) -> BaseBackend:
    """
    Quick plotting function for simple use cases.

    Args:
        tracks: Single track or list of tracks
        region: Genomic region string
        backend: Backend to use
        title: Figure title
        save_as: Output filename

    Returns:
        Backend with the created figure
    """
    # Convert single track to list
    if isinstance(tracks, BaseTrack):
        tracks = [tracks]

    return plot_genomic_tracks(
        tracks=tracks,
        genomic_region=region,
        backend=backend,
        figure_title=title,
        output_file=save_as,
    )
