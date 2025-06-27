"""
Genomic coordinate utilities for jdgenometracks V3.

Handles parsing, formatting, and intelligent tick placement for genomic coordinates.
Consolidates genomic_ticks.py functionality with improved organization.
"""

import math
import re
from typing import List, Optional, Tuple

from .constants import GenomicConstants


class GenomicRegion:
    """Represents a genomic region with chromosome, start, and end coordinates."""

    def __init__(self, chromosome: str, start: int, end: int):
        self.chromosome = chromosome
        self.start = start
        self.end = end

    @property
    def size(self) -> int:
        """Get the size of the region in base pairs."""
        return self.end - self.start

    @classmethod
    def from_string(cls, region_str: str) -> "GenomicRegion":
        """
        Parse a genomic region string.

        Args:
            region_str: String like "chr1:1000-2000" or "chr1:1,000-2,000"

        Returns:
            GenomicRegion object

        Raises:
            ValueError: If string format is invalid
        """
        # Remove commas from coordinates
        region_str = region_str.replace(",", "")

        # Match patterns like "chr1:1000-2000"
        match = re.match(r"([^:]+):(\d+)-(\d+)", region_str)
        if not match:
            raise ValueError(f"Invalid region format: {region_str}")

        chromosome = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))

        if start >= end:
            raise ValueError(
                f"Start coordinate must be less than end: {start} >= {end}"
            )

        return cls(chromosome, start, end)

    def __str__(self) -> str:
        """String representation of the region."""
        return f"{self.chromosome}:{self.start:,}-{self.end:,}"

    def __repr__(self) -> str:
        return f"GenomicRegion('{self.chromosome}', {self.start}, {self.end})"


class GenomicTickCalculator:
    """
    Calculate optimal tick marks for genomic coordinates based on the visible range.

    Provides intelligent tick placement that scales appropriately:
    - Base pair level for very small ranges (< 1kb)
    - 10s-100s of bp for small ranges (1kb - 10kb)
    - Hundreds of bp to kb for medium ranges (10kb - 1Mb)
    - kb to Mb for large ranges (> 1Mb)
    """

    # Tick spacing preferences for different scales
    BASE_SPACINGS = [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500]  # Base pairs
    KILO_SPACINGS = [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500]  # Kilobases
    MEGA_SPACINGS = [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500]  # Megabases

    # Target number of ticks (will try to get close to this)
    TARGET_TICK_COUNT = 8
    MIN_TICK_COUNT = 4
    MAX_TICK_COUNT = 12

    @classmethod
    def calculate_optimal_ticks(
        cls, range_start: int, range_end: int
    ) -> Tuple[List[int], str, int]:
        """
        Calculate optimal tick locations and formatting for a genomic range.

        Args:
            range_start: Start of the genomic range
            range_end: End of the genomic range

        Returns:
            Tuple of (tick_locations, unit_label, scale_factor)
            - tick_locations: List of genomic coordinates for tick marks
            - unit_label: String label for the unit ("b", "Kb", "Mb")
            - scale_factor: Factor to divide coordinates by for display
        """
        range_size = range_end - range_start

        if range_size <= 0:
            return [range_start, range_end], "b", 1

        # Determine the appropriate scale
        if range_size < 1000:  # < 1kb - use base pairs
            unit_label = GenomicConstants.BASE_UNIT
            scale_factor = 1
            spacings = cls.BASE_SPACINGS
        elif range_size < 1_000_000:  # < 1Mb - use kilobases
            unit_label = GenomicConstants.KILOBASE_UNIT
            scale_factor = 1000
            spacings = [s * 1000 for s in cls.KILO_SPACINGS]  # Convert to bp
        else:  # >= 1Mb - use megabases
            unit_label = GenomicConstants.MEGABASE_UNIT
            scale_factor = 1_000_000
            spacings = [s * 1_000_000 for s in cls.MEGA_SPACINGS]  # Convert to bp

        # Find the best spacing
        best_spacing = cls._find_best_spacing(range_size, spacings)

        # Calculate tick locations
        ticks = cls._generate_ticks(range_start, range_end, best_spacing)

        return ticks, unit_label, scale_factor

    @classmethod
    def _find_best_spacing(cls, range_size: int, spacings: List[int]) -> int:
        """Find the spacing that gives the most appropriate number of ticks."""
        best_spacing = spacings[0]
        best_score = float("inf")

        for spacing in spacings:
            if spacing > range_size:
                continue

            tick_count = range_size / spacing

            # Score based on how close we are to the target tick count
            if tick_count < cls.MIN_TICK_COUNT:
                score = cls.TARGET_TICK_COUNT - tick_count + 10  # Penalty for too few
            elif tick_count > cls.MAX_TICK_COUNT:
                score = tick_count - cls.TARGET_TICK_COUNT + 5  # Penalty for too many
            else:
                score = abs(tick_count - cls.TARGET_TICK_COUNT)

            if score < best_score:
                best_score = score
                best_spacing = spacing

        return best_spacing

    @classmethod
    def _generate_ticks(
        cls, range_start: int, range_end: int, spacing: int
    ) -> List[int]:
        """Generate tick marks at regular intervals."""
        # Find the first tick position (round up to nearest spacing multiple)
        first_tick = math.ceil(range_start / spacing) * spacing

        # Generate all tick positions
        ticks = []
        current_tick = first_tick
        while current_tick <= range_end:
            ticks.append(current_tick)
            current_tick += spacing

        # Ensure we have at least start and end
        if not ticks or ticks[0] > range_start:
            ticks.insert(0, range_start)
        if not ticks or ticks[-1] < range_end:
            ticks.append(range_end)

        return ticks


def format_genomic_coordinate(
    coord: int, unit_label: str, scale_factor: int, precision: int = 1
) -> str:
    """
    Format a genomic coordinate with appropriate units.

    Args:
        coord: Genomic coordinate
        unit_label: Unit to display ("b", "Kb", "Mb")
        scale_factor: Factor to divide coordinate by
        precision: Number of decimal places

    Returns:
        Formatted coordinate string

    Examples:
        >>> format_genomic_coordinate(1500, "Kb", 1000)
        "1.5 Kb"
        >>> format_genomic_coordinate(2000000, "Mb", 1000000)
        "2.0 Mb"
    """
    scaled_coord = coord / scale_factor

    # Use minimal precision if coordinate is a round number
    if scaled_coord == int(scaled_coord):
        return f"{int(scaled_coord)} {unit_label}"
    else:
        return f"{scaled_coord:.{precision}f} {unit_label}"


def format_genomic_ticks(
    tick_locations: List[int],
    unit_label: str,
    scale_factor: int,
    show_unit_on_last: bool = True,
    precision: int = 1,
) -> List[str]:
    """
    Format a list of genomic tick locations.

    Args:
        tick_locations: List of genomic coordinates
        unit_label: Unit to display
        scale_factor: Factor to divide coordinates by
        show_unit_on_last: Whether to show unit only on last tick
        precision: Number of decimal places

    Returns:
        List of formatted tick labels
    """
    labels = []

    for i, coord in enumerate(tick_locations):
        scaled_coord = coord / scale_factor

        # Format the number
        if scaled_coord == int(scaled_coord):
            label = str(int(scaled_coord))
        else:
            label = f"{scaled_coord:.{precision}f}"

        # Add unit to last tick if requested, or to all if not
        if show_unit_on_last:
            if i == len(tick_locations) - 1:
                label += f" {unit_label}"
        else:
            label += f" {unit_label}"

        labels.append(label)

    return labels


def calculate_axis_shift(
    regions: List[GenomicRegion], relative_x_axis: bool = False
) -> int:
    """
    Calculate the axis shift for relative positioning.

    Args:
        regions: List of genomic regions
        relative_x_axis: Whether to use relative positioning

    Returns:
        Axis shift value
    """
    if not relative_x_axis or not regions:
        return 0

    # Use the minimum start coordinate across all regions
    return min(region.start for region in regions)


def parse_region_string(region_str: Optional[str]) -> Optional[GenomicRegion]:
    """
    Parse a region string, returning None if invalid or None.

    Args:
        region_str: Region string or None

    Returns:
        GenomicRegion object or None
    """
    if not region_str:
        return None

    try:
        return GenomicRegion.from_string(region_str)
    except ValueError:
        return None
