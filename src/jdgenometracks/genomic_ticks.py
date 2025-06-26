"""
Enhanced genomic axis utilities for intelligent tick mark placement.

This module provides utilities for calculating optimal tick marks and formatting
based on the genomic coordinate range being displayed.
"""

import math
from typing import List, Tuple


class GenomicTickCalculator:
    """
    Calculate optimal tick marks for genomic coordinates based on the visible range.

    This class provides intelligent tick placement that scales appropriately:
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
            unit_label = "b"
            scale_factor = 1
            spacings = cls.BASE_SPACINGS
        elif range_size < 1_000_000:  # < 1Mb - use kilobases
            unit_label = "Kb"
            scale_factor = 1000
            spacings = [s * 1000 for s in cls.KILO_SPACINGS]  # Convert to bp
        else:  # >= 1Mb - use megabases
            unit_label = "Mb"
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
        # Start from the first tick that's >= range_start
        first_tick = math.ceil(range_start / spacing) * spacing

        ticks = []
        current = first_tick
        while current <= range_end:
            ticks.append(current)
            current += spacing

        # Avoid crowding at boundaries - be aggressive about preventing visual overlap
        if not ticks:
            ticks = [range_start, range_end]
        else:
            # Use 70% of spacing as minimum distance to prevent label overlap
            # This ensures sufficient visual space between tick labels
            min_boundary_distance = spacing * 0.7

            # Add start point only if it's far enough from the first tick
            if (
                ticks[0] - range_start > min_boundary_distance
                and range_start not in ticks
            ):
                ticks.insert(0, range_start)

            # Add end point only if it's far enough from the last tick
            if range_end - ticks[-1] > min_boundary_distance and range_end not in ticks:
                ticks.append(range_end)

        return ticks


def format_genomic_coordinate(coord: int, unit_label: str, scale_factor: int) -> str:
    """
    Format a genomic coordinate for display.

    Args:
        coord: The genomic coordinate
        unit_label: Unit label ("b", "Kb", "Mb")
        scale_factor: Factor to divide coordinate by

    Returns:
        Formatted coordinate string
    """
    scaled_coord = coord / scale_factor

    # Choose appropriate number of decimal places
    if scale_factor == 1:  # Base pairs
        return f"{coord:,d}"
    elif scaled_coord == int(scaled_coord):  # Whole number
        return f"{int(scaled_coord):,d}"
    elif scaled_coord < 10:  # Show 1 decimal for small numbers
        return f"{scaled_coord:.1f}"
    else:  # No decimals for larger numbers
        return f"{scaled_coord:.0f}"


def format_genomic_ticks(
    ticks: List[int], unit_label: str, scale_factor: int, show_unit_on_last: bool = True
) -> List[str]:
    """
    Format a list of tick coordinates for display.

    Args:
        ticks: List of genomic coordinates
        unit_label: Unit label ("b", "Kb", "Mb")
        scale_factor: Factor to divide coordinates by
        show_unit_on_last: Whether to show unit label on the last tick

    Returns:
        List of formatted tick labels
    """
    if not ticks:
        return []

    labels = []
    for i, tick in enumerate(ticks):
        formatted = format_genomic_coordinate(tick, unit_label, scale_factor)

        # Add unit to the last tick (or all ticks if requested)
        if show_unit_on_last and i == len(ticks) - 1:
            formatted += f" {unit_label}"

        labels.append(formatted)

    return labels
