"""
Layout optimization utilities for genomic track plotting.

This module provides efficient algorithms for assigning y-coordinates to genomic regions
to prevent overlaps, replacing the memory-intensive sparse matrix approach.
"""

import heapq
from typing import List

import pandas as pd


def assign_y_levels_sweep_line(
    regions_df: pd.DataFrame, rect_height: float = 1.0, rect_padding: float = 0.0
) -> List[float]:
    """
    Assign y-levels to BED regions using a sweep line algorithm to prevent overlaps.

    This algorithm is more efficient than the sparse matrix approach:
    - Time complexity: O(n log n) where n is the number of regions
    - Space complexity: O(n) instead of O(genomic_range_width)
    - Optimal y-level assignment: Uses minimal number of y-levels needed

    Args:
        regions_df: DataFrame with 'chromStart' and 'chromEnd' columns
        rect_height: Height of each rectangle
        rect_padding: Padding between rectangles

    Returns:
        List of y-values for each region in the same order as input DataFrame
    """
    if regions_df.empty:
        return []

    # Create events for start and end of each region
    events = []
    for idx, region in regions_df.iterrows():
        events.append((region["chromStart"], "start", idx))
        events.append((region["chromEnd"], "end", idx))

    # Sort events by position, with 'end' events before 'start' events at same position
    # This ensures that when a region ends at position X and another starts at X,
    # the y-level is freed before being potentially reused
    events.sort(key=lambda x: (x[0], x[1] == "start"))

    # Track active regions and their y-levels
    active_regions = {}  # region_idx -> y_level
    available_y_levels = []  # min-heap of available y-levels
    next_y_level = 0.0
    region_y_levels = {}  # region_idx -> y_level

    for pos, event_type, region_idx in events:
        if event_type == "start":
            # Assign y-level to starting region
            if available_y_levels:
                y_level = heapq.heappop(available_y_levels)
            else:
                y_level = next_y_level
                next_y_level += rect_height + rect_padding

            active_regions[region_idx] = y_level
            region_y_levels[region_idx] = y_level

        else:  # event_type == 'end'
            # Free up y-level from ending region
            y_level = active_regions.pop(region_idx)
            heapq.heappush(available_y_levels, y_level)

    # Return y-levels in original DataFrame order
    return [region_y_levels[i] for i in range(len(regions_df))]


def assign_y_levels_binned(
    regions_df: pd.DataFrame,
    genomic_start: int,
    genomic_end: int,
    bin_size: int = 1000,
    rect_height: float = 1.0,
    rect_padding: float = 0.0,
) -> List[float]:
    """
    Alternative approach using binning for very large genomic regions.

    This approach divides the genomic region into bins and tracks which y-levels
    are occupied in each bin. Less optimal than sweep line but can be useful
    for specific use cases.

    Args:
        regions_df: DataFrame with 'chromStart' and 'chromEnd' columns
        genomic_start: Start of the genomic region being plotted
        genomic_end: End of the genomic region being plotted
        bin_size: Size of each bin in base pairs
        rect_height: Height of each rectangle
        rect_padding: Padding between rectangles

    Returns:
        List of y-values for each region in the same order as input DataFrame
    """
    if regions_df.empty:
        return []

    # Calculate number of bins needed
    num_bins = (genomic_end - genomic_start + bin_size - 1) // bin_size

    # Track occupied y-levels for each bin
    bin_occupancy = [set() for _ in range(num_bins)]

    region_y_levels = []

    for idx, region in regions_df.iterrows():
        # Find which bins this region overlaps
        start_bin = max(0, (region["chromStart"] - genomic_start) // bin_size)
        end_bin = min(num_bins - 1, (region["chromEnd"] - genomic_start) // bin_size)

        # Find the lowest available y-level across all overlapping bins
        occupied_levels = set()
        for bin_idx in range(start_bin, end_bin + 1):
            occupied_levels.update(bin_occupancy[bin_idx])

        # Find the first available y-level
        y_level = 0.0
        while y_level in occupied_levels:
            y_level += rect_height + rect_padding

        # Mark this y-level as occupied in all overlapping bins
        for bin_idx in range(start_bin, end_bin + 1):
            bin_occupancy[bin_idx].add(y_level)

        region_y_levels.append(y_level)

    return region_y_levels


def get_max_y_level(y_levels: List[float], rect_height: float = 1.0) -> float:
    """
    Get the maximum y-coordinate that will be used for plotting.

    Args:
        y_levels: List of y-levels assigned to regions
        rect_height: Height of each rectangle

    Returns:
        Maximum y-coordinate (bottom of topmost rectangle + height)
    """
    if not y_levels:
        return 0.0
    return max(y_levels) + rect_height
