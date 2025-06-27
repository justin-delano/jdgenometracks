"""
Track implementations for jdgenometracks V3.

This module provides all track classes for genomic data visualization:
- BaseTrack: Abstract base class for all tracks (NEW V3 ARCHITECTURE)
- XAxisTrack: For genomic coordinate axis display (REFACTORED)
- SpacerTrack: For empty space between tracks (REFACTORED)
- BedTrack: For genomic regions (BED format) (REFACTORED)
- BedGraphTrack: For continuous genomic data (bedGraph format) (REFACTORED)

All tracks now inherit from BaseTrack and use shared functionality for:
- File reading and validation
- Data filtering and coordinate processing
- Horizontal line plotting
- Common plotting utilities

Also provides factory functions for creating tracks from configuration.
"""

from .axis_track import XAxisTrack
from .base_track import BaseTrack

# Fully refactored tracks (BaseTrack-compatible)
from .bed_track import BedTrack
from .bedgraph_track import BedGraphTrack
from .factory import (
    TrackFactory,
    TrackRegistry,
    create_track,
    get_available_track_types,
    register_track_type,
)
from .spacer_track import SpacerTrack

__all__ = [
    # New V3 architecture
    "BaseTrack",
    # Fully refactored tracks (BaseTrack-compatible)
    "XAxisTrack",
    "SpacerTrack",
    "BedTrack",
    "BedGraphTrack",
    # Factory classes and functions
    "TrackRegistry",
    "TrackFactory",
    "create_track",
    "register_track_type",
    "get_available_track_types",
]
