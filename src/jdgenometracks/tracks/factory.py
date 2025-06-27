"""
Track factory and registry for jdgenometracks V3.

This module provides a factory system for creating track instances
based on track type and configuration. It handles track registration
and creation with proper validation.

Author: Assistant
Date: 2024
"""

from typing import Any, Dict, Optional, Type

from .axis_track import XAxisTrack
from .base_track import BaseTrack
from .bed_track import BedTrack
from .bedgraph_track import BedGraphTrack
from .spacer_track import SpacerTrack


class TrackRegistry:
    """
    Registry for track classes.

    Maintains a mapping of track types to their corresponding classes
    and provides methods for registration and retrieval.
    """

    def __init__(self):
        self._tracks: Dict[str, Type[BaseTrack]] = {}
        self._register_default_tracks()

    def _register_default_tracks(self) -> None:
        """Register the default track types."""
        # Register tracks that have been refactored to inherit from BaseTrack
        self.register("spacer", SpacerTrack)
        self.register("x_axis", XAxisTrack)
        self.register("axis", XAxisTrack)  # Alias

        # Register the refactored genomic data tracks
        self.register("bed", BedTrack)
        self.register("bedgraph", BedGraphTrack)

    def register(self, track_type: str, track_class: Type[BaseTrack]) -> None:
        """
        Register a track class for a given track type.

        Args:
            track_type: The string identifier for the track type
            track_class: The track class to register

        Raises:
            ValueError: If track_class is not a subclass of BaseTrack
        """
        if not issubclass(track_class, BaseTrack):
            raise ValueError(f"Track class {track_class} must inherit from BaseTrack")

        self._tracks[track_type.lower()] = track_class

    def get_track_class(self, track_type: str) -> Type[BaseTrack]:
        """
        Get the track class for a given track type.

        Args:
            track_type: The track type identifier

        Returns:
            The track class

        Raises:
            ValueError: If track type is not registered
        """
        track_type = track_type.lower()
        if track_type not in self._tracks:
            available_types = list(self._tracks.keys())
            raise ValueError(
                f"Unknown track type: {track_type}. Available types: {available_types}"
            )

        return self._tracks[track_type]

    def list_track_types(self) -> list[str]:
        """
        Get a list of all registered track types.

        Returns:
            List of track type identifiers
        """
        return list(self._tracks.keys())

    def is_registered(self, track_type: str) -> bool:
        """
        Check if a track type is registered.

        Args:
            track_type: The track type to check

        Returns:
            True if the track type is registered, False otherwise
        """
        return track_type.lower() in self._tracks


class TrackFactory:
    """
    Factory for creating track instances.

    Uses the TrackRegistry to create track instances based on
    configuration dictionaries. Handles validation and setup.
    """

    def __init__(self, registry: Optional[TrackRegistry] = None):
        self.registry = registry or TrackRegistry()

    def create_track(self, track_config: Dict[str, Any]) -> BaseTrack:
        """
        Create a track instance from a configuration dictionary.

        Args:
            track_config: Dictionary containing track configuration
                Required keys:
                - track_type: String identifier for track type
                Optional keys:
                - file_path: Path to data file
                - track_name: Human-readable name
                - track_options: Styling and display options
                - show_legend: Whether to show legend
                - hlines: List of horizontal reference lines

        Returns:
            An instance of the appropriate track class

        Raises:
            ValueError: If required configuration is missing or invalid
        """
        # Validate required fields
        if "track_type" not in track_config:
            raise ValueError("track_config must contain 'track_type'")

        track_type = track_config["track_type"]

        # Get the track class
        track_class = self.registry.get_track_class(track_type)

        # Extract configuration parameters for BaseTrack
        track_kwargs = {
            "track_type": track_type,
            "file_path": track_config.get("file_path"),
            "track_name": track_config.get("track_name"),
            "track_options": track_config.get(
                "track_options", track_config.get("options", {})
            ),
            "show_legend": track_config.get("show_legend", False),
            "hlines": track_config.get("hlines", []),
            "subplot_x": track_config.get("subplot_x", 0),  # Default to column 0
            "subplot_y": track_config.get("subplot_y", 0),  # Default to row 0
        }

        # Handle data if provided directly
        if "data" in track_config:
            track_kwargs["data"] = track_config["data"]

        # Create and return the track instance
        try:
            return track_class(**track_kwargs)
        except Exception as e:
            raise ValueError(f"Failed to create {track_type} track: {e}") from e

    def create_tracks_from_config(self, config: Dict[str, Any]) -> list[BaseTrack]:
        """
        Create multiple tracks from a complete configuration.

        Args:
            config: Configuration dictionary containing 'tracks' list

        Returns:
            List of track instances

        Raises:
            ValueError: If configuration is invalid
        """
        if "tracks" not in config:
            raise ValueError("Configuration must contain 'tracks' list")

        tracks = []
        for i, track_config in enumerate(config["tracks"]):
            try:
                track = self.create_track(track_config)
                tracks.append(track)
            except Exception as e:
                raise ValueError(f"Failed to create track {i}: {e}") from e

        return tracks


# Global instances for convenience
_default_registry = TrackRegistry()
_default_factory = TrackFactory(_default_registry)


def create_track(track_config: Dict[str, Any]) -> BaseTrack:
    """
    Convenience function to create a track using the default factory.

    Args:
        track_config: Track configuration dictionary

    Returns:
        Track instance
    """
    return _default_factory.create_track(track_config)


def register_track_type(track_type: str, track_class: Type[BaseTrack]) -> None:
    """
    Convenience function to register a track type with the default registry.

    Args:
        track_type: Track type identifier
        track_class: Track class to register
    """
    _default_registry.register(track_type, track_class)


def get_available_track_types() -> list[str]:
    """
    Get list of all available track types.

    Returns:
        List of track type identifiers
    """
    return _default_registry.list_track_types()
