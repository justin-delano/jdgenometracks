"""
Configuration parser for jdgenometracks V3.

This module provides comprehensive configuration parsing and validation for
genomic visualizations, integrating with the new backend architecture and
track factory system.

Author: Assistant
Date: 2024
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union

from ..backends import get_available_backends
from ..tracks.factory import TrackFactory
from ..utils.constants import PlotConstants
from ..utils.options import (
    get_internal_options,
    merge_options,
    normalize_options,
    translate,
)
from ..utils.units import convert_to_inches, convert_to_pixels, parse_size_with_units
from ..utils.validation import validate_genomic_region
from .plotter import GenomicPlotter, PlotConfiguration


@dataclass
class ConfigParser:
    """
    Parser for genomic visualization configurations.

    Handles loading, validation, and conversion of configuration dictionaries
    or JSON files into plot configurations that can be used with the plotting system.

    Attributes:
        track_factory: Factory for creating track instances
        default_backend: Default backend to use if not specified
        validate_files: Whether to validate that referenced files exist
    """

    # Config parser constants
    DEFAULT_BACKEND = "plotly"  # Default plotting backend
    SUPPORTED_CONFIG_FORMATS = [".json", ".yaml", ".yml"]  # Supported config file formats
    REQUIRED_CONFIG_KEYS = ["tracks"]  # Required top-level keys
    DEFAULT_GENOMIC_REGION = "chr1:1-1000000"  # Default region if not specified

    track_factory: TrackFactory = field(default_factory=TrackFactory)
    default_backend: str = field(default=DEFAULT_BACKEND)
    validate_files: bool = field(default=True)

    def __post_init__(self):
        """Validate configuration parser settings."""
        available_backends = get_available_backends()
        if self.default_backend not in available_backends:
            available = list(available_backends.keys())
            raise ValueError(f"Invalid default backend '{self.default_backend}'. Available: {available}")

    def load_config_file(self, config_path: str) -> Dict[str, Any]:
        """
        Load configuration from a file.

        Args:
            config_path: Path to configuration file (JSON or YAML)

        Returns:
            Configuration dictionary

        Raises:
            FileNotFoundError: If config file doesn't exist
            ValueError: If file format is unsupported or content is invalid
        """
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Configuration file not found: {config_path}")

        # Check file extension
        _, ext = os.path.splitext(config_path)
        if ext.lower() not in self.SUPPORTED_CONFIG_FORMATS:
            raise ValueError(
                f"Unsupported config format '{ext}'. Supported: {self.SUPPORTED_CONFIG_FORMATS}"
            )

        try:
            if ext.lower() == '.json':
                with open(config_path, 'r') as f:
                    config = json.load(f)
            else:
                # For YAML support (optional, requires PyYAML)
                try:
                    import yaml
                    with open(config_path, 'r') as f:
                        config = yaml.safe_load(f)
                except ImportError:
                    raise ValueError("YAML support requires PyYAML package. Install with: pip install PyYAML")
            
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in config file: {e}")
        except Exception as e:
            raise ValueError(f"Error reading config file: {e}")

        if not isinstance(config, dict):
            raise ValueError("Configuration must be a dictionary/object at top level")

        return config

    def parse_config(self, config: Union[str, Dict[str, Any]]) -> PlotConfiguration:
        """
        Parse a configuration into a PlotConfiguration object.

        Args:
            config: Configuration dictionary or path to config file

        Returns:
            PlotConfiguration object ready for plotting

        Raises:
            ValueError: If configuration is invalid
            FileNotFoundError: If config file doesn't exist
        """
        # Load from file if string path provided
        if isinstance(config, str):
            config_dict = self.load_config_file(config)
        elif isinstance(config, dict):
            config_dict = config
        else:
            raise ValueError("Config must be a dictionary or file path string")

        # Validate required keys
        self._validate_config_structure(config_dict)

        # Parse track configurations
        tracks = self._parse_tracks(config_dict["tracks"])

        # Parse figure options
        figure_options = self._parse_figure_options(config_dict)

        # Parse genomic region
        genomic_region = self._parse_genomic_region(config_dict)

        # Create plot configuration
        plot_config_kwargs = {
            "tracks": tracks,
            "genomic_region": genomic_region,
            "vertical_spacing": figure_options.get("vertical_spacing", PlotConstants.VERTICAL_SPACING),
            "horizontal_spacing": figure_options.get("horizontal_spacing", PlotConstants.HORIZONTAL_SPACING)
        }

        # Add optional parameters only if provided
        if "title" in figure_options:
            plot_config_kwargs["figure_title"] = figure_options["title"]
        
        if "width" in figure_options:
            plot_config_kwargs["figure_width"] = figure_options["width"]
        elif "total_width" in figure_options:
            plot_config_kwargs["figure_width"] = figure_options["total_width"]
        else:
            plot_config_kwargs["figure_width"] = PlotConfiguration.DEFAULT_FIGURE_WIDTH
            
        if "height" in figure_options:
            plot_config_kwargs["figure_height"] = figure_options["height"]
        elif "total_height" in figure_options:
            plot_config_kwargs["figure_height"] = figure_options["total_height"]
        else:
            plot_config_kwargs["figure_height"] = PlotConfiguration.DEFAULT_FIGURE_HEIGHT

        if "track_heights" in figure_options:
            plot_config_kwargs["track_heights"] = figure_options["track_heights"]
        
        if "height_props" in figure_options:
            # Map height_props to row_heights for subplot row proportions
            plot_config_kwargs["row_heights"] = figure_options["height_props"]
        
        if "track_titles" in figure_options:
            plot_config_kwargs["track_titles"] = figure_options["track_titles"]
            
        if "column_titles" in figure_options:
            plot_config_kwargs["column_titles"] = figure_options["column_titles"]

        plot_config = PlotConfiguration(**plot_config_kwargs)

        plot_config = PlotConfiguration(**plot_config_kwargs)

        return plot_config

    def _auto_detect_track_type(self, file_path: str) -> Optional[str]:
        """
        Auto-detect track type from file extension.
        
        Args:
            file_path: Path to the track data file
            
        Returns:
            Detected track type or None if cannot detect
        """
        if not file_path:
            return None
            
        import os
        ext = os.path.splitext(file_path)[1].lower()
        
        # Define extension to track type mapping
        extension_mapping = {
            '.bed': 'bed',
            '.bedgraph': 'bedgraph', 
            '.bg': 'bedgraph',
            '.bigwig': 'bigwig',
            '.bw': 'bigwig',
            '.gtf': 'gtf',
            '.gff': 'gff',
            '.gff3': 'gff'
        }
        
        return extension_mapping.get(ext)

    def _validate_config_structure(self, config: Dict[str, Any]) -> None:
        """
        Validate the basic structure of a configuration dictionary.

        Args:
            config: Configuration dictionary to validate

        Raises:
            ValueError: If configuration structure is invalid
        """
        # Check required keys
        for key in self.REQUIRED_CONFIG_KEYS:
            if key not in config:
                raise ValueError(f"Configuration missing required key: '{key}'")

        # Validate tracks list
        tracks = config["tracks"]
        if not isinstance(tracks, list):
            raise ValueError("'tracks' must be a list")
        
        if len(tracks) == 0:
            raise ValueError("'tracks' list cannot be empty")

        # Validate each track configuration
        for i, track_config in enumerate(tracks):
            if not isinstance(track_config, dict):
                raise ValueError(f"Track {i} must be a dictionary")
            if "track_type" not in track_config:
                # Try to auto-detect track type from file path
                file_path = track_config.get("file_path")
                if file_path:
                    detected_type = self._auto_detect_track_type(file_path)
                    if detected_type:
                        # Auto-detection successful - this will be applied in _parse_tracks
                        continue
                raise ValueError(f"Track {i} missing 'track_type' field and cannot auto-detect from file path: {file_path}")

    def _parse_tracks(self, track_configs: List[Dict[str, Any]]) -> List[Any]:
        """
        Parse track configurations into track instances.

        This method processes each track configuration, normalizes options,
        translates unified styling options to backend-specific format,
        and creates track instances using the factory.

        Args:
            track_configs: List of track configuration dictionaries

        Returns:
            List of track instances

        Raises:
            ValueError: If any track configuration is invalid
        """
        tracks = []
        
        for i, track_config in enumerate(track_configs):
            try:
                # Make a copy to avoid modifying original config
                processed_config = dict(track_config)
                
                # Auto-detect track type if missing
                if "track_type" not in processed_config:
                    file_path = processed_config.get("file_path")
                    if file_path:
                        detected_type = self._auto_detect_track_type(file_path)
                        if detected_type:
                            processed_config["track_type"] = detected_type
                        else:
                            # Default to bedgraph if can't detect
                            processed_config["track_type"] = "bedgraph"
                
                # Set default subplot positioning if not specified
                if "subplot_x" not in processed_config:
                    processed_config["subplot_x"] = 0  # Default to single column layout
                if "subplot_y" not in processed_config:
                    processed_config["subplot_y"] = i  # Each track gets its own row
                
                # Validate file paths if enabled
                if self.validate_files and "file_path" in processed_config:
                    file_path = processed_config["file_path"]
                    if file_path and not os.path.exists(file_path):
                        raise FileNotFoundError(f"Track data file not found: {file_path}")

                # Process unified options if present
                if "style" in processed_config or "options" in processed_config:
                    processed_config = self._process_track_options(processed_config)

                # Create track using factory
                track = self.track_factory.create_track(processed_config)
                tracks.append(track)

            except Exception as e:
                raise ValueError(f"Error creating track {i} (type: {track_config.get('track_type', 'unknown')}): {e}")

        return tracks

    def _process_track_options(self, track_config: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process and translate unified options for a track configuration.

        Args:
            track_config: Raw track configuration dictionary

        Returns:
            Processed track configuration with translated options
        """
        # Extract style and options
        style_options = track_config.pop("style", {})
        general_options = track_config.pop("options", {})
        
        # Normalize all options
        if style_options:
            style_options = normalize_options(style_options)
        if general_options:
            general_options = normalize_options(general_options)
        
        # Merge style and general options (general takes precedence)
        unified_options = merge_options(style_options, general_options)
        
        if unified_options:
            # Get the backend for this configuration or use default
            backend = track_config.get("backend", self.default_backend)
            
            # Map backend names to option translation targets
            backend_target = "plotly" if backend == "plotly" else "mpl"
            
            # Translate unified options to backend-specific format
            translated_options = translate(unified_options, backend_target)
            
            # Extract internal options that don't get translated
            internal_options = get_internal_options(unified_options)
            
            # Apply translated options to track config
            for element_type, element_options in translated_options.items():
                if element_type in track_config:
                    # Merge with existing options
                    if isinstance(track_config[element_type], dict):
                        track_config[element_type].update(element_options)
                    else:
                        # If not a dict, replace entirely
                        track_config[element_type] = element_options
                else:
                    # Add new element options
                    track_config[element_type] = element_options
            
            # Apply internal options directly to track config
            track_config.update(internal_options)
        
        return track_config

    def _parse_figure_options(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """
        Parse and validate figure options from configuration.

        This method processes figure-level options, handles unified styling,
        and prepares backend-specific configuration.

        Args:
            config: Full configuration dictionary

        Returns:
            Processed figure options dictionary
        """
        figure_options = config.get("figure", config.get("figure_options", {}))
        
        # Parse backend
        backend = figure_options.get("backend", self.default_backend)
        available_backends = get_available_backends()
        if backend not in available_backends:
            available = list(available_backends.keys())
            raise ValueError(f"Unknown backend '{backend}'. Available: {available}")

        # Parse dimensions with unit conversion
        processed_options = {"backend": backend}
        
        # Handle width and height with unit conversion
        for dimension in ["width", "height"]:
            if dimension in figure_options:
                value = figure_options[dimension]
                if isinstance(value, str):
                    # Parse size with units
                    numeric_value, unit = parse_size_with_units(value, "px")
                    
                    # Convert to appropriate units for backend
                    if backend == "plotly":
                        # Plotly uses pixels
                        processed_options[dimension] = convert_to_pixels(numeric_value, unit)
                    else:
                        # Matplotlib uses inches
                        processed_options[dimension] = convert_to_inches(numeric_value, unit)
                else:
                    processed_options[dimension] = value

        # Process unified styling options for figure-level elements
        if "style" in figure_options or "options" in figure_options:
            processed_options.update(self._process_figure_style_options(figure_options, backend))

        # Copy other figure options including total_width and total_height
        for key in ["title", "height_props", "width_props", "row_titles", "column_titles", "backend_options", "total_width", "total_height"]:
            if key in figure_options:
                processed_options[key] = figure_options[key]

        return processed_options

    def _process_figure_style_options(self, figure_options: Dict[str, Any], backend: str) -> Dict[str, Any]:
        """
        Process unified styling options for figure-level elements.

        Args:
            figure_options: Figure options dictionary
            backend: Target backend

        Returns:
            Processed style options
        """
        # Extract style and options
        style_options = figure_options.get("style", {})
        general_options = figure_options.get("options", {})
        
        # Normalize all options
        if style_options:
            style_options = normalize_options(style_options)
        if general_options:
            general_options = normalize_options(general_options)
        
        # Merge style and general options (general takes precedence)
        unified_options = merge_options(style_options, general_options)
        
        processed_style = {}
        
        if unified_options:
            # Map backend names to option translation targets
            backend_target = "plotly" if backend == "plotly" else "mpl"
            
            # Translate unified options to backend-specific format
            translated_options = translate(unified_options, backend_target)
            
            # Apply translated options to figure config
            for element_type, element_options in translated_options.items():
                if element_type in ["title", "legend", "layout"]:
                    # Figure-level styling
                    processed_style[element_type] = element_options
                elif element_type == "xaxis" or element_type == "yaxis":
                    # Axis styling - store for later application
                    processed_style[element_type] = element_options
                else:
                    # Other styling options
                    processed_style[element_type] = element_options
        
        return processed_style

    def _parse_genomic_region(self, config: Dict[str, Any]) -> Optional[str]:
        """
        Parse and validate genomic region from configuration.

        Args:
            config: Configuration dictionary

        Returns:
            Validated genomic region string, or None if not specified (to trigger auto-inference)

        Raises:
            ValueError: If genomic region is invalid
        """
        # Look for region in multiple possible locations
        region = config.get("region", config.get("genomic_region"))
        
        if region is None:
            # Return None to trigger automatic inference from track data
            return None

        # Validate region format if provided
        try:
            validate_genomic_region(region)
        except ValueError as e:
            raise ValueError(f"Invalid genomic region '{region}': {e}")

        return region

    def create_plot_from_config(
        self, 
        config: Union[str, Dict[str, Any]], 
        backend_name: Optional[str] = None,
        save_path: Optional[str] = None,
        show: bool = False
    ) -> Any:
        """
        Complete workflow: parse config and create plot.

        Args:
            config: Configuration dictionary or path to config file
            backend_name: Backend to use (overrides config if provided)
            save_path: Optional path to save the plot
            show: Whether to display the plot

        Returns:
            Backend instance with the created figure

        Raises:
            ValueError: If configuration is invalid
            FileNotFoundError: If files don't exist
        """
        # Create plotter
        plotter = self.create_plotter(config, backend_name=backend_name)
        
        # Generate plot
        backend = plotter.create_plot()

        # Save if requested
        if save_path:
            backend.save_figure(save_path)

        # Show if requested
        if show:
            backend.show_figure()

        return backend

    def _validate_tracks(self, config: Dict[str, Any]) -> None:
        """
        Validate track configurations in the given dictionary.

        This method checks that all tracks have valid configurations,
        including required fields and correct types.

        Args:
            config: Configuration dictionary containing track definitions

        Raises:
            ValueError: If any track configuration is invalid
        """
        tracks = config.get("tracks", [])
        
        for i, track_config in enumerate(tracks):
            if not isinstance(track_config, dict):
                raise ValueError(f"Track {i} must be a dictionary")
            
            # Check required fields for each track type
            track_type = track_config.get("track_type")
            if track_type == "bed":
                # BED tracks require 'file_path' and have optional 'name' and 'description'
                if "file_path" not in track_config:
                    raise ValueError(f"BED track {i} missing required 'file_path' field")
            elif track_type == "bigwig":
                # BigWig tracks require 'file_path' and have optional 'name' and 'description'
                if "file_path" not in track_config:
                    raise ValueError(f"BigWig track {i} missing required 'file_path' field")
            elif track_type == "gtf" or track_type == "gff":
                # GTF/GFF tracks require 'file_path' and have optional 'name' and 'description'
                if "file_path" not in track_config:
                    raise ValueError(f"GTF/GFF track {i} missing required 'file_path' field")
            else:
                raise ValueError(f"Track {i} has unknown track_type: {track_type}")

            # Auto-detect track type from file extension if not explicitly set
            if "file_path" in track_config and "track_type" not in track_config:
                detected_type = self._auto_detect_track_type(track_config["file_path"])
                if detected_type:
                    track_config["track_type"] = detected_type
                else:
                    raise ValueError(f"Cannot auto-detect track type for file: {track_config['file_path']}. Please specify 'track_type' explicitly.")

            # Validate file path existence if enabled
            if self.validate_files and "file_path" in track_config:
                file_path = track_config["file_path"]
                if file_path and not os.path.exists(file_path):
                    raise FileNotFoundError(f"Track data file not found: {file_path}")

    def get_backend_info(self, config: Union[str, Dict[str, Any]]) -> Tuple[str, Dict[str, Any]]:
        """
        Extract backend information from config.

        Args:
            config: Configuration dictionary or path to config file

        Returns:
            Tuple of (backend_name, backend_options)
        """
        # Load config if it's a file path
        if isinstance(config, str):
            config_dict = self.load_config_file(config)
        else:
            config_dict = config.copy()

        # Parse figure options to get backend info
        figure_options = self._parse_figure_options(config_dict)
        
        backend_name = figure_options.get("backend", self.default_backend)
        backend_options = figure_options.get("backend_options", {})
        
        return backend_name, backend_options

    def create_plotter(
        self, 
        config: Union[str, Dict[str, Any]], 
        backend_name: Optional[str] = None
    ) -> 'GenomicPlotter':
        """
        Parse config and create a GenomicPlotter instance.

        Args:
            config: Configuration dictionary or path to config file
            backend_name: Backend to use (overrides config if provided)

        Returns:
            GenomicPlotter instance ready to create plots

        Raises:
            ValueError: If configuration is invalid
            FileNotFoundError: If files don't exist
        """
        from .plotter import GenomicPlotter
        
        # Parse configuration
        plot_config = self.parse_config(config)
        
        # Get backend info
        if backend_name is None:
            backend_name, backend_options = self.get_backend_info(config)
        else:
            # Use provided backend_name and get options from config
            _, backend_options = self.get_backend_info(config)

        # Create plotter
        return GenomicPlotter(plot_config, backend_name=backend_name, backend_options=backend_options)
