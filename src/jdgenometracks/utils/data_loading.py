"""
Data loading and preprocessing utilities for jdgenometracks V3.

Handles reading and preprocessing of genomic data files including:
- BED format files
- BedGraph format files
- Data validation and type conversion
- Empty file handling
"""

import os
from typing import Any, Dict, List, Optional

import pandas as pd

from .constants import ErrorMessages, FileConstants


class DataLoader:
    """Handles loading and preprocessing of genomic data files."""

    @staticmethod
    def load_bed_file(file_path: str) -> pd.DataFrame:
        """
        Load a BED format file.

        Args:
            file_path: Path to the BED file

        Returns:
            DataFrame with standardized BED columns

        Raises:
            FileNotFoundError: If file doesn't exist
            ValueError: If file format is invalid
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

        try:
            # Read tab-separated file without headers
            data = pd.read_csv(file_path, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            # Return empty DataFrame with proper BED columns
            return DataLoader._create_empty_bed_dataframe()
        except Exception as e:
            raise ValueError(f"Error reading BED file {file_path}: {e}")

        # Assign column names based on number of columns
        num_cols = len(data.columns)
        if num_cols > len(FileConstants.BED_COLUMNS):
            raise ValueError(f"BED file has too many columns: {num_cols}")

        # Use appropriate number of BED columns
        data.columns = FileConstants.BED_COLUMNS[:num_cols]

        # Apply data type conversions
        data = DataLoader._apply_bed_dtypes(data)

        return data

    @staticmethod
    def load_bedgraph_file(file_path: str) -> pd.DataFrame:
        """
        Load a BedGraph format file.

        Args:
            file_path: Path to the BedGraph file

        Returns:
            DataFrame with standardized BedGraph columns

        Raises:
            FileNotFoundError: If file doesn't exist
            ValueError: If file format is invalid
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

        try:
            # Read tab-separated file without headers
            data = pd.read_csv(file_path, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            # Return empty DataFrame with proper BedGraph columns
            return DataLoader._create_empty_bedgraph_dataframe()
        except Exception as e:
            raise ValueError(f"Error reading BedGraph file {file_path}: {e}")

        # BedGraph should have exactly 4 columns, but we support 5 with optional name
        num_cols = len(data.columns)
        if num_cols < 4:
            raise ValueError(
                f"BedGraph file must have at least 4 columns, got {num_cols}"
            )
        elif num_cols > 5:
            raise ValueError(f"BedGraph file has too many columns: {num_cols}")

        # Use appropriate number of BedGraph columns
        data.columns = FileConstants.BEDGRAPH_COLUMNS[:num_cols]

        # Apply data type conversions
        data = DataLoader._apply_bedgraph_dtypes(data)

        return data

    @staticmethod
    def load_file_by_type(file_path: str, file_type: str) -> pd.DataFrame:
        """
        Load a file based on its type.

        Args:
            file_path: Path to the file
            file_type: Type of file ("bed" or "bedgraph")

        Returns:
            Loaded and preprocessed DataFrame
        """
        if file_type == "bed":
            return DataLoader.load_bed_file(file_path)
        elif file_type == "bedgraph":
            return DataLoader.load_bedgraph_file(file_path)
        else:
            raise ValueError(f"Unsupported file type: {file_type}")

    @staticmethod
    def infer_file_type(file_path: str) -> str:
        """
        Infer file type from file extension.

        Args:
            file_path: Path to the file

        Returns:
            Inferred file type ("bed" or "bedgraph")

        Raises:
            ValueError: If file type cannot be inferred
        """
        ext = os.path.splitext(file_path)[1].lower()

        if ext in FileConstants.SUPPORTED_EXTENSIONS:
            return FileConstants.SUPPORTED_EXTENSIONS[ext]
        else:
            raise ValueError(f"Cannot infer file type from extension: {ext}")

    @staticmethod
    def _create_empty_bed_dataframe() -> pd.DataFrame:
        """Create an empty DataFrame with BED columns."""
        return pd.DataFrame(
            columns=FileConstants.BED_COLUMNS[:5]
        )  # Minimum BED columns

    @staticmethod
    def _create_empty_bedgraph_dataframe() -> pd.DataFrame:
        """Create an empty DataFrame with BedGraph columns."""
        return pd.DataFrame(
            columns=FileConstants.BEDGRAPH_COLUMNS[:4]
        )  # Standard BedGraph columns

    @staticmethod
    def _apply_bed_dtypes(data: pd.DataFrame) -> pd.DataFrame:
        """Apply appropriate data types to BED DataFrame columns."""
        for col in data.columns:
            if col in FileConstants.BED_COLUMN_DTYPES:
                try:
                    # Handle missing values before type conversion
                    if col in [
                        "chromStart",
                        "chromEnd",
                        "thickStart",
                        "thickEnd",
                        "blockCount",
                    ]:
                        # Integer columns - fill missing with 0
                        data[col] = (
                            data[col]
                            .fillna(0)
                            .astype(FileConstants.BED_COLUMN_DTYPES[col])
                        )
                    elif col == "score":
                        # Float column - fill missing with 0.0
                        data[col] = (
                            data[col]
                            .fillna(0.0)
                            .astype(FileConstants.BED_COLUMN_DTYPES[col])
                        )
                    else:
                        # String columns - fill missing with empty string
                        data[col] = (
                            data[col]
                            .fillna("")
                            .astype(FileConstants.BED_COLUMN_DTYPES[col])
                        )
                except (ValueError, TypeError) as e:
                    # If conversion fails, keep as-is and warn
                    print(
                        f"Warning: Could not convert column {col} to {FileConstants.BED_COLUMN_DTYPES[col]}: {e}"
                    )

        return data

    @staticmethod
    def _apply_bedgraph_dtypes(data: pd.DataFrame) -> pd.DataFrame:
        """Apply appropriate data types to BedGraph DataFrame columns."""
        for col in data.columns:
            if col in FileConstants.BEDGRAPH_COLUMN_DTYPES:
                try:
                    if col in ["chromStart", "chromEnd"]:
                        # Integer columns - fill missing with 0
                        data[col] = (
                            data[col]
                            .fillna(0)
                            .astype(FileConstants.BEDGRAPH_COLUMN_DTYPES[col])
                        )
                    elif col == "value":
                        # Float column - fill missing with 0.0
                        data[col] = (
                            data[col]
                            .fillna(0.0)
                            .astype(FileConstants.BEDGRAPH_COLUMN_DTYPES[col])
                        )
                    else:
                        # String columns - fill missing with empty string
                        data[col] = (
                            data[col]
                            .fillna("")
                            .astype(FileConstants.BEDGRAPH_COLUMN_DTYPES[col])
                        )
                except (ValueError, TypeError) as e:
                    print(
                        f"Warning: Could not convert column {col} to {FileConstants.BEDGRAPH_COLUMN_DTYPES[col]}: {e}"
                    )

        return data


class DataPreprocessor:
    """Handles preprocessing of loaded genomic data."""

    @staticmethod
    def filter_by_region(
        data: pd.DataFrame, chromosome: str, start: int, end: int
    ) -> pd.DataFrame:
        """
        Filter data to a specific genomic region.

        Args:
            data: Input DataFrame
            chromosome: Target chromosome
            start: Start coordinate
            end: End coordinate

        Returns:
            Filtered DataFrame
        """
        if data.empty:
            return data

        # Filter by chromosome and coordinate overlap
        mask = (
            (data["chrom"] == chromosome)
            & (data["chromStart"] < end)
            & (data["chromEnd"] > start)
        )

        return data[mask].copy()

    @staticmethod
    def add_missing_columns(data: pd.DataFrame, data_type: str) -> pd.DataFrame:
        """
        Add missing columns with default values.

        Args:
            data: Input DataFrame
            data_type: Type of data ("bed" or "bedgraph")

        Returns:
            DataFrame with all expected columns
        """
        if data_type == "bed":
            expected_cols = FileConstants.BED_COLUMNS[:5]  # Minimum BED columns
            dtypes = FileConstants.BED_COLUMN_DTYPES
        elif data_type == "bedgraph":
            expected_cols = FileConstants.BEDGRAPH_COLUMNS
            dtypes = FileConstants.BEDGRAPH_COLUMN_DTYPES
        else:
            raise ValueError(f"Unsupported data type: {data_type}")

        # Add missing columns with appropriate defaults
        for col in expected_cols:
            if col not in data.columns:
                if col in dtypes:
                    dtype = dtypes[col]
                    if dtype == str:
                        default_val = ""
                    elif dtype == int:
                        default_val = 0
                    elif dtype == float:
                        default_val = 0.0
                    else:
                        default_val = None
                else:
                    default_val = None

                data[col] = default_val

        return data


def load_and_preprocess_file(
    file_path: str,
    file_type: Optional[str] = None,
    region_filter: Optional[Dict[str, Any]] = None,
) -> pd.DataFrame:
    """
    High-level function to load and preprocess a genomic data file.

    Args:
        file_path: Path to the file
        file_type: Type of file ("bed" or "bedgraph"), or None to infer
        region_filter: Optional dict with 'chromosome', 'start', 'end' for filtering

    Returns:
        Loaded and preprocessed DataFrame
    """
    # Infer file type if not provided
    if file_type is None:
        file_type = DataLoader.infer_file_type(file_path)

    # Load the file
    data = DataLoader.load_file_by_type(file_path, file_type)

    # Validate the data
    if file_type == "bed":
        from .validation import validate_bed_data

        validate_bed_data(data)
    elif file_type == "bedgraph":
        from .validation import validate_bedgraph_data

        validate_bedgraph_data(data)

    # Add missing columns if needed
    data = DataPreprocessor.add_missing_columns(data, file_type)

    # Apply region filter if provided
    if region_filter and not data.empty:
        data = DataPreprocessor.filter_by_region(
            data,
            region_filter["chromosome"],
            region_filter["start"],
            region_filter["end"],
        )

    return data
