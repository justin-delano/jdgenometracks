from __future__ import annotations

import pandas as pd

from .tracks.BedGraphTrack import BedGraphTrack
from .tracks.BedTrack import BedTrack
from .tracks.GenomeTrack import GenomeTrack
from .tracks.HelperTracks import SpacerTrack, XAxisTrack


def construct_track(file_path: str, **kwargs) -> GenomeTrack:
    """_summary_

    Args:
        file_path (str): Path to file containing data to be plotted in a track. If file path contains "spacer" or "axis", a SpacerTrack or XAxisTrack will be returned, respectively.

    Raises:
        NotImplementedError: Currently must be a .bed or .bedgraph file.

    Returns:
        GenomeTrack: A GenomeTrack object containing the data to be plotted.
    """
    dtypes = {
        "chrom": str,
        "chromStart": int,
        "chromEnd": int,
        "name": str,
        "score": float,
        "strand": str,
        "thickStart": int,
        "thickEnd": int,
        "itemRGB": str,
        "value": float,
    }

    track_name = kwargs.pop("track_name", file_path.split("/")[-1].split(".")[0])
    track_type = kwargs.pop("track_type", file_path.split(".")[-1])

    if "spacer" in file_path.lower():
        return SpacerTrack("Spacer", "Spacer", track_name, pd.DataFrame(), **kwargs)

    if "axis" in file_path.lower():
        return XAxisTrack("Axis", "Axis", track_name, pd.DataFrame(), **kwargs)

    data = pd.read_csv(file_path, sep="\t", header=None)

    if track_type == "bed":
        all_bed_columns = [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRGB",
        ]
        data.columns = all_bed_columns[: len(data.columns)]
        for column in data.columns:
            data[column].astype(dtypes[column])

        return BedTrack(file_path, track_type, track_name, data, **kwargs)
    elif track_type == "bedgraph":
        all_bedgraph_columns = ["chrom", "chromStart", "chromEnd", "value", "name"]
        data = data.iloc[:, : len(all_bedgraph_columns)]
        data.columns = all_bedgraph_columns[: len(data.columns)]
        for column in data.columns:
            data[column].astype(dtypes[column])

        return BedGraphTrack(file_path, track_type, track_name, data, **kwargs)
    else:
        raise NotImplementedError("Track type not yet supported")


def _get_height_props(tracks) -> list[float]:
    new_heights = []
    for track in tracks:
        if track.share_with_previous:
            continue
        elif track.height_prop is not None:
            new_heights.append(track.height_prop)
        else:
            new_heights.append(1)

    return new_heights


def _get_xlim_bedlim(tracks, column_region: str | None) -> tuple[int, int, int]:
    max_bed_regions = 0
    for idx, track in enumerate(tracks):
        if track.data is None or track.data.empty:
            continue

        formatted_data = track.format_data(subset_region=column_region)
        if formatted_data.empty:
            continue

        if idx == 0:
            xmin = formatted_data["chromStart"].min()
            xmax = formatted_data["chromEnd"].max()
        else:
            xmin = min(xmin, formatted_data["chromStart"].min())
            xmax = max(xmax, formatted_data["chromEnd"].max())

        if isinstance(track, BedTrack):
            max_bed_regions = max(max_bed_regions, formatted_data.shape[0])

    return xmin, xmax, max_bed_regions


def _get_col_limits(
    tracks_col, column_region: str | None, relative_x_axis: bool
) -> tuple[int, int, dict]:
    extra_options = {}
    if column_region is not None:
        chromosome = column_region.split(":")[0]
        extra_options["subset_region"] = column_region
    else:
        chromosome = tracks_col[0].data.iloc[0, 0]

    extra_options["chromosome"] = chromosome
    extra_options["max_regions"] = 0

    xmin, xmax, max_bed_regions = _get_xlim_bedlim(tracks_col, column_region)

    extra_options["max_regions"] = max_bed_regions
    if relative_x_axis and xmin and xmax:
        xmin -= xmin
        xmax -= xmin
        extra_options["axis_shift"] = xmin

    return xmin, xmax, extra_options
