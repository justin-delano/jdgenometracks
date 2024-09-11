# **jdgenometracks**

## Overview

**jdgenometracks** is a Python package designed to simplify the visualization of genomic data using **Plotly** and **Matplotlib**. The package supports visualizing different track types such as **BED**, **BEDGRAPH**, and **axis** tracks, making it ideal for managing and plotting large genomic datasets with customizable visual styles.

At the core of the package is the **TrackFactory**, which allows users to easily create and manage genomic tracks, enabling the use of both Plotly and Matplotlib with the same set of track definitions.

## Key Features

- **TrackFactory** for easy creation of tracks (e.g., BED, BEDGRAPH, axis, and spacer tracks).
- Supports **Plotly** and **Matplotlib** visualization.
- Customizable plotting options (e.g., colors, markers, lines).
- Easily handle large genomic datasets and visualize them in multiple columns or rows.

## Installation

```bash
pip install jdgenometracks
```

## Getting Started

### Example: Creating and Plotting Tracks with Plotly

```python
import jdgenometracks as jdg

# Define a set of tracks to visualize
tracks = [
    jdg.TrackFactory.create_track(
        file_path="example_data/track1.bedgraph",
        track_name="Track 1",
        track_type="bedgraph",
        plot_type="points",
        plotly_options={
            "marker_color": "blue",
            "marker_size": 4
        }
    ),
    jdg.TrackFactory.create_track(
        file_path="example_data/track2.bed",
        track_name="Track 2",
        track_type="bed",
        plotly_options={
            "fillcolor": "red",
            "line_color": "black"
        }
    ),
    jdg.TrackFactory.create_track(
        track_name="Bottom Axis",
        track_type="axis",
        axis_type="verbose"
    )
]

# Now, pass these tracks to a Plotly plotting utility (example):
plotter = jdg.PlotlyPlotter(tracks, total_height=800)
fig = plotter.plot_all_tracks(column_titles=["Sample Data"])
fig.show()
```

### Example: Creating and Plotting Tracks with Matplotlib

```python
import jdgenometracks as jdg

# Define a set of tracks for Matplotlib
tracks = [
    jdg.TrackFactory.create_track(
        file_path="example_data/track1.bedgraph",
        track_name="Track 1",
        track_type="bedgraph",
        plot_type="lines",
        mpl_plot_options={
            "color": "blue",
            "linewidth": 2
        }
    ),
    jdg.TrackFactory.create_track(
        file_path="example_data/track2.bed",
        track_name="Track 2",
        track_type="bed",
        mpl_rect_options={
            "color": "red",
            "linewidth": 2
        },
        mpl_text_options={
            "va": "center"
        },
        mpl_text_alignment="right"
    ),
    jdg.TrackFactory.create_track(
        track_name="Bottom Axis",
        track_type="axis",
        axis_type="verbose"
    )
]

# Use a Matplotlib plotter
plotter = jdg.MPLPlotter(tracks, total_height=8)
fig, axes = plotter.plot_all_tracks(plot_title="Genomic Data")
plt.show()
```

### Example: Displaying Two Columns of Tracks with Plotly

```python
import jdgenometracks as jdg
import numpy as np

# Define two sets of tracks to display in two columns
tracks = [
    # First Column (Column 1)
    jdg.TrackFactory.create_track(
        file_path="example_data/track1_column1.bedgraph",
        track_name="Track 1 - Column 1",
        track_type="bedgraph",
        plot_type="points",
        plotly_options={
            "marker_color": "blue",
            "marker_size": 4
        }
    ),
    jdg.TrackFactory.create_track(
        file_path="example_data/track2_column1.bed",
        track_name="Track 2 - Column 1",
        track_type="bed",
        plotly_options={
            "fillcolor": "red",
            "line_color": "black"
        }
    ),
    jdg.TrackFactory.create_track(
        track_name="Bottom Axis - Column 1",
        track_type="axis",
        axis_type="verbose"
    ),

    # Second Column (Column 2)
    jdg.TrackFactory.create_track(
        file_path="example_data/track1_column2.bedgraph",
        track_name="Track 1 - Column 2",
        track_type="bedgraph",
        plot_type="lines",
        plotly_options={
            "line_color": "green",
            "line_width": 2
        }
    ),
    jdg.TrackFactory.create_track(
        file_path="example_data/track2_column2.bed",
        track_name="Track 2 - Column 2",
        track_type="bed",
        plotly_options={
            "fillcolor": "purple",
            "line_color": "black"
        }
    ),
    jdg.TrackFactory.create_track(
        track_name="Bottom Axis - Column 2",
        track_type="axis",
        axis_type="verbose"
    )
]

# Reshape the tracks to 2 columns (first 3 tracks in column 1, next 3 in column 2)
tracks = np.array(tracks).reshape(-1, 2)

# Plot using PlotlyPlotter
plotter = jdg.PlotlyPlotter(tracks, total_height=800)
fig = plotter.plot_all_tracks(
    column_titles=["Column 1", "Column 2"],  # Titles for each column
    height_props=[1, 1, 0.1],  # Proportions for each row (tracks and axis)
    width_props=[1, 1]  # Equal width for both columns
)

# Display the figure
fig.show()
```

## TrackFactory Usage

The **TrackFactory** allows you to create different types of tracks. Here are some examples of how to use the factory to create tracks:

### BEDGRAPH Track

```python
jdg.TrackFactory.create_track(
    file_path="path/to/file.bedgraph",
    track_name="Example BEDGRAPH",
    track_type="bedgraph",
    plot_type="lines",
    plotly_options={
        "line_color": "green",
        "line_width": 2
    }
)
```

### BED Track

```python
jdg.TrackFactory.create_track(
    file_path="path/to/file.bed",
    track_name="Example BED",
    track_type="bed",
    plotly_options={
        "fillcolor": "purple",
        "line_width": 2
    }
)
```

### Axis Track

```python
jdg.TrackFactory.create_track(
    track_name="X-Axis",
    track_type="axis",
    axis_type="verbose"
)
```

## Customizing Plots

Each track can be customized using **Plotly** keyword arguments or **Matplotlib** keyword arguments. Example options include:

- **Colors**: Set line and marker colors.
- **Line Styles**: Adjust line widths and types.
- **Markers**: Control marker styles (e.g., size, symbol).
- **Axes**: Customize axis labels and ticks.

## Advanced Usage

For more advanced use cases, you can extend the package to support new track types or customize existing ones by modifying the `TrackFactory` to handle additional parameters or features.

