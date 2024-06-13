    # jdgenometracks

A plotting library for genomic data.

## Installation

'''bash
pip install jdgenometracks
'''

## Usage

Initialize track options:

```python
track_options = {
    "/path/to/data.bed": {
        "track_name": "BedData",
        "main_color": "red",
    },
    "/path/to/data.bedgraph": {
        "track_name": "BedGraphData",
        "main_color": "blue",
        "plot_type": "lines",
    },
}
```
(Possible options found within tracks/{filetype}Track.py)

Create track objects:

```python
import jdgenometracks as jdg
import numpy as np

track_datas = np.array([
    jdg.construct_track(track_path, **extra_options)
    for track_path, extra_options in track_options.items()
])

```

Plot tracks using matplotlib or plotly:

```python
tdp = jdg.MPLPlotter(track_datas, total_height=15)
fig, ax = tdp.plot_all_tracks()
fig.show()
```

```python
tdp = jdg.PlotlyPlotter(track_datas, total_height=15)
fig = tdp.plot_all_tracks()
fig.show()
```

