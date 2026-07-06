# FLATSTAT

ImageJ and Python scripts for automated measurement of sample planarity (tilt) from 3D fluorescence microscopy stacks, as described in Brocard, 2026 (https://doi.org/10.64898/2026.05.21.726891).

All raw data available on Zenodo (https://doi.org/10.5281/zenodo.20325392).

## Scripts

- **`FLATSTAT_v8.ijm`** — ImageJ/FIJI macro (Action Tool). Interactive analysis of a single 3D stack: estimates sample tilt (slope, µm/100µm) and direction (compass bearing, 0–360°) from the Z-position of maximum intensity.
- **`idr_flatstats_v6.py`** — Batch analysis of local TIFF stacks using the same algorithm as the macro, adapted for programmatic use.
- **`idr_scraping_full.py`** — Queries the IDR/OMERO API to identify projects and images matching eligibility criteria (imaging mode, stack size, pixel calibration).
- **`idr_flatstats_batch.py`** — Applies FlatStat to the IDR corpus retrieved by `idr_scraping_full.py`, producing the results table used in the manuscript.

## Requirements

- Python 3.9
- `omero-py`, `numpy`, `pandas`
- Conda environment example: `conda create -n omero_env python=3.9 omero-py numpy pandas`

## Usage

```bash
# 1. Scrape IDR/OMERO API and filter eligible projects/images
python idr_scraping_full.py --output idr_census.csv

# 2. Run FlatStat on the retrieved corpus
python idr_flatstats_batch.py --input idr_census.csv --output flatstat_results.csv

# 3. Run FlatStat on local stacks
python idr_flatstats_v6.py --input /path/to/stacks --output flatstat_results_local.csv
```

For the ImageJ macro: install `FLATSTAT_v8.ijm` as an Action Tool (drag into the ImageJ toolbar), open a stack, and click the tool icon.

