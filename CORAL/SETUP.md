# Setup and publishing

## Run the demo

Download the repository, extract it, and double-click `docs/index.html`. A browser runs the entire demo without Python or a server.

## Upload through GitHub in a browser

1. Sign in to GitHub and create a repository named `CORAL`.
2. Suggested description: “NASA-data terrain analog search and geometric relay-planning prototype comparing Earth corridors with a lunar reference.”
3. Choose Public if you want judges and portfolio reviewers to access it. Leave automatic README, license, and .gitignore creation off: these files are included here.
4. Choose “uploading an existing file” on the empty repository page (or Add file → Upload files).
5. Open the extracted CORAL folder. Drag its CONTENTS, including the folders, into the upload area. Upload the contents, not the ZIP and not an extra enclosing CORAL folder.
6. Commit with the message `Add CORAL terrain analysis and interactive demo`.
7. Verify that README.md and the docs, src, results, figures, and data folders appear at the repository root.

## Enable the live demo with GitHub Pages

Repository Settings → Pages → Deploy from a branch → main → /docs → Save. Use the published URL shown by GitHub when deployment completes. Site availability depends on your account and repository settings.

## Reproduce the Earth search

Use a Python environment that provides GDAL's `osgeo` bindings and NumPy. For figures, also install Matplotlib and Pillow. GDAL's Python bindings must match the installed GDAL native library; a GIS/conda environment is a convenient option. Exact package versions are not locked in this prototype.

Download the NASADEM n31w107 archive and place `n31w107.hgt`, `n31w107.num`, and `n31w107.swb` in `data/inputs/`. The supplied baseline geometry and fixed lunar reference are already included.

Run from the repository root:

```bash
python src/relay_earth_search.py --input-dir data/inputs --moon-results data/reference_moon --output-dir generated/earth
python src/relay_earth_present.py generated/earth
python src/relay_earth_plot.py generated/earth
```

The search can take time. Generated rasters are excluded from version control. The interactive demo embeds selected profile data and does not automatically rebuild itself when scripts run.

## Rebuild the lunar reference

Download the two native Site04 files named in DATA_SOURCES.md into `data/inputs/`. The lunar comparison also requires the exported `Kilbourne_profile01.csv` with its original field names. If absent from this package, obtain the original profile export before running this optional stage.

```bash
python src/relay_moon_comparison.py --input-dir data/inputs --output-dir generated/moon
python src/relay_comparison_plot.py generated/moon
```

The precomputed fixed lunar reference is included so this optional stage is not required for the Earth search.

## Improved Moon comparison commands

The Moon comparison scripts validate the native Site04 grid, lunar CRS, RMS values, and CSV columns before writing final outputs. Default heights reproduce the original trial results. Optional heights are in meters:

```bash
python src/relay_moon_comparison.py --input-dir data/inputs --output-dir generated/moon --base-height 2 --rover-height 2 --relay-height 10
python src/relay_comparison_plot.py generated/moon --dpi 200 --output generated/moon/comparison.png
```

These settings affect the single sampled-peak relay trial, not the Earth search or the web demo's multi-relay network. Use a new output folder for each experiment to retain previous results. Successful runs replace matching output files; invalid inputs leave existing final outputs untouched. File publication is atomic per file, not across the entire output bundle.

The CSV requires `layer`, `distance`, and `elevation` columns, one nonempty layer, finite elevations, and strictly increasing distances beginning at zero. UTF-8 files with a byte-order mark are accepted. An endpoint peak has no intermediate peak-relay trial; the corresponding metrics are JSON null. Profiles too short for the centered 90 m slope calculation have a null slope statistic.

### Run automated checks

The tests require Python and NumPy; the core tests do not require GDAL or Matplotlib:

```bash
python -m unittest discover -s tests -v
```

The updated analysis was also run on the original source rasters and CSV. All 2,176 existing numeric comparison values matched exactly. Both derived raster arrays matched the previous products. Custom antenna heights were exercised in a separate output folder. These checks establish implementation consistency, not independent terrain validation.

The web demo now includes all three candidates and the Kilbourne baseline. The optional `tests/test_demo.py` check uses Node.js to exercise 36 candidate/height/view combinations and compares relay counts and exported link clearances with a separate Python calculation. It is skipped if Node.js is unavailable. These are numeric and rendering-logic checks; a full browser visual check remains outstanding.
