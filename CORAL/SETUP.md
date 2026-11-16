# Setup and reproduction

## Open the demo

Visit **https://josebit408.github.io/CORAL/**. To run offline, extract a downloaded copy of this repository and open `docs/index.html` in a browser. The demo needs no Python installation, web server or account. Source-data links require internet access.

## Update the existing GitHub repository

1. Extract `CORAL_GitHub_Ready.zip` on your computer.
2. Open https://github.com/josebit408/CORAL and choose **Add file → Upload files**.
3. Open the extracted `CORAL` folder and drag its **contents** into the upload area. Keep the folder structure. Upload the files and folders, rather than the ZIP or an extra enclosing `CORAL` folder.
4. Commit to `main` with a message such as `Update CORAL documentation, checks, and research outputs`.
5. Verify that `README.md`, `docs`, `src`, `tests`, `results`, `figures` and `data` appear at the repository root.

Files with a leading dot may be hidden in your file browser. They include `.gitignore`, `.github/workflows/tests.yml` and `docs/.nojekyll`. Enable hidden files when uploading. If needed, add `.github/workflows/tests.yml` through **Add file → Create new file** using that exact path.

The repository uses GitHub Pages with branch `main` and folder `/docs`. Check **Settings → Pages** after uploading. A successful deployment publishes the demo at the address above. The repository's **About → Website** field can use the same URL.

## Install the analysis environment

The browser demo needs no dependencies. For the Python analysis, use a Conda-compatible installation with the `conda-forge` channel:

```bash
conda env create -f environment.yml
conda activate coral
```

The environment specification includes Python, NumPy, GDAL, Matplotlib, Pillow and Node.js. It specifies compatible version ranges rather than an exact lockfile. A fresh Conda environment solve has not been verified on every platform. GDAL's Python bindings must match its native library, so install them together through the GIS environment.

Confirm the raster tools are available:

```bash
python -c "from osgeo import gdal; import numpy, matplotlib, PIL; print(gdal.VersionInfo('--version'))"
```

## Run the automated checks

For a lightweight setup without GDAL:

```bash
python -m pip install -r requirements-test.txt
python -m unittest discover -s tests -v
```

Use Python 3.10 or newer. Install Node.js 20 or newer to include the demo checks. The tests explicitly skip demo checks if Node.js is unavailable. The GitHub workflow installs both dependencies and runs the same command.

These checks cover numerical consistency and input handling. They do not establish independent terrain accuracy or field suitability.

## Reproduce the Earth search

Run all commands below from the repository root. Download the NASADEM `n31w107` archive using [DATA_SOURCES.md](DATA_SOURCES.md), then put these extracted files in `data/inputs/`:

- `n31w107.hgt`
- `n31w107.num`
- `n31w107.swb`

The package already includes the Kilbourne baseline geometry and prepared lunar reference.

```bash
python src/relay_earth_search.py --input-dir data/inputs --moon-results data/reference_moon --output-dir generated/earth
python src/relay_earth_present.py generated/earth
python src/relay_earth_plot.py generated/earth
```

The search processes the full Earth tile and may take several minutes or longer, depending on the computer. It generates the 30 m Earth raster and candidate results. The presentation step creates a display-only hillshade and grid metadata. The plotting step needs those generated files to make maps. Bundled `results/` contains compact results, without the large Earth raster and hillshade.

To remake only the original Earth/Moon comparison figure from bundled data:

```bash
python src/relay_comparison_plot.py data/reference_moon --output generated/moon/comparison.png
```

This figure describes the original Kilbourne trial, which has 83 samples. The later Earth search and web demo use 75 samples. Their relay summaries describe different calculations. See [data/README.md](data/README.md).

## Rebuild the lunar reference

This optional stage requires the original Site04 elevation and RMS-error rasters listed in [DATA_SOURCES.md](DATA_SOURCES.md). Put them in `data/inputs/`. The original `Kilbourne_profile01.csv` is included.

```bash
python src/relay_moon_comparison.py --input-dir data/inputs --output-dir generated/moon
python src/relay_comparison_plot.py generated/moon
```

The scripts validate the native grid, lunar coordinate system, RMS values and CSV columns. The lunar surface stays in its original Moon projection.

Optional antenna heights, in meters:

```bash
python src/relay_moon_comparison.py --input-dir data/inputs --output-dir generated/moon-10m --base-height 2 --rover-height 2 --relay-height 10
```

These options change the original single-peak relay trial. They do not change the later Earth search or the demo's multiple-relay calculation. Use a separate output folder for each experiment. Successful runs may replace generated files in the chosen folder. Keep the packaged reference data unchanged when reproducing the published search.

The CSV must contain `layer`, `distance` and `elevation` columns, one nonempty layer, finite elevations, and strictly increasing distances beginning at zero. UTF-8 files with a byte-order mark are accepted. The scripts report undefined endpoint-peak trials and short-profile slope statistics as JSON `null`.

## Published data and demo updates

The HTML embeds its profile arrays so it works offline. Running Python does not rewrite `docs/index.html`. If the analysis changes, review the results and update the embedded profiles, evidence metrics and documentation together. Run the tests before publishing. Preserve the fixed reference when comparing versions.
