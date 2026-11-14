# CORAL

### Exploring Earth terrain to prepare for lunar communications challenges

CORAL is a geospatial research prototype that identifies Earth terrain corridors resembling a selected lunar elevation profile, then explores relay placement between a base and rover.

Developed for the **2026 NASA International Space Apps Challenge**, the project combines NASA elevation data, Python analysis, QGIS visualization, and an interactive web demonstration.

## The Problem

Terrain can obstruct communication between exploration vehicles and surface infrastructure. Earth locations with similar terrain profiles could provide useful environments for studying these challenges before testing equipment in space.

CORAL investigates a focused question:

**Where on Earth can we find terrain that resembles a selected lunar obstruction profile, and how does assumed antenna height affect geometric relay requirements?**

## What CORAL Does

- Processes Earth and lunar elevation datasets.
- Searches Earth terrain corridors for similarity to a fixed lunar reference.
- Compares three candidate corridors with a Kilbourne Hole baseline.
- Models geometric visibility between sampled antenna positions.
- Calculates relay requirements for different mast heights.
- Presents the results through GIS layers, scientific figures, and an interactive demo.

## Interactive Demo

The current demo compares Earth candidate **RS01** with the lunar reference.

Users can:

- Compare the two terrain profiles.
- Switch intermediate relay mast heights between 5, 10, and 20 meters.
- View the corresponding relay networks.
- Inspect elevation changes along the corridor.

To run locally:

1. Download and extract the repository.
2. Open `docs/index.html` in a web browser.
3. Select a mast height and chart view.

The demo runs without Python or a server. External data-source links require internet access.

For an online version, configure GitHub Pages to publish from the `main` branch and `/docs` folder.

## Methodology

1. Define an illustrative lunar corridor approximately 2.20 kilometers long in the Site04 Shackleton rim elevation model.
2. Aggregate the lunar elevation surface from 5-meter cells to a 30-meter grid.
3. Screen NASADEM source and water masks, then reproject the Earth elevation data into a 30-meter UTM grid.
4. Sample each terrain profile at 75 equally spaced positions.
5. Compare elevations relative to each corridor's starting elevation.
6. Rank Earth candidates using root mean squared error, without stretching terrain heights or corridor distances.
7. Evaluate forward relay links across the sampled terrain, including approximate planetary curvature.

## Preliminary Results

### Terrain Profile Comparison

| Corridor | Profile mismatch, RMSE |
| --- | ---: |
| RS01 | 13.98 m |
| RS02 | 23.08 m |
| RS03 | 23.47 m |
| Kilbourne baseline | 220.59 m |

Lower RMSE indicates a closer match to the selected lunar profile. It does not represent elevation accuracy, a similarity percentage, or mission success probability.

RS01's midpoint is approximately **31.020897° N, 106.542984° W**, in the Sierra La Candelaria area of Chihuahua, Mexico.

### Relay Mast Comparison

| Intermediate mast height | Earth RS01 relays | Lunar reference relays |
| --- | ---: | ---: |
| 5 m | 7 | 7 |
| 10 m | 6 | 5 |
| 20 m | 4 | 4 |

Base and rover antennas remain at 2 meters.

These are minimum intermediate relay counts among the 75 allowed positions, with links progressing forward along the corridor.

## Technology

| Tool | Role |
| --- | --- |
| Python | Terrain processing, candidate search, and relay calculations |
| GDAL | Raster processing and coordinate transformations |
| NumPy | Numerical analysis |
| Matplotlib | Scientific figures |
| QGIS | Spatial inspection and elevation-profile visualization |
| HTML, CSS, JavaScript | Interactive demonstration |

## Repository Layout

| Location | Contents |
| --- | --- |
| `docs/` | Interactive web demo |
| `src/` | Python analysis and visualization scripts |
| `results/` | Calculated metrics, profiles, and GIS layers |
| `figures/` | Generated charts and maps |
| `data/` | Reference data and baseline inputs |
| `SETUP.md` | Setup and reproduction instructions |
| `METHODS.md` | Assumptions and validation limits |
| `DATA_SOURCES.md` | Data provenance and acquisition details |
| `LICENSE` | MIT software license |

Large original elevation datasets must be downloaded separately. See `DATA_SOURCES.md`.

## Data Sources

**Earth:** NASADEM Merged DEM Global 1 arc second V001, tile `n31w107`.

**Moon:** NASA Goddard Planetary Geodynamics Data Archive, Site04 Shackleton rim elevation surface and associated RMS error estimates.

- [NASA Earthdata Search](https://search.earthdata.nasa.gov/)
- [NASADEM User Guide](https://lpdaac.usgs.gov/documents/2237/NASADEM_User_Guide_V13.pdf)
- [Lunar Elevation Data](https://pgda.gsfc.nasa.gov/products/78)

## Limitations

CORAL is an exploratory terrain-obstruction prototype.

- The Earth search currently covers one elevation tile.
- Candidates were selected and scored using the same elevation dataset.
- Independent elevation validation and field testing remain outstanding.
- Terrain uncertainty is not propagated into the displayed relay counts.
- The relay model excludes radio link budgets, Fresnel clearance, diffraction, vegetation, and buildings.
- Features between sampled positions may be missed.
- Geological equivalence and site access remain unverified.
- Profile lines represent antenna links, not rover travel routes.

Matching relay counts alone does not establish equivalent real-world communication performance.

## Planned Improvements

- Validate selected corridors against an independent elevation source.
- Add all three candidates and the baseline to the interactive demo.
- Document candidate geology, land cover, and access constraints.
- Display per-link terrain clearance.
- Evaluate explicit elevation-error scenarios.
- Expand the search to additional regions and lunar reference corridors.

## Author

**Jose Zaragoza**

Developed with AI-assisted coding, analysis, and documentation.

## License

Original software and documentation are available under the **MIT License**.

Third-party datasets retain their providers' applicable terms and attribution requirements. See `LICENSE` and `DATA_SOURCES.md`.
