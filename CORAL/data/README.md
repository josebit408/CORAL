# Packaged data

The original NASADEM tile and native 5 m lunar rasters require separate downloads. See [DATA_SOURCES.md](../DATA_SOURCES.md) for their providers and exact filenames.

| Path | Contents | Units and coordinates |
| --- | --- | --- |
| `inputs/Kilbourne_profile01.csv` | Original QGIS trial profile, 83 sample positions | Distance along the line and elevation in meters |
| `inputs/Kilbourne_profile01.gpkg` | Project-created baseline profile geometry | Inspect the layer CRS in QGIS |
| `reference_moon/Site04_30m_mean.tif` | Complete 6 × 6 block means of the 5 m Site04 surface | Meters, native lunar south polar stereographic projection |
| `reference_moon/Site04_30m_rms_upper.tif` | Maximum native error value in each complete block | Meters, same grid as the mean surface |
| `reference_moon/CORAL_comparison.json` | Original fixed-transect comparison with the 83-point Kilbourne export | Distances and elevations in meters |

The later Earth search resamples the fixed physical corridor at 75 equal-distance positions. Those profiles and candidate scores are in [`results/CORAL_Earth_search.json`](../results/CORAL_Earth_search.json). Use that file for the current candidate comparison and browser demo. Do not mix the original 83-point trial's single-peak relay diagnostics with the current sampled-network relay counts.

The lunar RMS grid is a conservative error proxy. It is neither a confidence interval nor an independent terrain check. A common 30 m grid spacing does not establish matching Earth and lunar data accuracy. Native data retain their source vertical references. The ranking compares changes relative to each corridor's starting elevation.

Keep original filenames and packaged reference files unchanged if reproducing the recorded search. Input hashes in the JSON results document provenance. These files contain derived data and project-created geometry, with source attribution retained separately from the MIT software license.
