# Published candidate results

`CORAL_Earth_search.json` stores the current 75-point Earth/Moon comparison. `CORAL_Earth_candidates.gpkg` stores candidate geometry for GIS inspection. The browser demo embeds copies of the profiles so it can run offline.

## JSON fields

| Field | Meaning |
| --- | --- |
| `input_sources`, `reference_result_sha256` | Input filenames and SHA-256 provenance hashes |
| `distance_m`, `length_m`, `sample_count` | The 75 equal-distance samples over a 2,202.181 m corridor |
| `moon_elevation_m` | Ground elevations for the fixed lunar reference |
| `baseline_elevation_m` | Kilbourne baseline elevations reevaluated on the current Earth grid |
| `candidates` | Three candidates, their profiles, coordinates, source checks and metrics |
| `metrics.rmse_m` | Root mean squared difference of start-relative elevations, in meters |
| `metrics.midpoint_lon_lat` | Candidate midpoint in longitude, latitude order, decimal degrees |
| `metrics.local_position_bearing_diagnostic_rmse_range_m` | Tested nearby position/heading range, not a confidence interval |
| `nominal_sampled_relay_network` | The stored 5 m intermediate-mast network for that profile |
| `moon_nominal_sampled_relay_network` | The stored lunar network with 5 m intermediate masts |

Earth profile coordinates use WGS 84 / UTM zone 13N (EPSG:32613), in meters. Lunar coordinates remain in the source Moon projection. Elevations are meters. Interpret source vertical references through the provider documentation.

## Relay interpretation

Stored networks use 2 m base/rover antennas and 5 m intermediate antennas at the sample positions. `relay_count` excludes the endpoints. Links progress along the corridor in increasing sample order. The browser demo also calculates 10 m and 20 m mast cases.

The weakest checked interior clearance excludes antenna endpoints. An adjacent-sample link has no interior terrain sample and reports `null` clearance. Such a link does not prove the terrain between samples is clear. Ground uncertainty, vegetation and radio propagation remain outside the nominal calculation.

Candidate selection and scoring use the same DEM. Read [METHODS.md](../METHODS.md) before interpreting these exploratory results.
