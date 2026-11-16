# Methods and evidence limits

Profiles use equal physical distances with start-relative elevations. Earth and lunar shape rankings have no vertical scaling or path warping. The Earth search evaluates a coarse grid followed by local refinement and retains three midpoint-separated candidates. Evaluation counts in results are not counts of unique sites.

Earth grids use EPSG:32613 with meters horizontally and vertically. Lunar grids retain their native Moon south polar stereographic parameters and spherical radius of 1,737,400 m; do not assign an Earth CRS. A 30 m cell spacing does not guarantee equal effective resolution or accuracy between source products.

Source-quality screening requires SRTM-origin NUM values 1–23, land in the water mask, valid elevations, a one-cell erosion, and passing bilinear support cells. It is not an independent elevation check or proof of bare-earth terrain.

The relay solver places 2 m endpoint antennas and 5, 10, or 20 m intermediate masts at the 75 sampled positions. Links go forward. Approximate curvature bulge is d(L−d)/(2R). Every sampled position must clear the modeled link. The solver minimizes hops; ties favor larger minimum sampled interior clearance. Adjacent-sample links have no interior samples. Subsample obstacles remain unresolved.

Lunar mean cells average full 6 × 6 source blocks; a trailing 10 m strip is excluded. RMS screening conservatively uses source error values without dividing by the square root of 36 because correlation is unknown. The earlier two-RMS obstruction envelope is a sensitivity scenario, not a confidence interval. This envelope is not propagated into nominal relay counts.

## Completed checks

Analysis work included independent arithmetic checks of profile mismatch, coordinate/support checks, GIS feature counts, and comparison of selected interpolation values with GDAL. The demo's chart calculations were checked for all three heights and views. A full browser visual check of the demo could not be completed in the build environment.

## Outstanding validation

A second suitably independent DEM has not checked the selected corridors. Earth elevation uncertainty is unquantified. Candidate selection and scoring use the same DEM. Nearby lunar lines, Earth line shifts, and mast trials are diagnostics only. Geological equivalence, access permissions, radio performance, and field suitability are unverified.

The demo candidate table reports existing nearby position/heading diagnostics, not confidence intervals. Sampled network clearances exclude antenna endpoints. Links between adjacent allowed positions have no interior sample; the demo reports their count explicitly and exports null for their minimum interior clearance. The interface retains the original fixed physical plot axes across candidates.
