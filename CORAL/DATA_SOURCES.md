# Data sources and provenance

## Earth

NASADEM Merged DEM Global 1 arc second V001 (NASADEM_HGT.001), tile n31w107. Source files: n31w107.hgt, n31w107.num, n31w107.swb. NASA / JPL / LP DAAC.

- Discovery: https://search.earthdata.nasa.gov/ (search NASADEM_HGT and n31w107; download may require Earthdata login)
- Dataset DOI: https://doi.org/10.5067/MEaSUREs/NASADEM/NASADEM_HGT.001
- User guide: https://lpdaac.usgs.gov/documents/2237/NASADEM_User_Guide_V13.pdf

NUM provides source/count metadata, not elevation uncertainty. SWB is a water mask. Original raster files are not bundled. The processed search JSON contains SHA-256 hashes for the search inputs. Baseline geometry is a project-created Kilbourne Hole trial profile.

Download and extract all three files from the same tile archive. The search expects the exact filenames above in `data/inputs/`. The analysis covers this single tile, bounded by 31–32° N and 107–106° W.

## Moon

NASA Goddard Planetary Geodynamics Data Archive, product 78, Site04 Shackleton rim elevation model:

https://pgda.gsfc.nasa.gov/products/78

Native inputs: Site04_final_adj_5mpp_surf.tif and Site04_final_adj_5mpp_toterr.tif. The latter supplies elevation RMS error estimates. Consult the provider page for its source publications, frame, and full attribution. Derived lunar mean and conservative RMS grids are included; the original 5 m rasters are not bundled.

The prepared lunar reference uses the native south polar stereographic grid and a Moon radius of 1,737,400 m. The workflow compares elevation changes from each corridor's start rather than treating Earth and lunar absolute elevation references as interchangeable. Packaged-data descriptions are in [data/README.md](data/README.md).

## Figures and imagery

Bundled figures are generated from the project's elevation analyses. The candidate map uses a DEM hillshade; it is not satellite photography. Satellite basemaps shown during QGIS exploration are not bundled. Any future exported basemap must carry its provider's required attribution.

## License scope

The MIT software license covers the project's original code and documentation. It does not replace third-party data terms or imply NASA endorsement. Retain source provenance when reusing derived data.
