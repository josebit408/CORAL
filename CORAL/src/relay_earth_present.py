"""Create a display-only hillshade and georeferencing metadata. Requires GDAL."""
from pathlib import Path
import argparse,json
from osgeo import gdal
ap=argparse.ArgumentParser();ap.add_argument('results_dir',nargs='?',default='results');args=ap.parse_args()
out=Path(args.results_dir);gdal.UseExceptions()
ds=gdal.Open(str(out/'Earth_tile_UTM30m.tif'));gt=ds.GetGeoTransform()
meta={'geotransform':gt,'width':ds.RasterXSize,'height':ds.RasterYSize,'epsg':32613,'hillshade':'Display only; azimuth 315 degrees, altitude 45 degrees, z factor 1, horizontal/vertical units meters.'}
(out/'display_grid.json').write_text(json.dumps(meta,indent=2))
hs=gdal.DEMProcessing(str(out/'Earth_tile_hillshade.tif'),ds,'hillshade',format='GTiff',computeEdges=True,azimuth=315,altitude=45,zFactor=1,scale=1,creationOptions=['COMPRESS=DEFLATE','TILED=YES']);hs=None
print('Prepared display-only hillshade and grid metadata.')
