"""Create a display-only hillshade and georeferencing metadata using GDAL."""
from pathlib import Path
import argparse
import json


def prepare_display(results_dir):
    out = Path(results_dir)
    source = out / 'Earth_tile_UTM30m.tif'
    if not source.is_file():
        raise ValueError(f'Missing input raster: {source}. Run relay_earth_search.py first.')
    try:
        from osgeo import gdal, osr
    except ImportError as error:
        raise ValueError('GDAL Python bindings are missing. Use the CORAL GIS environment.') from error
    gdal.UseExceptions()
    ds = gdal.Open(str(source))
    if ds is None or ds.RasterCount != 1:
        raise ValueError('Expected a single-band Earth elevation raster.')
    srs = osr.SpatialReference(wkt=ds.GetProjection())
    expected = osr.SpatialReference()
    expected.ImportFromEPSG(32613)
    if not srs.IsSame(expected):
        raise ValueError('Expected the Earth search raster in EPSG:32613.')
    gt = ds.GetGeoTransform()
    if gt[1] != 30 or gt[5] != -30 or gt[2] != 0 or gt[4] != 0:
        raise ValueError('Expected the aligned 30 m Earth search grid.')
    meta = {
        'geotransform': gt, 'width': ds.RasterXSize, 'height': ds.RasterYSize,
        'epsg': 32613,
        'hillshade': 'Display only; azimuth 315 degrees, altitude 45 degrees, '
                     'z factor 1, horizontal/vertical units meters.',
    }
    hs = gdal.DEMProcessing(
        str(out / 'Earth_tile_hillshade.tif'), ds, 'hillshade', format='GTiff',
        computeEdges=True, azimuth=315, altitude=45, zFactor=1, scale=1,
        creationOptions=['COMPRESS=DEFLATE', 'TILED=YES'])
    if hs is None:
        raise RuntimeError('GDAL did not produce a hillshade.')
    hs = None
    (out / 'display_grid.json').write_text(json.dumps(meta, indent=2), encoding='utf-8')
    print('Prepared display-only hillshade and grid metadata.')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results_dir', nargs='?', type=Path, default=Path('generated/earth'))
    args = parser.parse_args(argv)
    try:
        prepare_display(args.results_dir)
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(2, f'CORAL display error: {error}\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
