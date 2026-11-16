"""Compare fixed Site04 lunar transects with a supplied Earth elevation profile.

Requires NumPy and GDAL. Default settings reproduce the original CORAL trial.
Coordinates and ground heights are meters. Relay links are geometric trials;
RMS envelopes are sensitivity scenarios, not statistical confidence intervals.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import tempfile

import numpy as np

EARTH_RADIUS_M = 6_371_000.0
MOON_RADIUS_M = 1_737_400.0
NATIVE_TRANSFORM = (-9000.0, 5.0, 0.0, 1000.0, 0.0, -5.0)
COARSE_TRANSFORM = (-9000.0, 30.0, 0.0, 1000.0, 0.0, -30.0)
TRANSECT_Y_M = (-6000.0, -6500.0, -7000.0)
PRIMARY_INDEX = 1


def require(condition, message):
    """Validate even when Python is run with optimization (-O)."""
    if not condition:
        raise ValueError(message)


def finite_vector(values, name):
    array = np.asarray(values, dtype=np.float64)
    require(array.ndim == 1 and array.size >= 2,
            f'{name} must be a one-dimensional array with at least two samples.')
    require(np.isfinite(array).all(), f'{name} contains a missing or non-finite value.')
    return array


def validate_profile(distance, elevation):
    distance = finite_vector(distance, 'Distance')
    elevation = finite_vector(elevation, 'Elevation')
    require(distance.shape == elevation.shape, 'Distance and elevation lengths differ.')
    require(distance[0] == 0, 'The first distance must be 0 meters.')
    require((np.diff(distance) > 0).all(),
            'Distances must increase strictly; duplicate or reversed samples are not allowed.')
    return distance, elevation


def load_earth_profile(path):
    """Read the original single-layer QGIS export, including UTF-8 BOM files."""
    distance, elevation, layers = [], [], set()
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream)
        missing = {'distance', 'elevation', 'layer'} - set(reader.fieldnames or [])
        require(not missing, f'{path}: missing CSV columns: {", ".join(sorted(missing))}.')
        for row in reader:
            try:
                distance.append(float(row['distance']))
                elevation.append(float(row['elevation']))
                layer = (row['layer'] or '').strip()
                require(bool(layer), 'layer must not be empty')
                layers.add(layer)
            except (TypeError, ValueError) as error:
                raise ValueError(f'{path}: invalid CSV row {reader.line_num}: {error}') from error
    require(len(layers) == 1, f'{path}: expected exactly one nonempty profile layer.')
    return validate_profile(distance, elevation)


def load_gdal():
    """Keep numeric functions importable without the GIS dependencies."""
    try:
        from osgeo import gdal, ogr, osr
    except ImportError as error:
        raise ValueError('GDAL Python bindings are missing. Run this command in a '
                         'Python environment containing GDAL (osgeo) and NumPy.') from error
    gdal.UseExceptions()
    ogr.UseExceptions()
    osr.UseExceptions()
    return gdal, ogr, osr


def load_lunar_surface(surface_path, rms_path, gdal, osr):
    datasets = []
    for path in (surface_path, rms_path):
        require(path.is_file(), f'Input file does not exist: {path}')
        dataset = gdal.Open(str(path))
        require(dataset is not None, f'Cannot open raster: {path}')
        require(dataset.RasterCount == 1, f'{path.name}: expected a single raster band.')
        require((dataset.RasterXSize, dataset.RasterYSize) == (3200, 3200),
                f'{path.name}: expected the native 3200 × 3200 Site04 raster.')
        require(np.allclose(dataset.GetGeoTransform(), NATIVE_TRANSFORM, rtol=0, atol=1e-8),
                f'{path.name}: unexpected grid extent or pixel size; use native Site04 data.')
        datasets.append(dataset)
    surface, rms = datasets
    require(surface.GetGeoTransform() == rms.GetGeoTransform(),
            'Surface and RMS raster grids are not aligned.')
    wkt = surface.GetProjection()
    require(bool(wkt) and bool(rms.GetProjection()), 'Both rasters must include their lunar CRS.')
    srs = osr.SpatialReference(wkt=wkt)
    other_srs = osr.SpatialReference(wkt=rms.GetProjection())
    require(bool(srs.IsSame(other_srs)), 'Surface and RMS coordinate systems differ.')
    require(bool(srs.IsProjected()) and abs(srs.GetSemiMajor() - MOON_RADIUS_M) < 0.01
            and abs(srs.GetSemiMinor() - MOON_RADIUS_M) < 0.01
            and abs(srs.GetLinearUnits() - 1.0) < 1e-12
            and srs.GetAttrValue('PROJECTION') == 'Polar_Stereographic'
            and abs(srs.GetProjParm('latitude_of_origin') + 90) < 1e-8
            and abs(srs.GetProjParm('central_meridian')) < 1e-8
            and abs(srs.GetProjParm('scale_factor', 1) - 1) < 1e-8
            and abs(srs.GetProjParm('false_easting')) < 1e-8
            and abs(srs.GetProjParm('false_northing')) < 1e-8,
            'Expected the native Moon south polar stereographic CRS in meters.')
    arrays = []
    for dataset, path in zip(datasets, (surface_path, rms_path)):
        band = dataset.GetRasterBand(1)
        array = band.ReadAsArray()
        require(array is not None and np.isfinite(array).all(),
                f'{path.name}: raster contains missing or non-finite heights.')
        nodata = band.GetNoDataValue()
        if nodata is not None and math.isfinite(nodata):
            require(not np.any(array == nodata), f'{path.name}: raster contains NoData cells.')
        require(np.all(band.GetMaskBand().ReadAsArray() != 0),
                f'{path.name}: raster contains masked cells.')
        arrays.append(array)
    require((arrays[1] >= 0).all(), 'RMS values must be nonnegative.')
    return arrays[0], arrays[1], srs, wkt


def aggregate_surface(surface, rms, factor=6):
    """Mean complete blocks; preserve a conservative maximum RMS proxy."""
    require(isinstance(factor, int) and factor > 0, 'Aggregation factor must be a positive integer.')
    surface, rms = np.asarray(surface), np.asarray(rms)
    require(surface.ndim == 2 and surface.shape == rms.shape, 'Surface/RMS array shapes differ.')
    require(np.isfinite(surface).all() and np.isfinite(rms).all() and (rms >= 0).all(),
            'Surface/RMS arrays must be finite and RMS nonnegative.')
    rows, cols = (size // factor for size in surface.shape)
    require(rows >= 2 and cols >= 2, 'Raster is too small for aggregation and bilinear sampling.')
    shape = (rows, factor, cols, factor)
    cropped = np.s_[:rows * factor, :cols * factor]
    mean = surface[cropped].reshape(shape).mean(axis=(1, 3), dtype=np.float64).astype(np.float32)
    # Spatial independence is unknown: do not divide RMS by sqrt(block size).
    upper = rms[cropped].reshape(shape).max(axis=(1, 3)).astype(np.float32)
    return mean, upper


def sample_lunar_profile(surface, rms, x, y):
    """Bilinear heights with maximum RMS across the four supporting cells."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    require(x.ndim == 1 and x.shape == y.shape and np.isfinite(x).all()
            and np.isfinite(y).all(), 'Sample coordinates must be matching finite vectors.')
    require(surface.ndim == 2 and surface.shape == rms.shape, 'Surface/RMS grid shapes differ.')
    u = (x - COARSE_TRANSFORM[0]) / 30 - 0.5
    v = (COARSE_TRANSFORM[3] - y) / 30 - 0.5
    rows, cols = surface.shape
    require(((u >= 0) & (v >= 0) & (u < cols - 1) & (v < rows - 1)).all(),
            'Trial corridor extends beyond the lunar grid bilinear sampling bounds.')
    i, j = np.floor(u).astype(int), np.floor(v).astype(int)
    fx, fy = u - i, v - j
    z = ((1 - fx) * (1 - fy) * surface[j, i] + fx * (1 - fy) * surface[j, i + 1]
         + (1 - fx) * fy * surface[j + 1, i] + fx * fy * surface[j + 1, i + 1])
    sigma = np.maximum.reduce([rms[j, i], rms[j, i + 1], rms[j + 1, i], rms[j + 1, i + 1]])
    require(np.isfinite(z).all() and np.isfinite(sigma).all() and (sigma >= 0).all(),
            'Sampled heights or RMS values are invalid.')
    return z, sigma


def analyze_profile(distance, elevation, radius, sigma=None,
                    base_height=2.0, rover_height=2.0, relay_height=5.0):
    """Return metrics, endpoint line, obstruction and optional RMS envelope.

    Positive obstruction means ground is above the direct antenna path.
    A peak at an endpoint cannot serve as an intermediate peak relay.
    """
    d, z = validate_profile(distance, elevation)
    require(math.isfinite(radius) and radius > 0, 'Planetary radius must be positive and finite.')
    for name, value in [('Base', base_height), ('Rover', rover_height), ('Relay', relay_height)]:
        require(math.isfinite(value) and value > 0, f'{name} antenna height must be positive and finite.')
    length, peak = float(d[-1]), int(np.argmax(z))
    t = d / length
    line = (z[0] + base_height) * (1 - t) + (z[-1] + rover_height) * t
    obstruction = z + d * (length - d) / (2 * radius) - line
    worst = int(np.argmax(obstruction))
    centers = np.arange(45, length - 45, 30.0)
    slope = None
    if centers.size:
        grade = (np.interp(centers + 45, d, z) - np.interp(centers - 45, d, z)) / 90
        slope = float(np.percentile(np.degrees(np.arctan(abs(grade))), 95))
    left_min = right_min = both_clear = None
    status = 'not_applicable_peak_at_endpoint'
    if 0 < peak < len(z) - 1:
        p, top = d[peak], z[peak] + relay_height
        left = ((z[0] + base_height) + (top - z[0] - base_height) * d[:peak + 1] / p
                - z[:peak + 1] - d[:peak + 1] * (p - d[:peak + 1]) / (2 * radius))
        dr = d[peak:] - p
        right = (top + (z[-1] + rover_height - top) * dr / (length - p)
                 - z[peak:] - dr * (length - p - dr) / (2 * radius))
        left_min, right_min = float(left.min()), float(right.min())
        both_clear = bool(left_min > 0 and right_min > 0)
        status = 'both_legs_clear' if both_clear else 'one_or_both_legs_obstructed'
    metrics = {
        'length_m': length, 'sample_count': len(d), 'start_ground_m': float(z[0]),
        'end_ground_m': float(z[-1]), 'relief_m': float(np.ptp(z)),
        'max_direct_obstruction_m': float(obstruction[worst]),
        'max_obstruction_distance_m': float(d[worst]), 'peak_distance_m': float(d[peak]),
        'peak_ground_m': float(z[peak]), 'peak_index': peak, 'slope_90m_p95_abs_deg': slope,
        'relay_height_m': relay_height, 'base_height_m': base_height, 'rover_height_m': rover_height,
        'peak_relay_min_clearance_left_m': left_min, 'peak_relay_min_clearance_right_m': right_min,
        'peak_relay_nominal_both_legs_clear': both_clear, 'peak_relay_status': status,
        'curvature_radius_m': radius,
    }
    envelope = None
    if sigma is not None:
        sigma = finite_vector(sigma, 'RMS proxy')
        require(sigma.shape == z.shape and (sigma >= 0).all(),
                'RMS proxy must match the profile length and be nonnegative.')
        envelope = 2 * (sigma + (1 - t) * sigma[0] + t * sigma[-1])
        metrics.update({
            'rms_proxy_min_m': float(sigma.min()), 'rms_proxy_max_m': float(sigma.max()),
            'max_direct_obstruction_with_two_rms_favorable_perturbation_m': float(np.max(obstruction - envelope)),
            'two_rms_envelope_at_nominal_worst_point_m': float(envelope[worst]),
            'two_rms_is_probability_interval': False,
        })
    return metrics, line, obstruction, envelope


def write_raster(path, array, wkt, description, gdal):
    rows, cols = array.shape
    dataset = gdal.GetDriverByName('GTiff').Create(
        str(path), cols, rows, 1, gdal.GDT_Float32, options=['COMPRESS=DEFLATE', 'TILED=YES'])
    require(dataset is not None, f'Cannot create raster: {path}')
    dataset.SetProjection(wkt)
    dataset.SetGeoTransform(COARSE_TRANSFORM)
    band = dataset.GetRasterBand(1)
    band.WriteArray(array)
    band.SetNoDataValue(float('nan'))
    band.SetUnitType('m')
    band.SetDescription(description)
    dataset.FlushCache()
    band = dataset = None


def write_geopackage(path, profiles, distance, srs, ogr, osr):
    srs = srs.Clone()
    srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dataset = ogr.GetDriverByName('GPKG').CreateDataSource(str(path))
    require(dataset is not None, f'Cannot create GeoPackage: {path}')

    def layer(name, geometry_type, fields):
        result = dataset.CreateLayer(name, srs, geometry_type)
        require(result is not None, f'Cannot create GIS layer: {name}')
        for field, kind in fields:
            require(result.CreateField(ogr.FieldDefn(field, kind)) == 0, f'Cannot create field: {field}')
        return result

    def feature(target, values, points, geometry_type):
        result = ogr.Feature(target.GetLayerDefn())
        for key, value in values.items():
            result.SetField(key, value)
        geometry = ogr.Geometry(geometry_type)
        for x, y, z in points:
            geometry.AddPoint(float(x), float(y), float(z))
        result.SetGeometry(geometry)
        require(target.CreateFeature(result) == 0, 'Failed to write GIS feature.')

    lines = layer('lunar_reference_profiles', ogr.wkbLineString25D, [
        ('name', ogr.OFTString), ('role', ogr.OFTString), ('length_m', ogr.OFTReal),
        ('max_obstruction_m', ogr.OFTReal), ('relief_m', ogr.OFTReal)])
    for index, profile in enumerate(profiles):
        m = profile['metrics']
        feature(lines, {'name': m['id'], 'role': 'Primary trial' if index == PRIMARY_INDEX else 'Nearby parallel trial',
                        'length_m': float(distance[-1]), 'max_obstruction_m': m['max_direct_obstruction_m'],
                        'relief_m': m['relief_m']},
                zip(profile['x'], profile['y'], profile['elevation_m']), ogr.wkbLineString25D)
    points = layer('lunar_trial_points', ogr.wkbPoint25D, [
        ('name', ogr.OFTString), ('antenna_height_m', ogr.OFTReal), ('ground_elev_m', ogr.OFTReal),
        ('distance_m', ogr.OFTReal), ('status', ogr.OFTString)])
    profile = profiles[PRIMARY_INDEX]
    m = profile['metrics']
    trials = [('Base', 0, m['base_height_m']), ('Rover', len(distance) - 1, m['rover_height_m'])]
    if m['peak_relay_nominal_both_legs_clear'] is not None:
        trials.insert(1, ('Peak_relay_trial', m['peak_index'], m['relay_height_m']))
    for name, index, height in trials:
        feature(points, {'name': name, 'antenna_height_m': height, 'ground_elev_m': profile['elevation_m'][index],
                         'distance_m': float(distance[index]), 'status': 'Illustrative geometry: ' + m['peak_relay_status']},
                [(profile['x'][index], profile['y'][index], profile['elevation_m'][index])], ogr.wkbPoint25D)
    samples = layer('lunar_profile_checks', ogr.wkbPoint25D, [
        (name, ogr.OFTReal) for name in ['distance_m', 'ground_m', 'rms_upper_m', 'direct_line_m',
                                       'obstruction_m', 'sensitivity_2rms_m']])
    for index, d in enumerate(distance):
        values = {'distance_m': float(d)}
        for key, source in [('ground_m', 'elevation_m'), ('rms_upper_m', 'rms_proxy_m'),
                            ('direct_line_m', 'direct_line_m'), ('obstruction_m', 'obstruction_m'),
                            ('sensitivity_2rms_m', 'two_rms_envelope_m')]:
            values[key] = profile[source][index]
        feature(samples, values, [(profile['x'][index], profile['y'][index], profile['elevation_m'][index])], ogr.wkbPoint25D)
    lines = points = samples = None
    dataset = None
    check = ogr.Open(str(path))
    require(check is not None, 'Cannot reopen exported GeoPackage.')
    counts = [check.GetLayer(index).GetFeatureCount() for index in range(check.GetLayerCount())]
    require(counts == [len(profiles), len(trials), len(distance)], 'Exported GIS feature counts do not match.')


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def run(input_dir, output_dir, base_height=2.0, rover_height=2.0, relay_height=5.0):
    """Validate and calculate before publishing any final output files."""
    inp, out = Path(input_dir), Path(output_dir)
    surface_path = inp / 'Site04_final_adj_5mpp_surf.tif'
    rms_path = inp / 'Site04_final_adj_5mpp_toterr.tif'
    earth_path = inp / 'Kilbourne_profile01.csv'
    d, ze = load_earth_profile(earth_path)
    heights = dict(base_height=base_height, rover_height=rover_height, relay_height=relay_height)
    earth, earth_line, earth_obs, _ = analyze_profile(d, ze, EARTH_RADIUS_M, **heights)
    gdal, ogr, osr = load_gdal()
    surface, rms, srs, wkt = load_lunar_surface(surface_path, rms_path, gdal, osr)
    coarse, upper = aggregate_surface(surface, rms)
    profiles = []
    for y0 in TRANSECT_Y_M:
        x, y = -4000 + d, np.full_like(d, y0)
        z, sigma = sample_lunar_profile(coarse, upper, x, y)
        metrics, line, obs, envelope = analyze_profile(d, z, MOON_RADIUS_M, sigma, **heights)
        metrics.update({'id': f'Moon_y{int(y0)}', 'start_xy_m': [float(x[0]), y0],
                        'end_xy_m': [float(x[-1]), y0],
                        'relative_height_rmse_vs_earth_m': float(np.sqrt(np.mean(((z - z[0]) - (ze - ze[0])) ** 2)))})
        profiles.append({'metrics': metrics, 'x': x.tolist(), 'y': y.tolist(), 'elevation_m': z.tolist(),
                         'rms_proxy_m': sigma.tolist(), 'direct_line_m': line.tolist(),
                         'obstruction_m': obs.tolist(), 'two_rms_envelope_m': envelope.tolist()})
    summary = {
        'schema_version': 2,
        'purpose': 'Exploratory fixed-transect Earth/Moon terrain-obstruction comparison, not analog-site validation',
        'sources': [{'file': path.name, 'sha256': sha256(path)} for path in (surface_path, rms_path, earth_path)],
        'software': {'python': platform.python_version(), 'numpy': np.__version__, 'gdal': gdal.VersionInfo('--version')},
        'source_urls': ['https://pgda.gsfc.nasa.gov/products/78', 'https://doi.org/10.5067/MEaSUREs/NASADEM/NASADEM_HGT.001'],
        'sampling': f'Earth exported profile retained. Lunar 5m cells averaged in 6x6 blocks, bilinearly sampled at the same {len(d)} distances. Both elevation models have approximately 30m grid spacing; effective resolutions and errors differ. No vertical rescaling.',
        'uncertainty_method': 'Maximum native RMS in each 6x6 block, then maximum of four coarse neighboring cells at each profile sample. No independence assumption. Two-RMS envelope sums terrain and weighted endpoint RMS proxies; it is not a confidence interval or guaranteed error bound.',
        'selection': 'Fixed illustrative lines at y=-6000,-6500,-7000m, x=-4000m to -4000+L. Center line is primary. No optimization against Earth data. Not representative of the full Moon or all of Shackleton rim.',
        'physics': f'{base_height:g}m base, {rover_height:g}m rover, {relay_height:g}m trial relay at highest sampled ground point. Local spherical-curvature approximation. No atmospheric refraction, Fresnel clearance, RF propagation, buildings or subgrid obstacles.',
        'limits': ['Earth vertical uncertainty was not supplied and has not been propagated.',
                   'Lunar RMS does not capture every systematic or model error.',
                   'Equal pixel spacing does not equalize effective resolution.',
                   'Results apply only to these trial lines, not Kilbourne Hole as a whole.',
                   'Peak placement is a trial strategy, not an optimized network.',
                   'Keep Earth and Moon coordinates in separate QGIS projects.',
                   'If the highest sample is an endpoint, intermediate peak-relay metrics are null.',
                   'A 90m slope statistic is null when no full centered 90m window is available.'],
        'distance_m': d.tolist(), 'earth_elevation_m': ze.tolist(), 'earth_direct_line_m': earth_line.tolist(),
        'earth_obstruction_m': earth_obs.tolist(), 'earth_metrics': earth, 'lunar_profiles': profiles,
        'primary_lunar_profile_index': PRIMARY_INDEX, 'file_alignment_verified': True,
        'native_rms_percentiles_m': dict(zip(['min', 'median', 'p90', 'p95', 'p99', 'max'],
                                           map(float, np.percentile(rms, [0, 50, 90, 95, 99, 100])))),
    }
    encoded = json.dumps(summary, indent=2, allow_nan=False)
    out.mkdir(parents=True, exist_ok=True)
    # Stage all exports on the same filesystem; existing outputs survive validation/export failures.
    # Final replacements are atomic per file, not a transaction across the whole bundle.
    with tempfile.TemporaryDirectory(prefix='.coral-', dir=out) as temp:
        stage = Path(temp)
        write_raster(stage / 'Site04_30m_mean.tif', coarse, wkt, 'Mean height of 6x6 native 5m cells', gdal)
        write_raster(stage / 'Site04_30m_rms_upper.tif', upper, wkt, 'Maximum native RMS in each 6x6 block', gdal)
        write_geopackage(stage / 'CORAL_Moon_reference01.gpkg', profiles, d, srs, ogr, osr)
        (stage / 'CORAL_comparison.json').write_text(encoded, encoding='utf-8')
        for path in sorted(stage.iterdir()):
            path.replace(out / path.name)
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=Path('inputs'), help='Folder containing native Site04 rasters and the Earth CSV.')
    parser.add_argument('--output-dir', type=Path, default=Path('results'), help='Output folder; successful runs replace matching outputs.')
    parser.add_argument('--base-height', type=float, default=2.0, metavar='METERS')
    parser.add_argument('--rover-height', type=float, default=2.0, metavar='METERS')
    parser.add_argument('--relay-height', type=float, default=5.0, metavar='METERS')
    args = parser.parse_args(argv)
    try:
        summary = run(args.input_dir, args.output_dir, args.base_height, args.rover_height, args.relay_height)
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(2, f'CORAL error: {error}\n')
    print(f'CORAL comparison saved to {args.output_dir / "CORAL_comparison.json"}')
    print(f'Compared {len(summary["distance_m"])} samples across {summary["earth_metrics"]["length_m"]:.2f} m.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
