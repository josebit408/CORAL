"""Render a CORAL Earth/Moon comparison JSON using data-derived annotations.

Requires NumPy and Matplotlib. Importing this module does not parse arguments,
read data, or generate files. Older CORAL comparison JSON remains supported.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile

import numpy as np


def vector(data, key, size=None, nonnegative=False):
    try:
        result = np.asarray(data[key], dtype=float)
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f'Missing or invalid array: {key}') from error
    if result.ndim != 1 or len(result) < 2 or not np.isfinite(result).all():
        raise ValueError(f'{key} must contain at least two finite samples.')
    if size is not None and len(result) != size:
        raise ValueError(f'{key}: expected {size} samples, found {len(result)}.')
    if nonnegative and (result < 0).any():
        raise ValueError(f'{key} must be nonnegative.')
    return result


def relay_outcome(metrics):
    value = metrics.get('peak_relay_nominal_both_legs_clear')
    if value is True:
        return 'both sampled legs clear'
    if value is False:
        return 'at least one sampled leg obstructed'
    return 'no intermediate peak-relay result'


def load_comparison(path):
    with Path(path).open(encoding='utf-8-sig') as stream:
        result = json.load(stream)
    if not isinstance(result, dict):
        raise ValueError('Comparison JSON must be an object.')
    distance = vector(result, 'distance_m')
    if distance[0] != 0 or not (np.diff(distance) > 0).all():
        raise ValueError('Distances must start at 0 and increase strictly.')
    size = len(distance)
    vector(result, 'earth_elevation_m', size)
    vector(result, 'earth_obstruction_m', size)
    profiles = result.get('lunar_profiles')
    index = result.get('primary_lunar_profile_index')
    if not isinstance(profiles, list) or not profiles or type(index) is not int or not 0 <= index < len(profiles):
        raise ValueError('Missing lunar profiles or invalid primary profile index.')
    for profile in profiles:
        if not isinstance(profile, dict):
            raise ValueError('Each lunar profile must be an object.')
        for key in ['elevation_m', 'obstruction_m', 'two_rms_envelope_m']:
            vector(profile, key, size, nonnegative=key == 'two_rms_envelope_m')
    for metrics in [result.get('earth_metrics')] + [p.get('metrics') for p in profiles]:
        if not isinstance(metrics, dict):
            raise ValueError('Each profile must include its calculation metrics.')
        for key in ['base_height_m', 'rover_height_m', 'relay_height_m']:
            value = metrics.get(key)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value) or value <= 0:
                raise ValueError(f'Missing or invalid antenna height: {key}')
    return result


def render_comparison(result, output_path, dpi=170):
    """Render only validated results; close figures even when saving fails."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if type(dpi) is not int or not 50 <= dpi <= 600:
        raise ValueError('DPI must be an integer between 50 and 600.')
    distance = np.asarray(result['distance_m'])
    earth = np.asarray(result['earth_elevation_m'])
    profile = result['lunar_profiles'][result['primary_lunar_profile_index']]
    moon = np.asarray(profile['elevation_m'])
    moon_obs, earth_obs = np.asarray(profile['obstruction_m']), np.asarray(result['earth_obstruction_m'])
    envelope = np.asarray(profile['two_rms_envelope_m'])
    earth_metrics, moon_metrics = result['earth_metrics'], profile['metrics']
    maxima = [max(p['obstruction_m']) for p in result['lunar_profiles']]
    earth_color, moon_color = '#246b91', '#b66a20'
    output_path = Path(output_path)
    if output_path.suffix.lower() != '.png':
        raise ValueError('The output filename must end in .png.')
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({'font.family': 'DejaVu Sans', 'font.size': 11}):
        fig, (ax, bx) = plt.subplots(2, 1, figsize=(10.8, 9.4), sharex=True)
        try:
            fig.subplots_adjust(left=.1, right=.97, top=.85, bottom=.30, hspace=.42)
            fig.suptitle('CORAL | Earth–Moon trial comparison', x=.06, y=.98, ha='left', fontsize=18, fontweight='bold')
            fig.text(.06, .93, f'Same {distance[-1]:,.1f} m link length · {len(distance)} profile samples · No vertical rescaling.',
                     fontsize=11, color='#445469')
            ax.plot(distance, earth - earth[0], color=earth_color, linewidth=2.2, label='Earth: supplied trial profile')
            ax.plot(distance, moon - moon[0], color=moon_color, linewidth=2.2, label='Moon: primary Site04 trial')
            ax.axhline(0, color='#8a929a', linewidth=.7)
            ax.set(ylabel='Height relative to start (m)', title='Terrain profiles aligned at their starting ground elevation')
            ax.legend(frameon=False, loc='upper left', fontsize=9)
            bx.fill_between(distance, moon_obs - envelope, moon_obs + envelope, color=moon_color,
                            alpha=.22, label='Moon: two-RMS sensitivity envelope')
            bx.plot(distance, earth_obs, color=earth_color, linewidth=2.2, label='Earth obstruction')
            bx.plot(distance, moon_obs, color=moon_color, linewidth=2.2, label='Moon obstruction')
            bx.axhline(0, color='#5c6873', linewidth=1)
            bx.set(xlabel='Distance from trial base (m)', ylabel='Terrain above direct path (m)',
                   title='Positive values indicate direct-path obstruction')
            bx.legend(frameon=False, loc='upper left', fontsize=8)
            for label, obs, color, offset in [('Moon', moon_obs, moon_color, -28), ('Earth', earth_obs, earth_color, 28)]:
                index = int(np.argmax(obs))
                # Offset in screen points, not meters, so labels also work on short profiles.
                dx = -16 if distance[index] > distance[-1] * .65 else 16
                bx.annotate(f'{label}: {obs[index]:.1f} m', xy=(distance[index], obs[index]),
                            xytext=(dx, offset), textcoords='offset points',
                            ha='right' if dx < 0 else 'left', color=color, fontsize=10,
                            arrowprops={'arrowstyle': '->', 'color': color})
            for axis in (ax, bx):
                axis.grid(alpha=.2)
                axis.spines['top'].set_visible(False)
                axis.spines['right'].set_visible(False)
                axis.set_xlim(0, distance[-1])
            caption = (
                f'Lunar trials: {min(maxima):.1f}–{max(maxima):.1f} m maximum obstruction. Earth: {max(earth_obs):.1f} m.\n'
                f'Earth antennas: base {earth_metrics["base_height_m"]:g} m / rover {earth_metrics["rover_height_m"]:g} m. '
                f'Moon: base {moon_metrics["base_height_m"]:g} m / rover {moon_metrics["rover_height_m"]:g} m.\n'
                f'Earth peak-relay trial ({earth_metrics["relay_height_m"]:g} m): {relay_outcome(earth_metrics)}.\n'
                f'Moon peak-relay trial ({moon_metrics["relay_height_m"]:g} m): {relay_outcome(moon_metrics)}.'
            )
            fig.text(.06, .13, caption, fontsize=10, color='#293b4c', linespacing=1.5)
            fig.text(.06, .035,
                     'Exploratory examples, not analog-site validation. Curvature included. RMS sensitivity is not a confidence interval.\n'
                     'Earth elevation uncertainty and radio propagation are not modeled. See the comparison JSON for provenance.',
                     fontsize=9, color='#55616d')
            with tempfile.TemporaryDirectory(prefix='.coral-plot-', dir=output_path.parent) as temp:
                pending = Path(temp) / output_path.name
                fig.savefig(pending, dpi=dpi, bbox_inches='tight')
                pending.replace(output_path)
        finally:
            plt.close(fig)
    return output_path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results_dir', nargs='?', type=Path, default=Path('results'))
    parser.add_argument('--output', type=Path, help='PNG output path; defaults to the results folder.')
    parser.add_argument('--dpi', type=int, default=170, help='Resolution from 50 to 600 DPI (default: 170).')
    args = parser.parse_args(argv)
    try:
        result = load_comparison(args.results_dir / 'CORAL_comparison.json')
        output = args.output or args.results_dir / 'CORAL_EarthMoon_comparison.png'
        render_comparison(result, output, args.dpi)
    except (OSError, ValueError, RuntimeError, ImportError) as error:
        parser.exit(2, f'CORAL plot error: {error}\n')
    print(f'Chart saved to {output}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
