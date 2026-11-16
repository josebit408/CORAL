"""CLI failure behavior and sampled relay-network regression checks."""
import importlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from relay_earth_search import sampled_relay_network


class CommandLineTests(unittest.TestCase):
    def test_imports_do_not_parse_arguments_or_write_files(self):
        code = ('import sys; sys.path.insert(0, sys.argv[1]); '
                'import relay_earth_search, relay_earth_present, relay_earth_plot')
        with tempfile.TemporaryDirectory() as temp:
            result = subprocess.run([sys.executable, '-c', code, str(ROOT / 'src'), '--unknown'],
                                    cwd=temp, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout, '')
            self.assertEqual(list(Path(temp).iterdir()), [])

    def test_help_does_not_require_input_rasters(self):
        for script in sorted((ROOT / 'src').glob('*.py')):
            with self.subTest(script=script.name):
                result = subprocess.run([sys.executable, str(script), '--help'],
                                        cwd=tempfile.gettempdir(), capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('usage:', result.stdout)

    def test_missing_inputs_fail_cleanly_without_outputs(self):
        commands = [
            ('relay_earth_search.py', ['--input-dir', 'absent', '--moon-results', 'absent', '--output-dir', 'out']),
            ('relay_earth_present.py', ['absent']),
            ('relay_earth_plot.py', ['absent']),
        ]
        with tempfile.TemporaryDirectory() as temp:
            for script, args in commands:
                with self.subTest(script=script):
                    result = subprocess.run([sys.executable, str(ROOT / 'src' / script), *args],
                                            cwd=temp, capture_output=True, text=True)
                    self.assertEqual(result.returncode, 2, result.stderr)
                    self.assertIn('CORAL', result.stderr)
                    self.assertNotIn('Traceback', result.stderr)
                    self.assertEqual(list(Path(temp).iterdir()), [])


class EarthNetworkTests(unittest.TestCase):
    def test_adjacent_link_does_not_claim_interior_clearance(self):
        result = sampled_relay_network([0, 30], [0, 0], 6371000)
        self.assertEqual(result['relay_count'], 0)
        self.assertIsNone(result['links'][0]['min_sampled_interior_clearance_m'])
        self.assertIsNone(result['min_sampled_interior_clearance_m'])
        self.assertEqual(result['links_without_interior_samples'], 1)
        json.dumps(result, allow_nan=False)

    def test_checked_clearance_includes_curvature(self):
        result = sampled_relay_network([0, 50, 100], [0, 0, 0], 5000)
        self.assertEqual(result['relay_count'], 0)
        self.assertAlmostEqual(result['min_sampled_interior_clearance_m'], 1.75)
        self.assertEqual(result['links_without_interior_samples'], 0)

    def test_bundled_network_positions_are_preserved(self):
        data = json.loads((ROOT / 'results/CORAL_Earth_search.json').read_text())
        profiles = [(p['elevation_m'], p['nominal_sampled_relay_network'], 6371000)
                    for p in data['candidates']]
        profiles += [(data['moon_elevation_m'], data['moon_nominal_sampled_relay_network'], 1737400),
                     (data['baseline_elevation_m'], data['baseline_nominal_sampled_relay_network'], 6371000)]
        for elevation, stored, radius in profiles:
            actual = sampled_relay_network(data['distance_m'], elevation, radius)
            self.assertEqual(actual['sample_indices'], stored['sample_indices'])
            self.assertEqual(actual['relay_count'], stored['relay_count'])
            self.assertEqual(actual['links_without_interior_samples'],
                             stored['links_without_interior_samples'])
            for actual_link, stored_link in zip(actual['links'], stored['links']):
                expected = actual_link['min_sampled_interior_clearance_m']
                saved = stored_link['min_sampled_interior_clearance_m']
                if expected is None:
                    self.assertIsNone(saved)
                else:
                    self.assertAlmostEqual(expected, saved, places=8)

    def test_invalid_profile_or_radius(self):
        for d, z, radius in [([0, 0], [1, 2], 6371000), ([0, 1], [1, 2], 0),
                             ([0, 1], [1, float('nan')], 6371000)]:
            with self.subTest(distance=d, radius=radius), self.assertRaises(ValueError):
                sampled_relay_network(d, z, radius)


if __name__ == '__main__':
    unittest.main()
