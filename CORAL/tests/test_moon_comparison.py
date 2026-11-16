"""Run from the repository root: python -m unittest discover -s tests -v."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import relay_moon_comparison as moon
import relay_comparison_plot as plot


class ProfileTests(unittest.TestCase):
    def test_obstruction_and_clearance_hand_calculation(self):
        # 100 m span; 10 m central ridge; radius chosen for a 0.25 m central bulge.
        metrics, line, obstruction, envelope = moon.analyze_profile(
            [0, 50, 100], [0, 10, 0], 5000, sigma=[1, 1, 1])
        np.testing.assert_allclose(line, [2, 2, 2])
        np.testing.assert_allclose(obstruction, [-2, 8.25, -2])
        np.testing.assert_allclose(envelope, [4, 4, 4])
        self.assertEqual(metrics['max_obstruction_distance_m'], 50)
        self.assertTrue(metrics['peak_relay_nominal_both_legs_clear'])
        self.assertEqual(metrics['peak_relay_min_clearance_left_m'], 2)
        self.assertEqual(metrics['peak_relay_min_clearance_right_m'], 2)

    def test_endpoint_peak_is_not_an_intermediate_relay(self):
        metrics, *_ = moon.analyze_profile([0, 10, 20], [2, 1, 0], 6371000)
        self.assertIsNone(metrics['peak_relay_nominal_both_legs_clear'])
        self.assertIsNone(metrics['peak_relay_min_clearance_left_m'])
        self.assertIsNone(metrics['slope_90m_p95_abs_deg'])
        json.dumps(metrics, allow_nan=False)

    def test_changed_endpoint_heights(self):
        _, line, _, _ = moon.analyze_profile([0, 50, 100], [0, 10, 0], 5000,
                                            base_height=3, rover_height=7)
        np.testing.assert_allclose(line, [3, 5, 7])

    def test_invalid_profiles(self):
        for d, z in [([], []), ([0, 0], [1, 2]), ([1, 2], [1, 2]),
                     ([0, 2, 1], [1, 2, 3]), ([0, 1], [1, float('nan')]),
                     ([0, 1], [1, 2, 3])]:
            with self.subTest(distance=d), self.assertRaises(ValueError):
                moon.analyze_profile(d, z, 6371000)
        for kwargs in [{'relay_height': -1}, {'base_height': float('inf')}, {'sigma': [-1, 1]}]:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                moon.analyze_profile([0, 1], [1, 2], 6371000, **kwargs)

    def test_csv_errors_and_bom(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'profile.csv'
            for text in ['layer,distance\na,0\n', 'layer,distance,elevation\n',
                         'layer,distance,elevation\na,0,1\nb,10,2\n',
                         'layer,distance,elevation\na,0,bad\na,10,2\n']:
                path.write_text(text)
                with self.assertRaises(ValueError):
                    moon.load_earth_profile(path)
            path.write_text('\ufefflayer,distance,elevation\na,0,1\na,10,2\n', encoding='utf-8')
            d, z = moon.load_earth_profile(path)
            np.testing.assert_array_equal(d, [0, 10])
            np.testing.assert_array_equal(z, [1, 2])

    def test_incomplete_blocks_and_rms_maximum(self):
        surface = np.arange(25, dtype=float).reshape(5, 5)
        rms = np.ones((5, 5)); rms[1, 1] = 7; rms[-1, :] = 99
        mean, upper = moon.aggregate_surface(surface, rms, factor=2)
        np.testing.assert_array_equal(mean, [[3, 5], [13, 15]])
        np.testing.assert_array_equal(upper, [[7, 1], [1, 1]])

    def test_bilinear_center_and_bounds(self):
        surface = np.array([[0., 2.], [4., 6.]])
        rms = np.array([[1., 3.], [2., 4.]])
        z, sigma = moon.sample_lunar_profile(surface, rms, [-8970], [970])
        np.testing.assert_array_equal(z, [3])
        np.testing.assert_array_equal(sigma, [4])
        with self.assertRaisesRegex(ValueError, 'bounds'):
            moon.sample_lunar_profile(surface, rms, [-9000], [970])

    def test_optimized_python_keeps_input_checks(self):
        code = "import sys; sys.path.insert(0, sys.argv[1]); import relay_moon_comparison as m; m.validate_profile([0,0],[1,2])"
        result = subprocess.run([sys.executable, '-O', '-c', code, str(ROOT / 'src')], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Distances must increase strictly', result.stderr)

    def test_bad_input_does_not_touch_existing_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            output = folder / 'out'; output.mkdir()
            previous = output / 'CORAL_comparison.json'; previous.write_text('previous result')
            (folder / 'Kilbourne_profile01.csv').write_text('wrong,headers\n')
            with self.assertRaises(ValueError):
                moon.run(folder, output)
            self.assertEqual(previous.read_text(), 'previous result')


class PlotTests(unittest.TestCase):
    def test_old_reference_is_supported(self):
        result = plot.load_comparison(ROOT / 'data/reference_moon/CORAL_comparison.json')
        self.assertEqual(len(result['distance_m']), 83)

    def test_bad_json_schema(self):
        original = plot.load_comparison(ROOT / 'data/reference_moon/CORAL_comparison.json')
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'bad.json'
            for mutate in [lambda r: r.update(primary_lunar_profile_index=-1),
                           lambda r: r.update(earth_elevation_m=[1, 2]),
                           lambda r: r['lunar_profiles'][0].update(two_rms_envelope_m=[-1] * 83)]:
                result = copy.deepcopy(original); mutate(result)
                path.write_text(json.dumps(result))
                with self.assertRaises(ValueError):
                    plot.load_comparison(path)

    def test_import_has_no_cli_side_effects(self):
        code = "import sys; sys.path.insert(0,sys.argv[1]); import relay_moon_comparison,relay_comparison_plot"
        result = subprocess.run([sys.executable, '-c', code, str(ROOT / 'src'), '--unknown-option'],
                                capture_output=True, text=True, cwd=tempfile.gettempdir())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, '')


if __name__ == '__main__':
    unittest.main()
