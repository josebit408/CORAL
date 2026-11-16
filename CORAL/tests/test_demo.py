"""Exercise the demo calculations without a browser; Node.js is optional."""
import json
from pathlib import Path
import shutil
import subprocess
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which('node'), 'Node.js is required to check the web demo')
class DemoTests(unittest.TestCase):
    def test_all_candidates_against_independent_relay_calculation(self):
        script = r'''
const fs=require('fs'),vm=require('vm');
const html=fs.readFileSync(process.argv[1],'utf8'),elements={};
const context={document:{getElementById:id=>elements[id]??={value:'39',textContent:'',innerHTML:'',setAttribute(){},addEventListener(){}},querySelectorAll:()=>[]}};
vm.createContext(context);vm.runInContext(html.match(/<script>([\s\S]*?)<\/script>/)[1],context);
const results=[];
for(const id of ['RS01','RS02','RS03','baseline'])for(const h of [5,10,20]){
 for(const view of ['compare','earth','moon']){
  vm.runInContext(`selected='${id}';height=${h};view='${view}';draw()`,context);
  if(/NaN|undefined/.test(elements.chart.innerHTML))throw Error('Invalid chart');
  if(elements.score.textContent!==vm.runInContext('SITES[selected].rmse.toFixed(2)',context))throw Error('Stale score');
 }
 results.push({id,h,summary:vm.runInContext('currentSummary()',context)});
}
console.log(JSON.stringify({results, sites:vm.runInContext('SITES',context), moon:vm.runInContext('DATA.moon',context),length:vm.runInContext('L',context),samples:vm.runInContext('N',context)}));
'''
        completed = subprocess.run(['node', '-e', script, str(ROOT / 'docs/index.html')],
                                   check=True, capture_output=True, text=True)
        payload = json.loads(completed.stdout)
        records = payload['results']
        result = json.loads((ROOT / 'results/CORAL_Earth_search.json').read_text())
        earth = {c['metrics']['candidate_id']: c['elevation_m'] for c in result['candidates']}
        earth['baseline'] = result['baseline_elevation_m']
        distance = np.asarray(result['distance_m'])
        self.assertEqual(payload['samples'], result['sample_count'])
        self.assertAlmostEqual(payload['length'], result['length_m'])
        np.testing.assert_array_equal(payload['moon'], result['moon_elevation_m'])
        for site_id, z in earth.items():
            site = payload['sites'][site_id]
            np.testing.assert_array_equal(site['z'], z)
            relative = np.asarray(z) - z[0]
            moon_relative = np.asarray(payload['moon']) - payload['moon'][0]
            self.assertAlmostEqual(site['rmse'], float(np.sqrt(np.mean(
                (relative - moon_relative) ** 2))), places=9)
        for candidate in result['candidates']:
            metrics = candidate['metrics']
            site = payload['sites'][metrics['candidate_id']]
            self.assertEqual(site['coords'], metrics['midpoint_lon_lat'])
            self.assertEqual(site['sensitivity'],
                             metrics['local_position_bearing_diagnostic_rmse_range_m'])
        for record in records:
            for planet, z, radius in [('earth', earth[record['id']], 6371000),
                                      ('moon', result['moon_elevation_m'], 1737400)]:
                with self.subTest(site=record['id'], height=record['h'], planet=planet):
                    z = np.asarray(z)
                    heights = np.full(len(z), float(record['h'])); heights[[0, -1]] = 2
                    # Reverse recurrence, independent of the demo's forward solver.
                    hops = np.full(len(z), 999, dtype=int); hops[-1] = 0
                    def clearance(i, j):
                        t = (distance[i+1:j] - distance[i]) / (distance[j] - distance[i])
                        return ((z[i] + heights[i]) * (1 - t) + (z[j] + heights[j]) * t
                                - z[i+1:j] - (distance[i+1:j] - distance[i]) *
                                (distance[j] - distance[i+1:j]) / (2 * radius))
                    for i in range(len(z)-2, -1, -1):
                        hops[i] = min(1 + hops[j] for j in range(i+1, len(z))
                                      if (clearance(i, j) > 0).all())
                    network = record['summary'][planet]
                    self.assertEqual(network['count'], int(hops[0] - 1))
                    self.assertEqual(network['path'][0], 0)
                    self.assertEqual(network['path'][-1], len(z)-1)
                    self.assertEqual(len(network['links']), len(network['path']) - 1)
                    unchecked = 0
                    checked_clearances = []
                    for link in network['links']:
                        values = clearance(link['start_sample'], link['end_sample'])
                        if len(values):
                            self.assertGreater(values.min(), 0)
                            checked_clearances.append(float(values.min()))
                            self.assertAlmostEqual(link['minimum_interior_sample_clearance_m'],
                                                   float(values.min()), places=8)
                        else:
                            self.assertIsNone(link['minimum_interior_sample_clearance_m'])
                            unchecked += 1
                    self.assertEqual(network['links_without_interior_samples'], unchecked)
                    if checked_clearances:
                        self.assertAlmostEqual(network['minimum_interior_sample_clearance_m'],
                                               min(checked_clearances), places=8)
                    else:
                        self.assertIsNone(network['minimum_interior_sample_clearance_m'])
