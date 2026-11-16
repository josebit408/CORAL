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
console.log(JSON.stringify(results));
'''
        completed = subprocess.run(['node', '-e', script, str(ROOT / 'docs/index.html')],
                                   check=True, capture_output=True, text=True)
        records = json.loads(completed.stdout)
        result = json.loads((ROOT / 'results/CORAL_Earth_search.json').read_text())
        earth = {c['metrics']['candidate_id']: c['elevation_m'] for c in result['candidates']}
        earth['baseline'] = result['baseline_elevation_m']
        distance = np.asarray(result['distance_m'])
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
                    for link in network['links']:
                        values = clearance(link['start_sample'], link['end_sample'])
                        if len(values):
                            self.assertGreater(values.min(), 0)
                            self.assertAlmostEqual(link['minimum_interior_sample_clearance_m'],
                                                   float(values.min()), places=8)
                        else:
                            self.assertIsNone(link['minimum_interior_sample_clearance_m'])
