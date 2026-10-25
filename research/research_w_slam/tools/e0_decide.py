"""Parent's bounded forward qualification decision; no new scientific gates."""
from pathlib import Path
import json,gc
import numpy as np
from wslam.accounting import BudgetBook,json_default
from wslam.world import World
from wslam.physics import acquisition
from wslam.forward import Evaluator
root=Path(__file__).resolve().parents[1];out=root/'results/e0'
raw=json.loads((out/'qualification/RESULTS.json').read_text())['rows']
ref=json.loads((out/'refinement01_executable/RESULTS.json').read_text())['rows']
checks=[r for r in ref if 'point_grid_white' in r]
if len(checks)!=23:raise RuntimeError('Incomplete23-case comparison')
book=BudgetBook();tx,rx=acquisition(outer=True);extra=[]
for label,world,f in [('base_outer',World(np.ones((2,2)),np.array([[3.,4.],[3.,4.]]),.08,(.2,2.)),1.),
                      ('strong_outer',World(np.ones((2,2)),6.,.25,(1.,2.6)),1.15)]:
    a=Evaluator(book,'E0E1',tx=tx,rx=rx,h=.17,n=256,grid_refinement=2)
    b=Evaluator(book,'E0E1',tx=tx,rx=rx,h=.17,n=384,grid_refinement=3)
    u=a.evaluate(world,f);v=b.evaluate(world,f)
    difference=float(np.linalg.norm(u.data-v.data)/.01)
    extra.append(dict(case=label,point_grid_white=difference,residual=max(u.state.residual,v.state.residual)))
    np.savez_compressed(out/(label+'.npz'),primary=u.data,independent=v.data,tx=tx,rx=rx)
    a.clear();b.clear();del u,v;gc.collect()
pml=[r['point_difference_white'] for r in ref if r['case'].startswith('pml_')]
passed=all(r['point_grid_white']<=.1 for r in checks+extra) and all(v<=.1 for v in pml)
decision=dict(decision='PASS' if passed else 'FAIL',qualification='QUALIFIED_DISCRETE_FULLWAVE' if passed else 'UNQUALIFIED_NUMERICS',
 scope='new public grid,23 fixed cases, finite sources/frequencies/target cases; not uniform PDE certificate',
 numerical_threshold=.1,maximum_fixed_case_white=max(r['point_grid_white'] for r in checks),
 outer_aperture_checks=extra,primary=dict(h=.17,polygon_n=256,field_grid_refinement=2),
 independent=dict(h=.17,polygon_n=384,field_grid_refinement=3),
 initial_failed_tests_preserved=True,initial_maximum_white=max(r['point_grid_white'] for r in raw if 'point_grid_white' in r),
 geometry_contract='G_GRID_0125_V1; material support remains24x6 macrogrid, independent field subdivision does not change material class',
 shape_derivative_note='physical interface derivative is approximation to discrete remeshing; reported FD residuals, not machine precision; downstream rank/finite endpoint audits required',
 shape_FD_records=[r for r in ref if r['case']=='central_FD_new_primary'],
 continuous_certificate=False,legacy_witness='conditional fixed-support d_joint.798716/cal3.560512/dp.1 retained; not certified as member of new macrogrid prior',
 current_counts=book.totals('E0E1'))
(out/'gate.json').write_text(json.dumps(decision,indent=2,default=json_default));print(json.dumps(decision,default=json_default))
