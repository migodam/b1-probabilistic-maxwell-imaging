"""Cached-only numerical decision. No new Maxwell actions."""
from pathlib import Path
import json,csv,sys
import numpy as np
root=Path(__file__).resolve().parents[1]
def read(path):return json.loads((root/path).read_text())
wall=read('results/wall-qualification-resume/RESULTS.json')
pml=read('results/wall-pml/RESULTS.json')
grad=read('results/fem-derivatives/RESULTS.json')
cyl=read('results/independent/RESULTS.json')
fine=[r for r in wall if r['case']=='wall_complete' and r['h']==.145]
checks=[]
for r in fine:
    for key,value in r['mesh_white'].items():
        band,label=key.split('_');p=next(x for x in pml if x['case']=='wall_PML_complete' and x['seed']==r['seed'] and x['band']==band and x['label']==label)
        checks.append(dict(seed=r['seed'],band=band,state=label,mesh_white=value,PML_strength_white=p['PML_strength_white'],volume_probe_white=p['volume_probe_white'],sum_of_empirical_checks=value+p['PML_strength_white']+p['volume_probe_white'],is_rigorous_bound=False))
assert len(checks)==8
fd=min(r['relative_error'] for r in grad if r['case']=='FD')
adj=next(r['relative_defect'] for r in grad if r['case']=='real_adjoint')
rec=next(r['relative_error'] for r in grad if r['case']=='reciprocity')
passed=max(x['mesh_white'] for x in checks)<=.1 and max(x['sum_of_empirical_checks'] for x in checks)<=.1 and fd<1e-7 and adj<1e-10 and rec<1e-10
output=dict(stage='P0-Q',decision='PASS' if passed else 'FAIL',qualification='QUALIFIED_DISCRETE_FULLWAVE' if passed else 'UNQUALIFIED_NUMERICS',fixed_threshold=.1,checks=checks,adjoint_relative_defect=adj,finite_difference_best=fd,reciprocity_relative_defect=rec,reference_cylinder=[{k:r[k] for k in ('h','n','plane_relative','line_relative','line_white')} for r in cyl if r['case']=='cylinder'],PML_thickness_max_white=max(r['white_difference'] for r in pml if r['case']=='wall_PML_thickness'),continuum_error_bound_established=False,scope='Frozen two seeds, fixed calibration and probe8 target, fixed acquisition; new candidate endpoints require their own finer/independent check',P0_B_authorized=bool(passed),source_commit='db7146cae2d74ba231577d42f225dbfc9a389268')
(root/'P0_Q_QUALIFICATION.json').write_text(json.dumps(output,indent=2))
with (root/'results/P0_Q_ERROR_CHECKS.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=checks[0]);w.writeheader();w.writerows(checks)
# Difference-of-worlds is an added diagnostic, never a replacement gate.
a=np.load(root/'results/wall-qualification-resume/ARRAYS.npz');rows=[]
for band in ('wide','narrow'):
    for state in ('cal','task'):
        d15=a[f'seed2026102201_h0.15_{band}_{state}']-a[f'seed2026102202_h0.15_{band}_{state}']
        d145=a[f'seed2026102201_h0.145_{band}_{state}']-a[f'seed2026102202_h0.145_{band}_{state}']
        rows.append(dict(band=band,state=state,difference_response_mesh_white=float(np.linalg.norm(d145-d15)/.01),two_fixed_worlds_white=float(np.linalg.norm(d145)/.01),candidate_search=False))
(root/'results/DIFFERENCE_RESPONSE.json').write_text(json.dumps(rows,indent=2))
print(json.dumps(dict(decision=output['decision'],mesh_max=max(x['mesh_white'] for x in checks),empirical_sum_max=max(x['sum_of_empirical_checks'] for x in checks),scope=output['scope'])))
