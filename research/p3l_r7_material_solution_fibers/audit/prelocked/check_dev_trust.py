"""Reconstruct executed DEV correction steps from paid x; decoder only."""
import importlib.util
import json
import sys
from pathlib import Path
import numpy as np

P=Path(__file__).resolve().parents[2]
REPO=P.parents[1]
R4=REPO.parent/'b1-p3l-r4-conditional-transport'
sys.dont_write_bytecode=True
sys.path.insert(0,str(R4/'src'))
from p1_material.prior import MaterialPrior
spec=importlib.util.spec_from_file_location('r7_prelocked_jets_snapshot',P/'audit/prelocked/source_snapshot/material_jets.py')
jets=importlib.util.module_from_spec(spec);spec.loader.exec_module(jets)
prior=MaterialPrior(json.loads((R4/'research/p3l_r4_conditional_transport/configs/TARGET.json').read_text()),
    json.loads((R4/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points'])
cache={}
def endpoint_x(key):
    job=key.rsplit('_',1)[0]
    if not job.startswith('R6_CDEV_'):raise ValueError('audit restricted to exposed DEV endpoints')
    if key not in cache:
        path=P/'raw/jobs'/job/'endpoints'/f'{key}.npz'
        with np.load(path,allow_pickle=False) as z:cache[key]=z['x'].copy()
    return cache[key]

rows=[];mismatch=[];omitted=[]
for f in sorted((P/'E2/parents').glob('R6_CDEV_*/**/CORRECTION.json')):
    rel=f.relative_to(P/'E2/parents')
    radius=next((float(part[4:]) for part in rel.parts if part.startswith('RMS_')), .02)
    history=json.loads(f.read_text())['history']
    best=None;best_obj=float('inf')
    for i,row in enumerate(history):
        if row['feasible'] and row['objective']<best_obj:
            best=row['endpoint'];best_obj=row['objective']
        current=best if best is not None else row['endpoint']
        if 'step_norm' not in row:continue
        if i+1>=len(history):
            omitted.append(dict(file=str(f),iteration=i,reason='no later evaluated endpoint; proposal not charged/executed'))
            continue
        nextkey=history[i+1]['endpoint']
        x=endpoint_x(current);xn=endpoint_x(nextkey);step=xn-x
        latent=float(np.linalg.norm(step))
        if abs(latent-row['step_norm'])>1e-10:
            mismatch.append(dict(file=str(f),iteration=i,reconstructed=latent,logged=row['step_norm']))
        M=jets.real_material(jets.material_derivative(prior,x))
        linear=float(np.linalg.norm(M@step))
        finite=float(np.linalg.norm(jets.real_material(jets.decode(prior,xn)[0]-jets.decode(prior,x)[0])))
        trust=row['trust_radius'];limit=.5*radius
        rows.append(dict(parent=rel.parts[0],arm='/'.join(rel.parts[1:-1]),iteration=i,
            current_endpoint=current,next_endpoint=nextkey,distance=radius,
            latent_step=latent,latent_bound=trust,latent_ratio=latent/trust,
            linear_material_step=linear,finite_material_step=finite,material_bound=limit,
            linear_material_ratio=linear/limit,finite_material_ratio=finite/limit,
            latent_exceeds=latent>trust*(1+1e-8)+1e-12,
            linear_material_exceeds=linear>limit*(1+1e-8)+1e-12,
            finite_material_exceeds=finite>limit*(1+1e-8)+1e-12))
affected=[r for r in rows if r['latent_exceeds'] or r['linear_material_exceeds'] or r['finite_material_exceeds']]
report=dict(scope='exposed DEV correction metadata and saved paid endpoint x only; original decoder; no GT/Maxwell',
    status='EXCEEDANCES_FOUND' if affected else 'PASS',executed_steps_checked=len(rows),
    latent_exceedances=sum(r['latent_exceeds'] for r in rows),
    linear_material_exceedances=sum(r['linear_material_exceeds'] for r in rows),
    finite_material_exceedances=sum(r['finite_material_exceeds'] for r in rows),
    affected_steps=affected,reconstruction_mismatches=mismatch,omitted_unexecuted=omitted,
    new_Maxwell_RHS=0,GT_access=False)
assert not mismatch,mismatch
(P/'audit/prelocked/DEV_RECONSTRUCTED_TRUST_AUDIT.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
