"""Saved arrays and matrix products only; does not import a physics backend."""
import json,shutil
import numpy as np
from scipy.stats import chi2
from .common import ROOT,OUT,LXZ,dump,pack,rms
from prior_constrained import prior

def run():
    dest=OUT/'inputs';dest.mkdir(parents=True,exist_ok=True)
    predpath=prior.A26/'results/a26/evaluation/261009/PREDICTIONS.npz'
    valpath=prior.INPUTS/'validation.npz'
    with np.load(predpath) as p,np.load(valpath) as z:
        case='p250065_weak_a1_snr20';i=int(np.flatnonzero(p['case_id']==case)[0]);iv=int(np.flatnonzero(z['case_id']==case)[0])
        rows=[]
        for name,key,role in [('physics_R_base','BASE_N2','PHYSICS_ONLY_ONE_ROM_GRADIENT_STEP'),('N3','N3','LEARNED_POINT_ESTIMATE')]:
            chi=p[key][i];np.savez_compressed(dest/f'{name}.npz',chi=chi,truth=z['chi'][iv],truth_prediction=z['clean'][iv],
                observation=z['data'][iv],sigma=z['sigma'][iv],background=prior.background('weak'))
            rows.append(dict(id=name,role=role,case_id=case,background='weak',source=str(predpath),row=i,
                source_key=key,GT_initialization=False,GT_use='ORACLE_GEOMETRY_DIAGNOSTIC_ONLY',
                exact_anchor_prediction_cached=key=='N3',qualified_P0_anchor=False,
                prior_audit=prior.phase_audit(chi,prior.background('weak')),prior_chart=None))
    for name in ['250033_weak','250034_ring']:
        for ext in ['json','npz']:
            shutil.copyfile(LXZ/f'delegated/e0_2d/PUBLIC_TRAIN_CALIBRATION_{name}.{ext}',dest/f'{name}.{ext}')
        a=np.load(dest/f'{name}.npz');m=json.loads((dest/f'{name}.json').read_text())
        rows.append(dict(id=name,role='KNOWN_P0_TRAIN_CALIBRATION_NOT_INVERSE_RECOVERY',background=name.split('_')[-1],source=m['source'],
            row=m['source_index'],qualified_P0_anchor=True,prior_chart='same original eight uniform Gaussian coordinates as LXZ',
            cached_whitened_residual=float(np.sqrt(2)*np.linalg.norm(a['prediction']-a['observation'])/float(a['sigma'])),
            own_GT_difference='IDENTICALLY_ZERO_NOT_EVIDENCE_OF_RECOVERY',prior_audit=m['prior_certificate']))
    tau=float(np.sqrt(chi2.ppf(.95,512)))
    cache=prior.A30/'inputs/a29_release/inputs/a28/results/a28/exact_states'/f'{case}.npz'
    with np.load(cache) as z,np.load(dest/'N3.npz') as q:
        n3cached=dict(source=str(cache),prediction_distance_to_GT=float(np.sqrt(2)*np.linalg.norm(z['data'][0]-q['truth_prediction'])/float(q['sigma'])),
            residual=float(np.sqrt(2)*np.linalg.norm(z['data'][0]-q['observation'])/float(q['sigma'])),
            truth_residual=float(np.sqrt(2)*np.linalg.norm(q['truth_prediction']-q['observation'])/float(q['sigma'])),
            no_new_physics=True)
    report=dict(anchors=rows,selection='first pre-existing shared development case in frozen A17 order; same two frozen LXZ TRAIN smooth anchors; no outcome/GT-quality selection',
        new_physics=0,new_RHS=0,tau=tau,noise='complex circular sigma; real whitening sqrt(2)/sigma; 512 real observations; fixed gain 1',
        N3_cached=n3cached,independent_converged_physics_GN_or_TSOM='MISSING; R base is one update, not idealized sensitive-information recovery',
        warning='truth noise residual exceeds frozen 95% ball in first case: do not adjust threshold or infer model mismatch',
        no_locked_data=True,old_LXZ_budget_closed_and_unchanged=True)
    dump(OUT/'E0_COMPLETE.json',report)
    print(json.dumps(dict(new_physics=0,anchors=[r['id'] for r in rows],N3_cached=n3cached)))

if __name__=='__main__':run()
