"""A public-state witness: near-identical material shadows, different futures."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json,time,sys
import numpy as np


def main():
    root=Path(__file__).resolve().parents[1];sys.path[:0]=[str(root/'runtime')]
    from field_current import DenseFieldCurrentProjectors
    cid='B__r12_confirm24_0010'
    paths=[root/'RESULTS/development_v3_local'/f'{cid}__online_event',root/'RESULTS/development_v4_local'/f'{cid}__online_midpoint']
    reference=[]
    for p in paths:
        with np.load(p/'EVENT_GEOMETRY.npz') as z:reference.append({k:z[k].copy() for k in z.files})
    same=all(np.array_equal(reference[0][k],reference[1][k]) for k in ('cX','cQ','dX','dQ','reference_m'))
    with np.load(root/'data_public/development'/f'{cid}.npz') as z:pub={k:z[k].copy() for k in z.files}
    cpu,wall=time.process_time(),time.perf_counter()
    P=DenseFieldCurrentProjectors(pub['G'],pub['E0'],pub['S'],pub['Y'],float(pub['sigma']*pub['noise_radius']),scale=6,lower=.02,upper=20,grouped_ports=True)
    fields=[];metas=[]
    for p,z in zip(paths,reference):
        event=json.loads((p/'EVENTS.json').read_text())[0];t=event['target_t'];x=z['cX']+t*z['dX'];q=z['cQ']+t*z['dQ']
        cx,cq,m=P.material(x,q);fields.append(dict(X=x,Q=q,CX=cx,CQ=cq,m=m))
        metas.append(dict(run=str(p.relative_to(root)),target_t=t,ray_crossed_cells=int(np.sum(z['thresholds']<t)),
            future_primary=json.loads((p/'primary.json').read_text())))
    def distance(k):return float(np.linalg.norm(fields[0][k]-fields[1][k]))
    material=fields[0]['m']-fields[1]['m']
    result=dict(status='EMPIRICAL_PUBLIC_STATE_WITNESS',truth_access='NONE',exact_common_reference_arrays=same,
        state_distance=float(np.hypot(distance('X'),distance('Q'))),
        relative_state_distance=float(np.hypot(distance('X'),distance('Q'))/np.hypot(np.linalg.norm(fields[0]['X']),np.linalg.norm(fields[0]['Q']))),
        material_shadow_max_abs_difference=float(np.max(abs(material))),material_shadow_RMS_difference=float(np.sqrt(np.mean(material**2))),
        constitutive_shadow_field_distance=float(np.hypot(distance('CX'),distance('CQ'))),
        cases=metas,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        projector_cost=dict(P.costs),new_Maxwell_LU=0,
        interpretation='Same one-event interval and numerically equal material readout do not determine the future extended-state path. This is not a global instability/chaos/identifiability theorem.')
    out=root/'delegated/full_current/same_shadow_witness';out.mkdir(exist_ok=False)
    (out/'WITNESS.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(out/'STATES.npz',**{name+'_'+k:v for name,row in zip(('fractional','midpoint'),fields) for k,v in row.items()})
    print(json.dumps({k:result[k] for k in ['exact_common_reference_arrays','relative_state_distance','material_shadow_max_abs_difference','CPU_seconds']}))


if __name__=='__main__':main()
