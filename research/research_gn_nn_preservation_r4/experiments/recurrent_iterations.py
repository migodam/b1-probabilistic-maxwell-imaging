"""VAL public-input stability beyond the trained two steps; never reads GT."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import argparse,csv,json,sys,time
from pathlib import Path
STUDY=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(STUDY.parents[1]))
import numpy as np
import torch
from research.research_gn_nn_preservation_r4.src.networks import OperatorCache,normalized_features,choose_device
from research.research_gn_nn_preservation_r4.src.recurrent_feedback import RecurrentFeedback
from research.research_gn_nn_preservation_r4.src.predict import _tensor,_numpy,_sync
from research.research_z_alg_r4_diversity.src.physics import Model

def load(path):
    with np.load(path,allow_pickle=False)as z:return {k:z[k]for k in z.files}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--checkpoint',type=Path,required=True)
    p.add_argument('--device',default='mps')
    a=p.parse_args();contract=json.loads((STUDY/'audit/RECURRENT_ITERATION_DIAGNOSTIC.json').read_text())
    checkpoint=torch.load(a.checkpoint,map_location='cpu',weights_only=False)
    if checkpoint['format']!='gn_nn_r4_recurrent_feedback_v1':raise ValueError('R4 recurrent checkpoint required')
    cfg=checkpoint['config'];stats=checkpoint['feature_stats'];device=choose_device(a.device)
    operators=OperatorCache(STUDY/'data/operators.npz').to(device)
    network=RecurrentFeedback(cfg['mode'],cfg['rho']).to(device)
    network.load_state_dict(checkpoint['state_dict']);network.eval()
    network.passes=contract['iterations']
    public_ops=load(STUDY/'data/operators.npz');physics=Model(public_ops)
    out=STUDY/'results/recurrent_iterations'/a.checkpoint.parent.name;out.mkdir(parents=True,exist_ok=True)
    rows=[];start=time.perf_counter();inference_seconds=0.;operation_counts=[]
    for case in contract['cases']:
        obs=case['obs_id'];anchor=load(STUDY/'anchors'/f'{obs}.npz')
        x=normalized_features(anchor,stats,'j_full')[None].to(device)
        arguments=[_tensor(anchor[k][None].reshape(1,576)if k=='m0'else anchor[k][None],device,complex_value=k!='m0') for k in ('m0','j0','Y','Q0')]
        _sync(device);tick=time.perf_counter()
        with torch.no_grad():result=network(operators,x,*arguments[:3],stats,arguments[3],public_noise=_tensor(np.asarray(anchor['noise']).reshape(1),device))
        _sync(device);inference_seconds+=time.perf_counter()-tick
        materials=_numpy(result['stage_materials'])[0];currents=_numpy(result['stage_currents'],current=True)[0];fields=_numpy(result['stage_fields'],current=True)[0]
        operation_counts.append(dict(obs_id=obs,**result['preparation_counts']))
        previous_m=anchor['m0'];previous_j=anchor['j0'];y=anchor['Y'];yn=np.linalg.norm(y);jn=np.linalg.norm(anchor['j0'])
        for iteration,(m,j,q)in enumerate(zip(materials,currents,fields),1):
            real=physics.forward(np.asarray(m,float),label='VAL_recurrent_depth_diagnostic')
            rows.append(dict(**case,iteration=iteration,material_RMS_step=float(np.sqrt(np.mean((m-previous_m)**2))),
                current_relative_step=float(np.linalg.norm(j-previous_j)/jn),
                auxiliary_relative_data=float(np.linalg.norm(public_ops['S']@j-y)/yn),
                constitutive_relative_defect=float(np.linalg.norm(j-m[:,None]*q)/jn),
                physical_current_gap_relative=float(np.linalg.norm(real['j']-j)/jn),
                physical_internal_field_gap_relative=float(np.linalg.norm(real['Q']-q)/np.linalg.norm(anchor['Q0'])),
                true_FullVIE_relative_data=float(np.linalg.norm(real['y']-y)/yn)))
            previous_m,previous_j=m,j
        np.savez_compressed(out/f'{obs}.npz',stage_materials=materials,stage_currents=currents,stage_fields=fields)
    with(out/'trajectory.csv').open('w',newline='')as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (out/'RECEIPT.json').write_text(json.dumps(dict(status='COMPLETE',checkpoint=str(a.checkpoint),
        scope='VAL_ONLY_PUBLIC_RESIDUALS_NO_GT_SCORING',trained_iterations=2,diagnostic_iterations=network.passes,
        cases=len(contract['cases']),wall_seconds=time.perf_counter()-start,NN_inference_seconds=inference_seconds,
        physical_costcounts=dict(physics.counts),fixed_operator_counts=operation_counts,
        conclusion_limit='16-step stability is an empirical out-of-training-depth diagnostic, not proof of fixed-point or GT convergence'),indent=2)+'\n')
    print(out)

if __name__=='__main__':main()
