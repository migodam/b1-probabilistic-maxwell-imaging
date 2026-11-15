"""Postseal four-object state/readout diagnostics; never select an output by GT."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,csv,time,importlib.util

def main():
    root=Path(__file__).resolve().parents[1];w=time.perf_counter();c=time.process_time()
    seal=json.loads((root/'predictions/PREDICTION_SEAL_FORKS_AND_DIAGNOSTICS_V1.json').read_text())
    if not seal['authorized_scoring'] or len(seal['fork_jobs'])!=56:raise ValueError('Actual completed fork seal required')
    for item in seal['files']:
        f=root/item['path']
        if not f.is_file() or f.stat().st_size!=item['bytes']:raise ValueError('Sealed local material review changed: '+str(f))
    sys.path[:0]=[str(root/'runtime'),str(root/'src'),str(root.parent.parent)]
    import numpy as np
    import dualmetric_core as codec
    from queue_interventions import source_checkpoint
    helper=root.parent/'research_max_struct_r14/runtime/scoring_history/src/score_results.py';sys.path.insert(0,str(helper.parent))
    spec=importlib.util.spec_from_file_location('r15_fixed_diagnostic_metrics',helper);h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    out=root/'analysis/intervention_scoring_v1';out.mkdir(exist_ok=False,parents=True)
    rows=[];readouts=[]
    for i in (9,20,21,23):
        cid=f'r14_new_{i:04d}'
        truth=h.load_truth(dict(path=str(root.parent/'research_max_struct_r14/sealed_truth/independent24'/(cid+'.npz')),key='material_coarse',support_rule='fixed_chi_gt_1',host=.1))
        source=source_checkpoint(root,cid);state=codec._load_state_npz(source['path']);m0=state['m']
        for method in ('twofold_hop','twofold_dr','twofold_ap','twofold_full_dr','twofold_reset'):
            folder=root/'RESULTS/fullstate_interventions_v1'/(cid+'__'+method+'__base')
            for view,name in [('shared_initial',None),('paid_primary','primary.npz'),('state16','STATE_016.npz'),('state64','STATE_064.npz')]:
                if name is None:m=m0;reason=None
                elif (folder/name).is_file():
                    with np.load(folder/name,allow_pickle=False) as z:m=z['m'].copy()
                    reason=None
                else:m=None;reason='MISSING_RETAINED_STATE'
                metrics=h.material_metrics(m,truth,.02,20.)
                rows.append(dict(case_id=cid,method=method,view=view,**metrics,null_reason=reason,
                    role='OFFLINE_SAME_FULLSTATE_CAUSAL_DIAGNOSTIC; only paid_primary was selected online, other snapshots never GT-selected',
                    completed_updates=64,starting_global_iteration=source['chosen_iteration']))
        with np.load(root/'audit/readout_counterfactual_v1'/cid/'SOURCE_INPUT.npz',allow_pickle=False) as z:initial=z['m'].copy()
        initial_metrics=h.material_metrics(initial,truth,.02,20.)
        for arm in ('G19_nohop','S19_nohop','U11_nohop','U11_onehop'):
            folder=root/'audit/readout_counterfactual_v1'/cid
            with np.load(folder/(arm+'.npz'),allow_pickle=False) as z:m=z['m_new'].copy()
            metrics=h.material_metrics(m,truth,.02,20.);direction=m-initial;target=truth['m']-initial;denom=np.linalg.norm(direction)*np.linalg.norm(target)
            readouts.append(dict(case_id=cid,readout=arm,**metrics,initial_NRMSE=initial_metrics['NRMSE'],delta_NRMSE=metrics['NRMSE']-initial_metrics['NRMSE'],
                offline_cosine_to_true_material_direction=float(np.dot(direction,target)/denom) if denom>0 else None,
                role='OFFLINE_DIAGNOSTIC_ONLY; frozen16 directions, no candidate/parameter chosen by truth'))
    def write(name,records):
        (out/(name+'.json')).write_text(json.dumps(records,indent=2,allow_nan=False))
        keys=list(dict.fromkeys(k for row in records for k in row))
        with (out/(name+'.csv')).open('x',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(records)
    write('FORK_STATES80',rows);write('READOUT_MATERIAL16',readouts)
    receipt=dict(status='COMPLETE',material_objects=4,fork_state_rows=len(rows),readout_rows=len(readouts),new_LU=0,new_RHS=0,
        truth_access='Four original coarse truths after approved56+16 seal',CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w,
        interpretation='Early common-state continuation, not restart from a converged bad endpoint. No comparison with full120CPU inversion success, no new blind claim.')
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt))

if __name__=='__main__':main()
