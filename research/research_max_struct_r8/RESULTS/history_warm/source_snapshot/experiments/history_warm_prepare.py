"""Freeze exactly four historical public inputs and original10feature NN starts.

No new42 confirmation metadata/truth, no historical truth, no physics calls.
The fourth40ac topology control was selected by parent before any NN result.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,json
from pathlib import Path
import numpy as np
import torch
from research.research_max_struct_r8.experiments.run_case import load_public,REPO,ROOT,dump
from research.research_max_struct_r8.src.warm_features import FeatureOperator
from research.research_max_struct_r8.src.warm_network import load_predictor

IDS=['obs_7b8d5df517f1497d','obs_b20d5a55677d4eca','obs_42cb233e4f224bb3','obs_40ac5b79e8a74ac3']

def main():
    wall,cpu=time.perf_counter(),time.process_time();torch.set_num_threads(1);torch.set_num_interop_threads(1)
    folder=ROOT/'data_public/history_warm';folder.mkdir(parents=True,exist_ok=False)
    cases={r['case_id']:r for r in json.loads((REPO/'research/research_max_struct_r6/data_public/cases.json').read_text())['records']}
    weights=ROOT/'RESULTS/training/physics/weights.pt';t,c=time.perf_counter(),time.process_time();predict,checkpoint=load_predictor(weights)
    loadcost=dict(cpu_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t)
    rows=[];tasks=[]
    for ident in IDS:
        case=cases[ident];files=list(dict.fromkeys([case['operators_file'],case['observation_file']]))
        t,c=time.perf_counter(),time.process_time();public=load_public(files);operator=FeatureOperator(public);features,_,featurecost=operator(public['Y'])
        ft,fc=time.perf_counter(),time.process_time();m=np.clip(predict(features)[0].reshape(-1).astype(float),.02,20.);inference=dict(cpu_seconds=time.process_time()-fc,wall_seconds=time.perf_counter()-ft)
        init=folder/(ident+'__physics_nn.npz');np.savez_compressed(init,m_init=m)
        np.savez_compressed(folder/(ident+'__features.npz'),features=features)
        rows.append(dict(case_id=ident,cohort='EXPOSED_HISTORICAL_DEVELOPMENT_REGRESSION',public_files=files,initialization_file=str(init.relative_to(REPO)),feature_setup=operator.setup,feature_evaluation=featurecost,inference=inference,total_cpu_seconds=time.process_time()-c,total_wall_seconds=time.perf_counter()-t,material_min=float(m.min()),material_max=float(m.max()),sources=int(public['E0'].shape[1]),receivers=int(public['S'].shape[0]),amplitude=float(case['amplitude']),truth_access='NONE',physical_calls=0))
        for mode,arm in [('nn','cie_dual'),('nn','a_plus'),('nn','a_plus_sourcewise'),('nn','physical_lbfgs'),('uniform','a_plus_sourcewise')]:
            tag=f'{ident}__{mode}__{arm}';task=dict(task_id=tag,case_id=ident,cohort='EXPOSED_HISTORICAL_DEVELOPMENT_REGRESSION',arm=arm,initialization='original10feature_physics_NN' if mode=='nn' else 'uniform0p1',public_files=files,config='research/research_max_struct_r8/configs/solver.json',output_dir=str((ROOT/'RESULTS/history_warm'/tag).relative_to(REPO)),initialize_physical_current=False)
            if mode=='nn':task['initialization_file']=str(init.relative_to(REPO))
            tasks.append(task)
    base=next(t for t in tasks if t['case_id']=='obs_42cb233e4f224bb3' and t['arm']=='a_plus');extra=dict(base);extra.pop('initialization_file');extra.update(task_id='obs_42cb233e4f224bb3__uniform__a_plus',initialization='uniform0p1',output_dir='research/research_max_struct_r8/RESULTS/history_warm/obs_42cb233e4f224bb3__uniform__a_plus',reason='Singleparentrequested samebudgetuniformglobalAplus pairedcontrol; nohypersearch');tasks.append(extra)
    manifest=dict(status='FIXED21_HISTORICAL_TASKS_NO_NEW_CONFIRMATION',tasks=tasks,config='research/research_max_struct_r8/configs/solver.json',weights=str(weights.relative_to(REPO)),fourth_selection='40ac slit-shell different-topology disclosed control, parentdecision beforeNNresults; R6 fourth63ab not included',initial_Phi_rejection=False)
    dump(ROOT/'experiments/history_warm_manifest.json',manifest);dump(folder/'PREPROCESS_RECEIPT.json',dict(status='ORIGINAL_PHYSICS_NN_INITIALIZATIONS_FROZEN',weights=str(weights),training_epoch=checkpoint.get('selected_epoch',checkpoint.get('epoch')),weights_feature_mode=checkpoint['feature_mode'],model_load=loadcost,records=rows,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,truth_access='NONE',LU=0,complex_RHS=0,NN_initial_Phi_rejection=False,scope='Fourpublicdatasets, native tenfeatures exactly asoriginalphysics training; no paper-domain coordinate correction or extra channels'))
    print(json.dumps(dict(status='PREPARED',tasks=len(tasks),objects=len(rows),manifest=str(ROOT/'experiments/history_warm_manifest.json'))),flush=True)

if __name__=='__main__':main()
