"""Freeze observation-only initial materials; never reject by initial physical Phi."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[name]='1'
import argparse,json,time,sys
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter
import torch
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r8.src.warm_network import load_predictor

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cases',required=True);ap.add_argument('--weights',required=True);ap.add_argument('--name',required=True);args=ap.parse_args()
    begin,cpu=time.perf_counter(),time.process_time();torch.set_num_threads(1)
    rows=json.loads((REPO/args.cases).read_text())['records'];out=ROOT/'data_public/initializers'/args.name;out.mkdir(parents=True,exist_ok=False)
    predictor,checkpoint=load_predictor(REPO/args.weights);loadcpu=time.process_time()-cpu;loadwall=time.perf_counter()-begin;receipts=[];tasks=[]
    for i,row in enumerate(rows):
        with np.load(REPO/row['features_file']) as z:X=z['features'].copy();born=z['born'].copy();dpcs=z['dpcs'].copy()
        t,c=time.perf_counter(),time.process_time();nn= predictor(X)[0].ravel().astype(np.float64)
        rec=dict(case_id=row['case_id'],inference_cpu_seconds=time.process_time()-c,inference_wall_seconds=time.perf_counter()-t,truth_access='NONE',LU=0,complex_RHS=0)
        rng=np.random.default_rng(2026110610+i);random=.1+2*gaussian_filter(rng.random((24,24)),sigma=2).ravel()
        for mode,m in [('nn',nn),('born',born),('dpcs',dpcs),('random',random)]:
            path=out/(row['case_id']+'__'+mode+'.npz');np.savez_compressed(path,m_init=np.clip(m,.02,20.))
            tasks.append(dict(case_id=row['case_id'],mode=mode,path=str(path.relative_to(REPO))))
        receipts.append(rec)
    (out/'MANIFEST.json').write_text(json.dumps(dict(initializers=tasks,cases=args.cases,weights=args.weights),indent=2))
    (out/'RECEIPT.json').write_text(json.dumps(dict(status='OBSERVATION_ONLY_INITIALIZERS_FROZEN',objects=len(rows),weights=args.weights,selected_training_epoch=checkpoint['selected_epoch'],feature_mode=checkpoint['feature_mode'],records=receipts,model_load_cpu_seconds=loadcpu,model_load_wall_seconds=loadwall,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-begin,truth_access='NONE',physical_evaluations=0,initial_physical_Phi_rejection=False,shared_setup_amortization='one model load per deployed batch; report cold first use separately',feature_cost_source='generation/feature receipts; public feature preprocessing not free'),indent=2))
    print(json.dumps(dict(status='INITIALIZED',objects=len(rows),output=str(out))))

if __name__=='__main__':main()
