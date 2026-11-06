"""Measure frozen public feature/inference overhead on the solver hardware.

Does not replace initializers or run Maxwell; post-hoc deployment timing only.
"""
import time,os
BEGIN=time.perf_counter();CPU=time.process_time()
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json,sys,platform,argparse
import numpy as np
import torch
from threadpoolctl import threadpool_limits,threadpool_info
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r8.src.warm_features import FeatureOperator
from research.research_max_struct_r8.src.warm_network import load_predictor

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--manifest',default='warm_xinan_manifest.json');ap.add_argument('--output',default='XINAN_FEATURE_INFERENCE.json');args=ap.parse_args()
    torch.set_num_threads(1);imports_cpu=time.process_time()-CPU;imports_wall=time.perf_counter()-BEGIN
    tasks=json.loads((ROOT/'configs'/args.manifest).read_text())['tasks'];cases={t['case_id']:t for t in tasks if t.get('initialization')=='nn' and t['arm']=='a_plus'}
    records=[];loads={};predictors={};ops={};setups={}
    for mode in ('physics','born'):
        t,c=time.perf_counter(),time.process_time();predictors[mode],_=load_predictor(ROOT/f'RESULTS/training/{mode}/weights.pt');loads[mode]=dict(cpu_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t)
    for ident,t in sorted(cases.items()):
        public={}
        for filename in t['public_files']:
            with np.load(REPO/filename) as z:public.update({k:z[k].copy() for k in z.files})
        public.setdefault('E0',public.get('E'));key=t['public_files'][0]
        if key not in ops:ops[key]=FeatureOperator(public);setups[key]=ops[key].setup
        ft,fc=time.perf_counter(),time.process_time()
        for repeat in range(32):X,_,feature=ops[key](public['Y'])
        feature.update(cpu_seconds=(time.process_time()-fc)/32,wall_seconds=(time.perf_counter()-ft)/32,timing_repetitions=32)
        row=dict(case_id=ident,operator_key=key,feature=feature,models={})
        for mode,predict in predictors.items():
            start,cpu=time.perf_counter(),time.process_time()
            for repeat in range(64):m=predict(X)[0].ravel()
            measurement=dict(cpu_seconds=(time.process_time()-cpu)/64,wall_seconds=(time.perf_counter()-start)/64,timing_repetitions=64)
            if mode=='physics':
                with np.load(REPO/t['initialization_file']) as z:old=z['m_init'].ravel()
                measurement['maximum_absolute_difference_from_frozen_Mac_initializer']=float(np.max(abs(m-old)))
            row['models'][mode]=measurement
        records.append(row)
    rec=dict(status='INDEPENDENT_SAME_HARDWARE_DEPLOYMENT_OVERHEAD_MEASURED',host=platform.node(),system=platform.system(),torch=torch.__version__,records=records,operator_setup=setups,model_load=loads,imports_cpu_seconds=imports_cpu,imports_wall_seconds=imports_wall,cpu_seconds=time.process_time()-CPU,wall_seconds=time.perf_counter()-BEGIN,threads=threadpool_info(),CPU_threads=1,LU=0,complex_RHS=0,truth_access='NONE',initializers_replaced=False,scope='Timing counterfactual local deployment of identical frozen feature/network functions; actual primary initializers remain original Mac-generated files. Arithmetic backend differences recorded, not used for selection.')
    (ROOT/'audit'/args.output).write_text(json.dumps(rec,indent=2));print(json.dumps(dict(status=rec['status'],objects=len(records),CPU=rec['cpu_seconds'],wall=rec['wall_seconds'])))

if __name__=='__main__':
    with threadpool_limits(limits=1):main()
