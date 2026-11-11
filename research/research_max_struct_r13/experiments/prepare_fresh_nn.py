"""Observation-only frozen R8 NN inference for sixteen predeclared fresh objects.

Reads public observations/operators and one existing checkpoint only. No truth,
training, physical forward, objective filtering, inversion or scoring.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse,json,time,sys,traceback
from pathlib import Path


def main():
    cpu0,wall0=time.process_time(),time.perf_counter()
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);a=p.parse_args();root=a.root.resolve();repo=root.parents[1];sys.path.insert(0,str(repo))
    import numpy as np
    import torch
    torch.set_num_threads(1);torch.set_num_interop_threads(1)
    from research.research_max_struct_r8.src.warm_features import FeatureOperator,FEATURE_NAMES
    from research.research_max_struct_r8.src.warm_network import load_predictor
    out=root/'data_public/fresh16_nn';audit=root/'delegated/evidence/nnprep';audit.mkdir(parents=True,exist_ok=True)
    out.mkdir(parents=True,exist_ok=False)
    def dump(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    ids=[f'r13_fresh_{i:04d}' for i in range(0,48,3)]
    protocol=dict(status='PREDECLARED_BEFORE_INFERENCE',case_ids=ids,objects=16,weights='research/research_max_struct_r8/RESULTS/training/physics/weights.pt',initialization='original frozen R8 10-feature physics supervised CNN; material only',training_performed=False,initial_objective_rejection=False,material_bounds=[.02,20.],original_initial_material=.1,truth_access='NONE',CPU_threads=1,BLAS_threads=1)
    dump(audit/'PROTOCOL.json',protocol)
    weight=repo/protocol['weights'];stat=weight.stat();w,c=time.perf_counter(),time.process_time();predict,checkpoint=load_predictor(weight,device='cpu');load=dict(cpu_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w)
    count=sum(x.numel() for x in checkpoint['state_dict'].values())
    if count!=24737:raise ValueError(f'Expected frozen24737parameter architecture, found{count}')
    if checkpoint['feature_mode']!='physics':raise ValueError('Wrong feature mode checkpoint')
    metadata={k:checkpoint[k] for k in ['architecture','feature_mode','selected_epoch','selection','training_objects','validation_objects','seed','mean','std'] if k in checkpoint}
    dump(audit/'CHECKPOINT_METADATA.json',dict(exact_weights_path=str(weight),weights_filename=weight.name,weights_bytes=stat.st_size,weights_mtime_ns=stat.st_mtime_ns,parameter_count=count,checkpoint=metadata,torch_version=torch.__version__,new_training_cost=0,old_training_CPU_seconds=None,old_training_GPU_hours=None,old_training_physical_forward_calls=None,old_training_cost_scope='Inherited supervised prior is additional information; historic training cost not remeasured or treated as zero.',model_load=load))
    cache=[];records=[]
    try:
        for cid in ids:
            wholecpu,wholewall=time.process_time(),time.perf_counter();path=root/'data_public/fresh48'/f'{cid}.npz';w,c=time.perf_counter(),time.process_time()
            with np.load(path,allow_pickle=False) as f:public={k:f[k].copy() for k in f.files}
            inputcost=dict(cpu_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w,bytes=path.stat().st_size,input_file_loads=1)
            original=public['m_init'].copy()
            if original.shape!=(576,) or not np.array_equal(original,np.full(576,.1)):raise ValueError('Original public initialization contract differs')
            w,c=time.perf_counter(),time.process_time();match=None
            for item in cache:
                if item['sources']==public['E0'].shape[1] and item['receivers']==public['S'].shape[0] and all(np.array_equal(public[key],item['public'][key]) for key in ['G','E0','S','coords','grid_shape']):match=item;break
            checkcost=dict(cpu_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w)
            hit=match is not None
            if match is None:
                op=FeatureOperator(public);match=dict(setup_id=f'feature_setup_{len(cache):02d}',sources=int(public['E0'].shape[1]),receivers=int(public['S'].shape[0]),operator=op,public={k:public[k] for k in ['G','E0','S','coords','grid_shape']},first_case_id=cid);cache.append(match)
            op=match['operator'];features,unused_simple_baselines,featurecost=op(public['Y'])
            w,c=time.perf_counter(),time.process_time();raw=predict(features)[0].ravel().astype(float);m=np.clip(raw,.02,20.);infer=dict(cpu_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w)
            if m.shape!=(576,) or not np.isfinite(m).all():raise ValueError('Invalid NN output')
            public['m_init_original_uniform']=original;public['m_init']=m;public['initialization_id']=np.asarray('frozen_R8_physics_24737_CNN');public['initialization_prior']=np.asarray('existing supervised prior; no retraining or objective filtering')
            w,c=time.perf_counter(),time.process_time();dest=out/f'{cid}.npz';np.savez_compressed(dest,**public);writecost=dict(cpu_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w,bytes=dest.stat().st_size)
            record=dict(case_id=cid,input_file=str(path.relative_to(root)),output_file=str(dest.relative_to(root)),feature_setup_id=match['setup_id'],feature_setup_cache_hit=hit,exact_cache_check_keys=['G','E0','S','coords','grid_shape'],input_load=inputcost,cache_equality_check=checkcost,features=featurecost,NN_inference=infer,write=writecost,process_case_CPU_seconds=time.process_time()-wholecpu,process_case_wall_seconds=time.perf_counter()-wholewall,sources=match['sources'],receivers=match['receivers'],original_uniform_initial_value=.1,material_min=float(m.min()),material_max=float(m.max()),clip_count=int(np.count_nonzero(raw!=m)),initial_objective_evaluated=False,initial_objective_rejection=False,truth_access='NONE',LU=0,complex_solve_RHS=0)
            records.append(record);dump(audit/'CASE_COSTS.json',records)
            print(cid,'NN input prepared',flush=True)
        setups=[dict(setup_id=x['setup_id'],sources=x['sources'],receivers=x['receivers'],first_case_id=x['first_case_id'],**x['operator'].setup) for x in cache]
        totals=dict(input_file_loads=len(records),input_bytes=sum(r['input_load']['bytes'] for r in records),output_bytes=sum(r['write']['bytes'] for r in records),receiver_SVD=sum(r['receiver_SVD'] for r in setups),Green_action_columns=sum(r['features']['Green_actions_columns'] for r in records),receiver_adjoint_columns=sum(r['features']['receiver_adjoint_columns'] for r in records),DPCS_pseudoinverse_columns=sum(r['features']['DPCS_pseudoinverse_columns'] for r in records),LU=0,complex_solve_RHS=0,NN_inferences=len(records),model_loads=1)
        receipt=dict(status='READY',objects=16,case_ids=ids,output_directory=str(out),weights_path=str(weight),checkpoint_metadata='CHECKPOINT_METADATA.json',feature_names=FEATURE_NAMES,feature_setups=setups,model_load=load,totals=totals,process_CPU_seconds=time.process_time()-cpu0,process_wall_seconds=time.perf_counter()-wall0,training_performed=False,old_training_CPU_seconds=None,old_training_GPU_hours=None,extra_prior='Inherited R8 supervised training distribution/statistics/weights are additional prior information, not a zero-cost solver.',truth_access='NONE',new_Maxwell_forward=0,new_adjoint_solve=0,new_optimization=0,initial_objective_rejection=False,comparison_scope='Only preprocessing. Parent will independently run NN→DR/G19/direct after freezing its manifest.',cost_scope='Process interval includes imports, NN load, exact cache comparisons, feature setup/evaluation, inference and per-case writes; final receipt/public setup record/manifest writes and exit outside interval. Shared setup charged once here; independent deployment must charge its applicable share and disclose amortization.')
        dump(audit/'PREPROCESS_RECEIPT.json',receipt);dump(out/'PUBLIC_SETUP_RECORD.json',dict(protocol=protocol,feature_setups=setups,checkpoint_metadata=metadata,weights_path=str(weight),parameter_count=count,model_load=load,original_initialization_retained_in_each_NPZ='m_init_original_uniform',old_training_cost=None,records=records,totals=totals))
        dump(out/'MANIFEST.json',dict(status='READY_FOR_PARENT_FROZEN_DISPATCH',case_ids=ids,input_files=[r['output_file'] for r in records],truth_access='NONE',initialization='frozen_R8_physics_CNN',original_case_ids_unchanged=True,new_training=False))
        print(json.dumps(dict(status='READY',objects=16,totals=totals,CPU=receipt['process_CPU_seconds'],wall=receipt['process_wall_seconds'])),flush=True)
    except BaseException:
        dump(audit/'FAILURE.json',dict(error=traceback.format_exc(),completed_cases=len(records),costs=records,CPU_seconds=time.process_time()-cpu0,wall_seconds=time.perf_counter()-wall0,truth_access='NONE'));raise
if __name__=='__main__':main()
