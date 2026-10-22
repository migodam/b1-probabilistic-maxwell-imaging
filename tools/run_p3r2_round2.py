"""Continue three independently motivated lines after initial paid sequence."""
from pathlib import Path
import sys,time,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import OUT,write,save
from run_p3r2_physics_sequence import run

def main():
    while True:
        p=OUT/'CHECKPOINT.json'
        state=json.loads(p.read_text()) if p.exists() else {}
        if state.get('status')=='STOPPED_ERROR':raise RuntimeError('initial phase stopped; no automatic paid replay')
        if state.get('status')=='INITIAL_PHYSICS_SEQUENCE_COMPLETE':break
        time.sleep(10)
    old='../b1-p3l-local-posterior/research/p3l_local_posterior_imaging'
    run('C_policy_modelbased',['tools/train_p3r2_policy.py','--kind','model_based','--job','C_policy_modelbased','--trajectory-jobs',str(OUT/'C/jobs/C_train_r0')])
    run('C_eval_modelbased',['tools/run_p3r2_geometry.py','--job','C_exposed_eval_modelbased','--parents',*[f'DEV_{k:03d}' for k in range(24,32)],'--arms','dynamic_mbrl','--model-policy',str(OUT/'C/policy/FROZEN_MODEL_POLICY.json'),'--split','exposed_eval','--steps','8','--max-wall','600'])
    for parent in ('DEV_000','DEV_030'):
        run('B_partition_'+parent,['tools/run_p3r2_partitioned_slice.py','--job','B_partition_'+parent,'--parent',parent,'--anchor',old+f'/A/anchors/{parent}.npz'])
        run('B_source_'+parent,['tools/run_p3r2_source_smc.py','--job','B_source_'+parent,'--parent',parent])
    run('B_full_metric_pair',['tools/run_p3r2_metric_mala.py','--job','B_metric_pair_DEV030','--initial-cache',old+'/B/joint_gaussian/DEV_030/RAW_MODE_CACHE.npz','--parent','DEV_030','--arm','both','--comparison','equal-rhs'])
    write(OUT/'configs/A_DATA_SCALE_REVISION_002.json',dict(reason='Actual validation material MSE of mediumMLP falls from about.059 ROI at4096 to.038 at16384, while truth-forward is cheap; extend to32768 without changing prior or selection data',new_parents=16384,maximum_new_RHS=131072,master_seed_unchanged=2026102202,locked_data_generated=False,model_families=['mean','ridge','mlp_small','mlp_medium'],CNN_stopping_reason='Higher training cost and no validation superiority in prior scale; retained all completed CNN runs'))
    run('A_data_train_scale2',['tools/generate_p3r2_data.py','train','--start','16384','--count','16384'])
    import numpy as np
    blocks=[]
    for name in ('train_0_16384.npz','train_16384_16384.npz'):
        with np.load(OUT/'A/data'/name) as a:blocks.append({k:a[k].copy() for k in a.files})
    save(OUT/'A/data/train32768.npz',**{k:np.concatenate([v[k] for v in blocks]) for k in blocks[0]})
    write(OUT/'ROUND2_PHYSICS_COMPLETE.json',dict(status='COMPLETE',paid_algorithms='MBRL, partitioned reference, source-sequential assimilation, fullmetric MH, independent training extension; no locked observations yet'))
if __name__=='__main__':main()
