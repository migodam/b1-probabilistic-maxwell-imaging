"""One public observation -> strict GN -> trained R4 proposal -> physical top1."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import argparse,json,sys,time
from pathlib import Path
S=Path(__file__).resolve().parents[1];sys.path.insert(0,str(S.parents[1]))
import numpy as np
from research.research_gn_nn_preservation_r4.src.anchors import generate_anchor,prepare_reference
from research.research_gn_nn_preservation_r4.src.predict import predict_anchors
from research.research_gn_nn_preservation_r4.src.evaluate import evaluate_arm,save_result

def load(path):
    with np.load(path,allow_pickle=False)as z:return {k:z[k]for k in z.files}

def main():
    p=argparse.ArgumentParser();p.add_argument('--observation',type=Path,required=True)
    p.add_argument('--operators',type=Path,required=True);p.add_argument('--checkpoint',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--device',default='mps')
    a=p.parse_args();start=time.perf_counter();a.output.mkdir(parents=True,exist_ok=True)
    ref=a.output/'public_reference.npz';setup=prepare_reference(a.operators,ref)
    anchor_path=a.output/'public_anchor.npz';anchor_cost=generate_anchor(a.operators,a.observation,anchor_path,ref)
    pred_dir=a.output/'public_prediction'
    info=predict_anchors(a.checkpoint,[anchor_path],a.operators,pred_dir,a.device)[0]
    prediction=load(pred_dir/'public_anchor.npz');public={**load(a.operators),**load(a.observation)}
    receipt,arrays=evaluate_arm(public,load(anchor_path),prediction,arm=info['arm'],
        noise='clean'if float(public['noise'])==0 else'noisy',preparation_cost=info['preparation_cost'])
    receipt.update(public_observation=str(a.observation),checkpoint=str(a.checkpoint),
        shared_reference_setup=setup,anchor_construction=anchor_cost,
        cold_pipeline_wall_seconds=time.perf_counter()-start,GT_access='NONE')
    save_result(a.output,'public_result',receipt,arrays)
    print(json.dumps({'status':receipt['selection_status'],'data_residual':None if receipt['top1']is None else receipt['top1']['relative_residual'],
                      'cold_wall_seconds':receipt['cold_pipeline_wall_seconds'],'output':str(a.output),'GT_access':'NONE'}),flush=True)

if __name__=='__main__':main()
