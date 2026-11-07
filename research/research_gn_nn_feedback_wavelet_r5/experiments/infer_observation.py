"""Measured Y -> GN -> closed CNN -> original full-physics correction; no GT."""
from pathlib import Path
import sys,argparse,time,json
import numpy as np
import torch
from threadpoolctl import threadpool_limits
STUDY=Path(__file__).resolve().parents[1];sys.path.insert(0,str(STUDY.parents[1]))
from research.research_gn_nn_feedback_wavelet_r5.experiments.predict_public import infer
from research.research_gn_nn_preservation_r4.src.anchors import prepare_reference,generate_anchor
from research.research_gn_nn_preservation_r4.src.evaluate import evaluate_arm,save_result

def load(path):
    with np.load(path,allow_pickle=False)as z:return {k:np.array(z[k])for k in z.files}

def main():
    p=argparse.ArgumentParser();p.add_argument('--observation',type=Path,required=True);p.add_argument('--operators',type=Path,required=True);p.add_argument('--checkpoint',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--device',default='mps');a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True);started=time.perf_counter()
    header=torch.load(a.checkpoint,map_location='cpu',weights_only=False)
    mode=header['config'].get('anchor_mode','strict');del header
    reference=a.output/'public_reference.npz'
    setup=prepare_reference(a.operators,reference)
    anchor=a.output/'public_anchor.npz'
    config={'discrepancy_each_scale':True}if mode=='noiseaware'else{}
    anchor_cost=generate_anchor(a.operators,a.observation,anchor,reference,config=config)
    pred_path=a.output/'public_CNN_proposal.npz'
    nn_cost=infer(a.checkpoint,[anchor],a.operators,pred_path,steps=4,batch_size=1,device=a.device)
    prediction=load(pred_path)
    operators,observed=load(a.operators),load(a.observation)
    receipt,arrays=evaluate_arm({**operators,'Y':observed['Y']},load(anchor),{'m_proposal':prediction['m_final'][0]},arm='v_physics',noise='clean'if float(observed['noise'])==0 else'noisy',preparation_cost={'wall_seconds':nn_cost['wall_seconds'],'LU':nn_cost['physical_validation']['counts']['LU'],'actual_complex_RHS':nn_cost['physical_validation']['counts']['actual_complex_RHS'],'complete':True})
    receipt.update(method='GN-plus-closed-CNN',anchor_mode=mode,shared_reference=setup,actual_GN=anchor_cost,actual_NN=nn_cost,cold_total_wall_seconds=time.perf_counter()-started,GT_access='NONE')
    save_result(a.output,'public_final',receipt,arrays)
    print(json.dumps({'status':receipt['selection_status'],'GT_access':'NONE','output':str(a.output),'cold_wall_seconds':receipt['cold_total_wall_seconds']}))

if __name__=='__main__':
    torch.set_num_threads(1)
    with threadpool_limits(limits=1):main()
