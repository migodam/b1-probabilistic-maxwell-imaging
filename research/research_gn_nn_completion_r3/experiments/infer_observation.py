"""Same-acquisition unknown observation -> GN -> NN -> physical top1.

The observation NPZ contains only Y and noise. No labels, masks, amplitudes or
ground truth are accepted by this entry point.
"""
from pathlib import Path
import argparse,json,os,sys,time
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
S=Path(__file__).resolve().parents[1];sys.path.insert(0,str(S.parents[1]))
import numpy as np
from research.research_gn_nn_completion_r3.src.anchors import prepare_reference,generate_anchor
from research.research_gn_nn_completion_r3.src.predict import predict_anchors
from research.research_gn_nn_completion_r3.src.evaluate import evaluate_arm,save_result
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--observation',type=Path,required=True);p.add_argument('--operators',type=Path,required=True);p.add_argument('--checkpoint',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--device',default='auto');a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
reference=a.output/'reference_jacobian.npz';setup=prepare_reference(a.operators,reference)
anchor_path=a.output/'public_anchor.npz';anchor_receipt=generate_anchor(a.operators,a.observation,anchor_path,reference)
pred_dir=a.output/'prediction';predict_anchors(a.checkpoint,[anchor_path],a.operators,pred_dir,a.device)
def read_npz(path):
    with np.load(path,allow_pickle=False)as z:return {k:z[k]for k in z.files}
public={**read_npz(a.operators),**read_npz(a.observation)}
info=json.loads((pred_dir/'public_anchor.json').read_text());prediction=read_npz(pred_dir/'public_anchor.npz')
receipt,arrays=evaluate_arm(public,read_npz(anchor_path),prediction,arm=info['arm'],noise='clean'if float(public['noise'])==0 else'noisy',preparation_cost=info['preparation_cost'])
save_result(a.output/'physical','reconstruction',receipt,arrays)
summary=dict(status=receipt['selection_status'],truth_access='NONE',acquisition='Must match the checkpoint frozen 24x24/12x48 acquisition',whole_process_wall_seconds=time.perf_counter()-start,shared_reference_setup=setup,anchor=anchor_receipt,post=receipt,top1_path=str(a.output/'physical/reconstruction.npz'),material_error='UNKNOWN_WITHOUT_GT')
(a.output/'INFERENCE_SUMMARY.json').write_text(json.dumps(summary,indent=2));print(json.dumps({k:summary[k]for k in ('status','truth_access','whole_process_wall_seconds','top1_path','material_error')}))
