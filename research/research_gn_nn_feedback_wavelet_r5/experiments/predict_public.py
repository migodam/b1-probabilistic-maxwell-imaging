"""Frozen CNN inference on explicit public GN anchors; no truth or class read."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
import torch
from threadpoolctl import threadpool_limits

STUDY=Path(__file__).resolve().parents[1];sys.path.insert(0,str(STUDY.parents[1]))
from research.research_gn_nn_feedback_wavelet_r5.src.closed_loop import FeedbackCNN,FixedOperators
from research.research_gn_nn_preservation_r4.src.networks import normalized_features
from research.research_gn_nn_preservation_r4.src.train import FullVIEBridge,to_device

R4=STUDY.parent/'research_gn_nn_preservation_r4'

def sync(device):
    if device.type=='mps':torch.mps.synchronize()

def infer(checkpoint,anchor_paths,operators_path,output,steps=4,batch_size=8,device='mps'):
    start=time.perf_counter();device=torch.device(device)
    ckpt=torch.load(checkpoint,map_location='cpu',weights_only=False)
    cfg,stats=ckpt['config'],ckpt['stats']
    model=FeedbackCNN(cfg['representation'],cfg['closed'],noise_ball=cfg.get('noise_ball',False)).to(device)
    model.load_state_dict(ckpt['state_dict']);model.eval()
    ops=FixedOperators(operators_path).to(device)
    bridge=FullVIEBridge(operators_path)
    with np.load(operators_path)as z:snorm=float(z['S_singular_values'][0])
    with np.load(STUDY.parent/'research_gn_green_error_20261106/results/fixed_spectra.npz')as z:gnorm=float(z['in_s'][0])
    out={key:[]for key in ('obs_id','m_final','m_prior','stage_materials','stage_currents','data_residual','qualified','energy_trace','energy_before_trace','current_update_trace','constitutive_trace','gradient_trace','accepted_trace','inference_seconds_amortized')}
    loop_seconds=0.
    with torch.no_grad():
        for first in range(0,len(anchor_paths),batch_size):
            paths=anchor_paths[first:first+batch_size]
            anchors=[]
            for path in paths:
                with np.load(path,allow_pickle=False)as z:anchors.append({k:np.array(z[k])for k in ('m0','j0','Q0','Y','gb','noise')})
            features=torch.stack([normalized_features(a,stats,'v_physics')for a in anchors]).to(device)
            # Cast before device transfer: MPS does not acceptfloat64/complex128.
            # Final exact forward independently uses its originalf64operators.
            inputs={k:to_device(torch.from_numpy(np.stack([a[k]for a in anchors]).astype(
                np.float32 if k=='m0'else np.complex64)),device)for k in ('m0','j0','Q0','Y')}
            taus=[]
            for a in anchors:
                an=float(np.sum(np.abs(a['Y'])**2));bn=float(np.sum(np.abs(a['j0'])**2));dm=float(np.max(np.abs(a['m0'])))*gnorm
                taus.append(1/(1.5*(snorm**2/an+(1+dm)**2/bn)))
            sync(device);t=time.perf_counter()
            pred=model(ops,features,inputs['m0'],inputs['j0'],inputs['Q0'],inputs['Y'],stats,
                       torch.tensor(taus,device=device,dtype=torch.float32),steps=steps)
            sync(device);elapsed=time.perf_counter()-t;loop_seconds+=elapsed
            y=bridge(pred['m_proposal'])
            public_y=inputs['Y']if not inputs['Y'].is_complex()else torch.view_as_real(inputs['Y'])
            rel=torch.linalg.vector_norm((y-public_y).flatten(1),dim=1)/torch.linalg.vector_norm(public_y.flatten(1),dim=1)
            noise=np.array([float(a['noise'])for a in anchors]);tol=np.where(noise>0,.003,.001)
            def cpu(x):return x.detach().cpu().numpy()
            n=len(paths)
            values={'obs_id':np.array([p.stem for p in paths]),'m_final':cpu(pred['m_proposal']),
                    'm_prior':cpu(pred['prior']),'stage_materials':cpu(pred['stage_materials']),
                    'stage_currents':cpu(pred['stage_currents']),
                    'data_residual':cpu(rel),'qualified':cpu(rel)<=tol,
                    'energy_trace':cpu(pred['trace']['energy_after']),
                    'energy_before_trace':cpu(pred['trace']['energy_before_row']),
                    'current_update_trace':cpu(pred['trace']['relative_current_update']),
                    'constitutive_trace':cpu(pred['trace']['relative_constitutive_defect']),
                    'gradient_trace':cpu(pred['trace']['gradient_norm2']),
                    'accepted_trace':cpu(pred['trace']['accepted']),
                    'inference_seconds_amortized':np.full(n,elapsed/n)}
            for k,v in values.items():out[k].append(v)
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output,**{k:np.concatenate(v)for k,v in out.items()})
    receipt={'scope':'PUBLIC_INPUT_ONLY','checkpoint':str(checkpoint),'epoch':ckpt['epoch'],'config':cfg,'steps':steps,
             'cases':len(anchor_paths),'wall_seconds':time.perf_counter()-start,'NN_loop_seconds':loop_seconds,
             'input_GN_cost':'cached actual R4 input, must be paid in independent deployment; not included in inference-only seconds',
             'fixed_operator_calls':ops.call_counts,'physical_validation':bridge.receipt(),'truth_access':'NONE'}
    output.with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt

def main():
    p=argparse.ArgumentParser();p.add_argument('--checkpoint',type=Path,required=True);p.add_argument('--anchors',nargs='+',type=Path,required=True);p.add_argument('--operators',type=Path,default=R4/'data/operators.npz');p.add_argument('--output',type=Path,required=True);p.add_argument('--steps',type=int,default=4);p.add_argument('--batch',type=int,default=8);p.add_argument('--device',default='mps');a=p.parse_args()
    paths=[]
    for f in a.anchors:paths.extend(sorted(f.glob('*.npz'))if f.is_dir()else[f])
    print(json.dumps(infer(a.checkpoint,paths,a.operators,a.output,a.steps,a.batch,a.device),default=str),flush=True)

if __name__=='__main__':
    torch.set_num_threads(1)
    with threadpool_limits(limits=1):main()
