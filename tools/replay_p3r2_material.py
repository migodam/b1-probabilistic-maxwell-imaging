"""Portable deterministic material-risk replay: zero Maxwell queries.

Use extracted Release evidence plus the published source revision. This does
not reopen a sealed experiment or turn an exposed observation into new data.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.risk_models import load_predictor
from p3r2.residual_risk import load_residual_predictor
from p3r2.fusion import passive_box


def main():
    p=argparse.ArgumentParser();p.add_argument('--bundle-root',type=Path,required=True)
    p.add_argument('--parent',default='R2_locked_00000');p.add_argument('--method',default='fused',choices=['fused','direct','residual'])
    p.add_argument('--output',type=Path,required=True);p.add_argument('--compare-saved',action='store_true');a=p.parse_args()
    import torch
    torch.set_num_threads(1)
    if torch.get_num_interop_threads()!=1:torch.set_num_interop_threads(1)
    root=a.bundle_root.resolve()
    if a.output.exists():raise FileExistsError('Do not overwrite a prior replay')
    m=json.loads((root/'METHODS_FROZEN_FOR_LOCKED.json').read_text())
    if a.parent not in m['parents']:raise ValueError('Parent is not in the preserved cohort')
    def relocate(value):
        text=str(value);marker='/research/p3l_r2_risk_posterior_exploration/'
        result=root/text.split(marker,1)[1] if marker in text else root/text
        if not result.resolve().is_relative_to(root):raise ValueError('Replay artifact outside extracted bundle')
        return result
    start=time.perf_counter()
    with np.load(root/'L/contexts'/f'{a.parent}.npz',allow_pickle=False) as z:
        data={k:z[k][None].copy() for k in ('observation','q0','chi0','J','prediction','residual')}
    cfg=json.loads((root/'configs/TARGET.json').read_text())
    if a.method=='fused':
        action=json.loads(relocate(m['fusion_json']).read_text())['selected']
        path=relocate(action['checkpoint']);residual=action['name'].startswith('residual_')
    else:
        action=None;path=relocate(m[a.method+'_checkpoint']);residual=a.method=='residual'
    model=load_residual_predictor(path,device='cpu') if residual else load_predictor(path,device='cpu')
    chi=model.predict_full(data if residual else data['observation'])[0]
    if action is not None:
        if action.get('projection','none')=='passive_rectangle':
            roi=model.p['geometry']['roi'] if residual else model.p['roi_indices']
            chi=passive_box(chi,cfg,roi)
        chi=data['chi0'][0]+float(action['alpha'])*(chi-data['chi0'][0])
    metadata=dict(parent=a.parent,method=a.method,target_id=m['target_id'],source_revision=m['source_revision'],
        physical_queries=0,truth_read=False,wall_seconds=time.perf_counter()-start,
        interpretation='Supervised material-risk estimator replay, not a calibrated posterior mean/sample',
        material_image_probability=None,full_posterior_qualification=False)
    if a.compare_saved:
        with np.load(root/'L/estimates'/f'{a.parent}.npz',allow_pickle=False) as z:previous=z[a.method]
        metadata['RMS_difference_from_preserved_estimate']=float(np.linalg.norm(chi-previous)/24)
        metadata['maximum_difference_from_preserved_estimate']=float(np.max(abs(chi-previous)))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('xb') as f:np.savez_compressed(f,chi=chi,metadata_json=np.array(json.dumps(metadata)))
    a.output.with_suffix('.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(metadata),flush=True)


if __name__=='__main__':main()
