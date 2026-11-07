"""Prespecified exposed-VAL depth trajectory; never chooses deployment depth."""
from pathlib import Path
import sys,argparse,json,csv
import numpy as np
import torch
from threadpoolctl import threadpool_limits
STUDY=Path(__file__).resolve().parents[1];sys.path.insert(0,str(STUDY.parents[1]))
from research.research_gn_nn_feedback_wavelet_r5.experiments.predict_public import infer,R4
from research.research_gn_nn_preservation_r4.src.train import FullVIEBridge

ARMS=['pixel_closed_coupled','haar_closed_coupled','db2_closed_coupled','haar_closed_no_fullfield','haar_free_coupled']

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--anchor-mode',choices=['strict','noiseaware'],default='strict');parser.add_argument('--arms',nargs='+',default=ARMS);parser.add_argument('--tag',default='');args=parser.parse_args()
    root=STUDY/'results/noiseaware' if args.anchor_mode=='noiseaware' else STUDY/'results'
    modelroot=STUDY/'models_noiseaware' if args.anchor_mode=='noiseaware' else STUDY/'models'
    anchorroot=STUDY/'inputs/noiseaware_anchors' if args.anchor_mode=='noiseaware' else R4/'anchors'
    records=json.loads((root/'VAL_MANIFEST.json').read_text())
    parents=[]
    for family,digit in [('digits',2),('digits',8),('islands',None)]:
        r=next(r for r in records if r['family']==family and (digit is None or r.get('digit')==digit) and abs(r['amplitude']-.2)<1e-9)
        parents.append(r['parent_id'])
    chosen=[r for r in records if r['parent_id']in parents]
    folder=root/('depth_'+args.tag if args.tag else 'depth');folder.mkdir(exist_ok=True)
    (folder/'SELECTION.json').write_text(json.dumps({'rule':'firstA=.2 VAL digit2,digit8,islands mother; retainclean/noisy; no quality selection','mothers':parents,'deployment_depth_unchanged':4,'diagnostic_depths':[1,2,4,8,16],'scope':'DEVELOPMENT_ONLY'},indent=2)+'\n')
    paths=[anchorroot/(r['obs_id']+'.npz')for r in chosen]
    # Numerical model inference accepts only public anchor inputs.
    for arm in args.arms:infer(modelroot/arm/'best.pt',paths,R4/'data/operators.npz',folder/(arm+'.npz'),steps=16,batch_size=6)
    # Already disclosed VAL GT is read only after every public trajectory saved.
    truths=[];ys=[]
    for r in chosen:
        with np.load(R4/'data'/r['truth_path'])as z:truths.append(np.array(z['m']).reshape(576))
        with np.load(anchorroot/(r['obs_id']+'.npz'))as z:ys.append(np.array(z['Y']))
    y=np.array(ys);y_pair=torch.from_numpy(np.stack((y.real,y.imag),axis=-1).astype(np.float32))
    bridge=FullVIEBridge(R4/'data/operators.npz');rows=[]
    for arm in args.arms:
        with np.load(folder/(arm+'.npz'))as z:
            for depth in (1,2,4,8,16):
                m=z['stage_materials'][:,depth-1]
                if depth==16:residual=np.array(z['data_residual'])
                else:
                    with torch.no_grad():yp=bridge(torch.from_numpy(m))
                    residual=(torch.linalg.vector_norm((yp-y_pair).flatten(1),dim=1)/torch.linalg.vector_norm(y_pair.flatten(1),dim=1)).numpy()
                for i,r in enumerate(chosen):
                    err=float(np.sqrt(np.mean((m[i]-truths[i])**2)))/r['amplitude']
                    rows.append(dict(arm=arm,obs_id=r['obs_id'],parent_id=r['parent_id'],noise=r['noise'],depth=depth,material_NRMSE=err,data_residual=float(residual[i]),constitutive_defect=float(z['constitutive_trace'][i,depth-1]),current_update=float(z['current_update_trace'][i,depth-1]),energy=float(z['energy_trace'][i,depth-1])))
    csvname='DEPTH_DIAGNOSTICS_'+args.tag+'.csv' if args.tag else 'DEPTH_DIAGNOSTICS.csv'
    with (root/csvname).open('w')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (folder/'PHYSICAL_TRACE_COST.json').write_text(json.dumps(bridge.receipt(),indent=2)+'\n')
    print(json.dumps({'rows':len(rows),'scope':'DEVELOPMENT, notGTdepthselection','physics':bridge.counts}))

if __name__=='__main__':
    torch.set_num_threads(1)
    with threadpool_limits(limits=1):main()
