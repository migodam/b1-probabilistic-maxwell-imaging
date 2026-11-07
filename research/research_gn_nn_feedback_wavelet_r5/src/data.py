"""R5 development reuses R4 TRAIN/VAL only; public inference stays separate."""
from pathlib import Path
import json
import numpy as np
import torch
from research.research_gn_nn_preservation_r4.src.train import (
    PhysicsSupervisedAnchors, load_anchor, to_device,
)
from research.research_gn_nn_preservation_r4.src.networks import fit_feature_stats

STUDY=Path(__file__).resolve().parents[1]
R4=STUDY.parent/'research_gn_nn_preservation_r4'

class DevelopmentAnchors(PhysicsSupervisedAnchors):
    def _read(self,index):
        item=super()._read(index)
        path=self.anchor_root/(self.records[index]['obs_id']+'.npz')
        with np.load(path,allow_pickle=False)as z:
            q=np.asarray(z['Q0'],np.complex64)
            item['q0']=torch.from_numpy(q.copy())
            item['Y']=torch.from_numpy(np.asarray(z['Y'],np.complex64).copy())
            a=float(np.sum(np.abs(z['Y'])**2))
            b=float(np.sum(np.abs(z['j0'])**2))
            maxm=float(np.max(np.abs(z['m0'])))
        # Fixed public-anchor bound, reused at every macro step. Armijo checks
        # actual candidate energy if later material states differ from anchor.
        bound=self.snorm**2/a+(1+maxm*self.gnorm)**2/b
        item['step_size']=torch.tensor(1/(1.5*bound),dtype=torch.float32)
        return item

    def __init__(self,records,stats,gnorm,snorm,anchor_root=None):
        self.gnorm,self.snorm=float(gnorm),float(snorm)
        super().__init__(records,R4/'data',anchor_root or R4/'anchors',stats,'v_physics')

def cohort(anchor_root=None):
    records=json.loads((R4/'data/manifests/evaluation.json').read_text())['records']
    selected=json.loads((R4/'audit/OPTIONAL_FT_CONTRACT.json').read_text())['parent_ids']
    train=[r for r in records if r['split']=='train'and r['parent_id']in selected]
    val=[r for r in records if r['split']=='val']
    anchor_root=Path(anchor_root) if anchor_root else R4/'anchors'
    stats=fit_feature_stats(load_anchor(anchor_root/(r['obs_id']+'.npz'))for r in train)
    with np.load(STUDY.parent/'research_gn_green_error_20261106/results/fixed_spectra.npz')as z:
        gnorm,snorm=float(z['in_s'][0]),float(z['ob_s'][0])
    return train,val,stats,gnorm,snorm

def move(batch,device):
    return {key:to_device(value,device)for key,value in batch.items()}

def run_model(model,ops,batch,stats,steps):
    return model(ops,batch['x'],batch['m0'],batch['j0'],batch['q0'],batch['Y'],
                 stats,batch['step_size'],steps=steps)
