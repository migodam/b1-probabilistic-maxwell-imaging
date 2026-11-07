"""One synthetic end-to-end optimizer check, no cohort reads or fitted weights."""
import sys,json,time
from pathlib import Path
import torch,numpy as np
study=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(study/'src'))
from networks import *
torch.manual_seed(7);torch.set_num_threads(2)
device=choose_device('mps')
ops=OperatorCache(study/'data/operators.npz').to(device)
# Synthetic anchor obeys the constitutive relation for one fixed uniform m.
m0=torch.full((2,NC),.12,device=device)
e=torch.stack((ops.E_re,ops.E_im),dim=-1)
j0=m0[...,None,None]*e[None]
x=torch.randn(2,27,24,24,device=device)
x[:,0]=m0.reshape(2,24,24)/MATERIAL_SCALE
stats={'current_scale':.04}
target=m0+torch.linspace(0,.025,NC,device=device)[None]
qnorm=abs_square(ops.total_field(j0)).reshape(2,-1).sum(1)
rows=[]
for arm in ARMS:
    model=ResidualUNet(arm).to(device)
    optim=torch.optim.AdamW(model.parameters(),lr=.001)
    gradients=[]
    for step in range(2):
        optim.zero_grad(set_to_none=True)
        out=proposal(model,ops,x,m0,j0,stats,2,.3)
        loss,_=material_loss(out['m_proposal'],target)
        regs={}
        if arm!='v_physics':
            for name in ('const','internal'):
                r=public_regularizer(out,j0,qnorm,ops,name)
                regs[name]=float(r.detach().cpu())
                loss=loss+.01*r
        loss.backward()
        hg=float(torch.sqrt(sum(p.grad.square().sum() for p in model.head.parameters())).detach().cpu())
        sg=float(torch.sqrt(sum(p.grad.square().sum() for p in model.stem.parameters())).detach().cpu())
        if not np.isfinite(hg+sg) or hg<=0 or (step==1 and sg<=0): raise RuntimeError((arm,step,hg,sg))
        gradients.append({'step':step,'head_l2':hg,'stem_l2':sg,'loss':float(loss.detach().cpu())})
        optim.step()
    rows.append({'arm':arm,'stage_materials_shape':list(out['stage_materials'].shape),
                 'stage_currents_shape':list(out['stage_currents'].shape) if 'stage_currents'in out else None,
                 'gradients':gradients,'regularizers':regs,'parameter_count':sum(p.numel() for p in model.parameters())})
receipt={'status':'PASS','device':str(device),'cohort_or_truth_access':'NONE_SYNTHETIC_ONLY','old_weights_loaded':False,
         'model_weights_retained':False,'projectors':ops.projector_metadata,'arms':rows}
(study/'audit/training_wiring.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
