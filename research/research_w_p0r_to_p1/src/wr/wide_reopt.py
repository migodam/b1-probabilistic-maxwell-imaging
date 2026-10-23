"""Fixed-pair versus reoptimized-pair wide-band PHYSICAL existence control."""
import sys,json
import numpy as np
from .accounting import ROOT,Book,write,BudgetStop
from .phase_diagram import _validated_pair
from .native_world import LinearWorld
from .physical_witness import mesh_for
from .physics import Forward
from .optimization import trust_region

def main(manifest):
    cfg,q,_,gate,_=_validated_pair(manifest);book=Book('C')
    out=ROOT/'results/C/wide_reoptimized';out.mkdir(exist_ok=False)
    freq=json.loads((ROOT/'configs/ORIGINAL_P0.json').read_text())['frequencies']['wide']
    p0=((-0.9,2.15),.12,3.);p1=((-0.8,2.15),.12,3.)
    def fw(probe,label,coarse=False):
        h=.145 if coarse else .13;tol=5e-5 if coarse else 2.5e-5;res=2 if coarse else 4;n=256 if coarse else 512
        return Forward(LinearWorld(cfg['seed'],probe),mesh_for(cfg['seed'],probe,h=h,tol=tol,res=res,n=n),freq,book,
            f'C/wide_reoptimized/{label}',phase='C_reoptimized_wide_pair',identity='ORACLE_EXISTENCE_FIXED_TARGET_SHIFT',
            h=h,pml_strength=6,max_dofs=110000,max_factor_bytes=3500000000,retain_current=not coarse)
    a=fw(p0,'reference');y,_=a.evaluate(np.zeros(108))
    candidate=fw(p1,'candidate');raw,_=candidate.evaluate(q)
    history=[]
    def callback(r):
        history.append(r);write(out/'ITERATIONS.json',history)
        if r.get('event')=='anchor':np.savez_compressed(out/'ACCEPTED_LATENT.npz',x=r['x'])
    fit=trust_region(candidate.residual(y),q,np.full(108,-1.),np.full(108,1.),max_iter=18,max_calls=42,
          radius=.2,max_radius=.7,l2_radius=False,callback=callback,objective_floor=.05,fatal_exceptions=(Exception,))
    d,_=candidate.evaluate(fit.x)
    np.savez_compressed(out/'FINAL.npz',x=fit.x,initial_q=q,data=d,reference=y,initial_data=raw)
    # Independent grid control only for a candidate meeting the physical gate.
    numerical=None
    if np.linalg.norm(d-y)/.01<=1:
        yc,_=fw(p0,'reference_coarse',True).evaluate(np.zeros(108));dc,_=fw(p1,'candidate_coarse',True).evaluate(fit.x)
        numerical=dict(world0_grid_white=float(np.linalg.norm(y-yc)/.01),world1_grid_white=float(np.linalg.norm(d-dc)/.01),coarse_distance=float(np.linalg.norm(dc-yc)/.01))
    result=dict(fixed_pair_wide_distance=float(np.linalg.norm(raw-y)/.01),reoptimized_pair_wide_distance=float(np.linalg.norm(d-y)/.01),
        target_shift=.1,mass_relative_change=0,iterations=fit.nit,stopping=str(fit.message),projected_gradient=float(fit.projected_gradient),
        numerical=numerical,physical_scope='Oracle existence; updated wall not substituted into earlier narrow witness or controlled prior',
        no_global_lower_bound_or_nonexistence_claim=True,totals=book.totals())
    write(out/'RESULTS.json',result);print(json.dumps(result),flush=True)

if __name__=='__main__':main(sys.argv[1])
