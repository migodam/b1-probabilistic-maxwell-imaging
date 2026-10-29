from pathlib import Path
from dataclasses import replace
import json,numpy as np
from r2.accounting import BudgetBook,plain
from r2.physics import Eval,WallGrid,S1_TX
from compensation import base
ROOT=Path(__file__).resolve().parents[1]
def main():
    b=BudgetBook();out=ROOT/'results/final_pair_refinement';out.mkdir(exist_ok=True);w0,t0=base();q=np.load(ROOT/'results/compensation_final/ENDPOINT.npz');w=WallGrid(w0.x,w0.y,q['eps'],q['ell']);t=replace(t0,pose=t0.pose+[.05,0]);rows=[]
    for i,f in enumerate([.85,1.,1.15]):
        ev=Eval(b,'E2',f,S1_TX,w.x,w.y,wall_refine=3,h=.12,targetstep=.075,max_dofs=100000)
        z0=ev.evaluate(w0,t0,label='third-layer final compensation base').data.copy();s=ev.evaluate(w,t,label='third-layer final compensation endpoint');z1=s.data.copy();old=np.load(ROOT/f'results/compensation_final/INDEPENDENT_f{f}.npz')
        rows.append(dict(f=f,base_refinement_difference=float(np.linalg.norm(z0-old['base'])/.01),endpoint_refinement_difference=float(np.linalg.norm(z1-old['endpoint'])/.01),pair_refinement_difference=float(np.linalg.norm((z1-z0)-(old['endpoint']-old['base']))/.01),distance=float(np.linalg.norm(z1-z0)/.01),ndofs=ev.p.ndofs))
        np.savez_compressed(out/f'f{f}.npz',base=z0,endpoint=z1,mesh_p=ev.p.mesh.p,mesh_t=ev.p.mesh.t,field=s.u)
    result=dict(rows=rows,joint_distance=float(np.sqrt(sum(r['distance']**2 for r in rows))),pair_refinement_difference=float(np.sqrt(sum(r['pair_refinement_difference']**2 for r in rows))),dx=.05,status='NUMERICALLY_SUPPORTED_DISCRETE_WITNESS',continuum_error_bound=False)
    (out/'RESULTS.json').write_text(json.dumps(result,default=plain,indent=2));print(json.dumps(result))
if __name__=='__main__':main()
