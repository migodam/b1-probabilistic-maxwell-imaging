from pathlib import Path
from dataclasses import replace
import json,numpy as np
from r2.accounting import BudgetBook,plain
from r2.physics import Eval,WallGrid,S1_TX
from compensation import base
ROOT=Path(__file__).resolve().parents[1]
def main():
    b=BudgetBook();out=ROOT/'results/compensation_extended';w0,t0=base();t=replace(t0,pose=t0.pose+[.05,0]);coarse=np.load(ROOT/'results/compensation/BASE.npz')['data'];rows=[]
    for i,f in enumerate([.85,1.,1.15]):
        ev=Eval(b,'E2',f,S1_TX,w0.x,w0.y,wall_refine=2,h=.145,targetstep=.1,max_dofs=120000)
        basefile=ROOT/f'results/compensation/VALIDATED_f{f}.npz'
        basefine=np.load(basefile)['base'] if basefile.exists() else ev.evaluate(w0,t0,label='added-channel independent base').data.copy()
        for name in (['A2'] if f==.85 else ['A1','A2']):
            d=np.load(out/f'{name}_ENDPOINT.npz');w=WallGrid(w0.x,w0.y,d['eps'],d['ell']);s=ev.evaluate(w,t,label=f'{name} independently refined endpoint')
            rows.append(dict(optimized_for=name,f=f,error=float(np.linalg.norm(s.data-d['data'][i])/.01),base_error=float(np.linalg.norm(basefine-coarse[i])/.01),distance_full4=float(np.linalg.norm(s.data-basefine)/.01),distance_first2=float(np.linalg.norm((s.data-basefine)[:2])/.01)))
            np.savez_compressed(out/f'INDEPENDENT_{name}_f{f}.npz',base=basefine,endpoint=s.data)
    (out/'VALIDATION.json').write_text(json.dumps(rows,default=plain,indent=2));print(json.dumps(rows))
if __name__=='__main__':main()
