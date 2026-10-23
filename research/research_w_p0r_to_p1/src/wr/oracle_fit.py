"""ORACLE R1-C: inverse pixel expression + fully iterated fitting diagnostic."""
import sys,json
import numpy as np
from wq.inverse_support import PixelGrid,pixel_mesh
from .representation import projection
from .oracle_reference import reference
from .physics import Forward,PixelModel
from .optimization import trust_region,lbfgsb
from .accounting import ROOT,Book,write,BudgetStop

def main(seed,band='wide',method='TR',iterations=14):
    book=Book('R');g=PixelGrid();z=projection(seed,g,8)
    occupied=np.flatnonzero(z[:g.n_pixels]+z[g.n_pixels:]);indices=np.r_[occupied,g.n_pixels+occupied]
    model=PixelModel(g,indices);cfg=json.loads((ROOT/'configs/ORIGINAL_P0.json').read_text())
    folder=f'R3/oracle_{seed}_{band}_{method}';fw=Forward(model,pixel_mesh(g),cfg['frequencies'][band],book,folder,phase='R3_oracle_FWI',identity='ORACLE_DIAGNOSTIC')
    dest=ROOT/'results'/folder;y=reference(seed,band);x0=z[indices]*.95
    rows=[]
    def callback(row):rows.append(row);write(dest/'ITERATIONS.json',rows)
    solver=trust_region if method=='TR' else lbfgsb
    kwargs=dict(max_iter=iterations,max_calls=2*iterations+2,callback=callback,fatal_exceptions=(BudgetStop,))
    if method=='TR':kwargs.update(radius=.15,max_radius=.5,l2_radius=False,objective_floor=.5)
    fit=solver(fw.residual(y),x0,np.zeros(model.n),np.ones(model.n),prior_lambda=0.,**kwargs)
    d,_=fw.evaluate(fit.x)
    np.savez_compressed(dest/'FINAL.npz',x=fit.x,initial=x0,full_parameters=model.full(fit.x),data=d,indices=indices)
    row=dict(seed=seed,band=band,method=method,identity='ORACLE_DIAGNOSTIC',initialization='95% of oracle subcell projection, no online claim',d_mean=float(np.linalg.norm(d-y)/.01),objective=float(fit.fun),projected_gradient=float(fit.projected_gradient),success=bool(fit.success),stop=str(fit.message),iterations=fit.nit,minimum_feature='NOT_CERTIFIED, this is expression/optimizer diagnostic',support_area_bound=len(occupied)*g.cell_area,totals=book.totals())
    write(dest/'RESULTS.json',row);print(json.dumps(row),flush=True)
if __name__=='__main__':main(int(sys.argv[1]),*(sys.argv[2:4]),**({'iterations':int(sys.argv[4])} if len(sys.argv)>4 else {}))
