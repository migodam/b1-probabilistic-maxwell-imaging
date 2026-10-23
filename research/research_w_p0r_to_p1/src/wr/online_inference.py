"""Measurement-only inversion. This module imports no true-world generator.

All initialization/support scores depend on the data and public antenna/grid
information. Oracle-near diagnostics run through a separate module.
"""
import sys,json
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from scipy.ndimage import gaussian_filter,binary_closing,binary_opening,label
from wq.fem import line_green
from wq.inverse_support import PixelGrid,pixel_mesh
from .physics import Forward,PixelModel,pack,packj
from .accounting import ROOT,Book,write,BudgetStop
from .optimization import trust_region,lbfgsb

class GroupPixels(PixelModel):
    """Feature-sized square tile amplitudes, selected using measurement only.

Each 2x2 tile has side .16667/.1875, fully excluded if intersecting query T.
Its support is not asserted to cover every legal wall. Corner-only contacts
are removed from the data-selected support. It is a restricted inverse class.
"""
    def __init__(self,grid,groups):
        super().__init__(grid); self.groups=np.asarray(groups,int);self.n=2*len(groups)
        self.G=np.zeros((grid.n_pixels,len(groups)))
        for k,(iy,ix) in enumerate(groups):
            for y in [2*iy,2*iy+1]:
                for x in [2*ix,2*ix+1]:self.G[y*grid.nx+x,k]=1.
    def full(self,x):return np.r_[8*(self.G@x[:len(self.groups)]),.15*(self.G@x[len(self.groups):])]
    def jacobian(self,p,s,x):
        j=__import__('wq.inverse_support',fromlist=['volume_pixel_jacobian']).volume_pixel_jacobian(p,s,self.grid)
        n=self.grid.n_pixels
        return np.c_[8*(j[:,:n]@self.G),.15*(j[:,n:]@self.G)],0.
    def material_metadata(self,x):return dict(model='measurement_selected_feature_tiles',groups=self.groups,physical_parameters=self.full(x),prior_status='restricted prospective feature-tile class, no truth mask')

def born_initialize(y,frequencies,g):
    xx,yy=np.meshgrid((g.x_edges[:-1]+g.x_edges[1:])/2,(g.y_edges[:-1]+g.y_edges[1:])/2)
    pts=np.c_[xx.ravel(),yy.ravel()];tx=np.c_[np.linspace(-1.4,1.4,4),np.full(4,-.3)];rx=np.c_[np.linspace(-1.5,1.5,12),np.full(12,-.3)]
    scale=np.r_[np.full(g.n_pixels,8.),np.full(g.n_pixels,.15)];js=[];direct=[]
    for f in frequencies:
        k=2*np.pi*f;e=line_green(pts,tx,k);r=k*k*g.cell_area*line_green(rx,pts,k)
        a=np.einsum('rn,ns->srn',r,e).reshape(48,g.n_pixels);a[:,g.calibration_mask(pts)]=0
        js.append(np.c_[a,1j*a/f]*scale[None,:]);direct.append(line_green(rx,tx,k).T)
    j=packj(np.concatenate(js,axis=0));b=pack(y-np.stack(direct))
    z=j.T@cho_solve(cho_factor(j@j.T+.02*np.eye(len(b)),check_finite=False),b,check_finite=False)
    z=np.maximum(z,0);z*=min(1.,.35/max(float(z.max()),1e-30));return z

def selected_groups(z,g):
    energy=(z[:g.n_pixels]**2+z[g.n_pixels:]**2).reshape(g.ny,g.nx)
    scores=energy.reshape(g.ny//2,2,g.nx//2,2).sum((1,3));ny,nx=scores.shape
    yy,xx=np.indices(scores.shape)
    xmin=-2+xx*8/g.nx;xmax=xmin+8/g.nx;ymin=yy*6/g.ny;ymax=ymin+6/g.ny
    forbidden=(xmax>-1.2)&(xmin<1.2)&(ymax>1.6)&(ymin<2.8)
    scores[forbidden]=-np.inf; maxgroups=int(3/(4*g.cell_area));ids=np.argsort(-scores.ravel(),kind='stable')[:maxgroups]
    mask=np.zeros(scores.shape,bool);mask.ravel()[ids]=True
    # Isolated tiles and diagonal contacts cannot be read as wall topology.
    # Keep connected 4-neighbor components; do not silently add pixels/area.
    labels,n=label(mask); counts=np.bincount(labels.ravel());keep=counts>=2;keep[0]=False
    if np.any(keep):mask=keep[labels]
    return np.argwhere(mask)

def run_one(seed,band,method='TR',iterations=16,max_calls_override=None):
    if not (ROOT/'results/R2/COMPLETE.json').exists():raise RuntimeError('derivative audit missing')
    cfg=json.loads((ROOT/'configs/ORIGINAL_P0.json').read_text());freq=cfg['frequencies'][band]
    # This file contains measurement arrays only, not material/GT coordinates.
    with np.load(ROOT/f'inputs/mean_{seed}_{band}.npz') as z:y=z['data'].copy()
    book=Book('R');g=PixelGrid();z=born_initialize(y,freq,g)
    groups=selected_groups(z,g);model=GroupPixels(g,groups);m=len(groups)
    x0=np.r_[model.G.T@z[:g.n_pixels]/4,model.G.T@z[g.n_pixels:]/4]
    folder=f'R3/online_{seed}_{band}_{method}';fw=Forward(model,pixel_mesh(g),freq,book,folder,phase='R3_online_FWI',identity='ONLINE_INFERENCE')
    dest=ROOT/'results'/folder;hist=[]
    def callback(row):hist.append(row);write(dest/'ITERATIONS.json',hist)
    limit=2*iterations+10 if max_calls_override is None else min(2*iterations+10,int(max_calls_override))
    if method=='TR':
        result=trust_region(fw.residual(y),x0,np.zeros(model.n),np.ones(model.n),prior_lambda=.002,max_iter=iterations,radius=.15,l2_radius=False,max_radius=.5,max_calls=limit,callback=callback,objective_floor=.5,fatal_exceptions=(BudgetStop,))
    else:
        result=lbfgsb(fw.residual(y),x0,np.zeros(model.n),np.ones(model.n),prior_lambda=.002,max_iter=iterations,max_calls=limit,callback=callback,fatal_exceptions=(BudgetStop,))
    d,_=fw.evaluate(result.x,jac=False)
    np.savez_compressed(dest/'FINAL.npz',x=result.x,full_parameters=model.full(result.x),data=d,initial=x0,groups=groups)
    summary=dict(seed=seed,band=band,method=method,identity='ONLINE_INFERENCE',d_mean=float(np.linalg.norm(d-y)/.01),objective=float(result.fun),projected_gradient=float(result.projected_gradient),iterations=result.nit,success=bool(result.success),stopping=str(result.message),group_count=m,occupied_area_max=float(4*m*g.cell_area),geometry_status='prospective feature-tile validator still required',rhs_totals=book.totals())
    write(dest/'RESULTS.json',summary);print(json.dumps(summary),flush=True)

if __name__=='__main__':run_one(int(sys.argv[1]),sys.argv[2],*(sys.argv[3:4]),**({'iterations':int(sys.argv[4])} if len(sys.argv)>4 else {}))
