"""Necessary real Full-J, adjoint and boundary-direction diagnostic."""
import json
import numpy as np
from wq.fem import FEM
from wq.inverse_support import PixelGrid,pixel_mesh,pixel_material,volume_pixel_jacobian
from .accounting import ROOT,Book,write
from .physics import pack,packj

def main():
    book=Book('R');g=PixelGrid();mesh=pixel_mesh(g)
    xx,yy=np.meshgrid((g.x_edges[:-1]+g.x_edges[1:])/2,(g.y_edges[:-1]+g.y_edges[1:])/2)
    region=(abs(xx)<.75)&(yy>.4)&(yy<1.15)
    ring=region&~((abs(xx)<.25)&(yy>.65)&(yy<.9))
    cases=[('air',np.zeros_like(region),0.,0.,'one-sided'),('weak',region,.2,.008,'central'),('high',region,7.,.1,'central'),('box_boundary',region,8.,.15,'one-sided'),('area_boundary',(yy<.75),3.,.03,'central'),('cavity',ring,3.5,.025,'central')]
    rows=[];out=ROOT/'results/R2';out.mkdir(parents=True,exist_ok=True)
    tx=np.c_[np.linspace(-1.4,1.4,4),np.full(4,-.3)];rx=np.c_[np.linspace(-1.5,1.5,12),np.full(12,-.3)]
    rng=np.random.default_rng(2026102310)
    for name,mask,c,si,mode in cases:
        theta=g.pack(c*mask,si*mask)
        if name=='air':v=g.pack(.1*region,.002*region)
        elif mode=='one-sided':v=g.pack(-.2*mask,-.003*mask)
        else:v=g.pack(.2*mask*np.cos(3*xx),.003*mask*np.sin(2*yy))
        p=FEM(1.,tx,rx,book,h=.18,order=3,mesh=mesh,phase='R2_derivative',receiver_mode='volume')
        s=p.solve(pixel_material(theta,1.,g));j=volume_pixel_jacobian(p,s,g)
        pred=j@v;cot=rng.normal(size=s.data.shape)+1j*rng.normal(size=s.data.shape)
        adj=p.vjp(s,cot)
        def direction(pts):
            ids=g.indices(pts);vv=v[:g.n_pixels]+1j*v[g.n_pixels:];z=vv[ids];z[g.calibration_mask(pts)]=0;return z
        lhs=np.real(np.vdot(cot,pred.reshape(s.data.shape)));rhs=np.real(np.vdot(adj,direction(s.material_points)))
        base=s.data.copy();arrays={'theta':theta,'direction':v,'base':base,'J':j,'Jv':pred,'adjoint':adj}
        for h in [1e-3,1e-4,1e-5]:
            a=p.solve(pixel_material(theta+h*v,1.,g));arrays[f'plus_{h}']=a.data.copy()
            if mode=='central':
                b=p.solve(pixel_material(theta-h*v,1.,g));arrays[f'minus_{h}']=b.data.copy();fd=(a.data-b.data)/(2*h)
            else:fd=(a.data-base)/h
            err=float(np.linalg.norm(fd.ravel()-pred)/np.linalg.norm(pred))
            row=dict(case=name,step=h,scheme=mode,relative_error=err,adjoint_relative_defect=float(abs(lhs-rhs)/max(abs(lhs),abs(rhs),1e-30)),status='discrete derivative test only, not topology or projection differentiability')
            rows.append(row);write(out/'RESULTS.json',rows);print(json.dumps(row),flush=True)
            np.savez_compressed(out/f'{name}.npz',**arrays)
        del p,s
    write(out/'COMPLETE.json',dict(totals=book.totals(),central_best_max=max(min(r['relative_error'] for r in rows if r['case']==c) for c in ['weak','high','area_boundary','cavity']),boundary_one_sided_convergence='reported separately'))
if __name__=='__main__':main()
