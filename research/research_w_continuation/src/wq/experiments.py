"""Controlled qualification jobs; raw complex arrays precede interpretation."""
from pathlib import Path
import json,sys,time,gc,os
import numpy as np
from scipy.special import hankel1
from p0.vie import VIE,source_integral,square_self_green
from p0.geometry import antennas,material_grid,wall_values,in_query,probes
from p0.qualification import disk_fraction,cylinder_reference
from .accounting import Book,encode

ROOT=Path(__file__).resolve().parents[2]
def save(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False,default=encode))
def continuous(seed,f,probe=None):
    def fun(points):
        er,loss=wall_values(points,seed);m=in_query(points)
        er=np.where(m,1.,er);loss=np.where(m,0.,loss)
        if probe is not None:
            c,r,v=probe;m=np.linalg.norm(points-np.asarray(c),axis=1)<=r
            er=np.where(m,1+v,er);loss=np.where(m,.02,loss)
        return er+1j*loss/f
    return fun
def smooth(points):
    q=((points[:,0]/1.2)**2+((points[:,1]-1.2)/.6)**2)
    v=np.zeros(len(points));m=q<1;v[m]=np.exp(1-1/(1-q[m]))
    return 1+(1.5+.025j)*v
def main(mode):
    out=ROOT/'results'/mode
    if out.exists():out=ROOT/'results'/(mode+'_'+os.environ.get('WQ_JOB_ID','repeat'))
    out.mkdir(parents=True,exist_ok=False)
    book=Book(ROOT/'results');tx,rx=antennas();rows=[];arrays={}
    def row(**kw):
        rows.append(kw);save(out/'RESULTS.json',rows)
        np.savez_compressed(out/'ARRAYS.npz',**arrays)
        print(json.dumps(kw,default=encode),flush=True)
    if mode=='attribution':
        # Every material solve, even diagnostic fixtures, uses the one new Book.
        angles=np.linspace(1.3,1.8,4);center=np.array([0,1.5]);radius=.3;ep=4+.03j
        ref=cylinder_reference(rx,2*np.pi,center,radius,ep,angles);arrays['MIE']=ref
        for nx,q,sampling in [(64,4,'fraction'),(96,4,'center'),(96,4,'fraction'),(96,8,'fraction'),(96,16,'fraction'),(128,4,'fraction')]:
            p=VIE(nx,1.,tx,rx,book,quadrature=q,source_cell_average=True)
            frac=disk_fraction(p.points,p.h,center,radius) if sampling=='fraction' else (np.linalg.norm(p.points-center,axis=1)<radius).astype(float)
            p.incident=np.stack([np.exp(1j*p.k*((p.points-center)@np.array([np.cos(a),np.sin(a)]))) for a in angles],axis=1)
            p.incident*=np.array([np.sinc(p.k*np.cos(a)*p.h/(2*np.pi))*np.sinc(p.k*np.sin(a)*p.h/(2*np.pi)) for a in angles])[None,:]
            p.direct*=0;s=p.solve(1+(ep-1)*frac)
            key=f'cylinder_nx{nx}_q{q}_{sampling}';arrays[key]=s.scattered
            row(case='cylinder',nx=nx,q=q,sampling=sampling,relative_error=float(np.linalg.norm(s.scattered-ref)/np.linalg.norm(ref)),white_error=float(np.linalg.norm(s.scattered-ref)/.01),residual=s.residual,seconds=s.timings['state_seconds']);del p,s;gc.collect()
        p=VIE(96,1.,tx,rx,book,source_cell_average=True)
        for seed in (2026102201,2026102202):
            last=None
            for sub in (4,8,16,32):
                ep=material_grid(p,seed,subsamples=sub);s=p.solve(ep)
                arrays[f'wall_{seed}_sub{sub}']=s.data
                row(case='wall_coefficient',seed=seed,nx=96,sub=sub,white_change=None if last is None else float(np.linalg.norm(s.data-last)/.01),residual=s.residual,seconds=s.timings['state_seconds']);last=s.data.copy();del s
        del p;gc.collect()
        for nx in (48,64,96,128):
            p=VIE(nx,1.,tx,rx,book,source_cell_average=True);s=p.solve(smooth(p.points))
            arrays[f'smooth_{nx}']=s.data;row(case='smooth',nx=nx,relative_to_previous=None if nx==48 else float(np.linalg.norm(s.data-last)/np.linalg.norm(s.data)),white_change=None if nx==48 else float(np.linalg.norm(s.data-last)/.01),residual=s.residual,seconds=s.timings['state_seconds']);last=s.data.copy();del p,s;gc.collect()
    elif mode=='kernels':
        h=4/96;k=2*np.pi
        for name,offset in [('self',(0.,0.)),('nearest',(h,0.)),('diagonal',(h,h)),('far',(5*h,3*h))]:
            dst=np.asarray([offset]);prev=None
            exact=square_self_green(k,h) if name=='self' else source_integral(dst,np.zeros((1,2)),k,h,96)[0,0]
            for q in (4,8,16,32,64):
                z=source_integral(dst,np.zeros((1,2)),k,h,q)[0,0]
                row(cell=name,q=q,real=float(z.real),imag=float(z.imag),relative_error=float(abs(z-exact)/abs(exact)),method_reference='analytic radial square integral' if name=='self' else 'independent higher-order GL96')
    elif mode=='fem-cylinder':
        from .fem import FEM,PMLConfig
        angles=np.linspace(1.3,1.8,4);center=np.array([0.,1.5]);radius=.3;ep=4+.03j
        ref=cylinder_reference(rx,2*np.pi,center,radius,ep,angles);arrays['MIE']=ref
        def cyl(pts):return 1+(ep-1)*(np.linalg.norm(pts-center,axis=1)<radius)
        for h,order in [(0.12,2),(0.08,2),(0.06,2)]:
            try:p=FEM(1.,tx,rx,book,h=h,order=order)
            except MemoryError as e:
                row(case='cylinder',route='FEM',h=h,order=order,status='NOT_RUN_MEMORY_PREFLIGHT',reason=str(e));continue
            s=p.solve(cyl,angles=angles,plane_wave_origin=center)
            arrays[f'data_h{h}_p{order}']=s.scattered
            row(case='cylinder',route='FEM',h=h,order=order,relative_error=float(np.linalg.norm(s.scattered-ref)/np.linalg.norm(ref)),white_error=float(np.linalg.norm(s.scattered-ref)/.01),residual=s.residual,timings=s.timings);del p,s;gc.collect()
    elif mode=='galerkin-cylinder':
        from .vie_upgrade import VIEUpgrade
        angles=np.linspace(1.3,1.8,4);center=np.array([0.,1.5]);radius=.3;ep=4+.03j
        ref=cylinder_reference(rx,2*np.pi,center,radius,ep,angles);arrays['MIE']=ref
        for nx in (64,96,128):
            p=VIEUpgrade(nx,1.,tx,rx,book,field_mode='active',source_cell_average=True,field_average_axes=(True,True),field_quadrature=4)
            frac=disk_fraction(p.points,p.h,center,radius)
            p.incident=np.stack([np.exp(1j*p.k*((p.points-center)@np.array([np.cos(a),np.sin(a)]))) for a in angles],axis=1)
            p.incident*=np.array([np.sinc(p.k*np.cos(a)*p.h/(2*np.pi))*np.sinc(p.k*np.sin(a)*p.h/(2*np.pi)) for a in angles])[None,:]
            p.direct*=0;s=p.solve(1+(ep-1)*frac)
            arrays[f'data_nx{nx}']=s.scattered
            row(case='cylinder',route='double_average_VIE',nx=nx,relative_error=float(np.linalg.norm(s.scattered-ref)/np.linalg.norm(ref)),white_error=float(np.linalg.norm(s.scattered-ref)/.01),residual=s.residual,seconds=s.timings);del p,s;gc.collect()
    elif mode=='fitted-cylinder':
        from .fem import FEM,PMLConfig
        from .meshes import interface_mesh,circle_loop
        angles=np.linspace(1.3,1.8,4);center=np.array([0.,1.5]);radius=.3;ep=4+.03j
        ref=cylinder_reference(rx,2*np.pi,center,radius,ep,angles);arrays['MIE']=ref
        def cyl(pts):return 1+(ep-1)*(np.linalg.norm(pts-center,axis=1)<radius)
        for h,order,n in [(0.16,2,128),(0.16,3,256),(0.12,3,512)]:
            mesh=interface_mesh(h=h,loops=[circle_loop(center,radius,n)])
            p=FEM(1.,tx,rx,book,h=h,order=order,mesh=mesh)
            s=p.solve(cyl,angles=angles,plane_wave_origin=center)
            arrays[f'data_h{h}_p{order}_n{n}']=s.scattered
            row(case='cylinder',route='fitted_FEM',h=h,order=order,polygon_nodes=n,ndofs=p.ndofs,relative_error=float(np.linalg.norm(s.scattered-ref)/np.linalg.norm(ref)),white_error=float(np.linalg.norm(s.scattered-ref)/.01),residual=s.residual,timings=s.timings);del p,s;gc.collect()
    else:raise ValueError(mode)
    np.savez_compressed(out/'ARRAYS.npz',**arrays)
    save(out/'TOTALS.json',book.totals())
if __name__=='__main__':main(sys.argv[1])
