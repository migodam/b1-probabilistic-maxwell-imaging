"""Independent cylindrical Maxwell reference and full-data grid qualification."""
from pathlib import Path
import json,time,sys,gc
import numpy as np
from scipy.integrate import quad
from scipy.special import jv,jvp,hankel1,h1vp
from .accounting import Book,BudgetStop
from .vie import VIE
from .geometry import antennas,material_grid,probes

def disk_fraction(points,h,center,radius):
    result=np.zeros(len(points));x,y=(points-np.asarray(center)).T
    far=np.hypot(x,y)>radius+np.sqrt(2)*h/2
    inside=np.hypot(x,y)<radius-np.sqrt(2)*h/2
    result[inside]=1.
    for i in np.flatnonzero(~far&~inside):
        left,right=max(x[i]-h/2,-radius),min(x[i]+h/2,radius)
        low,high=y[i]-h/2,y[i]+h/2
        def length(t):
            v=np.sqrt(max(radius**2-t*t,0))
            return max(0.,min(high,v)-max(low,-v))
        splits=[left,right]
        for yy in (low,high):
            if abs(yy)<radius:
                xx=np.sqrt(radius**2-yy**2)
                splits += [s for s in [-xx,xx] if left<s<right]
        splits=sorted(set(splits))
        result[i]=sum(quad(length,a,b,epsabs=1e-12)[0] for a,b in zip(splits[:-1],splits[1:]))/h**2
    return np.clip(result,0,1)

def cylinder_reference(rx,k,center,radius,epsilon,angles):
    rel=rx-np.asarray(center);rr=np.linalg.norm(rel,axis=1);theta=np.arctan2(rel[:,1],rel[:,0])
    modes=np.arange(-28,29);index=np.sqrt(epsilon);z=k*radius;zn=index*z
    am=(index*jvp(modes,zn)*jv(modes,z)-jv(modes,zn)*jvp(modes,z))/(jv(modes,zn)*h1vp(modes,z)-index*jvp(modes,zn)*hankel1(modes,z))
    outgoing=hankel1(modes[:,None],k*rr[None,:])*np.exp(1j*modes[:,None]*theta[None,:])
    return np.array([((1j**modes)*np.exp(-1j*modes*angle)*am)@outgoing for angle in angles])

def save_json(path,data):path.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False))

def main(root,mode):
    root=Path(root)
    averaging=mode.endswith('-averaged')
    if averaging:mode=mode[:-9]
    out=root/('results/qualification_source_average' if averaging else 'results/qualification');out.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((root/'w_config.json').read_text());book=Book(root/'results')
    tx,rx=antennas();rows=[];arrays={}
    if mode=='benchmark':
        p=VIE(32,1.,tx,rx,book);eps=1+2*disk_fraction(p.points,p.h,[0,1.5],.3)+.01j*disk_fraction(p.points,p.h,[0,1.5],.3)
        t=time.perf_counter();s=p.solve(eps,jacobian=True)
        payload=dict(nx=32,n=p.N,active=len(s.active),wall_seconds=time.perf_counter()-t,setup_seconds=p.setup_seconds,residual=s.residual,bytes_material_factor=s.factor[0].nbytes,
                     numerical_output_used_for_selection=False,decision='Original hard budget retained; active-only exact LU and matrix blocks keep memory bounded')
        save_json(out/'BENCHMARK.json',payload);print(json.dumps(payload));return
    if mode=='cylinder':
        angles=np.linspace(1.3,1.8,4);epsc=4+.03j;center=[0,1.5];radius=.3
        reference=cylinder_reference(rx,2*np.pi,center,radius,epsc,angles)
        arrays['mie_reference']=reference
        for nx in ([96,128] if averaging else [32,64,96,128]):
            p=VIE(nx,1.,tx,rx,book,phase='repair' if averaging else 'validation',source_cell_average=averaging);fraction=disk_fraction(p.points,p.h,center,radius)
            p.incident=np.stack([np.exp(1j*p.k*((p.points-np.asarray(center))@np.array([np.cos(a),np.sin(a)]))) for a in angles],axis=1)
            if averaging:
                p.incident*=np.array([np.sinc(p.k*np.cos(a)*p.h/(2*np.pi))*np.sinc(p.k*np.sin(a)*p.h/(2*np.pi)) for a in angles])[None,:]
            p.direct=np.zeros_like(p.direct)
            s=p.solve(1+(epsc-1)*fraction)
            diff=s.scattered-reference
            row=dict(nx=nx,N=p.N,active=len(s.active),relative_error=float(np.linalg.norm(diff)/np.linalg.norm(reference)),whitened_error=float(np.linalg.norm(diff)/.01),residual=s.residual,seconds=s.timings['state_seconds'])
            rows.append(row);arrays[f'cylinder_{nx}']=s.scattered
            print(json.dumps(row),flush=True);del p,s;gc.collect()
        save_json(out/'CYLINDER.json',rows);np.savez_compressed(out/'CYLINDER_ARRAYS.npz',**arrays);return
    if mode=='derivatives':
        p=VIE(32,1.,tx,rx,book,phase='repair' if averaging else 'validation',source_cell_average=averaging)
        xy=p.points-np.array([0,1.35]);g=np.exp(-np.sum((xy/.6)**2,axis=1))
        eps=1+(.8+.02j)*g;s=p.solve(eps,jacobian=True)
        rng=np.random.default_rng(2026102205);direction=g*rng.normal(size=p.N);pred=s.jacobian[:,:p.N]@direction
        for h in [1e-3,1e-4,1e-5,1e-6]:
            plus=p.solve(eps+h*direction);minus=p.solve(eps-h*direction)
            fd=(plus.data-minus.data).ravel()/(2*h)
            rows.append(dict(step=h,relative_error=float(np.linalg.norm(fd-pred)/np.linalg.norm(pred))))
        a=rng.normal(size=48)+1j*rng.normal(size=48);v=rng.normal(size=2*p.N)
        left=float(np.real(np.vdot(a,s.jacobian@v)));right=float(np.dot(np.real(s.jacobian.conj().T@a),v))
        result=dict(finite_differences=rows,real_adjoint_relative_error=abs(left-right)/max(abs(left),abs(right),1e-15),maximum_backward_residual=s.residual)
        save_json(out/'DERIVATIVES.json',result);print(json.dumps(result));return
    if mode=='reciprocity':
        # Original center-source and matched cell-average-source versions are
        # compared explicitly; all four native solves are charged to repair.
        swapped_tx=rx[:4].copy();swapped_rx=np.vstack([tx,rx[4:]])
        report=[]
        for averaging in [False,True]:
            p=VIE(64,1.,tx,rx,book,phase='repair',source_cell_average=averaging)
            eps=material_grid(p,cfg['truth_seeds'][0]);s=p.solve(eps)
            q=VIE(64,1.,swapped_tx,swapped_rx,book,phase='repair',source_cell_average=averaging)
            t=q.solve(eps)
            diff=s.data[:,:4]-t.data[:,:4].T
            report.append(dict(source_cell_average=averaging,relative_error=float(np.linalg.norm(diff)/np.linalg.norm(s.data[:,:4])),absolute_error=float(np.linalg.norm(diff)),
                               epsilon_identity=bool(np.array_equal(s.epsilon,t.epsilon))))
        save_json(out/'RECIPROCITY_REPAIR.json',report);print(json.dumps(report));return
    if mode=='assembly-parity':
        old=np.load(out/'WALL_seed2026102201_nx64.npz')
        p=VIE(64,1.,tx,rx,book,phase='repair')
        state=p.solve(old['epsilon_cal'].ravel())
        reference=old['data_1.00000000_cal']
        result=dict(relative_error=float(np.linalg.norm(state.data-reference)/np.linalg.norm(reference)),whitened_difference=float(np.linalg.norm(state.data-reference)/.01),residual=state.residual,
                    change='chunked Fortran matrix, in-place LU, operator residual; algebra unchanged')
        save_json(out/'ASSEMBLY_PARITY.json',result);print(json.dumps(result));return
    if mode=='walls':
        # Only physics qualification, not any candidate/proposal search. Results
        # are saved even if the pre-registered convergence target fails.
        for seed in cfg['truth_seeds']:
            previous={};last={}
            for nx in ([128,192] if averaging else [64,96,128,192]):
                data={};materials={};fields={};timers=[]
                for f in sorted(set(cfg['frequencies']['wide']+cfg['frequencies']['narrow'])):
                    p=VIE(nx,f,tx,rx,book,phase='repair' if averaging else 'validation',source_cell_average=averaging)
                    for label,probe in [('cal',None),('task',probes()[8])]:
                        eps=material_grid(p,seed,probe=probe)
                        s=p.solve(eps)
                        key=f'{f:.8f}_{label}'
                        data[key]=s.data.copy();timers.append(s.timings)
                        if f==1.:
                            materials[label]=eps.reshape(p.ny,p.nx);fields[label]=s.field.reshape(p.ny,p.nx,4)
                        del s
                    del p;gc.collect()
                report={}
                for band,ff in cfg['frequencies'].items():
                    for label in ['cal','task']:
                        current=np.stack([data[f'{f:.8f}_{label}'] for f in ff])
                        if previous:
                            old=np.stack([previous[f'{f:.8f}_{label}'] for f in ff])
                            report[f'{band}_{label}_mesh_whitened']=float(np.linalg.norm(current-old)/.01)
                row=dict(seed=seed,nx=nx,previous_nx=last.get('nx'),mesh=report,state_wall_seconds=sum(t['state_seconds'] for t in timers),status='QUALIFICATION_ONLY_NOT_A_WITNESS')
                rows.append(row);print(json.dumps(row),flush=True)
                np.savez_compressed(out/f'WALL_seed{seed}_nx{nx}.npz',**{f'data_{k}':v for k,v in data.items()},**{f'epsilon_{k}':v for k,v in materials.items()},**{f'field_{k}':v for k,v in fields.items()})
                save_json(out/'WALL_GRIDS.json',rows)
                previous=data;last=row
                # Once a consistent pair is found, no extra refinement is paid.
                if report and max(report.values())<=cfg['numerics']['mesh_whitened_tolerance']:break
        return
    raise ValueError(mode)

if __name__=='__main__':main(Path(__file__).resolve().parents[2],sys.argv[1])
