"""Fitted FEM independent references and full original wall qualification."""
from pathlib import Path
import json,time,sys,gc,os
import numpy as np
from scipy.special import hankel1,jv,jvp,h1vp
from .accounting import Book,encode
from .fem import FEM,PMLConfig
from .meshes import interface_mesh,circle_loop,wall_loops
from .experiments import continuous,smooth
from p0.geometry import antennas,probes
from p0.qualification import cylinder_reference
ROOT=Path(__file__).resolve().parents[2]

def line_mie(tx,rx,k,center,radius,ep):
    # Graf addition for exterior line sources, independent of FEM/VIE.
    m=np.arange(-32,33);n=np.sqrt(ep);z=k*radius
    a=(n*jvp(m,n*z)*jv(m,z)-jv(m,n*z)*jvp(m,z))/(jv(m,n*z)*h1vp(m,z)-n*jvp(m,n*z)*hankel1(m,z))
    rs=tx-np.asarray(center);rr=rx-np.asarray(center)
    incoming=.25j*hankel1(m[:,None],k*np.linalg.norm(rs,axis=1)[None,:])*np.exp(-1j*m[:,None]*np.arctan2(rs[:,1],rs[:,0])[None,:])
    outgoing=hankel1(m[:,None],k*np.linalg.norm(rr,axis=1)[None,:])*np.exp(1j*m[:,None]*np.arctan2(rr[:,1],rr[:,0])[None,:])
    return (outgoing.T@(a[:,None]*incoming)).T

def main(mode):
    out=ROOT/'results'/mode
    if out.exists():out=ROOT/'results'/(mode+'_'+os.environ.get('WQ_JOB_ID','repeat'))
    out.mkdir(parents=True,exist_ok=False);book=Book(ROOT/'results');tx,rx=antennas();rows=[];arrays={}
    def row(**kw):
        rows.append(kw);(out/'RESULTS.json').write_text(json.dumps(rows,indent=2,default=encode,allow_nan=False))
        np.savez_compressed(out/'ARRAYS.npz',**arrays);print(json.dumps(kw,default=encode),flush=True)
    def fitted(h,loops,pml=None,q=8,order=3):
        pm=PMLConfig() if pml is None else pml
        mesh=interface_mesh(h=h,outer_box=pm.outer_box,physical_box=pm.physical_box,loops=loops)
        return FEM(1.,tx,rx,book,h=h,order=order,mesh=mesh,pml=pm,quadrature_order=q)
    if mode=='independent':
        center=np.array([0.,1.5]);radius=.3;ep=4+.03j;angles=np.linspace(1.3,1.8,4)
        def cyl(p):return 1+(ep-1)*(np.linalg.norm(p-center,axis=1)<radius)
        ref=cylinder_reference(rx,2*np.pi,center,radius,ep,angles);lr=line_mie(tx,rx,2*np.pi,center,radius,ep)
        arrays['plane_Mie']=ref;arrays['line_Mie']=lr
        for h,n in [(.16,512),(.13,256),(.13,512)]:
            try:p=fitted(h,[circle_loop(center,radius,n)])
            except MemoryError as e:row(case='cylinder',h=h,n=n,status='NOT_RUN_MEMORY_PREFLIGHT',reason=str(e));continue
            s=p.solve(cyl,angles=angles,plane_wave_origin=center);l=p.solve(cyl)
            key=f'h{h}_n{n}';arrays[key+'_plane']=s.scattered;arrays[key+'_line']=l.scattered
            row(case='cylinder',h=h,n=n,ndofs=p.ndofs,plane_relative=float(np.linalg.norm(s.scattered-ref)/np.linalg.norm(ref)),plane_white=float(np.linalg.norm(s.scattered-ref)/.01),line_relative=float(np.linalg.norm(l.scattered-lr)/np.linalg.norm(lr)),line_white=float(np.linalg.norm(l.scattered-lr)/.01),metadata=l.metadata,timings=l.timings,residual=max(s.residual,l.residual));del p,s,l;gc.collect()
        # At fixed fitted mesh size compare PML and exterior expansion separately.
        for name,pm in [('strength6',PMLConfig(strength=6)),('thickness1.3',PMLConfig(thickness=1.3)),('physical_expanded',PMLConfig(physical_box=(-2.8,2.8,-1.1,3.8)))]:
            p=fitted(.16,[circle_loop(center,radius,256)],pml=pm);s=p.solve(cyl,angles=angles,plane_wave_origin=center)
            arrays['pml_'+name]=s.scattered;row(case='cylinder_PML',variant=name,h=.16,plane_relative=float(np.linalg.norm(s.scattered-ref)/np.linalg.norm(ref)),metadata=s.metadata,residual=s.residual);del p,s;gc.collect()
        # Smooth field: independent FE route on same smooth coefficient as VIE.
        for h in (.18,.15,.13):
            p=fitted(h,[]);s=p.solve(smooth);arrays[f'smooth_h{h}']=s.data
            row(case='smooth',h=h,ndofs=p.ndofs,white_change=None if h==.18 else float(np.linalg.norm(s.data-last)/.01),metadata=s.metadata,residual=s.residual);last=s.data.copy();del p,s;gc.collect()
        p=fitted(.18,[]);s=p.solve(lambda x:np.ones(len(x)));arrays['vacuum']=s.data
        row(case='vacuum',scattered_norm=float(np.linalg.norm(s.scattered)),direct_identity=float(np.linalg.norm(s.data-s.direct)),residual=s.residual);del p,s;gc.collect()
    elif mode in ('wall-qualification','wall-qualification-resume'):
        cfg=json.loads((ROOT/'docs/P0_FROZEN_CONFIG.json').read_text());ff=sorted(set(sum(cfg['frequencies'].values(),[])))
        probe=probes()[8]
        for seed in cfg['truth_seeds']:
            last={}
            for h,tol,resolution in [(.18,1e-4,1),(.15,5e-5,1),(.145,5e-5,2)]:
                cached=ROOT/'results/wall-qualification'/f'WORLD_{seed}_h{h}.npz'
                if mode=='wall-qualification-resume' and cached.exists():
                    with np.load(cached) as z:current={k[5:]:z[k].copy() for k in z.files if k.startswith('data_')}
                    metrics={}
                    for band,frequencies in cfg['frequencies'].items():
                        for label in ('cal','task'):
                            a=np.stack([current[f'{f:.8f}_{label}'] for f in frequencies]);arrays[f'seed{seed}_h{h}_{band}_{label}']=a
                            if last:metrics[f'{band}_{label}']=float(np.linalg.norm(a-np.stack([last[f'{f:.8f}_{label}'] for f in frequencies]))/.01)
                    row(case='wall_complete',seed=seed,h=h,tolerance=tol,geometry_resolution=resolution,mesh_white=metrics,status='REUSED_EXISTING_Q_STATE',cache=str(cached.relative_to(ROOT)))
                    last=current;continue
                loops=wall_loops(seed,tol,resolution)+[circle_loop(probe[0],probe[1],256)]
                mesh=interface_mesh(h=h,loops=loops);current={};savefields={}
                for f in ff:
                    p=FEM(f,tx,rx,book,h=h,order=3,mesh=mesh)
                    for label,pr in [('cal',None),('task',probe)]:
                        s=p.solve(continuous(seed,f,pr));key=f'{f:.8f}_{label}';current[key]=s.data.copy()
                        if f==1.:
                            savefields.update({f'{label}_mesh_p':p.mesh.p,f'{label}_mesh_t':p.mesh.t,f'{label}_points':s.material_points,f'{label}_epsilon':s.epsilon,f'{label}_field':s.field,f'{label}_dofs':s.scattered_dofs})
                        row(case='wall_state',seed=seed,h=h,frequency=f,label=label,ndofs=p.ndofs,residual=s.residual,state_seconds=s.timings['state_seconds'],setup_seconds=s.timings['setup_seconds']);del s
                    del p;gc.collect()
                metrics={}
                for band,frequencies in cfg['frequencies'].items():
                    for label in ('cal','task'):
                        a=np.stack([current[f'{f:.8f}_{label}'] for f in frequencies]);arrays[f'seed{seed}_h{h}_{band}_{label}']=a
                        if last:metrics[f'{band}_{label}']=float(np.linalg.norm(a-np.stack([last[f'{f:.8f}_{label}'] for f in frequencies]))/.01)
                np.savez_compressed(out/f'WORLD_{seed}_h{h}.npz',**{f'data_{k}':v for k,v in current.items()},**savefields)
                row(case='wall_complete',seed=seed,h=h,tolerance=tol,geometry_resolution=resolution,mesh_white=metrics,status='QUALIFICATION_ONLY_NOT_WITNESS')
                last=current
    elif mode=='wall-pml':
        cfg=json.loads((ROOT/'docs/P0_FROZEN_CONFIG.json').read_text());ff=sorted(set(sum(cfg['frequencies'].values(),[])));probe=probes()[8]
        for seed in cfg['truth_seeds']:
            loops=wall_loops(seed,5e-5,2)+[circle_loop(probe[0],probe[1],256)]
            pm=PMLConfig(strength=6);mesh=interface_mesh(h=.145,outer_box=pm.outer_box,physical_box=pm.physical_box,loops=loops)
            data={};point={}
            for f in ff:
                p=FEM(f,tx,rx,book,h=.145,order=3,mesh=mesh,pml=pm,receiver_mode='volume')
                for label,pr in [('cal',None),('task',probe)]:
                    s=p.solve(continuous(seed,f,pr));key=f'{f:.8f}_{label}'
                    data[key]=s.data.copy();point[key]=s.probe_data.copy()
                    row(case='wall_PML_state',seed=seed,frequency=f,label=label,volume_probe_white=float(np.linalg.norm(s.data-s.probe_data)/.01),residual=s.residual,ndofs=p.ndofs)
                del p,s;gc.collect()
            refpath=ROOT/'results/wall-qualification-resume'/f'WORLD_{seed}_h0.145.npz'
            with np.load(refpath) as z:reference={k[5:]:z[k].copy() for k in z.files if k.startswith('data_')}
            for band,frequencies in cfg['frequencies'].items():
                for label in ('cal','task'):
                    v=np.stack([data[f'{f:.8f}_{label}'] for f in frequencies]);b=np.stack([point[f'{f:.8f}_{label}'] for f in frequencies]);r=np.stack([reference[f'{f:.8f}_{label}'] for f in frequencies])
                    arrays[f'seed{seed}_{band}_{label}_volume']=v;arrays[f'seed{seed}_{band}_{label}_probe']=b
                    row(case='wall_PML_complete',seed=seed,band=band,label=label,PML_strength_white=float(np.linalg.norm(b-r)/.01),volume_probe_white=float(np.linalg.norm(v-b)/.01),combined_white=float(np.linalg.norm(v-r)/.01))
            # Same continuous material/mesh-size, expanded PML at edge frequencies.
            for f in (.8,1.2):
                compare={}
                for thickness in (1.,1.3):
                    pm=PMLConfig(thickness=thickness);mesh=interface_mesh(h=.18,outer_box=pm.outer_box,physical_box=pm.physical_box,loops=loops)
                    p=FEM(f,tx,rx,book,h=.18,order=3,mesh=mesh,pml=pm,receiver_mode='volume')
                    for label,pr in [('cal',None),('task',probe)]:
                        s=p.solve(continuous(seed,f,pr));compare[thickness,label]=s.data.copy();arrays[f'seed{seed}_f{f}_thick{thickness}_{label}']=s.data
                    del p,s;gc.collect()
                for label in ('cal','task'):
                    row(case='wall_PML_thickness',seed=seed,frequency=f,label=label,white_difference=float(np.linalg.norm(compare[1.,label]-compare[1.3,label])/.01),note='Same h/geometry tolerance, remeshed expanded exterior; field-mesh and PML effects conservatively combined')
    elif mode=='fem-derivatives':
        p=fitted(.18,[]);p.receiver_mode='volume';s=p.solve(smooth);rng=np.random.default_rng(2026102205)
        def direction(pts):return (smooth(pts)-1)*(.2+.1j)*np.cos(2*pts[:,0])
        j=p.jvp(s,direction);r=rng.normal(size=s.data.shape)+1j*rng.normal(size=s.data.shape);g=p.vjp(s,r)
        lhs=float(np.real(np.vdot(r,j)));rhs=float(np.real(np.vdot(g,direction(s.material_points))))
        row(case='real_adjoint',relative_defect=abs(lhs-rhs)/max(abs(lhs),abs(rhs),1e-30))
        for step in (1e-3,1e-4,1e-5,1e-6):
            a=p.solve(lambda pts:smooth(pts)+step*direction(pts));b=p.solve(lambda pts:smooth(pts)-step*direction(pts));fd=(a.data-b.data)/(2*step)
            row(case='FD',step=step,relative_error=float(np.linalg.norm(fd-j)/np.linalg.norm(j)),residual=max(a.residual,b.residual));del a,b
        del p,s;gc.collect()
        pm=PMLConfig();mesh=interface_mesh(h=.18,loops=[])
        p=FEM(1.,tx,rx,book,h=.18,order=3,mesh=mesh,receiver_mode='volume');s=p.solve(smooth)
        q=FEM(1.,rx[:4],np.vstack([tx,rx[4:]]),book,h=.18,order=3,mesh=mesh,receiver_mode='volume');t=q.solve(smooth)
        row(case='reciprocity',relative_error=float(np.linalg.norm(s.data[:,:4]-t.data[:,:4].T)/np.linalg.norm(s.data[:,:4])),probe_relative_error=float(np.linalg.norm(s.probe_data[:,:4]-t.probe_data[:,:4].T)/np.linalg.norm(s.probe_data[:,:4])),note='volume and FE-point extraction independently retained')
    else:raise ValueError(mode)
    (out/'TOTALS.json').write_text(json.dumps(book.totals(),indent=2))
if __name__=='__main__':main(sys.argv[1])
