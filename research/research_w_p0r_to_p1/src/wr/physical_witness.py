"""PHYSICAL EXISTENCE SEARCH. The fixed reference world is explicitly public.

Not an online inverse algorithm and never imported by online_inference.
Coefficient edits preserve the reference support exactly; no new pixel holes,
gain fitting, noise enlargement or frequency-specific material is permitted.
"""
import json,sys,gc,time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from p0.geometry import wall_values,probes
from wq.meshes import wall_loops,circle_loop,interface_mesh
from .physics import Forward,coefficient_jacobian,pack,packj
from .accounting import ROOT,Book,ScopedBook,write,BudgetStop
from .oracle_reference import reference
from .optimization import trust_region

class ModulatedWorld:
    def __init__(self,seed,probe=None,shape=(6,4),cavities=()):
        self.seed=seed;self.probe=probe;self.shape=shape;self.m=shape[0]*shape[1];self.n=2*self.m;self.cavities=cavities
    def basis(self,pts):
        return np.column_stack([np.cos(i*np.pi*(pts[:,0]+2)/4)*np.cos(j*np.pi*pts[:,1]/3) for i in range(self.shape[0]) for j in range(self.shape[1])])
    def components(self,x,pts):
        er,si=wall_values(pts,self.seed)
        for center,radius in self.cavities:
            active=np.linalg.norm(pts-np.asarray(center),axis=1)<radius
            er[active]=1.;si[active]=0.
        b=self.basis(pts)
        a=b@x[:self.m];c=b@x[self.m:]
        re=1+(er-1)*(1+.6*np.tanh(a));loss=si*(1+.75*np.tanh(c))
        return re,loss,b,er-1,si,a,c
    def coefficient(self,x,f):
        def eps(pts):
            re,loss,*_=self.components(x,pts)
            if self.probe is not None:
                center,radius,contrast=self.probe;active=np.linalg.norm(pts-np.asarray(center),axis=1)<radius
                re[active]=1+contrast;loss[active]=.02
            return re+1j*loss/f
        return eps
    def derivative(self,x,pts,f):
        re,loss,b,co,si,a,c=self.components(x,pts)
        d=np.c_[co[:,None]*.6/np.cosh(a[:,None])**2*b,1j/f*si[:,None]*.75/np.cosh(c[:,None])**2*b]
        if self.probe is not None:
            center,radius,_=self.probe;d[np.linalg.norm(pts-np.asarray(center),axis=1)<radius]=0
        return d
    def jacobian(self,p,s,x):return coefficient_jacobian(p,s,lambda pts,f:self.derivative(x,pts,f))
    def material_metadata(self,x):return dict(model='smooth_fixed_support_DCT_modulation',seed=self.seed,basis=self.shape,probe=self.probe,cavities=self.cavities,identity='ORACLE_CONSTRUCTED_PHYSICAL_WITNESS_SEARCH',same_support_as_reference=not bool(self.cavities),material_bounds='contrast<=6.4, sigma<=.063, passive globally',new_features='none' if not self.cavities else 'explicit circular cavity, diameter .12; clearance audited separately')

def mesh_for(seed,probe=None,h=.145,tol=5e-5,res=2,n=256,cavities=()):
    loops=wall_loops(seed,tol,res)
    if probe is not None:loops+=[circle_loop(probe[0],probe[1],n)]
    loops += [circle_loop(c,r,n) for c,r in cavities]
    return interface_mesh(h=h,loops=loops)

def full_task_direction(j,k):
    c=j.T@j+1e-2*np.eye(j.shape[1]);a=k.T@k
    vals,vec=eigh(a,c);v=vec[:,-1];v/=np.linalg.norm(v)
    if v[np.argmax(abs(v))]<0:v=-v
    return v,float(vals[-1])

def main(seed,band='narrow',kind='cal'):
    book=ScopedBook('B2',f'cal_{seed}_{band}');cfg=json.loads((ROOT/'configs/ORIGINAL_P0.json').read_text());freq=cfg['frequencies'][band];p=probes()[8]
    out=ROOT/f'results/B2/{kind}_{seed}_{band}';out.mkdir(parents=True,exist_ok=False)
    # Calibration and task share the exact geometric integration mesh.
    model=ModulatedWorld(seed);fw=Forward(model,mesh_for(seed,p),freq,book,str(out.relative_to(ROOT/'results')/'cal'),phase='B2_oracle_cal',identity='PHYSICAL_EXISTENCE_SEARCH',h=.145,pml_strength=6,retain_current=True)
    task=ModulatedWorld(seed,p);ft=Forward(task,mesh_for(seed,p),freq,book,str(out.relative_to(ROOT/'results')/'task'),phase='B2_oracle_task',identity='PHYSICAL_EXISTENCE_SEARCH',h=.145,pml_strength=6)
    zero=np.zeros(model.n);d0,j0=fw.evaluate(zero,jac=True);t0,k0=ft.evaluate(zero,jac=True)
    j=packj(j0);k=packj(k0);v,ratio=full_task_direction(j,k)
    rng=np.random.default_rng(2026102315+seed);rand=rng.normal(size=model.n);rand/=np.linalg.norm(rand)
    grid=np.stack(np.meshgrid(np.linspace(-2,2,201),np.linspace(0,3,151)),axis=-1).reshape(-1,2)
    material0=model.coefficient(zero,1.)(grid)
    dm=model.derivative(zero,grid,1.);metric=np.real(dm.conj().T@dm)/len(grid)
    mv,mq=np.linalg.eigh(metric);valid=mv>mv.max()*1e-8
    T=mq[:,valid]/np.sqrt(mv[valid])[None,:]
    vv,qq=np.linalg.eigh(T.T@(j.T@j)@T);weak=T@qq[:,0];weak/=np.linalg.norm(weak)
    rows=[];np.savez_compressed(out/'BASE.npz',x=zero,cal=d0,task=t0,J=j0,K=k0,direction=v,random_direction=rand,task_ratio=ratio)
    derivatives=[]
    for h in [1e-3,1e-4,1e-5]:
        plus,_=fw.evaluate(h*v);minus,_=fw.evaluate(-h*v)
        err=float(np.linalg.norm((plus-minus).ravel()/(2*h)-j0@v)/np.linalg.norm(j0@v))
        derivatives.append(dict(step=h,relative_error=err));write(out/'DERIVATIVES.json',derivatives)
    if min(r['relative_error'] for r in derivatives)>1e-5:raise RuntimeError('continuous material Full-J mismatch')
    write(out/'SETUP_COST.json',dict(common_cal_forward_J_rhs=80,common_task_forward_J_rhs=80,derivative_audit_rhs=120,task_aware_additional_receiver_J_rhs_vs_common_task_forward=60,all_actual_reservations_in_ledger=True,standalone_cost_note='ordinary controls need baseline task F20 but not its receiver J60; setup must not be given free to task-aware'))
    for distance in [.25,.5,1.0]:
        for method,direction in [('TASK_AWARE',v),('RANDOM_MATCHED',rand),('FULL_J_WEAK',weak)]:
            wall_start=time.perf_counter();cpu_start=time.process_time();rhs_start=book.totals()['B2']['rhs']
            wanted=float(np.sqrt(np.mean(abs(model.coefficient(distance*v,1.)(grid)-material0)**2)))
            if method=='FULL_J_WEAK':
                # Select the least observable direction that can actually
                # reach the paired finite material radius inside the box.
                # This avoids making a saturated/numerical-kernel direction
                # an artificially weak ordinary Full-J control.
                direction=None
                for col in range(qq.shape[1]):
                    candidate=T@qq[:,col];candidate/=np.linalg.norm(candidate)
                    if candidate[np.argmax(abs(candidate))]<0:candidate=-candidate
                    limit=2./np.max(abs(candidate))
                    reach=float(np.sqrt(np.mean(abs(model.coefficient(limit*candidate,1.)(grid)-material0)**2)))
                    if reach>=wanted:direction=candidate;break
                if direction is None:
                    rows.append(dict(method=method,coefficient_progress=distance,status='NOT_RUN_NO_FINITE_RADIUS_DIRECTION',wanted_material_RMS=wanted));write(out/'RESULTS.json',rows);continue
            if method=='TASK_AWARE':target=distance*direction
            else:
                upper=2./np.max(abs(direction))
                def displacement(a):return float(np.sqrt(np.mean(abs(model.coefficient(a*direction,1.)(grid)-material0)**2)))
                if displacement(upper)<wanted:
                    rows.append(dict(method=method,coefficient_progress=distance,status='NOT_RUN_DECODER_RADIUS_UNREACHABLE',wanted_material_RMS=wanted));write(out/'RESULTS.json',rows);continue
                target=brentq(lambda a:displacement(a)-wanted,0,upper,xtol=1e-12)*direction
            target_material=model.coefficient(target,1.)(grid)-material0
            material_norm2=float(np.vdot(target_material,target_material).real)
            raw,_=fw.evaluate(target);rawtask,_=ft.evaluate(target)
            # Maintain explicit signed branch progress; data corrector cannot
            # erase the edit by returning x=0. A linear progress row is not
            # interpreted as likelihood or posterior density.
            weight=100.;rows_iter=[]
            def ev(x,jac=False):
                d,jc=fw.evaluate(x,jac)
                progress=float(np.vdot(target_material,model.coefficient(x,1.)(grid)-material0).real/material_norm2)
                r=np.r_[pack(d-d0),weight*(progress-1.)]
                pg=np.real(target_material.conj()@model.derivative(x,grid,1.))/material_norm2 if jc is not None else None
                jj=None if jc is None else np.vstack([packj(jc),weight*pg])
                return r,jj
            def cb(row):rows_iter.append(row);write(out/f'{method}_{distance}_ITERATIONS.json',rows_iter)
            fit=trust_region(ev,target,np.full(model.n,-2.),np.full(model.n,2.),max_iter=7,max_calls=18,radius=.2,max_radius=.7,l2_radius=False,prior_lambda=0.,callback=cb,objective_floor=.05,fatal_exceptions=(BudgetStop,))
            d,_=fw.evaluate(fit.x);td,_=ft.evaluate(fit.x)
            dc=float(np.linalg.norm(d-d0)/.01);dt=float(np.linalg.norm(td-t0)/.01)
            field=fw.model.coefficient(fit.x,1.)(grid);basefield=fw.model.coefficient(zero,1.)(grid)
            progress=float(np.vdot(target_material,field-basefield).real/material_norm2)
            row=dict(seed=seed,band=band,method=method,coefficient_progress=distance,actual_signed_material_progress=progress,matched_raw_material_RMS=wanted,d_cal=dc,d_T=dt,raw_d_cal=float(np.linalg.norm(raw-d0)/.01),raw_d_T=float(np.linalg.norm(rawtask-t0)/.01),material_RMS=float(np.sqrt(np.mean(abs(field-basefield)**2))),corrector_material_RMS=float(np.sqrt(np.mean(abs(field-model.coefficient(target,1.)(grid))**2))),secant_linear_error_white=float(np.linalg.norm((raw-d0).ravel()-j0@target)/.01),incremental_rhs=book.totals()['B2']['rhs']-rhs_start,incremental_wall_seconds=time.perf_counter()-wall_start,incremental_cpu_seconds=time.process_time()-cpu_start,status='CANDIDATE_REQUIRES_INDEPENDENT_VALIDATION',cal_threshold=bool(dc<=1),task_threshold=bool(dt>=3),identity='ORACLE_CONSTRUCTED_PHYSICAL_WITNESS',stop=str(fit.message))
            np.savez_compressed(out/f'{method}_{distance}_FINAL.npz',x=fit.x,proposal_x=target,cal=d,task=td,cal_base=d0,task_base=t0,direction=direction,material_grid=grid,material=field,reference_material=basefield)
            rows.append(row);write(out/'RESULTS.json',rows);print(json.dumps(row),flush=True)
        if any(r.get('cal_threshold') and r.get('task_threshold') and r.get('coefficient_progress')==distance for r in rows):
            write(out/'NOMINATION.json',dict(reason='first radius completed for all three methods with a provisional physical witness; independent validation next',not_yet_GO=True));break
    write(out/'COMPLETE.json',dict(totals=book.totals(),probabilities=None))

if __name__=='__main__':main(int(sys.argv[1]),*sys.argv[2:])
