"""PHYSICAL EXISTENCE SEARCH. The fixed reference world is explicitly public.

Not an online inverse algorithm and never imported by online_inference.
Coefficient edits preserve the reference support exactly; no new pixel holes,
gain fitting, noise enlargement or frequency-specific material is permitted.
"""
import json,sys,gc
import numpy as np
from scipy.linalg import eigh
from p0.geometry import wall_values,probes
from wq.meshes import wall_loops,circle_loop,interface_mesh
from .physics import Forward,coefficient_jacobian,pack,packj
from .accounting import ROOT,Book,write,BudgetStop
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
    book=Book('B2');cfg=json.loads((ROOT/'configs/ORIGINAL_P0.json').read_text());freq=cfg['frequencies'][band];p=probes()[8]
    out=ROOT/f'results/B2/{kind}_{seed}_{band}';out.mkdir(parents=True,exist_ok=False)
    # Calibration and task share the exact geometric integration mesh.
    model=ModulatedWorld(seed);fw=Forward(model,mesh_for(seed,p),freq,book,str(out.relative_to(ROOT/'results')/'cal'),phase='B2_oracle_cal',identity='PHYSICAL_EXISTENCE_SEARCH',h=.145,pml_strength=6)
    task=ModulatedWorld(seed,p);ft=Forward(task,mesh_for(seed,p),freq,book,str(out.relative_to(ROOT/'results')/'task'),phase='B2_oracle_task',identity='PHYSICAL_EXISTENCE_SEARCH',h=.145,pml_strength=6)
    zero=np.zeros(model.n);d0,j0=fw.evaluate(zero,jac=True);t0,k0=ft.evaluate(zero,jac=True)
    j=packj(j0);k=packj(k0);v,ratio=full_task_direction(j,k)
    rng=np.random.default_rng(2026102315+seed);rand=rng.normal(size=model.n);rand/=np.linalg.norm(rand)
    rows=[];np.savez_compressed(out/'BASE.npz',x=zero,cal=d0,task=t0,J=j0,K=k0,direction=v,random_direction=rand,task_ratio=ratio)
    derivatives=[]
    for h in [1e-3,1e-4,1e-5]:
        plus,_=fw.evaluate(h*v);minus,_=fw.evaluate(-h*v)
        err=float(np.linalg.norm((plus-minus).ravel()/(2*h)-j0@v)/np.linalg.norm(j0@v))
        derivatives.append(dict(step=h,relative_error=err));write(out/'DERIVATIVES.json',derivatives)
    if min(r['relative_error'] for r in derivatives)>1e-5:raise RuntimeError('continuous material Full-J mismatch')
    for method,direction in [('TASK_AWARE',v),('RANDOM_MATCHED',rand),('FULL_J_WEAK',np.linalg.svd(j,full_matrices=False)[2][-1])]:
        for distance in [.25,.5,1.0]:
            target=distance*direction
            raw,_=fw.evaluate(target);rawtask,_=ft.evaluate(target)
            # Maintain explicit signed branch progress; data corrector cannot
            # erase the edit by returning x=0. A linear progress row is not
            # interpreted as likelihood or posterior density.
            weight=100.;rows_iter=[]
            def ev(x,jac=False):
                d,jc=fw.evaluate(x,jac);r=np.r_[pack(d-d0),weight*(direction@x-distance)]
                jj=None if jc is None else np.vstack([packj(jc),weight*direction])
                return r,jj
            def cb(row):rows_iter.append(row);write(out/f'{method}_{distance}_ITERATIONS.json',rows_iter)
            fit=trust_region(ev,target,np.full(model.n,-2.),np.full(model.n,2.),max_iter=7,max_calls=18,radius=.2,max_radius=.7,l2_radius=False,prior_lambda=0.,callback=cb,objective_floor=.05,fatal_exceptions=(BudgetStop,))
            d,_=fw.evaluate(fit.x);td,_=ft.evaluate(fit.x)
            dc=float(np.linalg.norm(d-d0)/.01);dt=float(np.linalg.norm(td-t0)/.01)
            grid=np.stack(np.meshgrid(np.linspace(-2,2,401),np.linspace(0,3,301)),axis=-1).reshape(-1,2)
            field=fw.model.coefficient(fit.x,1.)(grid);basefield=fw.model.coefficient(zero,1.)(grid)
            row=dict(seed=seed,band=band,method=method,coefficient_progress=distance,actual_progress=float(direction@fit.x),d_cal=dc,d_T=dt,raw_d_cal=float(np.linalg.norm(raw-d0)/.01),raw_d_T=float(np.linalg.norm(rawtask-t0)/.01),material_RMS=float(np.sqrt(np.mean(abs(field-basefield)**2))),status='CANDIDATE_REQUIRES_INDEPENDENT_VALIDATION',cal_threshold=bool(dc<=1),task_threshold=bool(dt>=3),identity='ORACLE_CONSTRUCTED_PHYSICAL_WITNESS',stop=str(fit.message))
            np.savez_compressed(out/f'{method}_{distance}_FINAL.npz',x=fit.x,cal=d,task=td,cal_base=d0,task_base=t0,direction=direction)
            rows.append(row);write(out/'RESULTS.json',rows);print(json.dumps(row),flush=True)
    write(out/'COMPLETE.json',dict(totals=book.totals(),probabilities=None))

if __name__=='__main__':main(int(sys.argv[1]),*sys.argv[2:])
