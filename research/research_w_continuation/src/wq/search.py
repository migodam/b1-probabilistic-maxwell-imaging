"""Fixed-budget full-domain nonlinear wall search, conditional on Q PASS.

Online code below never imports wall_values/continuous truth. Reference loader
reads only four complex measurement arrays. Truth scoring is a separate job.
"""
from pathlib import Path
import os,json,time,gc,sys
import numpy as np
from scipy.linalg import cho_factor,cho_solve,eigh
from scipy.ndimage import gaussian_filter
from .accounting import SearchBook,BudgetStop,encode
from .fem import FEM,PMLConfig,line_green
from .inverse_support import PixelGrid,pixel_mesh,pixel_material,project_parameters,occupied_area,volume_pixel_jacobian
ROOT=Path(__file__).resolve().parents[2]
SCALE=None

def real(a):return np.r_[np.asarray(a).real.ravel(),np.asarray(a).imag.ravel()]/.01
def realj(a):return np.vstack([a.real,a.imag])/.01
def legal(z,g):return project_parameters(z*np.r_[np.full(g.n_pixels,8.),np.full(g.n_pixels,.15)],g)
def normed(theta,g):return theta/np.r_[np.full(g.n_pixels,8.),np.full(g.n_pixels,.15)]
def ridge_step(j,r,z,lam=.02):
    c=cho_factor(j@j.T+lam*np.eye(len(r)),check_finite=False)
    new=j.T@cho_solve(c,j@z-r,check_finite=False)
    step=new-z
    return step*min(1.,.2/max(np.max(abs(step)),1e-30))
def task_direction(j,k):
    lam=.01;c=cho_factor(j@j.T+lam*np.eye(j.shape[0]),check_finite=False)
    def rinv(x):return (x-j.T@cho_solve(c,j@x,check_finite=False))/lam
    action=rinv(k.T);a=k@action;a=(a+a.T)/2
    vals,vecs=eigh(a,subset_by_index=[len(a)-1,len(a)-1],check_finite=False)
    v=action@vecs[:,0];return v/max(np.max(abs(v)),1e-30)*.2,float(vals[0])
def random_smooth(rng,g,amplitude=.2):
    a=np.r_[gaussian_filter(rng.normal(size=(g.ny,g.nx)),1.5).ravel(),gaussian_filter(rng.normal(size=(g.ny,g.nx)),1.5).ravel()]
    return amplitude*a/max(np.max(abs(a)),1e-30)

class Forward:
    def __init__(self,frequencies,book,phase='search'):
        self.freq=frequencies;self.book=book;self.phase=phase;self.grid=PixelGrid();self.mesh=pixel_mesh(self.grid)
        self.tx=np.c_[np.linspace(-1.4,1.4,4),np.full(4,-.3)];self.rx=np.c_[np.linspace(-1.5,1.5,12),np.full(12,-.3)]
        self.scale=np.r_[np.full(self.grid.n_pixels,8.),np.full(self.grid.n_pixels,.15)]
        self.backends={};self.theta=None;self.data=None;self.states=[];self.jac=None
    def clear(self):self.backends.clear();self.states=[];self.theta=None;self.data=None;self.jac=None;gc.collect()
    def evaluate(self,theta,task=None,jac=False):
        if task is None and self.theta is not None and np.array_equal(theta,self.theta):
            if not jac or self.jac is not None:return self.data,self.jac
        # Stream frequencies instead of retaining five sparse factors/bases.
        # All calls needing J request it with F in the same frequency pass.
        self.clear();cached=False
        outputs=[];js=[];states=[]
        for i,f in enumerate(self.freq):
            p=self.backends.get(f)
            if p is None:
                p=FEM(f,self.tx,self.rx,self.book,h=.18,order=3,mesh=self.mesh,phase=self.phase,receiver_mode='volume');self.backends[f]=p
            if cached:s=self.states[i]
            else:
                base=pixel_material(theta,f,self.grid)
                if task is None:coef=base
                else:
                    center,radius,contrast=task
                    def coef(pts,base=base,center=center,radius=radius,contrast=contrast,f=f):
                        a=base(pts);m=np.linalg.norm(pts-np.asarray(center),axis=1)<=radius
                        a[m]=1+contrast+1j*.02/f;return a
                s=p.solve(coef)
            outputs.append(s.data.ravel());states.append(s)
            if jac:js.append(volume_pixel_jacobian(p,s,self.grid)*self.scale[None,:])
            self.backends.pop(f,None);del p,s;states=[];gc.collect()
        d=np.concatenate(outputs);j=np.concatenate(js,axis=0) if jac else None
        if task is None:self.theta=theta.copy();self.data=d.copy();self.states=[];self.jac=j
        else:self.clear()
        return d,j
    def bp(self,data):
        g=self.grid;x=(g.x_edges[:-1]+g.x_edges[1:])/2;y=(g.y_edges[:-1]+g.y_edges[1:])/2
        xx,yy=np.meshgrid(x,y);pts=np.c_[xx.ravel(),yy.ravel()];J=[];direct=[]
        for f in self.freq:
            k=2*np.pi*f;E=line_green(pts,self.tx,k);S=k*k*g.cell_area*line_green(self.rx,pts,k)
            a=np.einsum('rn,ns->srn',S,E).reshape(48,g.n_pixels);a[:,g.calibration_mask(pts)]=0
            J.append(np.c_[a,1j*a/f]*self.scale[None,:]);direct.append(line_green(self.rx,self.tx,k).T.ravel())
        j=realj(np.concatenate(J));b=real(data-np.concatenate(direct));c=cho_factor(j@j.T+.02*np.eye(len(b)),check_finite=False)
        z=j.T@cho_solve(c,b,check_finite=False);z*=min(1.,.35/max(np.max(abs(z)),1e-30))
        return legal(z,g)

def reference(seed,band,label):
    with np.load(ROOT/'results/wall-pml/ARRAYS.npz') as z:return z[f'seed{seed}_{band}_{label}_volume'].ravel().copy()
def counts(book):return book.totals()['rhs']
def objective(d,y,z):return .5*np.linalg.norm(real(d-y))**2+.01*np.dot(z,z)

def main(mode):
    q=json.loads((ROOT/'P0_Q_QUALIFICATION.json').read_text())
    if q['decision']!='PASS':raise RuntimeError('Q failed; search not authorized')
    out=ROOT/'results/p0b'/mode;out.mkdir(parents=True,exist_ok=False)
    book=SearchBook(ROOT/'results/p0b');g=PixelGrid();scale=np.r_[np.full(g.n_pixels,8.),np.full(g.n_pixels,.15)]
    rows=[];arr={}
    def row(**kw):
        rows.append(kw);(out/'RESULTS.json').write_text(json.dumps(rows,indent=2,default=encode,allow_nan=False));np.savez_compressed(out/'ARRAYS.npz',**arr);print(json.dumps(kw,default=encode),flush=True)
    if mode=='preflight':
        fw=Forward([1.],book,phase='repair');rng=np.random.default_rng(2026102301)
        raw=np.r_[np.full(g.n_pixels,.3),np.full(g.n_pixels,.01)];theta=project_parameters(raw,g)
        d,j=fw.evaluate(theta,jac=True);v=theta*rng.normal(size=len(theta))*.1
        # Stored Jacobian includes normalized-coordinate scales.
        pred=j@(v/scale);step=1e-4
        a,_=fw.evaluate(theta+step*v);b,_=fw.evaluate(theta-step*v)
        err=float(np.linalg.norm((a-b)/(2*step)-pred)/np.linalg.norm(pred))
        row(case='pixel_J_FD',relative_error=err,rhs=book.totals()['rhs'],occupied_area=occupied_area(theta,g),parameter_count=len(theta),mesh_vertices=fw.mesh.nvertices)
        arr['theta']=theta;arr['direction']=v;arr['Jdirection']=pred;arr['FD']=(a-b)/(2*step)
        np.savez_compressed(out/'ARRAYS.npz',**arr)
        if err>1e-6:raise RuntimeError('Pixel derivative inconsistent; stop search')
    elif mode in ('search','search-corrected'):
        cfg=json.loads((ROOT/'docs/P0_FROZEN_CONFIG.json').read_text());probe=cfg['probes'][8]
        for seed in cfg['truth_seeds']:
            for band,frequencies in cfg['frequencies'].items():
                search_used=sum(r.get('rhs',0) for r in [json.loads(x) for x in book.path.read_text().splitlines()] if r.get('event')=='reserve' and r.get('phase')=='search')
                if 3200-search_used<800:
                    row(case='context_not_run',seed=seed,band=band,reason='Insufficient remaining search quota for all four full method allowances after retained repair costs',remaining_rhs=3200-search_used);continue
                context_start=counts(book);fw=Forward(frequencies,book);y=reference(seed,band,'cal');common=fw.bp(y);z=normed(common,g)
                data,jc=fw.evaluate(common,jac=True);j=realj(jc);r=real(data-y);initial_obj=objective(data,y,z)
                shared_data=data.copy();shared_j=j.copy();arr[f'{seed}_{band}_COMMON']=common
                taskdata,kc=fw.evaluate(common,task=probe,jac=True);k=realj(kc)
                v,eigen=task_direction(j,k);s=ridge_step(j,r,z)
                def predicted(a):return np.linalg.norm(r+j@a)**2+.02*np.linalg.norm(z+a)**2
                if predicted(-v)<predicted(v):v=-v
                rng=np.random.default_rng(seed+(0 if band=='wide' else 100))
                proposals={'FULL_J':s,'RANDOM_MATCHED':s+random_smooth(rng,g),'MULTISTART_FWI':random_smooth(rng,g,.35),'TASK_AWARE':s+v}
                row(case='context_initial',seed=seed,band=band,d_cal=float(np.linalg.norm(r)),shared_rhs=counts(book)-context_start,task_ratio=eigen,representation='Full-domain48x32 passive pixels; no wall geometry input')
                for method,proposal in proposals.items():
                    started=time.perf_counter();mstart=counts(book);proposal*=min(1.,.2/max(np.max(abs(proposal)),1e-30))
                    theta=legal(z+proposal,g);candidate_d,cj=fw.evaluate(theta,jac=True)
                    best_theta=common.copy();best_data=shared_data.copy();best_obj=initial_obj;accepted=0
                    first_obj=objective(candidate_d,y,normed(theta,g))
                    if first_obj<best_obj:best_theta=theta.copy();best_data=candidate_d.copy();best_obj=first_obj;accepted+=1
                    # Same one full-J corrector and four standard trial scales.
                    # Rejected proposals never become the accepted corrector anchor.
                    # This repair changes no physical gate/threshold or trial budget.
                    if mode=='search-corrected' and not accepted:
                        theta=common.copy();candidate_d=shared_data.copy();cj_complex=None;cz=z.copy();cr=r.copy();corrector_j=shared_j
                    else:
                        candidate_d,cj=fw.evaluate(theta,jac=True);cz=normed(theta,g);cr=real(candidate_d-y);corrector_j=realj(cj)
                    correction=ridge_step(corrector_j,cr,cz)
                    for alpha in (1.,.5,.25,.125):
                        if counts(book)-mstart+4*len(frequencies)>160:break
                        trial=legal(cz+alpha*correction,g);td,_=fw.evaluate(trial);obj=objective(td,y,normed(trial,g))
                        row(case='corrector_trial',seed=seed,band=band,method=method,alpha=alpha,d_cal=float(np.linalg.norm(real(td-y))),objective=float(obj),accepted=bool(obj<best_obj))
                        if obj<best_obj:
                            best_theta=trial.copy();best_data=td.copy();best_obj=obj;accepted+=1;break
                    key=f'{seed}_{band}_{method}';arr[key]=best_theta;arr[key+'_data']=best_data.reshape(len(frequencies),4,12)
                    row(case='calibration_candidate',seed=seed,band=band,method=method,d_cal=float(np.linalg.norm(real(best_data-y))),objective=float(best_obj),accepted_updates=accepted,method_rhs=counts(book)-mstart,method_wall_seconds=time.perf_counter()-started,occupied_area=occupied_area(best_theta,g),material_step_from_BP=float(np.linalg.norm((best_theta-common)*np.sqrt(g.cell_area))),source='ONLINE_MEASUREMENT_AND_PRIOR_ONLY',independent_endpoint_validation=False)
                    fw.clear()
                fw.clear();del fw;gc.collect()
    elif mode=='task-endpoints':
        cfg=json.loads((ROOT/'docs/P0_FROZEN_CONFIG.json').read_text());probe=cfg['probes'][8]
        searchdir=ROOT/'results/p0b/search-corrected'
        candidates=json.loads((searchdir/'RESULTS.json').read_text())
        with np.load(searchdir/'ARRAYS.npz') as z:saved={k:z[k].copy() for k in z.files if not k.endswith('_data')}
        for r in candidates:
            if r['case']!='calibration_candidate':continue
            seed,band,method=r['seed'],r['band'],r['method'];key=f'{seed}_{band}_{method}'
            fw=Forward(cfg['frequencies'][band],book,phase='endpoint');d,_=fw.evaluate(saved[key],task=probe);ref=reference(seed,band,'task')
            arr[key+'_task_data']=d.reshape(len(fw.freq),4,12)
            row(case='fixed_task_endpoint',seed=seed,band=band,method=method,d_cal=r['d_cal'],task8_response_distance=float(np.linalg.norm(real(d-ref))),task_class_coverage='1/24; lower bound only',joint_centroid_difference=0.,status='CANDIDATE_ONLY_NOT_GO',rhs=4*len(fw.freq))
            fw.clear();del fw;gc.collect()
    else:raise ValueError(mode)
    (out/'TOTALS.json').write_text(json.dumps(book.totals(),indent=2))
if __name__=='__main__':main(sys.argv[1])
