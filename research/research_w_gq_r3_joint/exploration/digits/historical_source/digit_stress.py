"""Independent nonparametric material-grid stress; no digit/template imports online."""
from pathlib import Path
import os,sys,json,time,fcntl,argparse,gc
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import numpy as np
from scipy.optimize import minimize
from skfem import asm
ROOT=Path(__file__).resolve().parents[1];EX=ROOT/'exploration/digits';R2=ROOT.parent/'research_w_gq_r2'
for p in reversed([ROOT/'frozen/src',R2/'src',R2/'reuse/previous_gq/src',R2/'reuse/fixed_vendor',ROOT.parents[2]/'research-w-continuation/.deps']):sys.path.insert(0,str(p))
from r3.accounting import BudgetBook,BudgetExceeded,fem_lock,plain
from r3.physics import PublicWall
from r2.physics import Eval,S1_TX,line_green
from wslam.fem_backend import _contrast_mass,_forcing
WX=np.array([-1.5,-.5,.5,1.5]);WY=np.array([.5,.875,1.25]);TX=np.linspace(-.7,.7,8);TY=np.linspace(1.4,2.8,8);FREQ=[.85,1.,1.15];SIGMA=.01
LOW=np.r_[np.ones(6),np.zeros(6),np.zeros(98)];HIGH=np.r_[np.full(6,6.),np.full(6,.25),np.full(49,8.),np.full(49,.25)];SCALE=np.r_[np.ones(6),np.full(6,.1),np.full(49,2.),np.full(49,.1)]
def dump(p,a):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(a,indent=2,default=plain,allow_nan=False)+'\n')
class Book(BudgetBook):
 def check(self,new_job=False):
  super().check(new_job=False)
  if time.time()>1793298600:raise BudgetExceeded('digit exploration hard stop 18:30 UTC')
  marker=EX/'FIRST_HEAVY_STARTED.json'
  if new_job and not marker.exists():
   if time.time()>self.heavy_after:raise BudgetExceeded('first digit job missed original new-heavy cutoff')
   dump(marker,dict(time_unix=time.time(),pid=os.getpid(),scope='one authorized multi-process digit8 experiment'))
 def reserve(self,phase,label,rhs=0,factorizations=0,**meta):
  self.check();p=EX/'PHYSICS_LEDGER.jsonl'
  with (EX/'.budget.lock').open('a') as lock:
   fcntl.flock(lock,fcntl.LOCK_EX);rows=[json.loads(l) for l in p.read_text().splitlines()] if p.exists() else [];nr=sum(r['rhs'] for r in rows);nf=sum(r['factorizations'] for r in rows)
   if nr+rhs>3600 or nf+factorizations>450:raise BudgetExceeded('digit cap 450FE/3600RHS')
   super().reserve(phase,label,rhs=rhs,factorizations=factorizations,**meta)
   with p.open('a') as f:f.write(json.dumps(dict(time_unix=time.time(),pid=os.getpid(),rhs=rhs,factorizations=factorizations,label=label,phase=phase))+'\n')
class GridEval:
 def __init__(self):
  self.book=Book(root=ROOT);self.book.check(new_job=True);self.e=[]
  with fem_lock(ROOT):
   for f in FREQ:
    ev=Eval(self.book,'DIGIT_EXPLORE',f,S1_TX,WX,WY,h=.17,targetstep=.1,max_factor_bytes=1500000000);p=ev.p;pts=p.material_points
    inside=(pts[:,0]>=TX[0])&(pts[:,0]<=TX[-1])&(pts[:,1]>=TY[0])&(pts[:,1]<=TY[-1]);ix=np.clip(np.searchsorted(TX,pts[:,0],side='right')-1,0,6);iy=np.clip(np.searchsorted(TY,pts[:,1],side='right')-1,0,6)
    self.e.append((ev,inside,iy*7+ix))
 def run(self,z,y=None,gradient=False):
  z=np.asarray(z);wall=PublicWall(WX,WY,z[:6].reshape(2,3),z[6:12].reshape(2,3));data=[];g=np.zeros(110);diag=[]
  with fem_lock(ROOT):
   for k,(ev,inside,cells) in enumerate(self.e):
    self.book.check();p=ev.p;ct=np.zeros(len(cells),complex);ct[inside]=z[12:61][cells[inside]]+1j*z[61:110][cells[inside]]/p.freq
    A,b=ev.background(wall);M=asm(_contrast_mass,p.material_basis,contrast=ct.reshape(p._material_shape));A=(A-p.k**2*M[p._free_dofs][:,p._free_dofs]).tocsc();inc=line_green(p.material_points,p.tx,p.k)
    for j in range(4):b[:,j]+=asm(_forcing,p.material_basis,load=(p.k**2*ct*inc[:,j]).reshape(p._material_shape))[p._free_dofs]
    fac,A,storage=ev.factor(A,label='digit material LU');u=fac.solve(b,label='digit primal4');d=ev.direct+np.asarray(ev.R@u).T;data.append(d.copy());diag.append(dict(freq=p.freq,ndofs=p.ndofs,factor_bytes=storage,backward=fac.last_defect))
    if gradient:
     adj=fac.solve(np.asarray(ev.R.conj().T@((d-y[k])/SIGMA**2).T),trans='H',label='digit weighted adjoint4');gw=np.zeros(6,complex);gt=np.zeros(49,complex)
     for sl,(scattered,av) in ev._material_chunks(u,adj):
      total=scattered+inc[sl];prod=p.k**2*p.material_weights[sl]*np.einsum('qs,qs->q',np.conj(av),total);wm=ev._inside_wall[sl];tm=inside[sl];np.add.at(gw,ev._cells[sl][wm],prod[wm]);np.add.at(gt,cells[sl][tm],prod[tm])
     g[:6]+=gw.real;g[6:12]+=(1j/p.freq*gw).real;g[12:61]+=gt.real;g[61:110]+=(1j/p.freq*gt).real
  data=np.asarray(data)
  if y is None:return data,diag
  pe=(z[:6]-3.)/1.5;pl=(z[6:12]-.06)/.15;value=.5*float(np.sum(abs(data-y)**2))/SIGMA**2+.5*float(pe@pe+pl@pl)
  if gradient:g[:6]+=pe/1.5;g[6:12]+=pl/.15
  return value,g,data,diag
 def maps(self,z):return dict(wall_eps=z[:6].reshape(2,3),wall_loss=z[6:12].reshape(2,3),target_contrast=z[12:61].reshape(7,7),target_loss=z[61:].reshape(7,7),target_x_edges=TX,target_y_edges=TY)
class CallCap(Exception):pass
def initial(start=0):return np.r_[np.full(6,3. if start==0 else 3.8),np.full(6,.06 if start==0 else .08),np.full(49,.1),np.full(49,.03)]
def online(measurement=None,output=None,config=None):
 if config is not None:
  public=json.loads(Path(config).read_text())
  if not np.allclose(public['target_x_edges'],TX) or not np.allclose(public['target_y_edges'],TY):raise ValueError('public grid mismatch')
 out=Path(output) if output is not None else EX/'results/digit8';out.mkdir(parents=True,exist_ok=True)
 if (out/'ONLINE_SEALED.json').exists():raise RuntimeError('digit8 already sealed')
 y=np.load(Path(measurement) if measurement is not None else EX/'online/observation.npz')['y'];ev=GridEval();records=[];totalcalls=0;best=[float('inf'),None,None];audit=initial();rng=np.random.default_rng(33105008);v=rng.normal(size=110)*SCALE;v/=np.linalg.norm(v/SCALE);h=1e-5
 base,grad,_,_=ev.run(audit,y,True);plus=ev.run(audit+h*v,y,False)[0];minus=ev.run(audit-h*v,y,False)[0];fd=(plus-minus)/(2*h);analytic=float(grad@v);dump(EX/'WEIGHTED_GRADIENT_AUDIT.json',dict(relative_error=abs(fd-analytic)/max(abs(fd),1e-12),finite_difference=fd,adjoint_directional=analytic,scope='one legal public mixed wall/target direction; not reconstruction evidence'));totalcalls=3
 if abs(fd-analytic)>1e-3*max(abs(fd),1e-10):raise RuntimeError('mixed weighted gradient failed')
 for start in range(2):
  calls=0;last=[None,None]
  def fun(q):
   nonlocal calls,totalcalls
   z=q*SCALE
   if last[0] is not None and np.array_equal(q,last[0]):return last[1]
   if calls>=65 or totalcalls>=140:raise CallCap('frozen actual call cap')
   value,g,pred,diag=ev.run(z,y,True);calls+=1;totalcalls+=1
   np.savez_compressed(out/f's{start}_call{calls:03d}.npz',z=z,data=pred)
   with (out/'TRACE.jsonl').open('a') as f:f.write(json.dumps(dict(time_unix=time.time(),start=start,call=calls,totalcalls=totalcalls,objective=value,gradient_norm=float(np.linalg.norm(g*SCALE))))+'\n')
   if value<best[0]:best[:]=[value,z.copy(),pred.copy()]
   last[:]=[q.copy(),(value,g*SCALE)];return last[1]
  success=False;nit=None
  try:r=minimize(fun,initial(start)/SCALE,jac=True,method='L-BFGS-B',bounds=list(zip(LOW/SCALE,HIGH/SCALE)),options=dict(maxiter=40,maxls=8,ftol=1e-9,gtol=1e-4,maxfun=65));stop=str(r.message);success=bool(r.success);nit=int(r.nit)
  except (CallCap,BudgetExceeded) as e:stop=str(e)
  records.append(dict(start=start,calls=calls,iterations=nit,optimizer_success=success,stop=stop));dump(out/'STARTS.json',records)
  if best[1] is not None:np.savez_compressed(out/'endpoint.npz',z=best[1],data=best[2],**ev.maps(best[1]))
  if time.time()>1793298570:break
 if best[1] is None:raise RuntimeError('No actual grid reconstruction')
 dump(out/'ONLINE_SEALED.json',dict(time_unix=time.time(),objective=best[0],selection='minimum full measured objective plus frozen wall prior',truth_access=False,ellipse_used=False,template_used=False,target_free_cells=49,target_real_parameters=98,starts=records,total_objective_calls=totalcalls,scope='NEW_EXPLORATORY_NONPARAMETRIC_GRID same-discrete truth/inverse; not a primary R3 blind PASS'))
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--measurement',type=Path);parser.add_argument('--config',type=Path);parser.add_argument('--output',type=Path);args=parser.parse_args();online(args.measurement,args.output,args.config)
