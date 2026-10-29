from core import *
import argparse,traceback,resource
from scipy.special import ndtr
from scipy.linalg import solve

CFG=json.loads((P/'INVERSE_CONFIG.json').read_text())

def edges(idx):
 lookup={int(k):i for i,k in enumerate(idx)};ed=[]
 for k,i in lookup.items():
  for kk in (k+1 if k%24<23 else -1,k+24):
   if kk in lookup:ed.append((i,lookup[kk]))
 D=np.zeros((len(ed),len(idx)))
 for r,(i,j) in enumerate(ed):D[r,i]=1;D[r,j]=-1
 return D

class Problem:
 def __init__(self,ds,s):self.ds=ds;self.pr=prior();self.idx=s['ROI'];self.D=edges(self.idx)
 def decode(self,x):
  u=ndtr(x[-2:]);delta=complex(.8+u[0],.03+.15*u[1])
  if self.ds=='N':return x[:256],delta
  z=x[:64];ss=(self.pr.config['field_mean']+self.pr.A[self.idx]@z)/self.pr.temperature
  return saturation(ss),delta
 def jac(self,x,a,j):
  uder=np.exp(-x[-2:]**2/2)/np.sqrt(2*np.pi)
  if self.ds=='N':J=j['J_occ']
  else:
   ss=(self.pr.config['field_mean']+self.pr.A[self.idx]@x[:64])/self.pr.temperature
   Da=saturation_prime(ss)[:,None]*self.pr.A[self.idx]/self.pr.temperature;J=j['J_occ']@Da
  return white_jac(np.column_stack((J,j['J_delta']*uder[0],1j*.15*j['J_delta']*uder[1])))
 def reg(self,x,hess=False):
  n=len(x)
  if self.ds=='P':return (.5*np.dot(x,x),x.copy(),np.eye(n))
  c=CFG['N']['regularization'];a=x[:256];v=self.D@a;w=np.sqrt(v*v+c['TV_epsilon']**2)
  val=c['TV_Huber_lambda']*np.sum(w-c['TV_epsilon'])+.5*c['occupancy_ridge']*np.dot(a,a)+c['occupancy_L1']*a.sum()+.5*np.dot(x[-2:],x[-2:])
  g=np.r_[c['TV_Huber_lambda']*self.D.T@(v/w)+c['occupancy_ridge']*a+c['occupancy_L1'],x[-2:]]
  H=np.zeros((n,n));H[:256,:256]=c['TV_Huber_lambda']*(self.D.T/w)@self.D+c['occupancy_ridge']*np.eye(256);H[-2:,-2:]=np.eye(2)
  return val,g,H
 def project(self,x):
  x=x.copy()
  if self.ds=='N':x[:256]=np.clip(x[:256],0,1)
  return x

def make_solver(method,s,led,a0):
 if method=='FULL256':return ExactSolver(s['C'],s['O'],s['E'],s['Fh'],ledger=led)
 if method=='AGENT_A':
  from research.p3l_r11_stable_structure_digits.agent_a.solver import Solver
  return Solver(s['C'],s['O'],s['E'],s['Fh'],SIGMA,ledger=led)
 if method=='AGENT_B':
  from research.p3l_r11_stable_structure_digits.agent_b.fixed_total_field import TotalFieldSolver
  return TotalFieldSolver(s['C'],s['O'],s['E'],s['Fh'],a_anchor=a0,ranks=(32,64,96),basis_kind='anchor',ledger=led)
 if method=='AGENT_C':
  from research.p3l_r11_stable_structure_digits.agent_c.controlled_solver import GlobalControlledSolver
  return GlobalControlledSolver(s['C'],s['O'],s['E'],s['Fh'],SIGMA,ledger=led,a_anchor=a0,ranks=(32,64,96,128),basis_kind='green',rhs_cap=330000)
 raise ValueError(method)

def run_one(ds,oid,method):
 out=P/'results'/ds/method/oid
 if (out/'PREDICTION_SEALED.json').exists():print('EXIST',oid,method,flush=True);return
 out.mkdir(parents=True,exist_ok=True);led=Ledger(out/'COST_LEDGER.jsonl',dict(parent=oid,method=method),cap=12000)
 start=time.perf_counter();cpu=time.process_time();y=np.load(P/'data/online'/f'{oid}.npz')['y'];s=load_setup();problem=Problem(ds,s);cfg=CFG[ds];sc=CFG['shared'];n=cfg['dimension'];x=np.zeros(n)
 if ds=='N':x[:256]=cfg['initial_occupancy']
 a0,_=problem.decode(x);ts=time.perf_counter();solver=make_solver(method,s,led,a0);setup=time.perf_counter()-ts
 history=[];nf=nj=0;trust=sc['trust_radius'];damp=sc['initial_damping'];pat=0;status='MAX_J';err=None;res=None
 try:
  for it in range(cfg['max_J']):
   if time.perf_counter()-start>cfg['max_wall_seconds']:status='WALL_CAP';break
   a,d=problem.decode(x);res=solver.evaluate(a,d,need_jac=True,y=y);nf+=1;nj+=1
   r=white(res['F']-y);J=problem.jac(x,a,res);rv,rg,rH=problem.reg(x);obj=.5*r@r+rv;g=J.T@r+rg;H=J.T@J+rH
   scale=np.maximum(np.diag(H),1.);step=solve(H+damp*np.diag(scale),-g,assume_a='pos',check_finite=False)
   sn=np.linalg.norm(step)
   if sn>trust:step*=trust/sn
   accepted=False;newobj=obj;dx=0.
   for alpha in sc['line_search']:
    if nf>=cfg['max_forward']:break
    q=problem.project(x+alpha*step);dx=np.linalg.norm(q-x)
    if dx<sc['projected_step_stop']:break
    aq,dq=problem.decode(q);trial=solver.evaluate(aq,dq,need_jac=False,y=y);nf+=1;rr=white(trial['F']-y);ro=problem.reg(q)[0];newobj=.5*rr@rr+ro
    if newobj<obj:
     old=x;x=q;accepted=True;res=trial;break
   history.append(dict(iteration=it,objective=float(obj),new_objective=float(newobj),Phi=float(.5*r@r),accepted=accepted,step_norm=float(dx),trust=trust,damping=damp,support_size=int(np.count_nonzero(a>0)),forward_calls=nf,J_calls=nj,physics=res.get('diagnostics',{})))
   if accepted:
    rel=(obj-newobj)/max(abs(obj),1.);pat=pat+1 if rel<sc['relative_objective_stop'] else 0;trust=min(sc['max_trust_radius'],trust*1.2);damp=max(damp*.5,1e-8)
   else:
    damp=min(1e12,damp*10);trust=max(.001,trust*.5);pat+=1
   if pat>=sc['relative_objective_patience']:status='OBJECTIVE_STOP';break
   if nf>=cfg['max_forward']:status='FORWARD_CAP';break
  a,d=problem.decode(x)
  # Independent final full-Maxwell check, charged to samejob; noGT.
  reference=ExactSolver(s['C'],s['O'],s['E'],s['Fh'],ledger=led).evaluate(a,d,False,y)
  chi=np.full(576,HOST);chi[s['ROI']]+=d*a;af=np.zeros(576);af[s['ROI']]=a
  np.savez_compressed(out/'PREDICTION.npz',x=x,occupancy=af,material=chi,F=reference['F'],delta=np.array(d),approx_F=res['F'],initial_occupancy=a0,mask=af>=.5)
  write(out/'HISTORY.json',history)
  end=dict(status=status,parent=oid,dataset=ds,method=method,GT_used=False,template_used=False,digit_label_used=False,model='native258new' if ds=='N' else 'originalGaussianized66',source_blake2b=blake(Path(__file__)),inverseconfig_blake2b=blake(P/'INVERSE_CONFIG.json'),prediction_blake2b=blake(out/'PREDICTION.npz'),Phi=reference['phi'],noise_ball_pass=bool(reference['phi']<=sc['observed_noise_ball']),actual_final_white_approx_error=float(np.linalg.norm(white(reference['F']-res['F']))),forward_calls=nf+1,J_calls=nj,physics_RHS=led.rhs,physics_LU=led.lu,physics_kernel_wall=led.wall,setup_seconds=setup,wall_seconds=time.perf_counter()-start,CPU_seconds=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=err)
  write(out/'PREDICTION_SEALED.json',end);print('SEALED',oid,method,status,'Phi',round(end['Phi'],3),'J',nj,'RHS',led.rhs,'wall',round(end['wall_seconds'],3),flush=True)
 except BaseException as e:
  write(out/'FAILURE.json',dict(status='FAILED',error=repr(e),traceback=traceback.format_exc(),newRHS=led.rhs,LU=led.lu,history=history));print('FAIL',oid,method,repr(e),flush=True);raise
 finally:led.close()

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--dataset',choices=['N','P'],required=True);ap.add_argument('--method',required=True);ap.add_argument('--objects',default='0,1,2,3,4,5,6,7,8');ns=ap.parse_args()
 for i in map(int,ns.objects.split(',')):run_one(ns.dataset,f'{ns.dataset}_{i:03d}',ns.method)
