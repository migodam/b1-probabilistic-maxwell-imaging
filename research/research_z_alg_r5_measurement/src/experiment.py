from common import *
import subprocess,signal,copy
from scipy.optimize import least_squares,minimize
from certification import certify

def solve_trf(p,A,y,method='TRF64',budget=30.,initials=None,bounds=None):
 begin=time.monotonic();cpu=time.process_time();rows=[];unique=[];n=p.n;scale=max(np.linalg.norm(y),1e-12);bounds=BOX if bounds is None else bounds
 seeds=starts(n) if initials is None else initials
 reciprocal=method=='reciprocal64'
 def check():
  if time.monotonic()-begin>=budget:raise BudgetExceeded('discovery_wall_limit')
 def fun(x):
  check();m=1/x if reciprocal else x;return pack(p.response(m,A)-y)/scale
 def jac(x):
  check();m=1/x if reciprocal else x;_,J=p.response(m,A,True);return J*(-1/x**2)[None,:]/scale if reciprocal else J/scale
 for i,m0 in enumerate(seeds):
  if time.monotonic()-begin>=budget:
   rows.append(dict(start=i,status='NOT_RUN_BUDGET'));continue
  try:
   if reciprocal:b=(1/np.asarray(bounds[1]),1/np.asarray(bounds[0]));x0=1/np.clip(m0,*bounds)
   else:b=bounds;x0=np.clip(m0,*bounds)
   res=least_squares(fun,x0,jac=jac,bounds=b,max_nfev=120,ftol=1e-12,xtol=1e-12,gtol=1e-12)
   m=1/res.x if reciprocal else res.x;r=np.linalg.norm(res.fun);row=dict(start=i,m=m,residual=r,nfev=res.nfev,njev=res.njev,status=int(res.status))
   if r<=1e-8 and not any(np.linalg.norm(m-z)/np.sqrt(n)<=1e-6 for z in unique):unique.append(m)
  except BudgetExceeded as ex:row=dict(start=i,status='TIMEOUT',error=str(ex))
  except Exception as ex:row=dict(start=i,status='ERROR',error=repr(ex))
  rows.append(row)
 candidates=[r for r in rows if 'm' in r]
 best=min(candidates,key=lambda r:r['residual'])['m'] if candidates else np.full(n,15.25)
 return dict(method=method,rows=rows,roots=unique,best=best,wall_seconds=time.monotonic()-begin,CPU_seconds=time.process_time()-cpu,starts_completed=len(candidates),starts_requested=len(seeds))

def quotient(p,y):
 ker=p.old_kernel();t=y.reshape(1,-1);begin=time.monotonic();cpu=time.process_time()
 try:
  result=enumerate_quotient(ker,t) if p.n==4 else enumerate_general_quotient(ker,t)
  if p.n==4:result['system'].s=2
  filt=filter_general_roots(result,box=BOX)
  roots=[r['m'] for r in filt['roots'] if r['status']=='CANDIDATE']
  return dict(method='Q6' if p.n==4 else 'Q20',roots=roots,raw=result,filter=filt,wall_seconds=time.monotonic()-begin,CPU_seconds=time.process_time()-cpu)
 except Exception as ex:return dict(method='Q',roots=[],error=repr(ex),wall_seconds=time.monotonic()-begin,CPU_seconds=time.process_time()-cpu)

def validate(p,m,y):
 A=np.eye(p.L)[:,:p.s];st=p.model.forward(m,independent=True);pred=st['y'].reshape(-1)[:p.s];mat=np.eye(p.n)-m[:,None]*p.H
 e=np.linalg.norm(pred-y)/max(np.linalg.norm(y),1e-300);sig=np.linalg.svd(mat,compute_uv=False)[-1]
 return dict(response_relative=e,state_relative=st['state_relative_residual'],sigma_min=sig,accepted=bool(st['state_relative_residual']<=1e-10 and e<=1e-8 and sig>=1e-6 and np.all(m>=.5-1e-8) and np.all(m<=30+1e-8)))

def acquire(case,arm,label,A,folder):
 req=folder/(label+'_request.npz');out=folder/(label+'_observed.npz');np.savez(req,a=A)
 subprocess.run([sys.executable,str(ROOT/'src/broker.py'),case,arm,label,str(req),str(out)],check=True,stdout=subprocess.DEVNULL)
 return np.load(out)['y'].reshape(-1)

def select_design(p,Q,roots,gn,menu,source_pool=None):
 arrays=np.column_stack([v['a'] for v in menu]);pred=np.array([p.dual(m)@p.B@arrays for m in roots]).reshape(len(roots),len(menu))
 gaps=[]
 for j in range(len(menu)):
  gaps.append(min((abs(pred[i,j]-pred[k,j]) for i in range(len(roots)) for k in range(i)),default=0.))
 admissible=[j for j,v in enumerate(menu) if not v['label'].startswith('old')]
 ai=max(admissible,key=lambda j:gaps[j]) if len(roots)>1 else next(j for j,v in enumerate(menu) if v['label']=='primitive_00')
 old=np.eye(p.L)[:,:p.s];_,J0=p.response(gn,old,True);scores=[]
 for a in arrays.T:
  _,J=p.response(gn,a,True);sign,logdet=np.linalg.slogdet(J0.T@J0+J.T@J+1e-14*np.eye(p.n));scores.append(logdet if sign>0 else -np.inf)
 fi=max(admissible,key=lambda j:scores[j])
 choices=dict(fixed=next(j for j,v in enumerate(menu) if v['label']=='primitive_00'),random=next(j for j,v in enumerate(menu) if v['label']=='random_00'),maxmin=ai,fisher=fi)
 return dict(choices=choices,predictions=pred,min_pair_gaps=gaps,fisher_logdet=scores,menu=menu)

def completion(p):
 inds=list(range(p.s));remaining=list(range(p.s,p.L));history=[]
 while len(inds)<p.n:
  Q,_=np.linalg.qr(p.B[:,inds]);scores=[np.linalg.norm((np.eye(p.n)-Q@Q.conj().T)@p.B[:,j]) for j in remaining]
  j=remaining.pop(int(np.argmax(scores)));inds.append(j);history.append(dict(added=j,singular_values=np.linalg.svd(p.B[:,inds],compute_uv=False)))
 return inds,history

def local_tubes(p,roots,A,y,old_y):
 output=[];eps=EPS+MODEL_ALLOWANCE
 for m in roots:
  _,J=p.response(m,np.eye(p.L)[:,:p.s],True)
  rad=float(np.clip(3*np.linalg.norm(np.linalg.pinv(J),2)*eps,.02,2.))
  lo=np.maximum(.5,m-rad);hi=np.minimum(30.,m+rad)
  fit=solve_trf(p,A,y,method='local_tube_TRF',budget=3.,initials=[m],bounds=(lo,hi))
  candidate=fit['best'];res=abs(p.response(candidate,A)-y)
  # Least squares not a feasibility solver: a failed test means SEARCH_INCONCLUSIVE.
  feasible=bool(np.all(res<=eps))
  output.append(dict(center=m,radius=rad,fit=candidate,residual_components=res,compatible_found=feasible,status='LOCAL_NUMERIC_COMPATIBLE' if feasible else 'LOCAL_SEARCH_INCONCLUSIVE',coverage='NOT_CERTIFIED',fit_record=fit))
 return output

def run_case(case):
 book=Ledger(case);g=geometry(case);pos,Q=dictionary(case,g);p=Physics(g,pos,book);obs=np.load(ROOT/'inputs'/f'{case}.npz');old=np.eye(p.L)[:,:p.s];folder=ROOT/'results'/case;folder.mkdir(exist_ok=False);caseout={}
 cert_budget=600.;cert_spent=0.
 for arm in ('clean','noisy'):
  y=obs['old_'+arm];af=folder/arm;af.mkdir();armout={};methods={}
  for method in ('Q','TRF64','reciprocal64'):
   before=Counter(book.counts);a=time.process_time()
   result=quotient(p,y) if method=='Q' else solve_trf(p,old,y,method)
   result['operation_counts']=dict(book.counts-before);result['discovery_CPU_total']=time.process_time()-a
   valid=[];certs=[]
   for m in result['roots']:
    val=validate(p,m,y)
    if not val['accepted']:certs.append(dict(m=m,validation=val,status='REJECTED_FULL_FORWARD'));continue
    valid.append(np.array(m))
    if cert_spent<cert_budget:
     c=certify(p.H,p.B[:,:p.s],p.R,y.reshape(1,-1),m,timeout=min(30.,cert_budget-cert_spent));cert_spent+=c['wall_seconds']
    else:c=dict(status='UNKNOWN',reason='TARGET_CERTIFICATION_BUDGET')
    certs.append(dict(m=m,validation=val,certificate=c))
   result['accepted_roots']=valid;result['certificates']=certs;methods[method]=result
   save(af/(method+'.json'),result)
  roots=methods['Q']['accepted_roots'];gn=methods['TRF64']['best']
  certified=[np.array(x['m']) for x in methods['Q']['certificates'] if x.get('certificate',{}).get('status')=='CERTIFIED_LOCAL_REAL_ROOT']
  design_roots=certified if len(certified)==len(roots) else roots
  with book.op('SOURCE_DESIGN',design_tasks=1):
   menu=p.menu(Q,design_roots,gn);design=select_design(p,Q,design_roots,gn,menu);inds,hist=completion(p)
   design.update(completion_indices=inds,completion_history=hist,certified_design=bool(len(certified)==len(roots) and len(roots)>0),candidate_coverage='UNKNOWN',noise_feasible_coverage='UNKNOWN')
   # Multistart attribution: exactly same active design rule, its own found set.
   for baseline in ('TRF64','reciprocal64'):
    br=methods[baseline]['accepted_roots'];bm=p.menu(Q,br,methods[baseline]['best']);bd=select_design(p,Q,br,methods[baseline]['best'],bm)
    design[baseline+'_design']=bd
  save(af/'design_freeze.json',design)
  traces={}
  for strategy,j in design['choices'].items():
   a=menu[j]['a'];yp=acquire(case,arm,strategy,a,af);A=np.column_stack([old,a]);joint_y=np.r_[y,yp]
   predictions=np.array([p.dual(m)@p.B@a for m in roots]);hard=[i for i,pr in enumerate(predictions) if abs(pr-yp[0])<=(1e-8 if arm=='clean' else EPS+MODEL_ALLOWANCE)]
   tubes=local_tubes(p,roots,A,joint_y,y) if arm=='noisy' else []
   corrected=solve_trf(p,A,joint_y,method='joint_candidate_TRF',budget=10.,initials=[np.full(p.n,15.25)]+roots)
   traces[strategy]=dict(source_label=menu[j]['label'],source_a=a,observed=yp,predictions=predictions,hard_retained=hard,tubes=tubes,corrected=corrected,new_measurements=1,energy=1.,duration=1.,candidate_status='CANDIDATE_ONLY' if hard else 'NO_FEASIBLE_CANDIDATE',global_noise_status='UNKNOWN')
   save(af/(strategy+'_trace.json'),traces[strategy])
  for baseline in ('TRF64','reciprocal64'):
   bd=design[baseline+'_design'];j=bd['choices']['maxmin'];a=bd['menu'][j]['a'];yp=acquire(case,arm,baseline+'_active',a,af);A=np.column_stack([old,a]);yy=np.r_[y,yp]
   corr=solve_trf(p,A,yy,method=baseline+'_active_joint',budget=10.,initials=[methods[baseline]['best']]+methods[baseline]['accepted_roots'])
   traces[baseline+'_active']=dict(source_label=bd['menu'][j]['label'],source_a=a,observed=yp,corrected=corr,new_measurements=1,energy=1.,duration=1.)
   save(af/(baseline+'_active_trace.json'),traces[baseline+'_active'])
  # Strong direct completion is mandatory even if root coverage fails.
  Ac=np.eye(p.L)[:,inds[p.s:]];yc=acquire(case,arm,'direct_completion',Ac,af);yy=np.r_[y,yc];Bfull=p.B[:,inds]
  with book.op('DIRECT_COMPLETION',direct_linear_solve=1):
   w=np.linalg.solve(Bfull.T,yy);v=p.R.reshape(-1)+w@p.H;raw=w/v;projected=np.clip(raw.real,.5,30.)
  corr=solve_trf(p,np.eye(p.L)[:,inds],yy,method='direct_joint_TRF64',budget=30.)
  direct=dict(indices=inds,raw_material=raw,projected_material=projected,minimum_node=abs(v).min(),condition=np.linalg.cond(Bfull),rank=np.linalg.matrix_rank(Bfull,tol=1e-12),corrected=corr,new_measurements=p.n-p.s,energy=p.n-p.s,duration=p.n-p.s,observations=yy,coverage='UNKNOWN')
  save(af/'direct.json',direct)
  # Ordinary direct joint TRF on identical chosen active observation, full64 starts.
  active=traces['maxmin'];A=np.column_stack([old,active['source_a']]);yy=np.r_[y,active['observed']]
  joint=solve_trf(p,A,yy,method='same_active_joint_TRF64',budget=30.);save(af/'same_active_joint_TRF64.json',joint)
  cert_count=len(certified)
  # Disjointness of actual reported certificate boxes checked in summary later.
  armout=dict(root_counts={k:len(v['accepted_roots']) for k,v in methods.items()},local_certified_count=cert_count,global_root_coverage='UNKNOWN',noise_feasible_coverage='UNKNOWN',menu_count=len(menu),chosen_sources={k:menu[v]['label'] for k,v in design['choices'].items()},maxmin_gap=design['min_pair_gaps'][design['choices']['maxmin']],direct_condition=direct['condition'])
  caseout[arm]=armout
 save(folder/'summary.json',caseout);save(folder/'receipt.json',book.receipt());print(case,json.dumps(serial(caseout)),flush=True)

if __name__=='__main__':run_case(sys.argv[1])
