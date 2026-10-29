"""Separate frozen fresh10 mechanism check: generator, observation-only online, scorer."""
from __future__ import annotations
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[n]='1'
import argparse,copy,csv,fcntl,hashlib,json,sys,time,traceback
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent));import candidate_generation as cg
pv=cg.pv;P=cg.P;H=cg.H;OUT=P/'results/candidates_fresh10';CONFIG=P/'config/candidate_fresh10_protocol.json';ENABLE=P/'config/candidate_fresh10_enable.json';PRIVATE=P/'private/candidate_fresh10';cg.OUT=OUT
NEWARMS=cg.ARMS;OLDARMS=pv.ARMS;ALLARMS=NEWARMS+OLDARMS

def freeze():
 if CONFIG.exists():raise RuntimeError('no protocol refit')
 dev=json.loads(cg.CONFIG.read_text());old=json.loads((H/'config/possibility_volume_protocol.json').read_text());ns=[20261030,101009];entries=[]
 for i in range(10):entries.append(dict(parent=f'R10BI_FRESH_{i:03d}',index=i,prior_seed=ns+[101,i],noise_seed=ns+[102,i],target_seed=ns+[103,i],weak_seed_prefix=ns+[104,i],random_seed_prefix=ns+[105,i]))
 c=dict(id='R10_BILATERAL_FRESH10_MECHANISM_V1',approved_by_root=False,created_utc=time.time(),output=str(OUT),private_output=str(PRIVATE),parents=entries,arms=ALLARMS,K=8,seed_namespace=ns,exposure='Fresh10 bounded mechanism check, small N10; not paper qualification; all6arms fixed before any freshdata or DEVGT review by candidate worker',candidate_new=dev['candidate'],candidate=old['candidate'],anchor=old['anchor'],primary=old['primary'],model_prior_weight=1.,weak_prior_weight=1.,weak_rank=8,latent_radius=4.,radii=old['radii'],shared_setup=dev['shared_setup'],noise_ball=dev['noise_ball'],caps=dict(reduced_RHS256=155000,reduced_LU256=5000,full_RHS576=0,full_LU576=0,controlled_wall_seconds=1200),source_blake2b=hashlib.blake2b(Path(__file__).read_bytes(),digest_size=32).hexdigest(),dev_source_blake2b=cg.source_id(),frozen_dev_protocol_id=dev['id'],generation='original P1 x66Gaussian, probit amplitude q64:66; roles101/102 disjoint seeds; observation only online inputs',old_proposal_prefix='exact first8 original order pair0/radii .02,.05,.10,.15 each +/-; old6F4J and latentcap4',new_proposal_prefix='first8 independent paired real512 eta and Gaussian66 xi draws; exact frozen12F8J newarms; no latentcap',fairness='sameK8; optimization caps, prior-target distributions, radiusconstraints and correctorfreedoms differ across old/new arms and are explicitly structural; current3newarms matched; costprefix diagnostics, no costmatched whole-six claim',selection='no arm chosen or tuned using DEVGT; all3new+all3old predeclared',primary_endpoint='old primary_noise+radiusprogress; new legal truey noise gate; shared anchor when noisequalified',statistics='primary epsilon.075, qualification, offlineoracleMSE, actualuniform andexp(-Phi)heuristicactionMSE,rank,width; parentN10 statisticalunit; prefixK1/2/4/8 and actual cumulativeRHS',truth='online fresh executable never reads private truth or oldtruth; score only after complete480-row predictionseal',allocation='root supersedes workerinitialallocation: fresh155000RHS256/5000LU256; cumulativecandidate<=258824RHS256/9370LU256; global350000RHS256/10000LU256 unchanged',namespace_audit='new [20261030,101009], exact numeric role tuples not present in R9 frozen namespace[20261030,900901], R10DEV[20261030,1000701], R7 2026102807/09; pre-generation name+namespace audit in audit/candidate_fresh_namespace_audit.txt')
 pv.write(CONFIG,c);print('FRESH_FROZEN',CONFIG,flush=True)

def config():
 c=json.loads(CONFIG.read_text());en=json.loads(ENABLE.read_text())
 if not c['approved_by_root'] or en.get('protocol_id')!=c['id'] or en.get('enabled') is not True:raise RuntimeError('fresh root token missing')
 if c['source_blake2b']!=hashlib.blake2b(Path(__file__).read_bytes(),digest_size=32).hexdigest() or c['dev_source_blake2b']!=cg.source_id():raise RuntimeError('frozen source changed')
 return c

class FreshBook(cg.Book):
 def __init__(self,c,phase):
  self.c=c;self.job=phase;self.context={};self.jobdir=(PRIVATE if phase=='GENERATE' else OUT)/'physics'/phase;self.jobdir.mkdir(parents=True,exist_ok=False);self.started=time.perf_counter();self.cpu0=time.process_time();self.counts=tuple(k for k in set(cg.COUNTS+pv.Book.COUNTS) if k not in ('CPU_seconds','controlled_wall_seconds'));snap=OUT/'RESOURCE_SNAPSHOT.json';old=json.loads(snap.read_text()) if snap.exists() else {};self.total={k:old.get(k,0) for k in self.counts};self.priorwall=old.get('wall_seconds',0);self.priorcpu=old.get('CPU_seconds',0);self.seq=old.get('sequence',0);OUT.mkdir(parents=True,exist_ok=True);self.f=(OUT/'RESOURCE_LEDGER.jsonl').open('a');self.lock=(P/'candidate_physics_worker.lock').open('a');fcntl.flock(self.lock,fcntl.LOCK_EX|fcntl.LOCK_NB);self.event('JOB_BEGIN',pid=os.getpid(),phase=phase)
 def event(self,event,**kw):
  for k in self.counts:self.total[k]+=kw.get(k,0)
  self.seq+=1;self.f.write(json.dumps(pv.plain(dict(sequence=self.seq,utc=time.time(),event=event,phase=self.job,context=self.context.copy(),**kw,after=self.total.copy())),allow_nan=False)+'\n');self.f.flush()
 def check(self):
  if self.priorwall+time.perf_counter()-self.started>=self.c['caps']['controlled_wall_seconds']:raise cg.Stop('fresh controlled wallcap')
 def save(self,path,**data):
  st=time.perf_counter();self.event('ARRAY_SAVE_ATTEMPT',array_save_attempts=1)
  try:pv.save(path,**data)
  except Exception as e:self.event('ARRAY_SAVE_FAILURE',array_save_failures=1,error=repr(e));raise
  self.event('ARRAY_SAVE',array_saves=1,array_bytes=Path(path).stat().st_size,save_wall_seconds=time.perf_counter()-st)
 def close(self):
  wall=self.priorwall+time.perf_counter()-self.started;cpu=self.priorcpu+time.process_time()-self.cpu0;self.event('JOB_END',wall_seconds=wall,CPU_seconds=cpu);pv.write(OUT/'RESOURCE_SNAPSHOT.json',dict(**self.total,sequence=self.seq,wall_seconds=wall,CPU_seconds=cpu));self.f.close();self.lock.close()

def generate():
 c=config()
 if (OUT/'GENERATION_SEALED.json').exists():raise RuntimeError('no generation replay')
 PRIVATE.mkdir(parents=True,exist_ok=True);(PRIVATE/'.gitignore').write_text('*\n!.gitignore\n');book=FreshBook(c,'GENERATE');ph=pv.ExactROI(book)
 try:
  for entry in c['parents']:
   parent=entry['parent'];book.context=dict(parent=parent,role='PRIVATE_GENERATION');x=np.random.default_rng(np.random.SeedSequence(entry['prior_seed'])).normal(size=66);p=ph.evaluate(x,parent+'_truth_generation');rng=np.random.default_rng(np.random.SeedSequence(entry['noise_seed']));noise=pv.SIGMA/np.sqrt(2)*(rng.normal(size=(32,8))+1j*rng.normal(size=(32,8)));y=p.prediction+noise;book.save(PRIVATE/'truth'/f'{parent}.npz',x=x,q=p.q,material=p.material,prediction=p.prediction,noise=noise,observation=y,sigma=np.array(pv.SIGMA));book.save(OUT/'inputs/online'/f'{parent}.npz',observation=y,sigma=np.array(pv.SIGMA),parent=np.array(parent));print('FRESH_GENERATED',parent,flush=True)
  pv.write(OUT/'GENERATION_SEALED.json',dict(protocol_id=c['id'],utc=time.time(),parents=[e['parent'] for e in c['parents']],online_fields=['observation','sigma','parent'],truth_private=True))
 finally:ph.close();book.close()

def online():
 c=config()
 if not (OUT/'GENERATION_SEALED.json').exists() or (OUT/'PREDICTIONS_SEALED.json').exists():raise RuntimeError('generation missing or online alreadysealed')
 book=FreshBook(c,'ONLINE');ph=pv.ExactROI(book);rows=[];anchors={};cc=copy.deepcopy(c);cc['candidate']=c['candidate_new'];stopped=False
 try:
  L,meta=cg.bilateral(ph);pv.write(OUT/'BILATERAL_SETUP.json',meta)
  for entry in c['parents']:
   parent=entry['parent'];dest=OUT/'parents'/parent;y=pv.online_observation(c,parent);book.context=dict(parent=parent,role='ANCHOR');before=book.total.copy();ts=time.perf_counter();base=pv.anchor(ph,y,parent,dest/'anchor',c);anchors[parent]=dict(cost={k:book.total[k]-before[k] for k in book.counts},wall_seconds=time.perf_counter()-ts,Phi=ph.phi(base,y));cg.geometry_controls(base,parent);g0=pv.frame(ph.prior,base.x,base.J,8);book.save(dest/'ANCHOR_FRAME.npz',**g0);rng=np.random.default_rng(np.random.SeedSequence(entry['target_seed']));draws=[(rng.normal(size=512),rng.normal(size=66)) for _ in range(c['K'])]
   for k,(eta,rawxi) in enumerate(draws):
    ri=k//2;sign=1 if k%2==0 else -1;radius=c['radii'][ri];weak=np.random.default_rng(np.random.SeedSequence(entry['weak_seed_prefix']+[ri,0])).normal(size=8);fixed=None
    if g0['rank']>=8:fixed=g0['weaklift']@weak;fixed/=np.linalg.norm(g0['M']@fixed)
    rv=np.random.default_rng(np.random.SeedSequence(entry['random_seed_prefix']+[ri,0])).normal(size=66);norm=np.linalg.norm(g0['M']@rv);rv=rv/norm if norm>1e-12 else None
    # Rotate order per parent/slot, avoiding a fixed time-window arm advantage.
    order=ALLARMS[(entry['index']+k)%6:]+ALLARMS[:(entry['index']+k)%6]
    for arm in order:
     book.context=dict(parent=parent,arm=arm,k=k)
     if stopped:rows.append(dict(parent=parent,arm=arm,k=k,status='NOT_RUN_RIGHT_CENSORED',qualified=False,path=None,cost={},GT_used=False));continue
     if arm in NEWARMS:
      xi=rawxi.copy()
      if arm!='FULL66_RML':xi[64:]=0
      if arm=='BILATERAL_TARGET_RML_BLOCK':xi[:64]=L@xi[:64]
      row=cg.candidate(ph,base,y,eta,xi,cc,parent,arm,k)
     else:
      direction=rv if arm=='RANDOM66_FULL' else fixed;v=None if direction is None else sign*direction
      try:
       row=pv.candidate(ph,base,y,v,radius,arm,g0,dest/arm/f'{k:02d}',c);row.update(parent=parent,k=k,path=row.get('selected_path'),qualified=bool(row.get('primary_qualified')),GT_used=False)
      except (cg.Stop,pv.BudgetStop) as e:row=dict(parent=parent,arm=arm,k=k,status='RIGHT_CENSORED_GLOBAL_CAP',qualified=False,path=None,error=repr(e),cost={},GT_used=False)
      pv.write(dest/arm/f'{k:02d}'/'RECORD.json',row);print('FRESH_OLD',parent,arm,k,row['status'],row.get('Phi'),flush=True)
     rows.append(row)
     if 'CENSORED' in row['status']:stopped=True
   pv.write(dest/'COMPLETE.json',dict(parent=parent,rows=[r for r in rows if r['parent']==parent],anchor=anchors[parent],GT_used=False))
  pv.write(OUT/'PREDICTIONS_SEALED.json',dict(protocol_id=c['id'],utc=time.time(),GT_used=False,rows=rows,anchors=anchors,stopped=stopped,allmethods_frozen_before_DEVGT=True))
 finally:ph.close();book.close()

def score():
 # Separate executable phase: only here is private fresh truth opened.
 c=config();sealed=json.loads((OUT/'PREDICTIONS_SEALED.json').read_text())
 if len(sealed['rows'])!=480 or sealed['GT_used'] is not False:raise RuntimeError('full480predictionsealrequired')
 import candidate_statistics as cs
 pv.write(OUT/'GT_SCORING_OPENED.json',dict(utc=time.time(),prediction_seal_utc=sealed['utc'],fresh_parents=10,small_mechanism_check=True));rows=[]
 for entry in c['parents']:
  parent=entry['parent'];truth=pv.load_npz(PRIVATE/'truth'/f'{parent}.npz')['material'];anchor=pv.load_npz(OUT/'parents'/parent/'anchor/ANCHOR.npz');anchor['Phi']=sealed['anchors'][parent]['Phi'];gate=bool(anchor['Phi']<=c['noise_ball'])
  for arm in ALLARMS:
   records=sorted([r for r in sealed['rows'] if r['parent']==parent and r['arm']==arm],key=lambda r:r['k'])
   for K in (1,2,4,8):rows.append(dict(parent=parent,arm=arm,K=K,anchor_MSE=float(np.mean(abs(anchor['material']-truth)**2)),**cs.cloud_score(records[:K],anchor,truth,gate)))
 with (OUT/'METHOD_COMPARISON.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 summary={}
 for arm in ALLARMS:
  rs=[r for r in rows if r['arm']==arm and r['K']==8];summary[arm]=dict(parents=10,coverage_eps075_hits=sum(r['coverage_eps075'] for r in rs),mean_qualified_rate=float(np.mean([r['qualified_rate'] for r in rs])),mean_oracle_MSE=float(np.mean([r['oracle_MSE'] for r in rs if r['oracle_MSE'] is not None])),mean_uniform_action_MSE=float(np.mean([r['uniform_action_MSE'] for r in rs])),mean_likelihood_heuristic_action_MSE=float(np.mean([r['likelihood_heuristic_action_MSE'] for r in rs])),mean_pair_diameter=float(np.mean([r['pair_diameter'] for r in rs if r['pair_diameter'] is not None])),mean_cloud_rank=float(np.mean([r['cloud_rank'] for r in rs])),RHS256=sum(r['RHS256'] for r in rs),LU256=sum(r['LU256'] for r in rs),online_wall_seconds=sum(r['online_wall_seconds'] for r in rs),right_censored=sum(r['right_censored_count'] for r in rs))
 pv.write(OUT/'STATISTICS.json',dict(protocol_id=c['id'],summary=summary,rows=rows,anchors=sealed['anchors'],unit='parent N10',posterior_claim=False,small_mechanism_check=True,structural_differences=c['fairness'],allocation=c['allocation']))
 lines=['# Fresh10 candidate mechanism check','',c['exposure'], '', '| Arm | coverage .075 | qualification | oracle MSE | uniform action MSE | likelihood heuristic MSE | diameter | rank | RHS256 | LU256 |','|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
 for arm,s in summary.items():lines.append(f"| {arm} | {s['coverage_eps075_hits']}/10 | {s['mean_qualified_rate']:.3f} | {s['mean_oracle_MSE']:.6g} | {s['mean_uniform_action_MSE']:.6g} | {s['mean_likelihood_heuristic_action_MSE']:.6g} | {s['mean_pair_diameter']:.4f} | {s['mean_cloud_rank']:.2f} | {s['RHS256']} | {s['LU256']} |")
 lines+=['',c['fairness'], '', 'The shared anchor is a legal true-noise-qualified baseline whenever available. Oracle errors use private truth offline; uniform and exp(-Phi) actions use finite qualified clouds without posterior weighting claims. K1/2/4/8 cost-prefix rows are in METHOD_COMPARISON.csv. Fresh parents are independent, but N10 supports only a bounded mechanism check. All six methods and parameters were frozen before fresh generation; DEVGT did not select arms or retune them. All480onlineprediction records preceded fresh truth scoring.','', 'Amplitude block Schur solves are ordinary joint GN, not nonlinear profiling. Coupled nonlinear profile remains NOT_RUN. New/old latent caps and radius/progress qualification differ explicitly; equalK does not imply equal physics budget, geometry constraints, or posterior reference.']
 (P/'CANDIDATE_FRESH10_RESULTS.md').write_text('\n'.join(lines)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('action',choices=['freeze','generate','online','score']);a=p.parse_args();globals()[a.action]()
