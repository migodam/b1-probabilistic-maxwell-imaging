"""Truth scoring only after all online predictions are sealed; parent is statistical unit."""
from __future__ import annotations
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[n]='1'
import csv,json,sys,time
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P/'src'));import candidate_generation as cg
pv=cg.pv;H=cg.H;OUT=cg.OUT

def cloud_score(records,anchor,truth,gate):
 accepted=[];raw=[]
 for r in records:
  if r.get('path'):
   a=pv.load_npz(r['path']);raw.append(a['material'])
   if r['qualified']:accepted.append((a['material'],float(r['Phi'])))
 if gate:accepted.insert(0,(anchor['material'],float(anchor['Phi'])))
 unique=[]
 for m,phi in accepted:
  if not unique or min(np.linalg.norm(m-a[0])/24 for a in unique)>1e-9:unique.append((m,phi))
 result=dict(qualified_proposals=sum(bool(r.get('qualified')) for r in records),attempted_proposals=sum(bool(r.get('path')) for r in records),unique_qualified=len(unique),qualified_rate=sum(bool(r.get('qualified')) for r in records)/len(records),anchor_selected_proposals=sum(bool(r.get('anchor_selected')) for r in records),failure_count=sum(r['status']=='FAILED' for r in records),right_censored_count=sum('CENSORED' in r['status'] for r in records))
 if unique:
  X=np.stack([a[0] for a in unique]);phis=np.array([a[1] for a in unique]);err=np.mean(abs(X-truth)**2,axis=1);d=np.linalg.norm(X[:,None]-X[None,:],axis=2)/24;center=X.mean(axis=0);v=np.concatenate([(X-center).real,(X-center).imag],axis=1)/24;s=np.linalg.svd(v,compute_uv=False);keep=s>max((s[0] if len(s) else 0)*1e-8,1e-10);w=np.exp(-(phis-phis.min()));w/=w.sum();result.update(oracle_MSE=float(err.min()),oracle_RMS=float(np.sqrt(err.min())),coverage_eps075=bool(err.min()<=.075**2),pair_diameter=float(d.max()),mean_pair_distance=float(d[np.triu_indices(len(X),1)].mean()) if len(X)>1 else 0.,RMS_dispersion=float(np.sqrt(np.mean(abs(X-center)**2))),cloud_rank=int(keep.sum()),uniform_action_MSE=float(np.mean(abs(center-truth)**2)),likelihood_heuristic_action_MSE=float(np.mean(abs(w@X-truth)**2)),action_fallback=False)
 else:
  err=float(np.mean(abs(anchor['material']-truth)**2));result.update(oracle_MSE=None,oracle_RMS=None,coverage_eps075=False,pair_diameter=None,mean_pair_distance=None,RMS_dispersion=None,cloud_rank=0,uniform_action_MSE=err,likelihood_heuristic_action_MSE=err,action_fallback=True)
 result['all_selected_oracle_MSE_diagnostic']=None if not raw else float(np.min(np.mean(abs(np.stack(raw)-truth)**2,axis=1)))
 result['RHS256']=sum(r.get('cost',{}).get('reduced_RHS256',0) for r in records);result['LU256']=sum(r.get('cost',{}).get('reduced_LU256',0) for r in records);result['online_wall_seconds']=sum(r.get('wall_seconds',0) for r in records);result['CPU_seconds']=sum(r.get('CPU_seconds',0) for r in records);return result

def run():
 c=json.loads(cg.CONFIG.read_text());sealed=json.loads((OUT/'PREDICTIONS_SEALED.json').read_text())
 if sealed.get('GT_used') is not False:raise RuntimeError('invalid online seal')
 if len(sealed['rows'])!=len(c['parents'])*c['K']*len(c['arms']):raise RuntimeError('allproposalrecords must seal first')
 write=pv.write;write(OUT/'GT_SCORING_OPENED.json',dict(utc=time.time(),prediction_seal_utc=sealed['utc'],DEV_exposed=True,parents=c['parents'],no_confirmation_claim=True));rows=[]
 for parent in c['parents']:
  # Opening GT occurs here, after protocol and complete online seal.
  truth=pv.load_npz(H/'private/possibility_volume/truth'/f'{parent}.npz')['material'];anchor=pv.load_npz(H/'results/possibility_volume/parents'/parent/'anchor/ANCHOR.npz');y=pv.online_observation(dict(output=str(H/'results/possibility_volume')),parent);anchor['Phi']=float(np.sum(abs(anchor['prediction']-y)**2)/pv.SIGMA**2);gate=bool(anchor['Phi']<=c['noise_ball']);anchor_mse=float(np.mean(abs(anchor['material']-truth)**2))
  for arm in c['arms']:
   records=[r for r in sealed['rows'] if r['parent']==parent and r['arm']==arm];r=dict(parent=parent,arm=arm,role='matched current3arms',anchor_MSE=anchor_mse,**cloud_score(records,anchor,truth,gate));rows.append(r)
  old=json.loads((H/'results/possibility_volume/parents'/parent/'COMPLETE.json').read_text())
  for arm in pv.ARMS:
   records=[]
   for x in old['rows']:
    if x['arm']!=arm:continue
    records.append(dict(path=x.get('selected_path'),Phi=x.get('Phi'),qualified=bool(x.get('primary_qualified')),status=x['status'],cost=x.get('cost',{}),wall_seconds=x.get('wall_seconds',0),CPU_seconds=x.get('CPU_seconds',0)))
   rows.append(dict(parent=parent,arm='R9_'+arm,role='historical6F4J lowercap context only; not matched current cost',anchor_MSE=anchor_mse,**cloud_score(records,anchor,truth,gate)))
 keys=list(rows[0]);fp=P/'METHOD_COMPARISON.csv'
 with fp.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
 summary={}
 for arm in dict.fromkeys(r['arm'] for r in rows):
  rs=[r for r in rows if r['arm']==arm];summary[arm]=dict(parents=len(rs),coverage_eps075_hits=sum(r['coverage_eps075'] for r in rs),coverage_eps075_rate=float(np.mean([r['coverage_eps075'] for r in rs])),mean_qualified_rate=float(np.mean([r['qualified_rate'] for r in rs])),mean_unique_qualified=float(np.mean([r['unique_qualified'] for r in rs])),mean_oracle_MSE=float(np.mean([r['oracle_MSE'] for r in rs if r['oracle_MSE'] is not None])),mean_uniform_action_MSE=float(np.mean([r['uniform_action_MSE'] for r in rs])),mean_likelihood_heuristic_action_MSE=float(np.mean([r['likelihood_heuristic_action_MSE'] for r in rs])),mean_pair_diameter=float(np.mean([r['pair_diameter'] for r in rs if r['pair_diameter'] is not None])),mean_cloud_rank=float(np.mean([r['cloud_rank'] for r in rs])),RHS256=sum(r['RHS256'] for r in rs),LU256=sum(r['LU256'] for r in rs),online_wall_seconds=sum(r['online_wall_seconds'] for r in rs),right_censored_proposals=sum(r['right_censored_count'] for r in rs),failed_proposals=sum(r['failure_count'] for r in rs))
 comparisons={}
 for arm in c['arms'][1:]:
  base=[r for r in rows if r['arm']==c['arms'][0]];other=[r for r in rows if r['arm']==arm];comparisons[arm]={}
  for metric in ('oracle_MSE','uniform_action_MSE','likelihood_heuristic_action_MSE'):
   delta=[b[metric]-a[metric] for a,b in zip(base,other) if b[metric] is not None and a[metric] is not None];comparisons[arm][metric]=dict(paired_parent_delta=delta,mean_delta=float(np.mean(delta)),wins=sum(d<0 for d in delta),ties=sum(d==0 for d in delta),losses=sum(d>0 for d in delta),unit='independent parent, N8; proposals not replicate units')
 write(OUT/'STATISTICS.json',dict(protocol_id=c['id'],exposure=c['exposure'],epsilon=.075,gate=c['noise_ball'],summary=summary,paired_comparisons=comparisons,rows=rows,posterior_claim=False,weight_semantics='exp(-Phi) finite-cloud heuristic, no prior/proposal/Jacobian correction',historical_comparison='K16 equal; R9 old6F4J versus current12F8J: scientific descriptive context, no matched cost/time or mechanism attribution',global_ledger=json.loads((OUT/'RESOURCE_SNAPSHOT.json').read_text())))
 lines=['# Candidate results','',c['exposure']+'. Parent is the statistical unit (N=8). All online predictions were sealed before this scorer opened GT. These are heuristic finite material possibility clouds, not posterior samples.','', '| Arm | coverage .075 | qualified proposals | mean oracle MSE | uniform action MSE | likelihood heuristic MSE | mean diameter | mean rank | RHS256 | wall s |','|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
 for arm,s in summary.items():lines.append(f"| {arm} | {s['coverage_eps075_hits']}/8 | {s['mean_qualified_rate']:.3f} | {s['mean_oracle_MSE']:.6g} | {s['mean_uniform_action_MSE']:.6g} | {s['mean_likelihood_heuristic_action_MSE']:.6g} | {s['mean_pair_diameter']:.4f} | {s['mean_cloud_rank']:.2f} | {s['RHS256']} | {s['online_wall_seconds']:.2f} |")
 lines+=['','Current three arms have identical all66 corrector freedom, common random draws, and12F8J caps. Shape-target and bilateral-target differences are attributable to target distributions, not Schur elimination. Historical R9 controls are descriptive at16 proposals but have lower6F4J optimization caps and different construction; their timing is not a matched current-window comparison.','', 'Coupled nonlinear amplitude profile: NOT_RUN. Exact2x2 Schur elimination is ordinary joint GN and is not claimed as new physics acceleration. Fixed bilateral covariance is a normalized heuristic, not an inverse posterior covariance. Spectrum clipping is a declared design choice, not an invertibility theorem.','', 'Rejected/failed trials and right-censored proposals are retained in per-proposal records and resource ledger. Exact cached setup costs0 new full576RHS;264 historical setupRHS and reused old anchor/J costs must be included in standalone cost views.','', 'See METHOD_COMPARISON.csv for paired parent scores and results/candidates/STATISTICS.json for target-distribution comparisons and complete ledger.']
 (P/'CANDIDATE_RESULTS.md').write_text('\n'.join(lines)+'\n')
 import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 fig,ax=plt.subplots(1,2,figsize=(11,4));arms=c['arms']
 for i,arm in enumerate(arms):
  rs=[r for r in rows if r['arm']==arm];ax[0].plot(range(8),[r['oracle_MSE'] for r in rs],'o-',label=arm);ax[1].plot(range(8),[r['uniform_action_MSE'] for r in rs],'o-',label=arm)
 for a in ax:a.set_xlabel('exposed DEV parent');a.set_ylabel('material MSE');a.grid(alpha=.3)
 ax[0].set_title('Qualified cloud oracle (offline)');ax[1].set_title('Uniform finite-cloud action');ax[0].legend(fontsize=6);fig.tight_layout();fig.savefig(P/'figures/candidate_parent_MSE.png',dpi=160);plt.close(fig)
 print(json.dumps(summary,indent=2))
if __name__=='__main__':run()
