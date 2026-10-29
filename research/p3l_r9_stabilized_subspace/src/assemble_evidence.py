"""Read completed evidence only. No solver, model, GT, training or new experiment."""
from pathlib import Path
import json,csv,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[1]
def read(s):return json.loads((P/s).read_text())
def write(s,d):
 (P/s).parent.mkdir(parents=True,exist_ok=True);(P/s).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def csvout(s,rows):
 if not rows:return
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with (P/s).open('w') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
(P/'figures').mkdir(exist_ok=True)
e=read('results/exploration/SUMMARY.json');v=read('results/exploration_rev1/SUMMARY.json');r=read('results/stable_roi/RESULTS.json');pr=read('results/probability/aggregate_costs.json')
csvout('METHOD_COMPARISON.csv',[dict(version='original',**x) for x in e['table']]+[dict(version='radial_revision',**x) for x in v['table']])
s=read('audit/exploration/E1_HISTORICAL_SUBSPACE_STABILITY.json');ss=[]
fig,ax=plt.subplots(figsize=(7.5,4.5))
for row in s['rows']:
 ns=[n for n in row['nodes'] if 'weak_max_angle_degrees' in n]
 for n in ns:ss.append(dict(parent=row['parent'],**{k:v for k,v in n.items() if not isinstance(v,list)}))
 ax.plot([n['material_RMS_from_anchor'] for n in ns],[n['weak_max_angle_degrees'] for n in ns],marker='o',label=row['parent'].split('_')[-1])
ax.set(xlabel='Full-grid material RMS',ylabel='Maximum principal angle (degrees)',title='Reused R7 trajectories: rotation alone is not an outcome')
ax.legend(ncol=4,title='Exposed DEV parent',fontsize=8);ax.grid(alpha=.2);fig.tight_layout();fig.savefig(P/'figures/subspace_rotation.png',dpi=180);plt.close(fig)
csvout('STABILITY_METRICS.csv',ss)
fig,axs=plt.subplots(1,2,figsize=(12,4.8))
labels=['Native Full-J','Fixed-ray + Full-J','Dynamic-ray + Full-J','Frozen-J','Uncorrected']
for ax,radius in zip(axs,[.01,.02]):
 ar=['A_NATIVE_FULLJ','B_RAY','C_RAY','F_FROZENJ','U_UNCORRECTED'];allrows=e['table']+v['table'];rs=[next(x for x in allrows if x['arm']==a and x['radius']==radius) for a in ar]
 x=np.arange(len(rs));ax.bar(x-.17,[z['qualified'] for z in rs],.34,label='All physical gates');ax.bar(x+.17,[z['pair_le1'] for z in rs],.34,label='Pair distance <= 1')
 ax.set_xticks(x,labels,rotation=24,ha='right');ax.set(ylim=(0,9),title=f'Material RMS {radius:.2f}',ylabel='Count / 8 exposed DEV parents');ax.grid(axis='y',alpha=.2)
axs[0].legend(fontsize=8);fig.suptitle('Endpoint qualification; no posterior mass or fresh confirmation');fig.tight_layout();fig.savefig(P/'figures/endpoint_qualification.png',dpi=180);plt.close(fig)
rs=[]
for row in r['rows']:rs.append(dict(parent=row['parent'],role=row['role'],repeat=row['repeat'],direct_s=row['wall']['direct'],reduced_s=row['wall']['reduced'],scope='matched fixed-state forward+Jacobian kernel'))
csvout('FIXED_ROI_COST.csv',rs)
ids=sorted(set(x['parent'] for x in rs));full=[sum(x['direct_s'] for x in rs if x['parent']==p)/6*1000 for p in ids];red=[sum(x['reduced_s'] for x in rs if x['parent']==p)/6*1000 for p in ids]
fig,ax=plt.subplots(figsize=(8,4.5));x=np.arange(8);ax.bar(x-.18,full,.36,label='Original 576 LU');ax.bar(x+.18,red,.36,label='Exact ROI 256 LU');ax.set_xticks(x,[i.split('_')[-1] for i in ids]);ax.set(xlabel='Exposed DEV parent (2 states x 3 repeats)',ylabel='Mean kernel wall / ms',title='Ordinary exact support elimination: same full nonlinear model');ax.legend();ax.grid(axis='y',alpha=.2);fig.tight_layout();fig.savefig(P/'figures/exact_roi_cost.png',dpi=180);plt.close(fig)
# Actual material arrays, no truth selection: first fixed parent in manifest.
base=P/'results/exploration_rev1/parents/R7_LOCKED_000'
ms=[]
for arm in ['B_RAY','C_RAY']:
 with np.load(base/arm/'RMS_0.0200/ENDPOINT.npz') as z:ms.append(z['material'].copy())
with np.load(P/'results/stable_roi/R7_LOCKED_000_anchor_r0.npz') as z:m0=z['material'].copy()
fig,axs=plt.subplots(1,4,figsize=(12,3.4))
for ax,m,title in zip(axs,[m0,ms[0],ms[1],abs(ms[0]-ms[1])],['Physics anchor','Fixed weak ray','Dynamic weak space','Absolute B-C difference']):
 im=ax.imshow(np.real(m).reshape(24,24),origin='lower',cmap='viridis');ax.set_title(title,fontsize=10);ax.set_xticks([]);ax.set_yticks([]);fig.colorbar(im,ax=ax,fraction=.046)
fig.suptitle('Predetermined DEV 000; valid finite material exploration, no GT image selection',fontsize=11);fig.tight_layout();fig.savefig(P/'figures/material_endpoints.png',dpi=180);plt.close(fig)
# Import raw append-only rows once; original journals remain authoritative.
branches=[('exploration','results/exploration/RESOURCE_LEDGER.jsonl'),('exploration_rev1','results/exploration_rev1/RESOURCE_LEDGER.jsonl'),('probability','audit/probability/ledger.jsonl'),('stable_roi','results/stable_roi/COST_LEDGER.jsonl')]
ledger=P/'COST_LEDGER.jsonl';seen=set()
if ledger.exists():
 for line in ledger.read_text().splitlines():
  z=json.loads(line);seen.add((z.get('source'),z.get('source_line')))
accounts=[]
with ledger.open('a') as f:
 for branch,rel in branches:
  rows=[json.loads(l) for l in (P/rel).read_text().splitlines()];rhs=lu=0;wall=cpu=0
  for n,row in enumerate(rows,1):
   if branch.startswith('exploration'):a=row.get('RHS',0);b=row.get('LU',0)
   elif branch=='probability':a=row.get('full_complex_RHS',0);b=row.get('full_LU',0)
   else:a=row.get('full_RHS',0);b=row.get('full_LU',0)
   rhs+=a;lu+=b
   if row['event']=='JOB_END':wall+=row.get('wall_seconds',row.get('wall',0));cpu+=row.get('CPU_seconds',row.get('CPU',0))
   if (rel,n) not in seen:f.write(json.dumps(dict(source=rel,source_line=n,branch=branch,entry=row),ensure_ascii=False)+'\n')
  accounts.append(dict(branch=branch,full_RHS=rhs,full_LU=lu,job_wall_seconds_sum=wall,CPU_seconds_sum=cpu,source=rel))
resource=dict(branches=accounts,new_full_RHS=sum(x['full_RHS'] for x in accounts),full_LU=sum(x['full_LU'] for x in accounts),reduced_RHS256=pr['reduced_complex_RHS256']+r['reduced_RHS256'],reduced_LU256=r['reduced_LU256'],Schur256=1,initial_full_RHS_cap=28000,new_physics_after_formula_steering=0,CPU_seconds_sum=sum(x['CPU_seconds_sum'] for x in accounts),job_wall_seconds_sum=sum(x['job_wall_seconds_sum'] for x in accounts),note='Full and reduced RHS are different-sized systems, never convert one into the other; reused historical setup cost is not newly paid.')
write('RESOURCE_SUMMARY.json',resource)
write('RESULTS_SUMMARY.json',dict(status='SCOPED_DEVELOPMENT_COMPLETE_FORMULA_ATLAS_ADDED',fresh_confirmation='NOT_RUN_USER_FORMULA_FIRST_AND_NO_20_PERCENT_SUBSPACE_GAIN',full66_posterior='UNQUALIFIED',continuous_exact_fiber='NOT_ESTABLISHED',fixed_ray=dict(qualified_each_radius='8/8 DEV',pair01='8/8',pair02='7/8',cost_advantage='NOT_ESTABLISHED'),exact_roi=dict(kernel_setup_inclusive_saving=r['setup_inclusive_saving'],scope=r['scope'],full_system='NOT_RUN'),probability=read('results/probability/STATUS.json'),resource=resource))
print(json.dumps(resource,ensure_ascii=False))
