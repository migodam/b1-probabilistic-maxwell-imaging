"""Four parent-requested postseal figures, actual scores/arrays only, no solves."""
import os,time
START_CPU=time.process_time();START_WALL=time.perf_counter()
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
import argparse,csv,json,platform,sys
from pathlib import Path

def read(p):return json.loads(p.read_text())
def refresh_ecdf_inset(out):
 """Refresh only Figure1 from its existing complete288-row CSV."""
 import numpy as np
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 receipt_path=out/'PLOT_RECEIPT.json';receipt=read(receipt_path)
 with (out/'01_ALL48_PRIMARY.csv').open() as stream:rows=list(csv.DictReader(stream))
 methods=['cie','dual_cie','g19','dr','event_integer','direct'];labels=['CIE','Dual-CIE','G19','DR','Event integer','Direct'];colors=plt.get_cmap('tab10').colors
 plt.rcParams.update({'font.size':11,'pdf.fonttype':42})
 fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained');zoom=axes[0].inset_axes([.43,.12,.53,.39]);cpus=[]
 for i,(method,label) in enumerate(zip(methods,labels)):
  rr=[r for r in rows if r['method']==method]
  if len(rr)!=48:raise ValueError('Expected all48 rows for each method')
  vals=np.sort([float(r['NRMSE']) for r in rr if r['NRMSE'] not in ['', 'None']]);yy=np.arange(1,len(vals)+1)/48
  for ax in [axes[0],zoom]:ax.step(vals,yy,where='post',color=colors[i],label=f'{label} ({len(vals)}/48)')
  cpus.append(np.mean([float(r['external_child_CPU_seconds']) for r in rr if r['external_child_CPU_seconds'] not in ['', 'None']]))
 zoom.set(xlim=(0,.75),ylim=(.5,1));zoom.set_title('Zoom: same all48 ECDF',fontsize=9);zoom.tick_params(labelsize=8);zoom.grid(alpha=.2)
 axes[0].set(xlabel='Primary amplitude NRMSE (full recorded range)',ylabel='Fraction of all 48 objects',ylim=(0,1.02));axes[0].set_xlim(left=0);axes[0].legend(loc='center right',bbox_to_anchor=(.99,.75),fontsize=9);axes[0].grid(alpha=.2)
 axes[1].bar(labels,cpus,color=colors[:6]);axes[1].set(ylabel='Mean actual whole-child CPU (s), i9',xlabel='All 48 objects per method');axes[1].tick_params(axis='x',rotation=25);axes[1].grid(axis='y',alpha=.2)
 fig.suptitle('Fresh48: frozen primary estimates and actual run cost\nFull range and all48 denominator retained; inset magnifies the dense range',fontsize=13)
 for ext in ['png','pdf']:fig.savefig(out/'FIGURES'/('01_fresh48_ecdf_actual_cpu.'+ext),dpi=190,bbox_inches='tight')
 plt.close(fig);extra=time.process_time()-START_CPU;elapsed=time.perf_counter()-START_WALL
 # Existing total is authoritative; previous nested layout records stay refs.
 if receipt.get('ecdf_inset_refresh_attempt'):receipt.setdefault('earlier_ecdf_inset_refresh_attempts',[]).append(receipt['ecdf_inset_refresh_attempt'])
 receipt['ecdf_inset_refresh_attempt']=dict(CPU_seconds=extra,wall_seconds=elapsed,source='01_ALL48_PRIMARY.csv',only_figure='01_fresh48_ecdf_actual_cpu',new_Maxwell_LU=0)
 receipt['CPU_seconds']+=extra;receipt['wall_seconds']+=elapsed;receipt['timer_scope']='Sum of nonoverlapping original/layout plot attempts plus Figure1-only inset refresh; all attempt records nested references, wall active-stage sum'
 receipt.setdefault('input_source_files',[]).append(str((out/'01_ALL48_PRIMARY.csv').resolve()));write_text=json.dumps(receipt,indent=2)+'\n';receipt_path.write_text(write_text)
 print(json.dumps(dict(status='FIGURE1_ONLY_ECDF_INSET_REFRESHED',extra_CPU_seconds=extra,extra_wall_seconds=elapsed,total_CPU_seconds=receipt['CPU_seconds'],new_Maxwell_LU=0)))


def main():
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--authorization',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--refresh-layout',action='store_true');p.add_argument('--only-ecdf-inset',action='store_true');a=p.parse_args();root=a.root.resolve();auth=read(a.authorization)
 if auth.get('authorized_joint_offline_scoring') is not True or auth.get('authorized_aggregation') is not True:raise PermissionError('Joint parent authorization required')
 sys.path.insert(0,str(root/'src'));from seal_results import validate_seal
 for campaign in auth['required_campaigns']:
  seal=read(root/auth['seals'][campaign]);validate_seal(seal,require_authorization=True)
  if Path(seal['root']).resolve()!=root:raise ValueError('Seal root differs')
 if a.only_ecdf_inset:
  refresh_ecdf_inset(a.out.resolve());return
 import numpy as np
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.colors import Normalize
 out=a.out.resolve();previous=read(out/'PLOT_RECEIPT.json') if a.refresh_layout else None;out.mkdir(parents=True,exist_ok=a.refresh_layout);figdir=out/'FIGURES';figdir.mkdir(exist_ok=a.refresh_layout);sources=[str(a.authorization.resolve())];artifacts=[]
 def load(name):
  f=root/'analysis/combined'/(name+'.json');sources.append(str(f));return read(f)
 def table(name,rr):
  if not rr:return
  f=out/(name+'.csv');keys=list(dict.fromkeys(k for r in rr for k in r))
  with f.open('w',newline='') as stream:w=csv.DictWriter(stream,fieldnames=keys);w.writeheader();w.writerows(rr)
 def save(fig,name):
  for ext in ['png','pdf']:
   f=figdir/(name+'.'+ext);fig.savefig(f,dpi=190,bbox_inches='tight');artifacts.append(str(f.relative_to(out)))
  plt.close(fig)
 plt.rcParams.update({'font.size':11,'axes.titlesize':12,'axes.labelsize':11,'legend.fontsize':10,'pdf.fonttype':42,'ps.fonttype':42})
 records=load('ALL_ESTIMATOR_ROWS');deploy=load('CONDITIONAL_DEPLOYMENT_ALL48');costs=load('CONDITIONAL_DEPLOYMENT_COSTS');stats=load('PAIRED_OBJECT_BOOTSTRAP');core=[r for r in records if r['campaign']=='confirmation_core48' and r['view']=='primary'];methods=['cie','dual_cie','g19','dr','event_integer','direct'];labels=['CIE','Dual-CIE','G19','DR','Event integer','Direct'];colors=plt.get_cmap('tab10').colors
 cases=[f'r13_fresh_{i:04d}' for i in range(48)];lookup={(r['case_id'],r['method']):r for r in core};all48=[]
 fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained');cpus=[]
 for i,(method,label) in enumerate(zip(methods,labels)):
  rr=[lookup.get((cid,method),dict(case_id=cid,method=method,available=False)) for cid in cases];vals=np.array([r['NRMSE'] for r in rr if r.get('NRMSE') is not None]);vals.sort();axes[0].step(vals,np.arange(1,len(vals)+1)/48,where='post',label=f'{label} ({len(vals)}/48)',color=colors[i]);cpus.append(np.mean([r['external_child_CPU_seconds'] for r in rr if r.get('external_child_CPU_seconds') is not None]))
  all48 +=[{k:r.get(k) for k in ['case_id','method','available','availability_reason','NRMSE','IoU','external_child_CPU_seconds','noise_feasible']} for r in rr]
 axes[0].set(xlabel='Primary amplitude NRMSE (full recorded range)',ylabel='Fraction of all 48 objects',ylim=(0,1.02));axes[0].set_xlim(left=0);axes[0].legend();axes[0].grid(alpha=.2)
 axes[1].bar(labels,cpus,color=colors[:6]);axes[1].set(ylabel='Mean actual whole-child CPU (s), i9',xlabel='All 48 objects per method');axes[1].tick_params(axis='x',rotation=25);axes[1].grid(axis='y',alpha=.2)
 fig.suptitle('Fresh48: frozen primary estimates and actual run cost\nECDF denominator remains 48; no outlier trimming or new success threshold',fontsize=13);save(fig,'01_fresh48_ecdf_actual_cpu');table('01_ALL48_PRIMARY',all48)
 # Conditional all48 primary strategies, with9 actual trigger objects explicitly retained.
 strategies=sorted({r['method'] for r in deploy});base=[r for r in core if r['method']=='g19'];blook={r['case_id']:r for r in base};dlook={(r['case_id'],r['method']):r for r in deploy if r['view']=='primary'};triggers=sorted({r['case_id'] for r in costs if r['triggered']});bypass=sorted(set(cases)-set(triggers));means=[np.mean([r['NRMSE'] for r in base])];totals=[sum(r['external_child_CPU_seconds'] for r in base)];conditional=[];delta=np.full((len(triggers),len(strategies)),np.nan)
 for j,method in enumerate(strategies):
  rr=[dlook.get((cid,method),dict(case_id=cid,method=method,available=False)) for cid in cases];means.append(np.mean([r['NRMSE'] for r in rr if r.get('NRMSE') is not None]));totals.append(sum(r['deployment_external_CPU_seconds'] for r in costs if r['method']==method and r.get('deployment_external_CPU_seconds') is not None))
  for r in rr:
   cid=r['case_id'];change=r['NRMSE']-blook[cid]['NRMSE'] if r.get('NRMSE') is not None else None
   conditional.append(dict(case_id=cid,strategy=method,triggered=cid in triggers,bypass=cid in bypass,original_g19_NRMSE=blook[cid].get('NRMSE'),strategy_NRMSE=r.get('NRMSE'),delta_NRMSE=change,available=r.get('available')))
   if cid in triggers and change is not None:delta[triggers.index(cid),j]=change
 fig,axes=plt.subplots(1,3,figsize=(19,6),gridspec_kw={'width_ratios':[1,1,1.7]},layout='constrained');names=['Original G19']+[x.replace('_','\n') for x in strategies]
 axes[0].bar(range(len(names)),means,color=['#687c90']+list(colors[:len(strategies)]));axes[0].set_xticks(range(len(names)),names,rotation=30,ha='right');axes[0].set_ylabel('Mean primary amplitude NRMSE, all48');axes[1].bar(range(len(names)),totals,color=['#687c90']+list(colors[:len(strategies)]));axes[1].set_xticks(range(len(names)),names,rotation=30,ha='right');axes[1].set_ylabel('Total parent + actual rescue CPU (s), i9')
 bound=max(float(np.nanmax(abs(delta))),1e-12);im=axes[2].imshow(delta,cmap='coolwarm',vmin=-bound,vmax=bound,aspect='auto');axes[2].set_yticks(range(len(triggers)),triggers,fontsize=10);axes[2].set_xticks(range(len(strategies)),[x.replace('_','\n') for x in strategies],rotation=25,ha='right');axes[2].set_title(f'All {len(triggers)} trigger IDs: strategy minus original G19');fig.colorbar(im,ax=axes[2],label='Paired primary NRMSE change')
 for i in range(len(triggers)):
  for j in range(len(strategies)):axes[2].text(j,i,f'{delta[i,j]:+.3f}',ha='center',va='center',fontsize=8)
 fig.suptitle(f'Conditional deployment: all48 denominator; {len(triggers)} triggered + {len(bypass)} explicit G19 bypass objects\nCompeting strategy costs describe separate deployments and must not be summed as experimental consumption',fontsize=13);save(fig,'02_conditional_all48_cost_trigger_changes');table('02_ALL48_CONDITIONAL',conditional);table('02_ALL48_PARENT_RESCUE_COSTS',costs);table('02_BYPASS_IDS',[dict(case_id=cid,action='original G19 primary unchanged') for cid in bypass])
 # Fixed16 NN comparison uses same original uniform-start observations.
 fixed=cases[::3];nnrows=[];fig,axes=plt.subplots(2,3,figsize=(13,8),layout='constrained',gridspec_kw={'height_ratios':[2.4,1]})
 for j,method in enumerate(['dr','g19','direct']):
  original={r['case_id']:r for r in core if r['method']==method};nn={r['case_id']:r for r in records if r['campaign']=='confirmation_mechanism16' and r['method']=='nn_'+method and r['view']=='primary'};u=[];v=[]
  for i,cid in enumerate(fixed):
   aa=original.get(cid,{});bb=nn.get(cid,{});x=aa.get('NRMSE');y=bb.get('NRMSE');nnrows.append(dict(case_id=cid,method=method,uniform_NRMSE=x,NN_NRMSE=y,uniform_available=aa.get('available'),NN_available=bb.get('available'),uniform_CPU=aa.get('external_child_CPU_seconds'),NN_CPU=bb.get('external_child_CPU_seconds')))
   if x is not None and y is not None:axes[0,j].plot([0,1],[x,y],color=colors[i%10],alpha=.65,lw=1,marker='o',ms=3);u.append(x);v.append(y)
  axes[0,j].plot([0,1],[np.mean(u),np.mean(v)],color='black',lw=3,marker='D',label='Matched-object mean');axes[0,j].set_xticks([0,1],['Uniform 0.1','Frozen R8 NN']);axes[0,j].set_title(f'{method.upper()}: fixed16 ({len(u)}/16 pairs)');axes[0,j].set_ylabel('Primary amplitude NRMSE');axes[0,j].legend();axes[0,j].grid(alpha=.2)
  s=next(x for x in stats if x['block']=='NN_common_predeclared16' and x['view']=='primary' and x['method']=='nn_'+method and x['reference']==method);value=s['NRMSE_paired_mean'];lo=s['bootstrap95_low'];hi=s['bootstrap95_high'];axes[1,j].errorbar([value],[0],xerr=[[value-lo],[hi-value]],fmt='o',capsize=5,color='black');axes[1,j].axvline(0,color='gray',ls=':');axes[1,j].set_yticks([]);axes[1,j].xaxis.set_major_locator(plt.MaxNLocator(4));axes[1,j].tick_params(axis='x',labelsize=9);axes[1,j].set_xlabel('NN minus uniform paired mean NRMSE\nExisting object-bootstrap 95% interval');axes[1,j].set_title(f'{value:+.4f} [{lo:+.4f}, {hi:+.4f}]',fontsize=10)
 pre=root/'delegated/evidence/nnprep/PREPROCESS_RECEIPT.json';sources.append(str(pre));prep=read(pre);fig.suptitle(f'Frozen NN interfaces on the predeclared common16\nExtra preprocessing: Apple CPU {prep["process_CPU_seconds"]:.6f}s; historical supervised training cost unknown; hosts not pooled',fontsize=13);save(fig,'03_fixed16_nn_interfaces');table('03_FIXED16_NN_PAIRS',nnrows)
 # Only these parent-selected posthoc illustrations; every method shares row scale.
 illustrations=[cases[i] for i in [0,2,28,29,30,31]];fig,axes=plt.subplots(6,7,figsize=(20,17),layout='constrained');image_rows=[]
 for i,cid in enumerate(illustrations):
  gtpath=root/'sealed_truth/fresh48'/(cid+'.npz');sources.append(str(gtpath))
  with np.load(gtpath,allow_pickle=False) as z:gt=z['material_coarse'].reshape(24,24).copy()
  vmax=float(gt.max());arrays=[gt];counts=[0]
  for method in methods:
   row=lookup[(cid,method)];path=root/'RESULTS/confirmation_core48'/row['task_id']/'CANDIDATES.npz';sources.append(str(path))
   with np.load(path,allow_pickle=False) as z:m=z['m'][int(row['selection_candidate'])].reshape(24,24).copy()
   arrays.append(m);counts.append(int(np.sum(m>vmax)))
  norm=Normalize(vmin=0,vmax=vmax,clip=True)
  for j,m in enumerate(arrays):
   im=axes[i,j].imshow(m,cmap='viridis',norm=norm,origin='lower',interpolation='nearest');axes[i,j].set_xticks([]);axes[i,j].set_yticks([])
   if i==0:axes[i,j].set_title((['GT']+labels)[j],fontsize=12)
   if j==0:axes[i,j].set_ylabel(cid+f'\nrow max χ={vmax:.3f}',fontsize=10)
   if j:axes[i,j].set_xlabel(f'NRMSE={lookup[(cid,methods[j-1])]["NRMSE"]:.3f}\nabove-row-max pixels: {counts[j]}',fontsize=9)
   image_rows.append(dict(case_id=cid,column=(['GT']+methods)[j],row_vmin=0,row_vmax=vmax,saturated_above=counts[j],posthoc_illustration=True))
  fig.colorbar(im,ax=axes[i,:].tolist(),fraction=.012,pad=.012,label='Physical χ; common row scale')
 fig.suptitle('Six explicitly posthoc illustrative cases — not a selected confirmation subset\nAll methods use saved primary arrays; GT and all predictions share each row’s physical color scale; saturation counts shown',fontsize=14);save(fig,'04_posthoc_six_case_materials');table('04_ILLUSTRATION_SCALES_SATURATION',image_rows)
 receipt=dict(status='FOUR_REQUESTED_FIGURES_COMPLETE_PARENT_REVIEW_PENDING',host='Apple_local',platform=platform.platform(),python=sys.version,CPU_seconds=time.process_time()-START_CPU,wall_seconds=time.perf_counter()-START_WALL,input_source_files=sorted(set(sources)),artifacts=artifacts,new_Maxwell_LU=0,new_RHS=0,new_optimization=0,truth_access='Authorized frozen coarse truth for only six posthoc illustrations',timer_scope='Whole figure process from pre-import timer through CSV/figure serialization; final receipt write/exit outside interval; per-figure timers not added',qualification='Existing allcase galleries retained; no new success gate, GT selection or scientific winner interpretation')
 if previous:
  receipt['previous_layout_attempts']=previous.get('previous_layout_attempts',[])+[dict(CPU_seconds=previous['CPU_seconds'],wall_seconds=previous['wall_seconds'])];receipt['current_layout_CPU_seconds']=receipt['CPU_seconds'];receipt['current_layout_wall_seconds']=receipt['wall_seconds'];receipt['CPU_seconds']+=previous['CPU_seconds'];receipt['wall_seconds']+=previous['wall_seconds'];receipt['timer_scope']='Sum of two nonoverlapping whole plot attempts, including initial layout and tick-label refresh; stage refs not added again; wall is active-stage sum'
 (out/'PLOT_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:receipt[k] for k in ['status','CPU_seconds','wall_seconds','new_Maxwell_LU']}))
if __name__=='__main__':main()
