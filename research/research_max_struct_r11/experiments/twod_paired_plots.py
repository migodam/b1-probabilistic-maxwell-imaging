"""Complete root-authorized paired8 visual report; no solver/FW or rule changes."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import argparse,json,csv,time,html
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'delegated/twod/crossgrid8_paired';METHODS=['cie_single','cie_dual','global_dual','fullG19'];LABELS=['CIE','Dual CIE','Global Dual','Full G19'];GRIDS=['matched24','cross48to24'];COLORS=['#2166ac','#4d9221','#762a83','#d95f02']
def read(p):return json.loads(Path(p).read_text())
def csvrows(p):
 with Path(p).open(newline='') as f:return list(csv.DictReader(f))
def table(p,rows):
 with Path(p).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
def number(x):
 try:return float(x) if x not in (None,'') else None
 except (TypeError,ValueError):return None
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--scores',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--authorization',type=Path,required=True);a=ap.parse_args();wall=time.perf_counter()
 if read(a.authorization)['status']!='ROOT_EXPLICIT_POST64SEAL_SCORE' or read(a.scores/'RECEIPT.json')['status']!='SCORED_POST64_PUBLIC_SEAL':raise PermissionError('Completed authorized paired64 scores required')
 final=csvrows(a.scores/'FINAL64_SCORES.csv');allrows=csvrows(a.scores/'ALL_SELECTION_SCORES.csv');cases=sorted({r['case_id'] for r in final});lookup={(r['case_id'],r['grid_condition'],r['method']):r for r in final}
 if len(final)!=64 or len(cases)!=8 or len(lookup)!=64:raise ValueError('All8 objects,2 grids,4 methods required')
 a.output.mkdir(parents=True,exist_ok=False);(a.output/'PLOT_GATE_BEFORE_PRIVATE.json').write_text(json.dumps(dict(status='PASS_POSTSCORE_ALL64',root_authorization=str(a.authorization),GT_reopened_for_display_only=True,all_objects_included=True,new_FW=0),indent=2)+'\n');plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10});artifacts=[]
 def save(fig,name):
  for ext in ('png','svg'):
   p=a.output/(name+'.'+ext);fig.savefig(p,dpi=180,bbox_inches='tight');artifacts.append(dict(path=p.name,bytes=p.stat().st_size))
  plt.close(fig)
 selected=np.load(D/'public_seal/SELECTED_MATERIALS.npz',allow_pickle=False)['m'];selection=read(D/'public_seal/SELECTIONS.json');idx={(r['case_id'],r['grid_condition'],r['method']):r['selected']['material_index'] if r['selected'] else None for r in selection['rows'] if r['view']=='FINAL_CAP_ELIGIBLE'}
 # Authorized private reads only below this gate; material target is the frozen
 # projected coarse target, not a newly regenerated/discretized truth.
 truth={c:np.load(ROOT/'sealed_truth/crossgrid8_paired'/(c+'.npz'),allow_pickle=False)['m'].reshape(24,24) for c in cases};pages=[]
 for family in sorted({r['family'] for r in final}):
  cc=sorted({r['case_id'] for r in final if r['family']==family})
  if len(cc)!=2:raise ValueError('Each of4 families has high/mixed independent materials')
  fig,axs=plt.subplots(4,5,figsize=(13.8,10.6),squeeze=False);fig.subplots_adjust(left=.15,right=.94,top=.93,bottom=.05,wspace=.07,hspace=.2);image=None
  for i,(c,g) in enumerate((c,g) for c in cc for g in GRIDS):
   for j in range(5):
    ax=axs[i,j];ax.set_xticks([]);ax.set_yticks([]);key=(c,g,METHODS[j-1]) if j else None;k=idx.get(key) if j else None;value=truth[c] if j==0 else selected[k].reshape(24,24) if k is not None else None
    if value is None:ax.set_facecolor('#eee');ax.text(.5,.5,'FAILED / no paid state',ha='center',va='center',transform=ax.transAxes,fontsize=8)
    else:image=ax.imshow(value,origin='lower',vmin=0,vmax=20,cmap='viridis',interpolation='nearest')
    if i==0:ax.set_title('Projected GT' if not j else LABELS[j-1],fontsize=10)
    if j==0:ax.set_ylabel(c+'\n'+lookup[c,g,'cie_single']['regime']+' / '+g,rotation=0,ha='right',va='center',labelpad=9,fontsize=8)
    if j and lookup[key]['task_status'].startswith('FAILED'):
     for spine in ax.spines.values():spine.set_edgecolor('#b2182b');spine.set_linewidth(2)
     ax.text(.5,.03,'FAILED: paid state',transform=ax.transAxes,ha='center',fontsize=7,color='#b2182b')
  fig.colorbar(image,cax=fig.add_axes([.956,.16,.012,.65]),label='Material contrast');fig.suptitle(family+' — both independent materials, both grid conditions; common0..20 scale',fontsize=13);name='MATERIALS_'+family;save(fig,name);pages.append(name+'.png')
 fig,axs=plt.subplots(2,4,figsize=(14.5,7.8),squeeze=False)
 for j,m in enumerate(METHODS):
  for i,(metric,label) in enumerate([('NRMSE','NRMSE'),('fixed_support_IoU','Support IoU')]):
   ax=axs[i,j];points=[]
   for c in cases:
    l,r=lookup[c,'matched24',m],lookup[c,'cross48to24',m];x,y=number(l.get(metric)),number(r.get(metric))
    if x is not None and y is not None:points.append((x,y));ax.scatter(x,y,c='#2166ac' if l['regime']=='high' else '#d95f02',s=35);ax.annotate(c[-4:],(x,y),xytext=(3,3),textcoords='offset points',fontsize=7)
   upper=1 if metric.endswith('IoU') else max([v for p in points for v in p]+[.1])*1.1;ax.plot([0,upper],[0,upper],'--',c='.6',lw=.8);ax.set_xlim(0,upper);ax.set_ylim(0,upper);ax.set_xlabel('Matched24 '+label);ax.set_ylabel('Cross48to24 '+label);ax.set_title(LABELS[j]+'; paired8/8');ax.grid(alpha=.2)
 fig.suptitle('Same8 materials paired across conditions; blue=high, orange=mixed');fig.tight_layout(rect=[0,0,1,.95],h_pad=3);save(fig,'PAIRED_GRID_MATERIAL_METRICS')
 fig,axs=plt.subplots(2,2,figsize=(12,8.5));availability=[]
 for i,g in enumerate(GRIDS):
  for j,(metric,label) in enumerate([('NRMSE','Median NRMSE'),('fixed_support_IoU','Median support IoU')]):
   ax=axs[i,j]
   for k,m in enumerate(METHODS):
    ys=[]
    for cap in [60,120,240,480]:
     rr=[r for r in allrows if r['grid_condition']==g and r['method']==m and r['view']=='CPU_'+str(cap)];v=[number(r.get(metric)) for r in rr];v=[x for x in v if x is not None];ys.append(float(np.median(v)) if v else np.nan)
     if j==0:availability.append(dict(grid_condition=g,method=m,CPU_cap=cap,available=len(v),independent_materials=8))
    ax.plot([60,120,240,480],ys,'o-',c=COLORS[k],label=LABELS[k])
   ax.set_title(g);ax.set_xlabel('Allowed whole-child CPU cap (seconds)');ax.set_ylabel(label);ax.grid(alpha=.2)
 axs[0,0].legend(fontsize=8);fig.suptitle('Presealed paid original-Phi incumbents; no budget padding');fig.tight_layout(rect=[0,0,1,.96]);save(fig,'FIXED_CPU_CURVES');table(a.output/'CURVE_AVAILABILITY.csv',availability)
 fig,ax=plt.subplots(figsize=(10,6))
 for k,m in enumerate(METHODS):
  for g,marker in [('matched24','o'),('cross48to24','^')]:
   rr=[r for r in final if r['method']==m and r['grid_condition']==g];ax.scatter([number(r['actual_complete_child_CPU_seconds']) for r in rr],[number(r['NRMSE']) for r in rr],c=COLORS[k],marker=marker,label=LABELS[k]+' '+g,s=35,alpha=.75)
 ax.axvline(480,c='.5',ls='--',lw=.8);ax.set_xlabel('Actual complete valid-child CPU seconds');ax.set_ylabel('Final cap-eligible NRMSE');ax.set_title('Same480 allowance; actual consumption differs');ax.legend(fontsize=8,ncol=2);ax.grid(alpha=.2);fig.tight_layout();save(fig,'ACTUAL_CPU_COST')
 phase_names=['pure_low','pure_high','pure_background','partial_volume_boundary'];phase=[];heat=[];yl=[]
 for g in GRIDS:
  for m,label in zip(METHODS,LABELS):
   rr=[r for r in final if r['grid_condition']==g and r['method']==m];hr=[]
   for name in phase_names:
    v=[number(r.get(name+'_RMSE')) for r in rr];v=[x for x in v if x is not None];median=float(np.median(v)) if v else None;hr.append(median if median is not None else np.nan);phase.append(dict(grid_condition=g,method=m,phase=name,finite_regions=len(v),independent_materials=8,median_RMSE=median,mean_RMSE=float(np.mean(v)) if v else None,total_phase_cells=sum(int(r[name+'_cells']) for r in rr)))
   heat.append(hr);yl.append(g+' / '+label)
 fig,ax=plt.subplots(figsize=(10,6));im=ax.imshow(np.asarray(heat),aspect='auto',cmap='magma');ax.set_xticks(range(4),['Pure low','Pure high','Background','Partial-volume boundary']);ax.set_yticks(range(8),yl)
 for i,row in enumerate(heat):
  for j,v in enumerate(row):ax.text(j,i,'NA' if not np.isfinite(v) else f'{v:.3f}',ha='center',va='center',color='white' if not np.isfinite(v) or v>np.nanmax(heat)/2 else 'black',fontsize=9)
 fig.colorbar(im,ax=ax,label='Median phase RMSE over nonempty regions');ax.set_title('Frozen fine-defined phase partition; empty regions are not zero error');fig.tight_layout();save(fig,'PHASE_REGION_ERRORS');table(a.output/'PHASE_SUMMARY.csv',phase)
 model=csvrows(a.scores/'CACHED_MODEL_MISMATCH_8.csv');fig,ax=plt.subplots(figsize=(11,5));x=np.arange(len(model));ax.bar(x-.18,[float(r['over_matched_raw_sigma']) for r in model],width=.36,label='Cached model mismatch / matched sigma');ax.bar(x+.18,[float(r['over_cross_raw_sigma']) for r in model],width=.36,label='Cached model mismatch / cross sigma');ax.axhline(float(final[0]['CN99_radius']),ls='--',c='.4',label='CN99 raw-noise radius');ax.set_xticks(x,[r['case_id'][-4:] for r in model]);ax.set_ylabel('Frobenius norm / raw sigma');ax.set_xlabel('Independent material');ax.set_title('Cached coarse/fine cleanY difference; zero new forward calls');ax.legend(fontsize=8);ax.grid(axis='y',alpha=.2);fig.tight_layout();save(fig,'CACHED_MODEL_MISMATCH')
 summary=read(a.scores/'SUMMARY.json');lines=['# Paired-grid8 complete descriptive report','','Eight independent materials, each observed under matched24 and cross48to24 and solved by four fixed methods. All64 outcomes and failures are retained. Grid/method repeats are not16 or64 independent objects. Root-specific postseal authorization precedes these private reads; no new FW/adjoint/optimization or GT-selected state.','','Process deviation: the separate diagnostic private-text metadata search is disclosed in [INCIDENT_AUDIT.md](../../../diagnostics/private_metadata_search_incident/INCIDENT_AUDIT.md). It printed one private geometry-metadata line, not a GT array; broader private-text scan extent is unknown. No online solver/selection used it and no rule or draw changed. This report does not claim a perfectly blinded whole research process.','','|Condition|Method|Objects|Material PASS|Mean NRMSE|Median NRMSE|Actual valid-child CPU sum|CN99 pass|','|---|---|---:|---:|---:|---:|---:|---:|']
 for s in summary:
  rr=[r for r in final if r['grid_condition']==s['grid_condition'] and r['method']==s['method']];lines.append('|'+s['grid_condition']+'|'+s['method']+'|8|'+str(s['material_success'])+'|'+str(s['mean_NRMSE'])+'|'+str(s['median_NRMSE'])+'|'+str(sum(float(r['actual_complete_child_CPU_seconds']) for r in rr))+'|'+str(sum(r.get('CN99_pass')=='True' for r in rr))+'|')
 lines+=['','Success retains the frozen .15 NRMSE/.7 support-IoU/box criteria. Phase RMSEs use the frozen fine-defined pure-low/high/background/partial-volume partition; empty regions are explicitly missing. Cached cleanY mismatch is a discrete-model control, not a continuous-PDE error bound or an inferred recovery certificate. Same allowed CPU cap and actual paid CPU are distinct; commonB paired readouts do not erase either original full-child cost.','','Complete selected-state and paired-grid tables: ALL_SELECTION_SCORES.csv, FINAL64_SCORES.csv, FINAL_PAIRED_GRID_DIFFERENCES.csv, ACTUAL_SPEND_PAIRED_DIFFERENCES.csv and CACHED_MODEL_MISMATCH_8.csv in the sibling scored directory. Scientific interpretation remains root-owned.','']
 for name in pages+['PAIRED_GRID_MATERIAL_METRICS.png','FIXED_CPU_CURVES.png','ACTUAL_CPU_COST.png','PHASE_REGION_ERRORS.png','CACHED_MODEL_MISMATCH.png']:lines+=['!['+name+']('+name+')','']
 (a.output/'REPORT.md').write_text('\n'.join(lines));(a.output/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><style>body{max-width:1400px;margin:30px auto;font-family:system-ui}img{max-width:100%}p{line-height:1.6}</style><h1>Paired-grid8: all objects and both conditions</h1><p>Eight independent materials; all64 fixed outcomes. See REPORT.md for complete protocol/process-deviation disclosure and factual tables.</p>'+''.join('<figure><img src="'+html.escape(name)+'"><figcaption>'+html.escape(name)+'</figcaption></figure>' for name in pages+['PAIRED_GRID_MATERIAL_METRICS.png','FIXED_CPU_CURVES.png','ACTUAL_CPU_COST.png','PHASE_REGION_ERRORS.png','CACHED_MODEL_MISMATCH.png']))
 (a.output/'RECEIPT.json').write_text(json.dumps(dict(status='RENDERED_ALL64_PAIRED_GRID_POSTSCORE',independent_materials=8,grid_conditions=2,methods=4,material_pages=4,all64_outcomes_included=True,success_filtered=False,artifacts=artifacts,font='Bundled DejaVu Sans, English labels',CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-wall,new_FW=0,new_adjoint=0,process_deviation_audit='delegated/diagnostics/private_metadata_search_incident/INCIDENT_AUDIT.md',scientific_interpretation='ROOT_REVIEW_REQUIRED'),indent=2)+'\n');print('RENDERED_ALL64_PAIRED_GRID_POSTSCORE')
if __name__=='__main__':main()
