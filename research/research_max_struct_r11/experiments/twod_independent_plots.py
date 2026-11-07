"""Postscore all-object figures and factual report. No scoring or new physics."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import argparse,csv,json,time,html
from pathlib import Path
ENTRY_CPU=time.process_time()
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager,ft2font
from twod_recovery import active_task
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
METHODS=('cie_single','cie_dual','global_dual','fullG19','freezeA','contractedT')
LABELS=('CIE','Dual CIE','Global Dual','Full G19','Freeze A','Contracted T')
COLORS=plt.get_cmap('tab10').colors[:6]
def read(p):return json.loads(Path(p).read_text())
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def loadcsv(p):
 with Path(p).open(newline='') as f:return list(csv.DictReader(f))
def finite(v):
 try:return np.isfinite(float(v))
 except (TypeError,ValueError):return False
def number(v):return float(v) if finite(v) else None
def npz(p):
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def method(t):return t.get('material_readout',{}).get('label',t['arm'])
def resolve(p):
 p=Path(str(p).replace('\\','/'));return p if p.is_absolute() else REPO/p
def font_policy(language):
 english=Path(matplotlib.get_data_path())/'fonts/ttf/DejaVuSans.ttf';chosen=english;actual='en'
 if language=='zh':
  candidates=[Path('/System/Library/Fonts/PingFang.ttc'),Path('/System/Library/Fonts/Hiragino Sans GB.ttc'),Path('C:/Windows/Fonts/msyh.ttc'),Path('C:/Windows/Fonts/simhei.ttf')]
  candidates += [Path(p) for p in font_manager.findSystemFonts() if any(x in p.lower() for x in ('notosanscjk','sourcehansans'))]
  required={ord(c) for c in '真值失败全部对象材料对比度固定预算实际耗费未评分公开数据残差计算成本无状态统一色标缺失数另表列明已付重采样并非同蓝橙χε'}
  for candidate in candidates:
   if candidate.exists():
    try:
     if required.issubset(ft2font.FT2Font(str(candidate)).get_charmap()):chosen=candidate;actual='zh';break
    except (RuntimeError,OSError):pass
 font_manager.fontManager.addfont(str(chosen));prop=font_manager.FontProperties(fname=str(chosen));plt.rcParams.update({'font.family':prop.get_name(),'font.size':10,'axes.unicode_minus':False,'savefig.facecolor':'white'})
 return dict(requested=language,actual=actual,path=str(chosen),font=prop.get_name(),fallback=actual!=language)

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--scores',type=Path,required=True);ap.add_argument('--manifest',type=Path,required=True);ap.add_argument('--truth-dir',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--config',type=Path,default=ROOT/'configs/twod_independent48_plots.json');ap.add_argument('--language',choices=['en','zh'],default='en');ap.add_argument('--authorize-postscore',required=True);a=ap.parse_args();wall=time.perf_counter();body_cpu=time.process_time()
 if a.authorize_postscore!='ROOT_APPROVED_POSTSCORE_PLOTS':raise PermissionError('Root postscore plot authorization required; no GT opened')
 receipt=read(a.scores/'RECEIPT.json');gate=read(a.scores/'PUBLIC_GATE_BEFORE_TRUTH.json')
 if receipt.get('status')!='SCORED_POST_ALL288_TERMINAL' or gate.get('status')!='PASS_ALL288_TERMINAL' or receipt.get('tasks')!=288:raise RuntimeError('Completed root-authorized scoring required before plot GT access')
 rows=loadcsv(a.scores/'SCORES_AND_FIXED_BUDGETS.csv');final=[r for r in rows if r['view']=='FINAL_CAP_ELIGIBLE'];tasks=[active_task(t) for t in read(a.manifest)['tasks']];taskmap={(t['case_id'],method(t)):t for t in tasks};cases=sorted({r['case_id'] for r in final})
 if len(final)!=288 or len(cases)!=48 or len(taskmap)!=288:raise ValueError('All48 objects and all six methods required; no success filtering')
 for case in cases:
  if sorted(r['method'] for r in final if r['case_id']==case)!=sorted(METHODS):raise ValueError('Unbalanced scored objects')
 config=read(a.config);a.output.mkdir(parents=True,exist_ok=False);fonts=font_policy(a.language);zh=fonts['actual']=='zh';artifacts=[];availability=[]
 dump(a.output/'PLOT_GATE_BEFORE_TRUTH.json',dict(status='PASS_POSTSCORE_ALL48_ALL6',truth_opened=False,source_score_receipt=str(a.scores/'RECEIPT.json'),fonts=fonts,no_success_filter=True))
 def save(fig,name):
  fig.canvas.draw()
  for extension in ('png','svg'):
   p=a.output/(name+'.'+extension);fig.savefig(p,dpi=200,bbox_inches='tight');artifacts.append(dict(path=p.name,bytes=p.stat().st_size))
  plt.close(fig)
 lookup={(r['case_id'],r['method']):r for r in final};families=sorted({r['family'] for r in final});family_pages=[]
 for family in families:
  family_cases=sorted({r['case_id'] for r in final if r['family']==family})
  if len(family_cases)!=6:raise ValueError('Every family page must contain exactly6 objects')
  fig,axs=plt.subplots(6,7,figsize=(16.8,14),squeeze=False);fig.subplots_adjust(left=.105,right=.94,bottom=.06,top=.95,wspace=.05,hspace=.14)
  image=None
  for i,case in enumerate(family_cases):
   reference=lookup[(case,METHODS[0])]
   # First numeric GT access for plots is after the public gate above.
   truth=npz(a.truth_dir/(case+'.npz'))['m'].reshape(24,24)
   values=[truth]
   for m in METHODS:
    path=resolve(taskmap[(case,m)]['output_dir'])/'R11_CAP_PREDICTION.npz'
    values.append(npz(path)['m'].reshape(24,24) if path.exists() and lookup[(case,m)].get('available')=='True' else None)
   for j,value in enumerate(values):
    ax=axs[i,j];ax.set_xticks([]);ax.set_yticks([])
    if value is None:
     ax.set_facecolor('#eeeeee');ax.text(.5,.5,'失败 / 无状态' if zh else 'FAILED\nno paid state',ha='center',va='center',transform=ax.transAxes,color='#9c1720');ax.plot([0,1],[0,1],transform=ax.transAxes,color='#cc7777',lw=1)
    else:image=ax.imshow(value,origin='lower',vmin=0,vmax=20,cmap='viridis',interpolation='nearest')
    if j>0 and lookup[(case,METHODS[j-1])]['task_status'].startswith('FAILED'):
     for spine in ax.spines.values():spine.set_edgecolor('#b2182b');spine.set_linewidth(2)
     if value is not None:ax.text(.5,.04,'失败：已付状态' if zh else 'FAILED: paid state',ha='center',va='bottom',transform=ax.transAxes,color='#b2182b',fontsize=7,bbox=dict(facecolor='white',alpha=.9,edgecolor='none',pad=1))
    if j>0 and lookup[(case,METHODS[j-1])].get('engineering_retry')=='True':ax.text(.5,.91,'engineering retry',ha='center',va='top',transform=ax.transAxes,color='#885500',fontsize=7,bbox=dict(facecolor='white',alpha=.9,edgecolor='none',pad=1))
    if j>0 and lookup[(case,METHODS[j-1])]['task_status']=='CPU_BUDGET_REACHED':ax.text(.5,.04,'CPU cap',ha='center',va='bottom',transform=ax.transAxes,color='#444444',fontsize=7,bbox=dict(facecolor='white',alpha=.85,edgecolor='none',pad=1))
    if i==0:ax.set_title(('真值' if zh else 'GT') if j==0 else LABELS[j-1],fontsize=11)
    if j==0:ax.set_ylabel(case+'\n'+reference['regime']+' χ'+str(int(float(reference['high_chi']))),rotation=0,ha='right',va='center',labelpad=9,fontsize=9)
  cax=fig.add_axes([.955,.15,.012,.68]);fig.colorbar(image,cax=cax,label='材料对比度 χ' if zh else 'Material contrast χ = εr − 1')
  fig.suptitle(family+' — '+('全部6对象；统一色标0..20' if zh else 'all six objects; common 0..20 scale'),fontsize=15)
  name='MATERIALS_'+family;save(fig,name);family_pages.append(name+'.png')
 # Every available pair is shown, and unavailable counts remain explicit.
 fig,axs=plt.subplots(2,5,figsize=(18,9.2),squeeze=False)
 for j,m in enumerate(METHODS[1:]):
  for i,(field,label) in enumerate([('NRMSE','NRMSE'),('fixed_support_IoU','Support IoU')]):
   ax=axs[i,j];points=[]
   for case in cases:
    left=lookup[(case,METHODS[0])];right=lookup[(case,m)]
    if finite(left.get(field)) and finite(right.get(field)):
     x,y=float(left[field]),float(right[field]);points.append((x,y));ax.scatter(x,y,s=23,c='#2166ac' if left['regime']=='high' else '#d95f02',alpha=.72,edgecolors='white',linewidths=.3)
   limit=1 if field.endswith('IoU') else max([x for p in points for x in p]+[.1])*1.05
   ax.plot([0,limit],[0,limit],color='.5',lw=.8,ls='--');ax.set_xlim(0,limit);ax.set_ylim(0,limit);ax.set_aspect('equal',adjustable='box');ax.set_xlabel('CIE '+label);ax.set_ylabel(LABELS[j+1]+' '+label);ax.set_title('paired '+str(len(points))+'/48; missing '+str(48-len(points)));ax.grid(alpha=.15)
 fig.suptitle('全部对象配对；蓝=high 橙=mixed' if zh else 'All-object pairs: blue=high, orange=mixed; no success-threshold filtering',fontsize=14);fig.tight_layout(rect=[0,0,1,.95],h_pad=3.5);save(fig,'PAIRED_MATERIAL_METRICS')
 # Matched current/schedule control: isolate the Green-material contrasts.
 fig,axs=plt.subplots(2,3,figsize=(12,9.2),squeeze=False)
 for j,m in enumerate(METHODS[3:]):
  for i,(field,label) in enumerate([('NRMSE','NRMSE'),('fixed_support_IoU','Support IoU')]):
   ax=axs[i,j];points=[]
   for case in cases:
    left=lookup[(case,'global_dual')];right=lookup[(case,m)]
    if finite(left.get(field)) and finite(right.get(field)):
     x,y=float(left[field]),float(right[field]);points.append((x,y));ax.scatter(x,y,s=23,c='#2166ac' if left['regime']=='high' else '#d95f02',alpha=.72,edgecolors='white',linewidths=.3)
   limit=1 if field.endswith('IoU') else max([x for point in points for x in point]+[.1])*1.05
   ax.plot([0,limit],[0,limit],color='.5',lw=.8,ls='--');ax.set_xlim(0,limit);ax.set_ylim(0,limit);ax.set_aspect('equal',adjustable='box');ax.set_xlabel('Global Dual '+label);ax.set_ylabel(LABELS[j+3]+' '+label);ax.set_title('paired '+str(len(points))+'/48; missing '+str(48-len(points)));ax.grid(alpha=.15)
 fig.suptitle('Matched Global control versus Green material rules; all objects',fontsize=14);fig.tight_layout(rect=[0,0,1,.95],h_pad=3.5);save(fig,'MATCHED_GLOBAL_MATERIAL_METRICS')
 fig,axs=plt.subplots(1,2,figsize=(14,5.2));cost_lists=[]
 for k,m in enumerate(METHODS):
  selected=[r for r in final if r['method']==m];costs=[float(r['actual_child_CPU_seconds']) for r in selected if finite(r.get('actual_child_CPU_seconds'))];cost_lists.append(costs)
  for r in selected:
   if finite(r.get('actual_child_CPU_seconds')) and finite(r.get('NRMSE')):axs[0].scatter(float(r['actual_child_CPU_seconds']),float(r['NRMSE']),s=18,alpha=.65,c=[COLORS[k]],label=LABELS[k] if not any(x.get_label()==LABELS[k] for x in axs[0].collections) else None)
 axs[0].axvline(480,color='.4',lw=.8,ls='--',label='480 allowed CPU');axs[0].set_xlabel('Actual whole-child CPU seconds');axs[0].set_ylabel('Final cap-eligible NRMSE');axs[0].legend(fontsize=8);axs[0].grid(alpha=.15)
 boxes=axs[1].boxplot(cost_lists,patch_artist=True,showfliers=True)
 for patch,color in zip(boxes['boxes'],COLORS):patch.set_facecolor(color);patch.set_alpha(.5)
 axs[1].set_xticks(np.arange(1,7),LABELS,rotation=25,ha='right');axs[1].set_ylabel('Actual whole-child CPU seconds');axs[1].axhline(480,color='.4',lw=.8,ls='--');axs[1].grid(axis='y',alpha=.15)
 fig.suptitle('同允许CPU预算；实际消耗不同，不补空预算' if zh else 'Same allowed CPU cap; actual spent differs; no budget padding',fontsize=14);fig.tight_layout();save(fig,'ACTUAL_CPU_COST')
 fig,axs=plt.subplots(1,2,figsize=(13,5))
 for k,m in enumerate(METHODS):
  for field,ax in [('NRMSE',axs[0]),('fixed_support_IoU',axs[1])]:
   ys=[]
   for cap in config['CPU_checkpoints_seconds']:
    rr=[r for r in rows if r['method']==m and r['view']=='CPU_'+str(cap)];values=[float(r[field]) for r in rr if finite(r.get(field))];ys.append(float(np.median(values)) if values else np.nan)
    if field=='NRMSE':availability.append(dict(curve='fixed_CPU',method=m,budget=cap,available=len(values),total=48))
   ax.plot(config['CPU_checkpoints_seconds'],ys,'o-',label=LABELS[k],color=COLORS[k],lw=1.4)
 axs[0].set_ylabel('Median NRMSE (finite results)');axs[1].set_ylabel('Median support IoU (finite results)')
 for ax in axs:ax.set_xlabel('Allowed whole-child CPU cap (seconds)');ax.grid(alpha=.2)
 axs[0].legend(fontsize=8);fig.suptitle('固定CPU检查点；缺失数另表列明' if zh else 'Fixed CPU checkpoints; missingness shown in CURVE_AVAILABILITY.csv');fig.tight_layout();save(fig,'FIXED_CPU_CURVES')
 # Read scored frontiers only; choose by original Phi, never by GT metric.
 lu_grid=config['display_LU_grid'];best={};all_frontier=0
 with (a.scores/'ACTUAL_COST_FRONTIER.csv').open(newline='') as f:
  for r in csv.DictReader(f):
   all_frontier+=1
   if not finite(r.get('paid_LU')) or not finite(r.get('original_Phi')):continue
   for cap in lu_grid:
    if float(r['paid_LU'])<=cap:
     key=(r['case_id'],r['method'],cap)
     if key not in best or float(r['original_Phi'])<float(best[key]['original_Phi']):best[key]=r
 fig,axs=plt.subplots(1,2,figsize=(13,5))
 for k,m in enumerate(METHODS):
  for field,ax in [('NRMSE',axs[0]),('fixed_support_IoU',axs[1])]:
   ys=[]
   for cap in lu_grid:
    rr=[best[(case,m,cap)] for case in cases if (case,m,cap) in best];values=[float(r[field]) for r in rr if finite(r.get(field))];ys.append(float(np.median(values)) if values else np.nan)
    if field=='NRMSE':availability.append(dict(curve='fixed_physical_LU',method=m,budget=cap,available=len(values),total=48))
   ax.plot(lu_grid,ys,'o-',label=LABELS[k],color=COLORS[k],lw=1.4)
 axs[0].set_ylabel('Median NRMSE (finite results)');axs[1].set_ylabel('Median support IoU (finite results)')
 for ax in axs:ax.set_xlabel('Original physical LU allowance; auxiliary work is separate');ax.set_xscale('log');ax.grid(alpha=.2)
 axs[0].legend(fontsize=8);fig.suptitle('已付原Phi frontier重采样；并非同CPU' if zh else 'Resampled paid original-Phi frontier; equal physical LU is not equal CPU');fig.tight_layout();save(fig,'FIXED_PHYSICAL_LU_CURVES')
 with (a.output/'CURVE_AVAILABILITY.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=['curve','method','budget','available','total']);w.writeheader();w.writerows(availability)
 table=[]
 for label,m in zip(LABELS,METHODS):
  rr=[r for r in final if r['method']==m];ns=[float(r['NRMSE']) for r in rr if finite(r.get('NRMSE'))];ious=[float(r['fixed_support_IoU']) for r in rr if finite(r.get('fixed_support_IoU'))];cpu=[float(r['actual_child_CPU_seconds']) for r in rr if finite(r.get('actual_child_CPU_seconds'))]
  table.append(dict(method=label,finite=len(ns),success=sum(r.get('material_success')=='True' for r in rr),median_NRMSE=float(np.median(ns)) if ns else None,median_IoU=float(np.median(ious)) if ious else None,median_actual_CPU=float(np.median(cpu)) if cpu else None,total_actual_CPU=sum(cpu),CN99_pass=sum(r.get('CN99_pass')=='True' for r in rr),execution_failures=sum(r['task_status'].startswith('FAILED') for r in rr)))
 lines=['# Independent48: factual postscore figures','', 'All48 independent objects and all six frozen methods are included. Each family page contains six objects × seven columns (GT + six methods), with a common0..20 contrast scale. Failed or unavailable predictions keep their original cell. No object is selected by success thresholds.','', 'NRMSE is ||m−GT||₂ divided by24(max(GT)−.1), with fixed support defined by contrast>1. Material-success thresholds are inherited, while every continuous score and every object remains shown. Terminal/sealed does not mean numerical convergence; corrector stop reasons remain in the scored source table.','', 'Object0022 had controller-aborted attempts and explicit engineering retries. All48 objects remain primary; separate47-object paired sensitivity is labelled and never replaces primary. Valid-prefix childCPU and end-to-end cost including interruption/retry are distinct; exact end-to-end same-cost is indeterminate for affected comparisons. Original failure records/cost bounds remain visible.','', 'The four fixed CPU checkpoints are allowed budgets. Natural stopping consumes less; actual child expenditure comes from external GetProcessTimes and includes setup, auxiliary work, corrector, failures and finalization. The physical LU display grid resamples already-paid original-Phi frontiers; it does not equalize auxiliary work or CPU. Missingness is tabulated separately.','', '|Method|Finite/48|Material success/48|Median NRMSE|Median IoU|Median actual CPU s|Total actual CPU s|CN99 pass|Execution failures|','|---|---:|---:|---:|---:|---:|---:|---:|---:|']
 for r in table:lines.append('|'+str(r['method'])+'|'+str(r['finite'])+'|'+str(r['success'])+'|'+str(r['median_NRMSE'])+'|'+str(r['median_IoU'])+'|'+str(r['median_actual_CPU'])+'|'+str(r['total_actual_CPU'])+'|'+str(r['CN99_pass'])+'|'+str(r['execution_failures'])+'|')
 paired=[r for r in loadcsv(a.scores/'PAIRED_STRATIFIED_BOOTSTRAP.csv') if r['view']=='FINAL_CAP_ELIGIBLE' and r['left']=='cie_single']
 lines+=['','Paired stratified bootstrap, exactly as frozen/scored. Sign: CIE minus the comparison method; objects, not pixels, are the sampling unit.','','|Comparison|Finite pairs|Mean ΔNRMSE [95%]|Median ΔNRMSE [95%]|CIE-only success|Other-only success|','|---|---:|---|---|---:|---:|']
 for r in paired:
  lines.append('|'+r['right']+'|'+r['paired_finite']+'|'+r.get('mean_NRMSE_difference','')+' ['+r.get('mean95low','')+', '+r.get('mean95high','')+']|'+r.get('median_NRMSE_difference','')+' ['+r.get('median95low','')+', '+r.get('median95high','')+']|'+r['left_only_success']+'|'+r['right_only_success']+'|')
 global_pairs=[r for r in loadcsv(a.scores/'PAIRED_STRATIFIED_BOOTSTRAP.csv') if r['view']=='FINAL_CAP_ELIGIBLE' and r['left']=='global_dual' and r['right'] in METHODS[3:]]
 lines+=['','Matched-control contrasts must be read separately. CIE/Dual versus Global can explain changes already produced by the current/schedule control. Global versus each Green variant is the matched material-rule contrast; a benefit already present in Global is not Green-specific. No benefit is asserted by this template.','','|Global minus method|Finite pairs|Mean ΔNRMSE [95%]|Median ΔNRMSE [95%]|','|---|---:|---|---|']
 for r in global_pairs:lines.append('|'+r['right']+'|'+r['paired_finite']+'|'+r.get('mean_NRMSE_difference','')+' ['+r.get('mean95low','')+', '+r.get('mean95high','')+']|'+r.get('median_NRMSE_difference','')+' ['+r.get('median95low','')+', '+r.get('median95high','')+']|')
 lines+=['','These are factual summaries. Scientific interpretation, novelty judgment and final claims remain with the root reviewer.','']
 for image in family_pages+['PAIRED_MATERIAL_METRICS.png','MATCHED_GLOBAL_MATERIAL_METRICS.png','ACTUAL_CPU_COST.png','FIXED_CPU_CURVES.png','FIXED_PHYSICAL_LU_CURVES.png']:lines+=['!['+image+']('+image+')','']
 (a.output/'PLOT_REPORT.md').write_text('\n'.join(lines)+'\n')
 html_parts=['<!doctype html><meta charset="utf-8"><title>Independent48 postscore figures</title><style>body{max-width:1400px;margin:32px auto;padding:0 20px;font-family:system-ui,sans-serif;color:#222}img{max-width:100%;height:auto}table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:7px}p{line-height:1.6}</style><h1>Independent48: all objects and six methods</h1><p>Same allowed CPU cap; actual expenditure differs. No success-threshold filtering. All scientific interpretation remains with the root reviewer.</p>']
 for image in family_pages+['PAIRED_MATERIAL_METRICS.png','MATCHED_GLOBAL_MATERIAL_METRICS.png','ACTUAL_CPU_COST.png','FIXED_CPU_CURVES.png','FIXED_PHYSICAL_LU_CURVES.png']:html_parts.append('<figure><img src="'+html.escape(image)+'" alt="'+html.escape(image)+'"><figcaption>'+html.escape(image)+'</figcaption></figure>')
 (a.output/'PLOT_REPORT.html').write_text('\n'.join(html_parts))
 dump(a.output/'RECEIPT.json',dict(status='RENDERED_POSTSCORE_ALL_OBJECTS',material_pages=8,objects=48,methods=6,all_objects_included=True,success_filtered=False,physical_calls=0,GT_scoring='NO_NEW_SCORING; display of already-authorized GT only',frontier_rows_read=all_frontier,fonts=fonts,artifacts=artifacts,summary=table,CPU_seconds=time.process_time(),module_entry_CPU_seconds=ENTRY_CPU,body_entry_CPU_seconds=body_cpu,body_elapsed_CPU_seconds=time.process_time()-body_cpu,CPU_scope='Absolute process_time includes imports, reads, rendering and reporting through receipt preparation',wall_seconds=time.perf_counter()-wall,visual_review='ROOT_REVIEW_REQUIRED_AFTER_RENDER'))
 print(json.dumps(dict(status='RENDERED_POSTSCORE_ALL_OBJECTS',output=str(a.output),material_pages=8)))
if __name__=='__main__':main()
