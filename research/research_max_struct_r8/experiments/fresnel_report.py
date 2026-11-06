"""Aggregate actual receipts; draw derived-only scientific figures."""
from pathlib import Path
import os
os.environ['MPLBACKEND']='Agg'
import json,time,sys
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT.parents[1]));OUT=ROOT/'RESULTS/fresnel';FIG=ROOT/'FIGURES/fresnel'
ARMS=['cie_single','cie_dual','global_dual','a_material','a_plus']
LABELS={'cie_single':'CIE','cie_dual':'Dual-CIE','global_dual':'Global-Dual','a_material':'A','a_plus':'A+'}

def air_control_summary(air):
 """Posthoc distribution comparison of already paid endpoints; no physics rerun."""
 control=json.loads((OUT/'FoamTwinDielTM_multifrequency_control_8/RECEIPT.json').read_text())
 d=np.load(OUT/'FoamTwinDielTM_multifrequency_control_8/derived.npz');m=d['m'];h=.2/60
 background=np.max(np.abs(d['centers']),axis=1)>=.085
 def distribution(x):
  return dict(count=int(len(x)),mean=float(np.mean(x)),std=float(np.std(x)),quantiles=dict(zip(map(str,[0,.05,.25,.5,.75,.95,1]),map(float,np.quantile(x,[0,.05,.25,.5,.75,.95,1])))),fraction_at_declared_lower_bound=float(np.mean(np.abs(x-.02)<=1e-10)),fraction_below_old_0p02=float(np.mean(x<.02-1e-10)),fraction_near_zero=float(np.mean(x<=1e-8)))
 control['material_diagnostics'].update(whole_roi_distribution=distribution(m),outer_roi_band_distribution=distribution(m[background]),effective_area_raw_m2=float(h*h*m.sum()**2/np.sum(m*m)))
 summary=dict(status='DEV_MODEL_DIAGNOSTIC_POST_MEASUREMENT',comparison='One extra 8GHz continuation changes only lower bound .02 to 0; same paid anchor/prior/upper bound/model/300 accepted/1000 evaluations. Old control reused without rerun.',old_22_preserved=['FINAL_QA_BASE22.json','COST_AND_SELECTION_BASE22.json'],background_region='max(abs(x),abs(y))>=.085m; descriptive outer band only, not GT support or an optimization mask',common_foreground_area_definition='h2*(sum max(m-.02,0))2/sum max(m-.02,0)2 in both controls',old_floor_control=control,air0_control=air,no_pixel_GT=True,no_A_plus_advantage_claim=True)
 (OUT/'AIR0_CONTROL_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
 return summary

def main():
 cpu=time.process_time();rows=[]
 for p in sorted(OUT.glob('*/RECEIPT.json')):
  if p.parent.name.startswith('Foam') and not p.parent.name.endswith('_v3'):continue
  r=json.loads(p.read_text());r['result_directory']=str(p.parent.relative_to(ROOT))
  r['cost_scope']=('2GHz uncached A readout; no optimized-A speed claim' if r['case']!='FoamTwinDielTM' and r['arm']=='a_material' else 'Foam cached common protocol' if r['case']=='FoamTwinDielTM' else 'OOD NN regression; training separate' if r.get('initialization')=='neural' else '2GHz regression')
  rows.append(r)
 rows.sort(key=lambda r:({'FoamTwinDielTM':0,'dielTM_dec8f':1,'twodielTM_8f':2}[r['case']],r.get('initialization')=='neural',ARMS.index(r['arm'])))
 (OUT/'SUMMARY.json').write_text(json.dumps(rows,indent=2)+'\n')
 neural_differences=[]
 for case in ['dielTM_dec8f','twodielTM_8f']:
  for arm in ['cie_dual','a_plus']:
   u=[r for r in rows if r['case']==case and r['arm']==arm and r.get('initialization','uniform')=='uniform'];v=[r for r in rows if r['case']==case and r['arm']==arm and r.get('initialization')=='neural']
   if not u or not v:continue
   a=np.load(ROOT/u[0]['result_directory']/'derived.npz');b=np.load(ROOT/v[0]['result_directory']/'derived.npz');mask=a['train_mask']|a['withheld_mask']
   neural_differences.append(dict(case=case,arm=arm,estimated_material_relative_difference=float(np.linalg.norm(a['m']-b['m'])/np.linalg.norm(a['m'])),predicted_field_relative_difference_on_actual_channels=float(np.linalg.norm((a['prediction']-b['prediction'])[mask])/np.linalg.norm(a['prediction'][mask])),qualification='Differencebetweenestimates, notGTerror'))
 (OUT/'NN_ENDPOINT_DIFFERENCES.json').write_text(json.dumps(neural_differences,indent=2)+'\n')
 nie_path=OUT/'FoamTwinDielTM_8GHz_n60_independent_real_nie/RECEIPT.json';nie=json.loads(nie_path.read_text()) if nie_path.exists() else None
 if nie is not None:(OUT/'INDEPENDENT_NIE_SUMMARY.json').write_text(json.dumps(nie,indent=2)+'\n')
 mf=[]
 for mode in ['control_8','joint_6_8']:
  p=OUT/f'FoamTwinDielTM_multifrequency_{mode}'/'RECEIPT.json'
  if p.exists():mf.append(json.loads(p.read_text()))
 if mf:(OUT/'MULTIFREQUENCY_SUMMARY.json').write_text(json.dumps(mf,indent=2)+'\n')
 air_path=OUT/'FoamTwinDielTM_air0_control_8GHz_dev/RECEIPT.json';air=json.loads(air_path.read_text()) if air_path.exists() else None
 air_summary=air_control_summary(air) if air is not None else None
 corrected=[]
 for arm in ['cie_dual','a_plus']:
  p=OUT/f'FoamTwinDielTM_8GHz_n60_{arm}_dev_air0'/'RECEIPT.json'
  if p.exists():
   r=json.loads(p.read_text());r['result_directory']=str(p.parent.relative_to(ROOT));candidates=json.loads((p.parent/'auxiliary_candidates.json').read_text());r['selected_auxiliary_stage']=min(candidates,key=lambda v:v['objective'])['stage'];r['auxiliary_stage_stops']=[v for v in json.loads((p.parent/'auxiliary_history.json').read_text()) if v.get('event')=='stage_end'];corrected.append(r)
 if corrected:(OUT/'AIR0_CORRECTED_COMPARISON_SUMMARY.json').write_text(json.dumps(dict(status='DEV_MODEL_CORRECTION_POST_MEASUREMENT',original_23_preserved=['FINAL_QA_BASE23.json','COST_AND_SELECTION_BASE23.json'],runs=corrected,qualification='Freshinitial.1; onlylowerbound0; originalauxiliary/readout/corrector/prior/data; both6and8GHzholdoutpreviouslyexposed; noGTorA+victory'),indent=2)+'\n')
 failed=[]
 for p in sorted(OUT.glob('*/INTERRUPTED.json')):
  r=json.loads(p.read_text());r['result_directory']=str(p.parent.relative_to(ROOT));failed.append(r)
 for p in sorted(OUT.glob('*/FAILED.json')):
  r=json.loads(p.read_text());r['result_directory']=str(p.parent.relative_to(ROOT));failed.append(r)
 groups=[]
 for case in ['dielTM_dec8f','twodielTM_8f','FoamTwinDielTM']:
  group=[r for r in rows if r['case']==case and r.get('initialization','uniform')=='uniform'];mat=[]
  for r in group:mat.append(np.load(ROOT/r['result_directory']/'derived.npz')['m'])
  groups.append(dict(case=case,completed_arms=[r['arm'] for r in group],all_saved_m_bitwise_identical=all(np.array_equal(m,mat[0]) for m in mat) if mat else None,
   auxiliary_selected_stages={r['arm']:min(json.loads((ROOT/r['result_directory']/'auxiliary_candidates.json').read_text()),key=lambda v:v['objective'])['stage'] for r in group}))
 cost=dict(successful_inversion_CPU_seconds=sum(r['cpu_seconds'] for r in rows),successful_all_measured_tasks_CPU_seconds=sum(r['cpu_seconds'] for r in rows)+(nie['cpu_seconds'] if nie else 0.)+sum(r['cpu_seconds'] for r in mf)+(air['cpu_seconds'] if air else 0.),successful_all_measured_task_count=len(rows)+int(nie is not None)+len(mf)+int(air is not None),known_numerical_development_CPU_seconds_lower_bound=405.52,additional_first_layout_failure_CPU_seconds=None,additional_first_layout_failure_wall_seconds_lower_bound=302,development_failures=failed,reference_records=['REFERENCE_VALIDATION.json','REFERENCE_POSE_DIAGNOSTIC.json','ILLUMINATION_STABILITY_AND_REFERENCE_TRANSFER.json','FOAM_REFERENCE_VALIDATION.json'],
  independent_NIE_CPU_seconds=nie['cpu_seconds'] if nie else None,multifrequency_followup_CPU_seconds=sum(r['cpu_seconds'] for r in mf),scope='Successful total includes actual measured preprocessing/calibration/operator setup/auxiliary/corrector/prediction. In-loop checkpoint/history I/O is included; final serialization and plotting separate; interrupted costs are lower bounds or unknown explicitly, not zero. Retrieval receipts separate. Do not claim fully closed all-research CPU.',groups=groups)
 if air is not None:
  cost.update(original_22_completed_tasks_CPU_seconds=sum(r['cpu_seconds'] for r in rows)+nie['cpu_seconds']+sum(r['cpu_seconds'] for r in mf),air0_DEV_MODEL_DIAGNOSTIC_CPU_seconds=air['cpu_seconds'],air0_scope='Additional post-measurement feasible-set diagnostic; only lowerbound0; paid old control reused, not A+ advantage; both heldout frequencies previously exposed')
 if corrected:
  cost.update(original_23_completed_tasks_CPU_seconds=cost['successful_all_measured_tasks_CPU_seconds'],air0_fresh_start_DEV_MODEL_CORRECTION_CPU_seconds=sum(r['cpu_seconds'] for r in corrected),air0_fresh_start_completed_count=len(corrected),air0_fresh_start_scope='Two freshstart .1 correctedDual/A+; fulloriginal3000aux/1000corrector; onlylower0; previouslyexposedheldout; notblind/A+victory')
  cost['successful_all_measured_tasks_CPU_seconds']+=sum(r['cpu_seconds'] for r in corrected);cost['successful_all_measured_task_count']+=len(corrected)
 (OUT/'COST_AND_SELECTION_AUDIT.json').write_text(json.dumps(cost,indent=2)+'\n')
 table=['| 数据/网格 | 路径/初值 | train | 全留出 | 未见源 | 未见接收 | CPU s | cost scope | corrector accepted / stop |','|---|---|---:|---:|---:|---:|---:|---|---|']
 for r in rows:
  case=r['case'];name={'dielTM_dec8f':'2001 single2GHz/24²','twodielTM_8f':'2001 twin2GHz/24²','FoamTwinDielTM':'2005 Foam8GHz/60²'}[case];init=r.get('initialization','uniform');c=r['corrector'];stop='CONVERGED' if c['success'] else 'CAP/FAIL'
  table.append(f"| {name} | {LABELS[r['arm']]}/{init} | {r['train_relative_complex_residual']:.6f} | {r['withheld_all_relative_complex_residual']:.6f} | {r['withheld_sources_relative_complex_residual']:.6f} | {r['withheld_receivers_relative_complex_residual']:.6f} | {r['cpu_seconds']:.2f} | {r['cost_scope']} | {c['iterations']} / {stop} |")
 text='\n'.join(table)+'\n\n成功运行CPU之和为 **'+f"{cost['successful_inversion_CPU_seconds']:.2f}"+' s**，不包括上列未闭合的开发失败与检索/绘图成本。各阶段候选、真实物理修正器停止状态与逐步成本均在对应目录。\n\n'
 text+='两类2001回归与2005主实验采用不同频率/孔径/网格/停止合同，不能把跨任务残差当算法优势。NN初始化为独立合成网络的分布外试验。实测无逐像素GT；预测改善、参考几何相容性、材料准确性及算法优势必须分别判断。\n'
 if nie is not None:
  text+=f"\n独立实材料NIE方程适配另行完成：train {nie['train_relative_complex_residual']:.6f}，留出 {nie['withheld_all_relative_complex_residual']:.6f}，CPU {nie['cpu_seconds']:.2f}s；air0、box[0,20]、无prior/无corrector。该独立实现具有不同合同，非作者代码或成熟NIE复现，不混入上述五臂公平比较。\n"
 if mf:
  text+='\n多频有限跟进采用同一个8GHz Dual实材料锚点，6+8GHz等权归一化data loss加一次原prior，与仅8GHz继续优化对照。两支各最多300accepted/1000eval。6GHz是新增观测；8GHz留出已暴露，不称A+优势或相同信息下纯算法增益。\n\n| 跟进 | 6GHz train / heldout | 8GHz train / heldout | 两频平均heldout | peak contrast | CPU s | accepted / success |\n|---|---:|---:|---:|---:|---:|---|\n'
  for r in mf:
   v=r['metrics'];c=r['corrector'];text+=f"| {r['mode']} | {v['6']['train']:.6f} / {v['6']['heldout_all']:.6f} | {v['8']['train']:.6f} / {v['8']['heldout_all']:.6f} | {r['mean_6_8_heldout_relative']:.6f} | {r['material_diagnostics']['max_contrast']:.4f} | {r['cpu_seconds']:.2f} | {c['iterations']} / {c['success']} |\n"
  text+='\n跟进CPU为新增成本；共同8GHz锚点成本已计入主比较，不对两支重复计入。对照的6GHz“train”数值仅指既定训练空间mask上的事后预测，未用于对照反演。峰值/有效面积等只是材料分布统计，不是GT精度指标。各频校准、state残差、费用与停止状态单独保存。\n'
 for g in groups:
  text+=f"\n`{g['case']}` 已完成{len(g['completed_arms'])}/5个uniform路径；所有保存材料逐位相同：{g['all_saved_m_bitwise_identical']}；辅助选择阶段：`{g['auxiliary_selected_stages']}`。\n"
 if air_summary is not None:
  text+='\n仅下界改变的事后模型诊断（DEV_MODEL_DIAGNOSTIC）：同锚点/原prior/upper20/300accepted/1000eval；复用原floor.02对照，不重跑。两频留出此前已暴露，不是新测试或A+优势。外圈固定为max(|x|,|y|)≥85mm，仅描述分布，不作为GT支持或优化约束。\n\n| 8GHz物理对照 | train / heldout 8GHz | heldout 6GHz（事后） | peak | excess有效面积m² | 全ROI下界比例 | 外圈均值 / 中位数 / 95% | CPU s | accepted / success |\n|---|---:|---:|---:|---:|---:|---:|---:|---|\n'
  for label,r in [('原floor.02',air_summary['old_floor_control']),('air0',air)]:
   v=r['metrics'];md=r['material_diagnostics'];bg=md['outer_roi_band_distribution'];c=r['corrector'];text+=f"| {label} | {v['8']['train']:.6f} / {v['8']['heldout_all']:.6f} | {v['6']['heldout_all']:.6f} | {md['max_contrast']:.4f} | {md['foreground_effective_area_m2']:.8f} | {md['whole_roi_distribution']['fraction_at_declared_lower_bound']:.4f} | {bg['mean']:.5f} / {bg['quantiles']['0.5']:.5f} / {bg['quantiles']['0.95']:.5f} | {r['cpu_seconds']:.2f} | {c['iterations']} / {c['success']} |\n"
 if corrected:
  text+='\n随后仅追加两个fresh-start模型修正（DEV_MODEL_CORRECTION）：初值仍.1，原3000aux/1000corrector及prior/calibration/data，全部辅助/readout/corrector下界改0。原23条封存；两频留出已暴露，不称blind/A+胜利。Dual沿用sourcewise_cg=True，A+沿用sourcewise_cg=False（global PR），current beta固定6。\n\n| air0 fresh-start | 8GHz train / heldout | 6GHz heldout（事后） | peak | excess有效面积m² | aux steps / 选择stage | corrector accepted/eval/success | CPU s |\n|---|---:|---:|---:|---:|---:|---:|---:|\n'
  for r in corrected:
   md=r['material_diagnostics'];c=r['corrector'];counts=r['auxiliary']['operation_counts'];text+=f"| {LABELS[r['arm']]} | {r['train_relative_complex_residual']:.6f} / {r['withheld_all_relative_complex_residual']:.6f} | {r['posthoc_6GHz_metrics']['heldout_all']:.6f} | {md['max_contrast']:.4f} | {md['foreground_effective_area_m2']:.8f} | {counts['PR_CG_steps']} / {r['selected_auxiliary_stage']} | {c['iterations']} / {c['evaluations']} / {c['success']} | {r['cpu_seconds']:.2f} |\n"
 report=ROOT/'FRESNEL_REPORT.md';s=report.read_text();a=s.index('<!-- RESULTS_START -->');b=s.index('<!-- RESULTS_END -->');s=s[:a]+'<!-- RESULTS_START -->\n'+text+s[b:];report.write_text(s)
 FIG.mkdir(exist_ok=True,parents=True);fig,axs=plt.subplots(3,5,figsize=(16,10),constrained_layout=True)
 for i,case in enumerate(['dielTM_dec8f','twodielTM_8f','FoamTwinDielTM']):
  taskrows=[r for r in rows if r['case']==case and r.get('initialization','uniform')=='uniform'];rowmax=max([r['max_contrast'] for r in taskrows]+[1.])
  for j,arm in enumerate(ARMS):
   ax=axs[i,j];match=[r for r in rows if r['case']==case and r['arm']==arm and r.get('initialization','uniform')=='uniform']
   if not match:ax.set_axis_off();continue
   r=match[0];d=np.load(ROOT/r['result_directory']/'derived.npz');side=int(np.sqrt(len(d['m'])));width=100 if case=='FoamTwinDielTM' else 80
   im=ax.imshow(d['m'].reshape(side,side),origin='lower',extent=(-width,width,-width,width),vmin=0,vmax=rowmax,cmap='viridis');ax.set(title=f"{case.replace('TM_dec8f','').replace('TM_8f','')} / {LABELS[arm]}\nheld-out {r['withheld_all_relative_complex_residual']:.3f}; {r['cpu_seconds']:.0f} CPU s",xlabel='x [mm]',ylabel='y [mm]')
  fig.colorbar(im,ax=axs[i,:],label='Reconstructed real contrast',shrink=.8)
 fig.suptitle('Actual Fresnel measurements; per-task color scales shared across arms; no pixel GT')
 fig.savefig(FIG/'material_endpoints.png',dpi=160);plt.close(fig)
 fig,axs=plt.subplots(2,3,figsize=(12,8),constrained_layout=True)
 for i,case in enumerate(['dielTM_dec8f','twodielTM_8f']):
  options=[('cie_dual','uniform'),('cie_dual','neural'),('a_plus','neural')]
  for j,(arm,init) in enumerate(options):
   ax=axs[i,j];match=[r for r in rows if r['case']==case and r['arm']==arm and r.get('initialization','uniform')==init]
   if not match:ax.set_axis_off();continue
   r=match[0];d=np.load(ROOT/r['result_directory']/'derived.npz');im=ax.imshow(d['m'].reshape(24,24),origin='lower',extent=(-80,80,-80,80),vmin=0,vmax=20,cmap='viridis');ax.set(title=f"{case} / {LABELS[arm]} / {init}\nheld-out {r['withheld_all_relative_complex_residual']:.3f}; {r['cpu_seconds']:.0f} CPU s",xlabel='x [mm]',ylabel='y [mm]')
 fig.colorbar(im,ax=axs,label='Real reconstructed contrast; no pixel GT',shrink=.8);fig.suptitle('Frozen synthetic NN warm starts on actual 2 GHz measurements; out-of-distribution transfer')
 fig.savefig(FIG/'neural_warm_starts.png',dpi=160);plt.close(fig)
 from research.research_max_struct_r8.src.fresnel_adapter import load_file
 fig,axs=plt.subplots(1,3,figsize=(12,4),constrained_layout=True)
 for ax,case in zip(axs,['dielTM_dec8f','twodielTM_8f','FoamTwinDielTM']):
  match=[r for r in rows if r['case']==case and r['arm']=='cie_dual' and r.get('initialization','uniform')=='uniform']
  if not match:ax.set_axis_off();continue
  r=match[0];d=np.load(ROOT/r['result_directory']/'derived.npz');_,Y,_=load_file(ROOT/'data_local/fresnel'/f'{case}.exp',8 if case=='FoamTwinDielTM' else 2);mask=d['withheld_mask'];scale=np.linalg.norm(Y[mask])/np.sqrt(mask.sum());a=Y[mask]/scale;b=d['prediction'][mask]/scale
  ax.scatter(a.real,b.real,s=4,alpha=.25,label='Real');ax.scatter(a.imag,b.imag,s=4,alpha=.25,label='Imaginary');lim=max(abs(a).max(),abs(b).max());ax.plot([-lim,lim],[-lim,lim],'k--',lw=.8);ax.set_aspect('equal');ax.set(title=f"{case} / Dual-CIE\nheld-out relative error {r['withheld_all_relative_complex_residual']:.3f}",xlabel='Measured component / measured RMS',ylabel='Predicted component / measured RMS');ax.legend(fontsize=8)
 fig.suptitle('Actual held-out complex fields; common RMS scaling for display only')
 fig.savefig(FIG/'heldout_predictions.png',dpi=160);plt.close(fig)
 if len(mf)==2:
  anchor=np.load(OUT/'FoamTwinDielTM_8GHz_n60_cie_dual_v3/derived.npz')['m'];arrays=[anchor]+[np.load(OUT/f"FoamTwinDielTM_multifrequency_{r['mode']}"/'derived.npz')['m'] for r in mf]
  vmax=max(float(a.max()) for a in arrays);fig,axs=plt.subplots(1,3,figsize=(12,4),constrained_layout=True)
  for ax,a,title in zip(axs,arrays,['8 GHz Dual anchor']+[r['mode'] for r in mf]):
   im=ax.imshow(a.reshape(60,60),origin='lower',extent=(-100,100,-100,100),vmin=0,vmax=vmax,cmap='viridis');ax.set(title=f'{title}\npeak contrast {a.max():.3f}',xlabel='x [mm]',ylabel='y [mm]')
  fig.colorbar(im,ax=axs,label='Real reconstructed contrast; no pixel GT');fig.suptitle('Added 6 GHz observations versus 8 GHz continuation; shared material, prior and anchor')
  fig.savefig(FIG/'frequency_followup.png',dpi=160);plt.close(fig)
 if nie is not None:
  a=np.load(OUT/'FoamTwinDielTM_8GHz_n60_independent_real_nie/derived.npz')['m'];fig,ax=plt.subplots(figsize=(5,4),constrained_layout=True);im=ax.imshow(a.reshape(60,60),origin='lower',extent=(-100,100,-100,100),vmin=0,vmax=max(1.,a.max()),cmap='viridis');ax.set(title=f"Independent real NIE adaptation\nheld-out {nie['withheld_all_relative_complex_residual']:.3f}",xlabel='x [mm]',ylabel='y [mm]');fig.colorbar(im,ax=ax,label='Real contrast; no pixel GT');fig.suptitle('Air0 / no prior / no corrector; separate contract',fontsize=10);fig.savefig(FIG/'independent_nie.png',dpi=160);plt.close(fig)
 if air_summary is not None:
  arrays=[np.load(OUT/'FoamTwinDielTM_8GHz_n60_cie_dual_v3/derived.npz')['m'],np.load(OUT/'FoamTwinDielTM_multifrequency_control_8/derived.npz')['m'],np.load(air_path.parent/'derived.npz')['m']]
  vmax=max(float(a.max()) for a in arrays);fig,axs=plt.subplots(1,3,figsize=(12,4),constrained_layout=True)
  for ax,a,title in zip(axs,arrays,['Same paid 8 GHz Dual anchor','300-step control: lower .02','300-step diagnostic: lower 0']):
   im=ax.imshow(a.reshape(60,60),origin='lower',extent=(-100,100,-100,100),vmin=0,vmax=vmax,cmap='viridis');ax.set(title=f'{title}\npeak contrast {a.max():.3f}',xlabel='x [mm]',ylabel='y [mm]')
  fig.colorbar(im,ax=axs,label='Real contrast; no pixel GT');fig.suptitle('Post-measurement feasible-set diagnostic; same prior/model/anchor; no A+ claim')
  fig.savefig(FIG/'air0_control.png',dpi=160);plt.close(fig)
 if len(corrected)==2:
  arrays=[]
  for r in corrected:
   arrays.append([np.load(OUT/f"FoamTwinDielTM_8GHz_n60_{r['arm']}_v3/derived.npz")['m'],np.load(ROOT/r['result_directory']/'derived.npz')['m']])
  vmax=max(float(a.max()) for pair in arrays for a in pair);fig,axs=plt.subplots(2,2,figsize=(9,8),constrained_layout=True)
  for i,(r,pair) in enumerate(zip(corrected,arrays)):
   for j,(a,label) in enumerate(zip(pair,['Original lower .02','Corrected lower 0'])):
    ax=axs[i,j];im=ax.imshow(a.reshape(60,60),origin='lower',extent=(-100,100,-100,100),vmin=0,vmax=vmax,cmap='viridis');ax.set(title=f"{LABELS[r['arm']]} / {label}\npeak contrast {a.max():.3f}",xlabel='x [mm]',ylabel='y [mm]')
  fig.colorbar(im,ax=axs,label='Real contrast; no pixel GT');fig.suptitle('Fresh-start model correction; same original budgets/prior/calibration/data; post-measurement')
  fig.savefig(FIG/'air0_corrected_comparison.png',dpi=160);plt.close(fig)
 cost_path=OUT/'REPORT_COST.json';past=json.loads(cost_path.read_text()) if cost_path.exists() else []
 if isinstance(past,dict):past=[past]
 past.append(dict(cpu_seconds=time.process_time()-cpu,completed_main_pipelines=len(rows),multifrequency_followups=len(mf)));cost_path.write_text(json.dumps(past,indent=2))
 print(json.dumps(dict(completed=len(rows),successful_CPU_seconds=cost['successful_inversion_CPU_seconds'],groups=groups)))
if __name__=='__main__':main()
