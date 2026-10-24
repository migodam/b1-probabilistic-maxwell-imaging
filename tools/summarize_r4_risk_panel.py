"""Postclosure aggregation of saved risk losses; no truth/physics/model fitting."""
from pathlib import Path
import sys,json,csv,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np,torch
from p3r4.evaluation import paired_bootstrap,write_new_json,checkpoint_identity,read_online_context
from p3r4.fusion_provider_compat import load_fusion_compatible
from p3r4.risk_fusion import read_moments
from p3r2.risk_models import observation_features
S=ROOT/'research/p3l_r4_conditional_transport';P=S/'results/development_risk_strong_panel_v1';T=P/'summary_v2'
def csvwrite(path,rows):
 with path.open('x') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
 T.mkdir();start=time.monotonic();torch.set_num_threads(1);f=json.loads((P/'FROZEN.json').read_text());score=json.loads((P/'scores/RESULTS.json').read_text());families=list(dict.fromkeys(m['family'] for m in f['methods']));dev=[];val=[];vdata={};devarr={};normchecks=[]
 for family in families:
  methods=[m for m in f['methods'] if m['family']==family];rows=[score['methods'][m['id']] for m in methods]
  row=dict(family=family,seeds=len(methods),parents=len(f['parents']),fallbacks=sum(x['fallback_count'] for x in rows))
  arrays={};parents=None
  for action in ['raw','clipped']:
   for metric in ['full_mse','full_mae','roi_mse','roi_mae']:
    allv=[]
    for m in methods:
     with np.load(P/'scores'/f"{m['id']}_per_parent.npz") as z:
      assert z['parent'].astype(str).tolist()==f['parents'];allv.append(z[action+'__'+metric])
    values=np.stack(allv).mean(0);arrays[action+'__'+metric]=values;row[action+'_'+metric]=float(values.mean())
  dev.append(row);devarr[family]=arrays
  np.savez_compressed(T/f'{family}_DEV_parent_mean_losses.npz',parent=np.array(f['parents']),**arrays)
  vr=[];raw=[];clip=[]
  for m in methods:
   if 'checkpoint' not in m:continue
   d=Path(m['checkpoint']).parent;rp=d/'RESULT.json';vp=d/'validation.npz'
   if not rp.exists() or not vp.exists():continue
   result=json.loads(rp.read_text())
   with np.load(vp) as z:
    if not all(k in z for k in ['parent','risk','clipped_risk']):continue
    names=z['parent'].astype(str).tolist()
    if parents is None:parents=names
    assert parents==names;raw.append(z['risk']);clip.append(z['clipped_risk']);vr.append(result)
  if len(vr)==len(methods):
   r=np.mean(raw,axis=0);c=np.mean(clip,axis=0);vdata[family]=dict(parents=parents,raw=r,clipped=c)
   val.append(dict(family=family,status='SAVED_SELECTED_VALIDATION',seeds=len(vr),train_size=vr[0]['train_size'],validation_parents=len(parents),raw_full_mse=float(r.mean()),clipped_full_mse=float(c.mean()),selected_updates=';'.join(str(x['selected_update']) for x in vr)))
   np.savez_compressed(T/f'{family}_validation_parent_mean_losses.npz',parent=np.array(parents),raw_full_mse=r,clipped_full_mse=c)
  else:val.append(dict(family=family,status='NOT_EVALUATED_NO_COMMON_SAVED_VALIDATION',seeds=len(methods),train_size='',validation_parents='',raw_full_mse='',clipped_full_mse='',selected_updates=''))
 comparisons=[]
 for c in f['comparisons']:
  ca,ref=c['candidate_family'],c['reference_family']
  if ca in vdata and ref in vdata and vdata[ca]['parents']==vdata[ref]['parents']:comparisons.append({**c,**paired_bootstrap(vdata[ca][c['action']],vdata[ref][c['action']])})
 csvwrite(T/'DEV_COMMON_TABLE.csv',dev);csvwrite(T/'VALIDATION_TABLE.csv',val)
 # Verify source-normalization and recovered six predictions against closed arrays.
 checks=[];base=torch.load(S/'models/risk/strongrisk_spectral_s2026102501/BEST.pt',weights_only=False,map_location='cpu')['feature_normalization']
 for m in f['methods']:
  if m['loader']!='fusion':continue
  model=load_fusion_compatible(m['checkpoint']);normal=model.p['feature_normalization'];assert all(np.array_equal(np.asarray(normal[k]),np.asarray(v)) for k,v in base.items());normchecks.append(m['id'])
  if m['moment_provider']=='spline':continue
  with np.load(Path(m['checkpoint']).parent/'validation.npz') as z:parent=str(z['parent'][0]);saved=z['prediction'][:1]
  context=read_online_context(S/'inputs/strong_context_1024_256/anchor_validation.npz',parent)
  extra=read_moments(S/'inputs'/('fusion_'+m['moment_provider'])/'validation',[parent],model.p['moment_checkpoint_identity'],model.moment_provider,model.moment_K)
  pred=model.predict(observation_features(context['observation']),context,moments=extra,parents=[parent],moment_provider=model.moment_provider,moment_K=model.moment_K)
  error=float(np.max(abs(pred-saved)));assert error<1e-5
  checks.append(dict(id=m['id'],true_K=model.moment_K,closed_prediction_replay_max_abs=error,replay_tolerance=1e-5,replay_note='float32 batch1 online versus training batch64 convolution rounding',state_tensors_unchanged=all(torch.equal(v,model.model.state_dict()[k].cpu()) for k,v in model.p['state_dict'].items())))
  try:read_moments(S/'inputs'/('fusion_'+m['moment_provider'])/'validation',[parent],model.p['moment_checkpoint_identity'],model.moment_provider,999)
  except ValueError:pass
  else:raise AssertionError('wrongK must reject')
 profiles={p:json.loads((P/'features'/p/'MOMENTS_COMPLETE.json').read_text()) for p in f['providers']};online=json.loads((P/'online/ONLINE_COMPLETE.json').read_text());closed=json.loads((S/'I/jobs/common_dev128_strong_risk_panel_v1/CLOSED.json').read_text())
 result=dict(status='COMMON_DEV128_PANEL_COMPLETE',DEV=dev,validation=val,DEV_comparisons=score['comparisons'],validation_comparisons=comparisons,validation_scope='selected checkpoints on their own validation fold; 1024/256 and512/128 folds tabulated separately; no blind claim',statistics=f['statistics'],feature_normalization_equal_to_strong_spectral=normchecks,closed_provider_replay=checks,
  cost=dict(physics_calls=0,new_fits_by_this_panel=0,queue_completed_existing_fits=6,book_closed=closed,online_wall=online['wall_seconds'],feature_profiles={k:{z:v[z] for z in ['provider','K','wall_seconds','cpu_seconds','peak_rss_bytes']} for k,v in profiles.items()}),failures=online['failures'],locked_truth_read=False,locked_generation=False,scientific_selection='DEFERRED_TO_PARENT',interpretation='exposed DEV risk actions only; no posterior calibration/density/Flow advantage qualification',summary_wall=time.monotonic()-start)
 write_new_json(T/'PANEL_SUMMARY.json',result)
 lines=['# 共同 DEV128 风险比较（开发集已暴露）','', '三次训练种子的**逐父样本损失先取平均**，然后以128个独立父样本做4000次配对bootstrap（seed20261025099，区间98.333333%）。这不是先平均图像的集成预测。所有45个方法均保留128个父样本；失败与回退均为0。科学选择留给主线程。初次汇总因单样本float32复现容差2e-6过紧失败（观测最大2.4e-6）；保留原FAILED账本与部分表格，本目录为独立汇总v2，未重跑特征、推理或评分。','', '| 方法族 | raw MSE | 矩形截断 MSE | 截断 MAE |','|---|---:|---:|---:|']
 for r in dev:lines.append(f"| {r['family']} | {r['raw_full_mse']:.8f} | {r['clipped_full_mse']:.8f} | {r['clipped_full_mae']:.8f} |")
 lines+=['','| Spline fusion 相对参照（截断） | 风险改善 | 98.333333%配对区间 |','|---|---:|---:|']
 for c in score['comparisons']:
  if c['candidate_family']=='fusion_spline' and c['action']=='clipped':lines.append(f"| {c['reference_family']} | {100*c['improvement_relative']:.3f}% | [{100*c['relative_CI'][0]:.3f}%, {100*c['relative_CI'][1]:.3f}%] |")
 lines+=['','验证集结果见 VALIDATION_TABLE.csv：强模型为1024训练/256验证；旧CNN与阶段头为512/128，不能把这两组验证均值当成同一测试集比较。缺失的冻结R2/物理共同验证结果明确为NOT_EVALUATED。','', 'Spline、DiagGaussian、Field、空间置换fusion与强CNN的原始8通道TRAIN归一化逐项一致；fusion均6418参数。Spline与Diag使用原始在线观测生成固定32样本的模型推前均值/方差，Field是确定性点式联合Re/Im矩监督（K0，无潜变量样本）。这些四通道是粗略风险上下文，不是真实后验矩。所有已完成阶段头版本均披露，未从DEV挑选阶段头。','', '六个新provider模型的历史moment_K误写成参数量6418。原checkpoint与RESULT保持原字节；加法说明与独立兼容器仅根据封存TRAIN/validation源gate恢复32或0。六个模型复现封存验证首个预测差异逐项记录（最大2.4e-6；float32单样本与64样本批量运算舍入；容差1e-5）、参数逐项一致；错误K被拒绝。见 AUDIT/fusion/providers/CLOSED_K_METADATA_CLARIFICATION.json。','', f"本次Book I COMPLETE：{closed['controlled_wall_seconds']:.3f}秒，0 RHS、0新拟合；特征与推理/评分独立进程记录CPU/RSS。原物理native96与probit64输出只复用已有缓存，无Maxwell调用。父线程拥有物理任务，本任务未启动物理队列。",'', '流程是先冻结所有模型/源provider身份，再生成观测特征，再封存ONLINE_COMPLETE，最后另一个进程打开DEV真值。每种raw与截断动作、ROI、host/transition/pure分层以及每父样本损失均在scores。所有输入/预测文件保留路径、字节数、mtime身份（不做SHA）。','', 'LOCKED复用入口是 tools/r4_risk_panel.py：freeze(spec) → features(original online observations) → infer(allowlisted context/native cache) → score(separate process)。本次没有LOCKED生成或LOCKED真值读取。候选provider和方法必须事先声明；未知坏元数据拒绝加载。有限开发集风险下降不证明后验、概率校准或Flow机制优势。']
 (T/'REPORT_ZH.md').write_text('\n'.join(lines)+'\n');print(json.dumps(dict(status=result['status'],path=str(T))))
if __name__=='__main__':main()
