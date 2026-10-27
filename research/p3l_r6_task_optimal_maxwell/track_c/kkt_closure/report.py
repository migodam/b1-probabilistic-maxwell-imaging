"""Post-seal diagnostics and report. Zero physics/fit; CDEV development only."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
import json,csv,hashlib,datetime
from pathlib import Path
import numpy as np
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent;OUT=HERE.parents[1];TRACK=HERE.parent;REPO=OUT.parents[1];R5=REPO.parent/'b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk';R4=REPO.parent/'b1-p3l-r4-conditional-transport/research/p3l_r4_conditional_transport';torch.set_num_threads(1)
def read(p):return json.loads(Path(p).read_text())
def write(name,v):(HERE/name).write_text(json.dumps(v,indent=2,default=lambda v:v.item() if isinstance(v,np.generic) else v.tolist())+'\n')
def csvwrite(name,rows):
 with (HERE/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def ident(path):
 st=path.stat();return dict(path=str(path),bytes=st.st_size,mtime_ns=st.st_mtime_ns,blake2b512=hashlib.blake2b(path.read_bytes()).hexdigest())

def main():
 d=read(HERE/'RESULTS.json');diag=read(HERE/'FEATURE_DIAGNOSTICS.json');seal=read(HERE/'ALL_PREDICTIONS_SEALED.json');job=read(HERE/'JOB.json');parents=seal['parents'];arr=np.load(HERE/'PARENT_LOSSES.npz');configs=arr['configs'].tolist();modes=arr['modes'].tolist();risk=arr['meanseed_MSE'];baseline=arr['baseline_MSE'];orig={}
 for name in ['evaluation_CDEV','evaluation_CDEV_revision1']:
  with np.load(TRACK/name/'PARENT_LOSSES.npz') as z:
   mm=z['methods'].tolist().index('deterministic_matched')
   for i,c in enumerate(z['configs'].tolist()):orig[c]=z['MSE_meanseed'][:,i,mm]
 # Identical independent-parent resamples; same-config changes are exploratory.
 indices=np.random.default_rng(202610270608).integers(0,64,(20000,64));same=[]
 for i,config in enumerate(configs):
  q=risk[:,i*3];delta=q-orig[config];b=delta[indices].mean(1);same.append(dict(config=config,original_same_config_MSE=float(orig[config].mean()),q_MSE=float(q.mean()),q_minus_original_sameconfig_MSE=float(delta.mean()),paired_CI98p333_lower=float(np.quantile(b,.05/6)),paired_CI98p333_upper=float(np.quantile(b,1-.05/6)),interpretation='alreadyexposedCDEV exploratorydifference, notindependentKKT-specificgain'))
 csvwrite('SAME_CONFIG_ACTUAL_GRADIENT_COMPARISON.csv',same);write('SAME_CONFIG_COMPARISON.json',same);csvwrite('ALL_ARMS.csv',d['all_arms'])
 rr=[dict(parent=parents[r['parent_index']],**r) for r in diag['records']];csvwrite('FEATURE_KKT_DIAGNOSTICS.csv',rr);cfgsum=[]
 for c in configs:
  vv=[v for v in rr if v['config']==c];cfgsum.append(dict(config=c,median_KKT_error=float(np.median([v['KKT_g_minus_z_norm'] for v in vv])),p95_KKT_error=float(np.quantile([v['KKT_g_minus_z_norm'] for v in vv],.95)),median_cosine_g_z=float(np.median([v['cosine_g_z'] for v in vv])),median_normalized_back_difference=float(np.median([v['normalized_back_difference_norm'] for v in vv])),p95_normalized_back_difference=float(np.quantile([v['normalized_back_difference_norm'] for v in vv],.95)),max_raw_difference_to_bound_ratio=max(v['raw_back_difference_norm']/max(v['Lipschitz_bound'],1e-30) for v in vv)))
 csvwrite('CONFIG_FEATURE_SUMMARY.csv',cfgsum)
 identities=[]
 provider=R4/'models/density/extended_diag_gaussian_s2026102501/BEST.pt'
 for seed in [2026102601,2026102602,2026102603]:
  path=R5/f'experiments/deterministic_matched_s{seed}/BEST.pt';p=torch.load(path,map_location='cpu',weights_only=False);identities.append(dict(seed=seed,matched_head=ident(path),head_selected_update=p.get('selected_update'),head_parameter_count=sum(v.numel() for k,v in p['state_dict'].items() if k.endswith('weight') or k.endswith('bias')),original_feature_sigma=p['geometry']['config']['noise']['sigma'],j_is_whitened=p['geometry']['j_is_whitened'],feature_normalization_payload='checkpoint feature_payload; preserved unchanged'))
 write('SELECTED_WEIGHTS_IDENTITY.json',dict(created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),method='deterministic_matched',selected_development_config='native24_probit64',closure='q',seeds=[2026102601,2026102602,2026102603],provider=ident(provider),heads=identities,decoder_asset_inventory=ident(R5/'ASSET_INVENTORY.json'),weights_modified=False,new_training=False,checkpoint_loading='historical model_actions.Action→historical R5 inference.Action; no added trainable module'))
 write('CONTEXT_SOURCE.json',dict(selected_anchor='native24_probit64',selected_KKTclosure='q',baseline_reference='existing originaltruegradient native96_probit64 deterministic_matched',controls=['full native96_probit64 qclosure','ordinary native24_probit64 truegradient deterministic_matched','rawzero backclosure','negativeq closure'],source_context_registry=str(HERE/'SOURCE_CONTEXT_REGISTRY.json'),anchor_coordinates='nativeGaussianfield64 plusnativeuniformproperties2; chi0originaldecoderfield',API=dict(file=str(HERE/'closure_action.py'),class_name='Action',constructor="Action(method='deterministic_matched',seed=2026102601,closure_mode='q',noise_ratio=r)",predict='predict(yfloat32[N,512],online,parents=None,batch_size=64)',online_required=['q0[N,66]','chi0[N,576]'],unused_allowed=['J','residual','observation'],future_noise_ratio='r=sigma_new/sigma_featurecfg; rawg=r*z then rawg/max(norm(rawg),1)',sigma_featurecfg=.00119211798514519,supported_methods=['deterministic_matched','distilled_w128_n32768'],actual_CDEV_evaluated_methods=['deterministic_matched']),D_evaluation_not_run_by_this_worker=True,feature_closure_is_new_approximation=True,physicalJ_modified=False,physicalJ_identity_claim=False,new_physics_RHS=0))
 # Risk and calibration-control comparisons; pointmeans and fixedNI upper.
 fig,axs=plt.subplots(1,2,figsize=(12,4.5));x=np.arange(len(configs));labels=['96/64','24/16','48/32','72/48','Adapt72/48','24/64','48/64','96/32']
 for mode,mark in [('q','o'),('zero','s'),('negativeq','^')]:
  r=[v for v in d['all_arms'] if v['closure_mode']==mode];axs[0].plot(x,[v['MSE_ratio'] for v in r],mark+'-',label=mode);axs[1].plot(x,[v['MSE_ratio_upper98p333'] for v in r],mark+'-',label=mode)
 for ax,y in zip(axs,['Mean seed MSE / original fullDet','One-sided 98.333% upper risk ratio']):ax.axhline(1.01,color='black',ls=':',label='Unchanged NI1% gate');ax.set_xticks(x,labels,fontsize=8,rotation=25);ax.set_xlabel('Native / Probit allocation');ax.set_ylabel(y);ax.grid(alpha=.2);ax.legend(fontsize=8)
 fig.suptitle('KKT positive-sign closure versus zero/negative controls: exposed CDEV64 only');fig.tight_layout();fig.savefig(HERE/'KKT_CONTROL_RISK.png',dpi=180);plt.close(fig)
 fig,axs=plt.subplots(1,2,figsize=(11,4));axs[0].semilogy(x,[v['median_KKT_error'] for v in cfgsum],'o-',label='Median ||g-z||');axs[0].semilogy(x,[v['p95_KKT_error'] for v in cfgsum],'s--',label='p95 ||g-z||');axs[1].semilogy(x,[v['median_normalized_back_difference'] for v in cfgsum],'o-',label='Median feature difference');axs[1].semilogy(x,[v['p95_normalized_back_difference'] for v in cfgsum],'s--',label='p95 feature difference')
 for ax in axs:ax.set_xticks(x,labels,rotation=25,fontsize=8);ax.set_xlabel('Native / Probit allocation');ax.legend(fontsize=8);ax.grid(alpha=.2)
 axs[0].set_ylabel('Cached field stationarity residual');axs[1].set_ylabel('Normalized backchannel difference L2');fig.suptitle('Feature-space diagnostics only; no material-risk certificate');fig.tight_layout();fig.savefig(HERE/'KKT_FEATURE_DIAGNOSTICS.png',dpi=180);plt.close(fig)
 write('FIGURE_SOURCES.json',{'KKT_CONTROL_RISK.png':['ALL_ARMS.csv','RESULTS.json'],'KKT_FEATURE_DIAGNOSTICS.png':['FEATURE_KKT_DIAGNOSTICS.csv','CONFIG_FEATURE_SUMMARY.csv']})
 text=['# KKT closure：冻结功能模型的零新增物理特征算法','', '**开发结果：正号q closure在24/64与48/64配置通过原NI1%/wall20%数值门；zero-back和negativeq均未通过。** 这是已经暴露的CDEV64上的新近似算法选择，不是独立确认；单凭跨门前后不能称为KKT专属风险改善。','', '只替换六通道image中的复残差反投影方向。残差约定为y−F，因此无约束场Gaussian stationarity给出Jwhiteᵀrwhite=z，正号为+q0[:64]。新方向使用g/max(||g||,1)，保留原anchor、df/δ/A/weak/strong灵敏度、u2标准化、冻结Gaussian provider和三份匹配head。推断直接构造解析image，没有制造proxy J，没有修改物理J，也没有训练。','', '## 已冻结API与噪声约定','', '`closure_action.py: Action(method="deterministic_matched", seed=2026102601, closure_mode="q", noise_ratio=r).predict(yfloat32, online, parents, batch_size=64)`。online只需要原生q0和chi0；额外J/residual/observation被忽略。class实现也支持原蒸馏provider，但本次实际sealed pilot只执行deterministic_matched。','', '未来噪声比r=σnew/σfeaturecfg已在预测前声明：保存的Jwhite按新σ白化，而历史feature residual仍按旧σ白化，故raw closure用r*z，再按原式归一化。σfeaturecfg=.00119211798514519。只把它称为stationarity近似；早停、非收敛或其他偏离stationarity的对象仍可能不满足g=z。','', 'CONFIG_FROZEN.json先于任何预测保存；原5+revision3共8个context来源见SOURCE_CONTEXT_REGISTRY.json。72份预测（8配置×3模式×3seeds）全部BLAKE2b封存、重新验证后才打开本次CDEV truth。原fullDet预测直接使用原有seal，没有重算。无旧TEST/D truth访问。','', '## 固定门与实际结果','', '|配置/模式|MSE|ratio upper98.333%|wall saving lower98.333%|联合门|','|---|---:|---:|---:|---|']
 for v in d['all_arms']:
  if v['config'] in ['native24_probit64','native48_probit64','native96_probit64']:text.append(f"|{v['config']} / {v['closure_mode']}|{v['mean_seed_MSE']:.9f}|{v['MSE_ratio_upper98p333']:.6f}|{v['wall_saving_lower98p333']:.4f}|{'PASS' if v['development_joint_gate'] else 'FAIL'}|")
 text+=['', '基准原始fullDet MSE=.006353832110。固定20,000次独立parent重采样，先平均3seed loss，不是ensemble图像loss；one-sided98.333% upper≤1.01、wall saving lower≥.20未调整。whole wall使用真实anchor receipt加实测单parent推断；没有因closure不使用J而虚构扣除已经付费的J费用。计时仍受当前机器负载影响，需要后续独立/重复成本确认。','', '## 不把跨门当成KKT风险增益','', '24/64的原true-gradient Det上界1.015101，新q上界1.008195；但48/64原true-gradient本已通过，且q的MSE点估计比它略差。相同配置的配对差与98.333%双侧探索区间：']
 for v in same:
  if v['config'] in ['native24_probit64','native48_probit64','native96_probit64']:text.append(f"- {v['config']}: q−actual = {v['q_minus_original_sameconfig_MSE']:.9g}; [{v['paired_CI98p333_lower']:.9g}, {v['paired_CI98p333_upper']:.9g}]。")
 text+=['', '以上仅用于开发诊断；在新低成本q与普通低成本actual之间取得独立配对风险差证据前，不声称KKT专属风险改善。短probit配置的q仍未通过NI，不能说closure修复所有早停误差。','', '## 解析公式与特征界诊断','', '24个稀疏CPU重构检查将真实缓存gradient代回解析image，均满足预声明allclose(rtol=2e-5,atol=2e-6)，最大normalized image差2.8610e-6。这只是当前CPU公式实现资格；它与历史CPU/MPS预测比较的3e-6门是不同检查，没有覆盖或改变历史FAIL。','', '512个context的||g−z||、cosine、raw/normalized back差和decoder operator norm均保存。实际raw back差不超过测得的||δdiag(df)A||·||g−z||，但这些有限点数值界只在feature空间，不是风险证书。full配置的median ||g−z||仍为.006626，说明原全预算输出也不是每个parent都精确stationary。','', f'实际成功job CPU {job["CPU_seconds"]:.3f}s、wall {job["wall_seconds"]:.3f}s，CPU1，peakRSS {job["peak_RSS_bytes"]/1e9:.3f}GB；新增RHS=0、训练=0。一次启动import路径错误在任何预测前失败，startup_import_failure.log保留，随后修复依赖路径。','', 'SELECTED_WEIGHTS_IDENTITY.json记录原head/provider文件身份、selected update、feature噪声与未修改状态；CONTEXT_SOURCE.json可供独立QA/D scaffold加载。所有新文件仅在本kkt_closure目录。','', '![控制风险](KKT_CONTROL_RISK.png)','','![KKT特征诊断](KKT_FEATURE_DIAGNOSTICS.png)']
 (HERE/'REPORT_ZH.md').write_text('\n'.join(text)+'\n')
 write('PLAIN_CLAIMS.json',dict(development_positiveq_gate_pass=['native24_probit64','native48_probit64'],zero_or_negative_gate_pass=[],independent_KKT_specific_risk_gain='NOT_ESTABLISHED',lowcost24P64_original_actual_upper=1.0151012997534339,lowcost24P64_q_upper=1.0081950966260622,no_risk_certificate_from_feature_bound=True,no_fits=True,physics_RHS=0))
if __name__=='__main__':main()
