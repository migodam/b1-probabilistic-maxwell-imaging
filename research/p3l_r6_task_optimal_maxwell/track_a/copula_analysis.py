"""ROOT finalizer: read completed copula assets, render2figures and narrative.

No solver/training/GT/D-truth/source-model import. Matplotlib is imported ONLY
inside render_figures(), so importing this module is harmless to live timing.
Run after timing: python track_a/copula_analysis.py
Use --skip-figures for report/CSV-only work. Old raw/protocol/ledger stay read-only.
"""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
import argparse,json,csv,time
from pathlib import Path
import numpy as np
TRACK=Path(__file__).resolve().parent;OUT=TRACK/'copula_analysis';RAW=TRACK/'raw/copula_counterexample'
def read(p):return json.loads(Path(p).read_text())
def csvwrite(path,rows):
 with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def prepare():
 result=read(TRACK/'COPULA_COUNTEREXAMPLE_RESULTS.json');receipt=read(TRACK.parent/'raw/jobs/fixedmaterial_fibre_copula_counterexample/CLOSED.json');qa=read(TRACK/'COPULA_FIBRE_REFINEMENT_QA.json');OUT.mkdir(exist_ok=True);raw={};summary=[];points=[]
 if result['actual_new_RHS']!=448 or result['attempted_forward_calls']!=56:raise ValueError('actual result changed; revisit reportnumbers ratherthan silently reuse old statements')
 for row in result['rows']:
  case=row['case']
  with np.load(RAW/f'{case}.npz') as z:raw[case]={k:z[k] for k in z.files}
  z=raw[case];assert np.array_equal(z['Q_factorized'],z['Q_joint'][z['permutation']]);assert int(z['factorized_query_performed'].sum())==row['factorized_fresh_calls'];assert np.array_equal(z['prediction_difference'],z['factorized_prediction']-z['joint_reference_prediction'])
  v=z['saved_material_mode'][0]+1j*z['saved_material_mode'][1];fj=z['Q_joint']*z['D_joint_complex'];ff=z['Q_factorized']*z['D_joint_complex']
  for k in range(32):points.append(dict(case=case,index=k,Q_joint=float(z['Q_joint'][k]),Q_factorized=float(z['Q_factorized'][k]),deltaR=float(z['D_joint_complex'][k].real),deltaI=float(z['D_joint_complex'][k].imag),joint_V_real=float(fj[k].real),joint_V_imag=float(fj[k].imag),factorized_V_real=float(ff[k].real),factorized_V_imag=float(ff[k].imag),Full_complex_L2_over_sigma=float(z['complex_L2_over_sigma'][k]),Full_complex_RMS_over_sigma=float(z['complex_RMS_over_sigma'][k]),fresh_query=bool(z['factorized_query_performed'][k]),endpoint=str(z['factorized_endpoint'][k])))
  prior=row['conditional_prior_quadrature'];summary.append(dict(case=case,pairs=32,fresh_reference=1,fresh_factorized=row['factorized_fresh_calls'],exact_material_cache_hits=row['factorized_exact_material_cache_hits'],empirical_mean_shift_ROI_RMS=row['empirical_mean_shift_ROI_RMS'],empirical_material_variance_defect=row['empirical_material_variance_defect'],Full_response_complex_L2_over_sigma_mean=row['Full_response_complex_L2_over_sigma_mean'],Full_response_complex_L2_over_sigma_max=row['Full_response_complex_L2_over_sigma_max'],GL128_product_bias_real=prior['exact_product_mean_bias'][0],GL128_product_bias_imag=prior['exact_product_mean_bias'][1],GL128_product_material_variance=prior['exact_product_complex_material_variance'],one_permutation_mean_shift_real=row['empirical_mean_material_shift'][0],one_permutation_mean_shift_imag=row['empirical_mean_material_shift'][1],single_permutation_vs_empirical_marginal_product_real=row['single_permutation_mean_vs_empirical_product'][0]))
 csvwrite(OUT/'COPULA_SOURCE_SUMMARY.csv',summary);csvwrite(OUT/'COPULA_PLOT_SOURCE.csv',points)
 status=dict(status='COMPLETED_DATA_DESCRIPTIVE_COPULA_ANALYSIS',actual_new_RHS=result['actual_new_RHS'],fresh_forward_states=result['attempted_forward_calls'],exact_cache_hits=sum(r['factorized_exact_material_cache_hits'] for r in result['rows']),A_RHS_before=result['A_RHS_before'],A_RHS_after=result['A_RHS_after'],physics_job_wall_seconds=receipt['wall_seconds'],physics_job_CPU_seconds=receipt['CPU_seconds'],physics_job_peak_RSS_bytes=receipt['peak_RSS_bytes'],analysis_new_RHS=0,new_fits=0,scope='fixedmaterial conditional prior fibre; single32permutation; not iid productpopulation, posterior samples or Gamma calibration',old_protocol_raw_ledger_modified=False)
 (OUT/'ANALYSIS_SUMMARY.json').write_text(json.dumps(status,indent=2)+'\n');return result,receipt,qa,raw,summary

def render_figures(result,raw):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,axs=plt.subplots(2,2,figsize=(11,8))
 for col,row in enumerate(result['rows']):
  case=row['case'];z=raw[case];v=complex(*z['saved_material_mode']);Q=z['Q_joint'];Qf=z['Q_factorized'];D=z['D_joint_complex'];ax=axs[0,col];curve=np.linspace(min(Q.min(),Qf.min()),max(Q.max(),Qf.max()),200);ax.plot(curve,v.real/curve,color='gray',ls='--',label='Same material: deltaR=ReV/Q');ax.scatter(Q,D.real,s=32,marker='o',label='Joint fibre pairs');ax.scatter(Qf,D.real,s=38,marker='x',label='One Q permutation');ax.set_xlabel('Q');ax.set_ylabel('Native real increment deltaR');ax.set_title(case);ax.legend(fontsize=8);ax.grid(alpha=.2)
  ax=axs[1,col];Vf=Qf*D;ax.scatter(Vf.real,Vf.imag,s=30,alpha=.7,label='Permuted material increments');ax.scatter([v.real],[v.imag],s=100,marker='*',color='black',label='Fixed joint material V');emp=np.array(row['factorized_empirical_mean_increment']);pop=np.array(row['conditional_prior_quadrature']['exact_product_mean_increment']);ax.scatter([emp[0]],[emp[1]],s=70,marker='D',label='Single 32-pair permutation mean');ax.scatter([pop[0]],[pop[1]],s=80,marker='+',linewidths=2,label='GL128 exact-marginal product mean');ax.set_xlabel('Re material increment V');ax.set_ylabel('Im material increment V');ax.legend(fontsize=7);ax.grid(alpha=.2)
 fig.suptitle('Fixed-material fibre copula control: same empirical marginals, altered material coupling',fontsize=13);fig.tight_layout();fig.savefig(OUT/'COPULA_JOINT_AND_PRODUCT.png',dpi=180);plt.close(fig)
 fig,axs=plt.subplots(1,2,figsize=(11,4.2))
 for ax,row in zip(axs,result['rows']):
  case=row['case'];z=raw[case];norm=z['complex_L2_over_sigma'];order=np.argsort(norm);x=np.arange(1,33);ax.plot(x,norm[order],'o-',ms=4,label='Factorized - joint Full response');cached=~z['factorized_query_performed'][order];ax.scatter(x[cached],norm[order][cached],s=65,facecolors='none',edgecolors='black',label='Exact material cache reuse');ax.axhline(norm.mean(),ls='--',color='tab:red',label=f'32-pair mean={norm.mean():.3f}');ax.axhline(0,color='gray',lw=.8);ax.set_xlabel('Controlled pair sorted by response distance');ax.set_ylabel('Complex L2(Ffactor-Fjoint) / sigma');ax.set_title(case);ax.grid(alpha=.2);ax.legend(fontsize=8)
 fig.suptitle('Actual Full response distances; no posterior/Gamma calibration interpretation',fontsize=12);fig.tight_layout();fig.savefig(OUT/'COPULA_FULL_RESPONSE_DISTANCE.png',dpi=180);plt.close(fig)

def write_report(result,receipt,qa):
 rows=result['rows'];r,o=rows;ref=qa['comparisons'];maxbias=max(v['product_bias_real_abs_difference'] for v in ref);maxvar=max(v['variance_abs_difference'] for v in ref)
 text=f'''# Copula机制补充：精确边缘也不能代替正确联合

本次是**固定材料的条件prior-fibre copula控制**。两个V来自既有GH21 material likelihood mode，不是GT；每例按GL128 transition与pure质量取32个配对(Q,D=V/Q)，随后仅用冻结seed置乱一次Q相对于D。joint配对都解码为同一材料，置乱保留两个经验边缘但刻意改变联合关系；所有latent仍在原截断prior支持。它不是32个完整posterior样本，也不是32个独立物理世界的概率估计。条件推论由主线程定稿；这里不把标准Cauchy–Schwarz或decoder已有gauge称为新定理或Maxwell新零空间。[条件草稿](COPULA_COUNTEREXAMPLE_CONDITIONAL_PROOF_DRAFT.md)与[查询前协议](COPULA_COUNTEREXAMPLE_PROTOCOL.json)保持原版本。

真实root执行新增**448 RHS、56个forward states、0 adjoint**：两例各1个joint reference，factorized分别31与23个新状态，另有10个byte-identical decoded-material缓存复用，不能写成预算上限528已经全部发生。REGULAR单次32置乱材料均值偏移ROI RMS为{r['empirical_mean_shift_ROI_RMS']:.10g}、variance defect为{r['empirical_material_variance_defect']:.10g}，Full complex响应L2/σ均值{r['Full_response_complex_L2_over_sigma_mean']:.6g}；PURE_OVERLAP对应{o['empirical_mean_shift_ROI_RMS']:.10g}、{o['empirical_material_variance_defect']:.10g}、{o['Full_response_complex_L2_over_sigma_mean']:.6g}。这些响应差是无噪物理预测之间的控制距离，不是Gamma统计、PIT或校准通过率。[实际结果](COPULA_COUNTEREXAMPLE_RESULTS.json)、[receipt](../raw/jobs/fixedmaterial_fibre_copula_counterexample/CLOSED.json)及[图源CSV](copula_analysis/COPULA_SOURCE_SUMMARY.csv)保留完整费用和指标。

GL128条件纤维给出的精确边缘乘积实材料均值偏移是REGULAR .00263820、overlap .01802785；它与单次有限32置乱的均值不同，不能把一个经验置乱当成总体iid product distribution。独立GL64/128/256检查中，实均值偏移相对GL256的最大差约{maxbias:.3g}、product variance最大差约{maxvar:.3g}。这些是固定V上的float64数值refinement，仍只有Q1/Q2诊断，不提供global tail或Full66后验认证。[独立QA JSON](COPULA_FIBRE_REFINEMENT_QA.json)与[CSV](COPULA_FIBRE_REFINEMENT_QA.csv)由纯函数AST提取计算，没有import物理Book/Target。

A当前实际总费用为**158144 RHS**：原分层/局部参考157656、明确向量复测40、此copula补充448。旧analysis中的157656、后续posthoc时点的157696保留为各自阶段范围，不回写成仿佛当时已经发生后续查询。此补充实际wall {receipt['wall_seconds']:.6f}s、CPU {receipt['CPU_seconds']:.6f}s，旧raw、协议、主阈值与模型不变。图由下面的分析脚本在计时空闲后统一生成；source manifest同时列出NPZ、Full reference current/field cache及CSV，不新增物理或fit。

![joint纤维、一次置乱与材料增量](copula_analysis/COPULA_JOINT_AND_PRODUCT.png)

![真实Full响应控制距离](copula_analysis/COPULA_FULL_RESPONSE_DISTANCE.png)
'''
 (TRACK/'COPULA_COUNTEREXAMPLE_REPORT_ZH.md').write_text(text)
 manifest={'producer':'track_a/copula_analysis.py','render_command':'python track_a/copula_analysis.py','figures':[
  {'path':'track_a/copula_analysis/COPULA_JOINT_AND_PRODUCT.png','identity':'lowdim_local_probability','sources':['track_a/copula_analysis/COPULA_PLOT_SOURCE.csv','track_a/COPULA_COUNTEREXAMPLE_RESULTS.json','track_a/raw/copula_counterexample/REGULAR.npz','track_a/raw/copula_counterexample/PURE_OVERLAP.npz'],'scope':'controlledfixedmaterial priorfibre, notposterior sample gallery'},
  {'path':'track_a/copula_analysis/COPULA_FULL_RESPONSE_DISTANCE.png','identity':'lowdim_local_probability','sources':['track_a/copula_analysis/COPULA_PLOT_SOURCE.csv','track_a/COPULA_COUNTEREXAMPLE_RESULTS.json','track_a/raw/copula_counterexample/REGULAR_REFERENCE_CACHE.npz','track_a/raw/copula_counterexample/PURE_OVERLAP_REFERENCE_CACHE.npz'],'scope':'actualFullresponse difference, notGamma calibration'}],
  'independent_refinement':'track_a/COPULA_FIBRE_REFINEMENT_QA.json','actual_new_RHS':448,'actual_forward_states':56,'cache_hits':10,'A_total_RHS':158144,'postprocessing_new_RHS':0,'new_fits':0}
 (TRACK/'COPULA_FIGURE_SOURCES.json').write_text(json.dumps(manifest,indent=2)+'\n')

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--skip-figures',action='store_true');args=parser.parse_args();tick=time.process_time();result,receipt,qa,raw,_=prepare();write_report(result,receipt,qa)
 if not args.skip_figures:render_figures(result,raw)
 print(json.dumps(dict(status='COPULA_DESCRIPTIVE_ARTIFACTS_READY',figures_rendered=not args.skip_figures,new_physics_RHS=0,new_fits=0,CPU_seconds=time.process_time()-tick)))
if __name__=='__main__':main()
