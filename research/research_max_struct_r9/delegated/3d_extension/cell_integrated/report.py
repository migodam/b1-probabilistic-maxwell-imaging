"""Postprocess the sole fixed IGT control; no solves or inverse GT reads."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import time
CPU,WALL=time.process_time(),time.perf_counter()
from pathlib import Path
import json,datetime
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'RESULTS/vector_extension/cell_integrated_mie16';D=Path(__file__).resolve().parent
r=json.loads((OUT/'RECEIPT.json').read_text());checks=r['quadrature_checks']
with np.load(OUT/'fields.npz',allow_pickle=False) as z:
    theta=z['theta'];integrated=z['far'];reference=z['Mie_far'];old=z['retained_old_CMRR_far']
etheta=np.c_[np.cos(theta),theta*0,-np.sin(theta)]
fig,axes=plt.subplots(1,2,figsize=(11,4.2))
curves=[]
for name,field,color,style in [('analytic Mie',reference,'black','--'),('retained CM+RR16',old,'tab:blue','-'),('cell-integrated16',integrated,'tab:orange','-')]:
    amplitude=np.sum(field*etheta,axis=1);curves.append((name,amplitude));axes[0].plot(theta*180/np.pi,amplitude.real,label=name,color=color,ls=style);axes[1].plot(theta*180/np.pi,amplitude.imag,label=name,color=color,ls=style)
for ax,title in zip(axes,['Real far-theta amplitude','Imaginary far-theta amplitude']):
    ax.set_title(title);ax.set_xlabel('Scattering angle (degrees)');ax.set_ylabel('Far amplitude / lambda0');ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('epsilon9 sphere: SAME MG16 geometry, unchanged physical units\nCM+RR error43.057% / integrated error43.393%; no gain or frequency fit')
fig.tight_layout();figure=ROOT/'FIGURES/vector_extension/cell_integrated_mie16_angles.png';fig.savefig(figure,dpi=160);plt.close(fig)
with (OUT/'angular_components.csv').open('w') as file:
    file.write('theta_deg,analytic_real,analytic_imag,CMRR_real,CMRR_imag,IGT_real,IGT_imag\n')
    for i,angle in enumerate(theta):file.write(','.join(map(str,[angle*180/np.pi,*[v for _,a in curves for v in (a[i].real,a[i].imag)]]))+'\n')
mg=[]
for n,path in [(8,ROOT/'RESULTS/vector_extension/validation/grid8.npz'),(10,ROOT/'RESULTS/vector_extension/validation/grid10.npz'),(16,ROOT/'RESULTS/vector_extension/validation_refinement16/grid16.npz')]:
    with np.load(path,allow_pickle=False) as z:f=z['fraction'];x=z['chi']
    mg.append(dict(nside=n,boundary_cells=int(np.sum((f>0)&(f<1))),positive_cells=int(np.sum(f>0)),MG_sumchi_ratio_to_linear_volume=float(x.sum()/(8*f.sum())),
        MG_static_polarizability_sum_ratio=float(np.sum(x/(3+x))/np.sum(f*8/11))))
(OUT/'MG_DIAGNOSTIC.json').write_text(json.dumps(dict(rows=mg,scope='Descriptive totals of preexisting sphere fractions; not proof of a physical mean-epsilon error or uniquecause ofMieerror'),indent=2)+'\n')
lines=['# 单次体积分Green控制：Mie16阴性结果','',
    '**固定cell-integrated pulse/collocation控制未修复这个MG边界球的Mie大误差。** 新误差43.3934%，保留的CM+RR16为43.0572%；截面误差分别30.1015%/30.9126%。没有重跑旧CM+RR、改变频率/材料/半径、增加网格或运行新反演。旧29结果及旧Mie负证据保持不变。','',
    '两者同为εr9/χ8、r=.145λ、h=.025λ，完整复矢量散射；复用完全相同的16³/1136 active-cell MG24 fractions、中心入射与解析参考。新算子使用root预先指定的标准脉冲基/中心配点体积分：近场Chebyshev lag≤2采用Gauss8，其他非自项Gauss3；自项Duffy Gauss16并包含全部动态实/虚部。新远场严格乘立方体均匀电流的sinc乘积，不用拟合归一化。','',
    '|量|保留CM+RR16|新体积分16|','|---|---:|---:|',
    f'|复远场相对L2误差|{100*r["retained_old_CMRR_far_error"]:.4f}%|{100*r["integrated_far_complex_vector_relative_error"]:.4f}%|',
    f'|截面相对误差|{100*r["retained_old_CMRR_cross_section_error"]:.4f}%|{100*r["scattering_cross_section_relative_error"]:.4f}%|',
    f'|散射截面λ²|.1329469|{r["scattering_cross_section"]:.7f}|',
    f'|新线性残差|不重求|{r["linear_relative_residual"]:.3g}|','',
    '![三曲线复幅度](../../../FIGURES/vector_extension/cell_integrated_mie16_angles.png)','',
    '## 必要实现检查','',
    f'- 自项s={r["operator"]["self_density"][0]:.15g} + i{r["operator"]["self_density"][1]:.15g}；Duffy16/20绝对差{checks["self16_vs20_absolute"]:.3g}。',
    f'- 独立Cartesian regular-part Gauss24加已核准静态cube积分常数，与主自项差{checks["self16_vs_independent_absolute"]:.3g}；不是把高阶值替换主协议。',
    f'- 单个固定近lag(1,0,0)：Gauss8/12相对差{checks["near8_vs12_relative"]:.3g}；固定远lag(3,2,1)：Gauss3/5差{checks["far3_vs5_relative"]:.3g}。',
    f'- 单lag反演/轴置换对称差{checks["inversion_symmetry_relative"]:.3g}/{checks["axis_permutation_relative"]:.3g}，block转置互易差{checks["block_reciprocity_relative"]:.3g}。',
    '这些是固定必要spot checks，没有按结果调阶，也不能凭它们给出连续被动性或整个VIE的误差证书。','',
    '## 同MG控制的解释边界','',
    '此配对同时保留了原MG边界本构，**不能隔开界面本构近似与传播离散误差**。它说明本次积分核/动态self替换没有修复这个完整控制，不证明体积分理论错误，也不证明MG是唯一原因。','',
    '|Grid|切边cells/positive cells|Σχ/(8Σf)|MG静态极化率总和比|','|---|---:|---:|---:|']
for v in mg:lines.append(f'|{v["nside"]}|{v["boundary_cells"]}/{v["positive_cells"]}|{v["MG_sumchi_ratio_to_linear_volume"]:.6f}|{v["MG_static_polarizability_sum_ratio"]:.6f}|')
lines+=['','MG满足χeff=3fχ/[3+(1−f)χ]，静态局部极化率求和按构造保持；Σχ低于线性χ体积分数只是这两种边界模型的描述性差异。不能将其直接等同“真实体平均ε误差”，更不能由此反推出材料或频率来追求Mie吻合。','',
    '## 可复用接口与成本','',
    '[cell_integrated_green.py](../../../src/cell_integrated_green.py)独立于冻结旧backend。`CellIntegratedGreen(h,k).assemble_dense(integer_coordinates)`返回已积分G，无额外h³；`self_density`用于一致地构造局部H=Glocal−sI。`external_electric_operator(points,receivers)`也做source-cell Gauss3积分，其精确离散伴随为同一S.conj().T；`far_field`包含h³和cube formfactor。可供root后续OriginalModel/FFT kernel与A+自项一致接入，本轮未接入或重跑旧反演。','',
    f'本次XINAN P-core[8]、一BLAS线程、connected SSH+kill-on-close JobObject，无GPU；恰好1新LU/1源RHS。Whole-process CPU={r["wholeprocess_CPU_seconds"]:.5f}s，script wall={r["script_wall_seconds"]:.5f}s，peak working set={r["process_peak_RSS_bytes"]/1e6:.2f}MB。',
    f'LUT共{r["operator"]["unique_integer_lags"]}唯一lags、{r["operator"]["dyadic_quadrature_evaluations"]}dyadic quadrature节点；LUT CPU={r["operator"]["LUT_CPU_seconds"]:.5f}s，含LUT的dense组装CPU={r["operator"]["assemble_total_CPU_seconds"]:.5f}s（均嵌套在whole-process中，不重复相加）。',
    '新费用在本诊断receipt/COST.json单列，不回写旧29的消费账；小型传输/打包helper CPU未独立计量，未知不写为0。root拥有后续方向/物理资格判断；本worker没有自行扩展grid、check或inverse。']
(D/'REPORT.md').write_text('\n'.join(lines)+'\n')
cpu=time.process_time()-CPU
cost=dict(status='COMPLETED_NEGATIVE_PHYSICAL_DIAGNOSTIC',remote_wholeprocess_CPU_seconds=r['wholeprocess_CPU_seconds'],local_postprocess_CPU_seconds=cpu,
    measured_new_CPU_seconds=r['wholeprocess_CPU_seconds']+cpu,postprocess_wall_seconds=time.perf_counter()-WALL,new_forward_LU=1,new_source_RHS=1,new_inverse_runs=0,
    remote_script_wall_seconds=r['script_wall_seconds'],remote_peak_RSS_bytes=r['process_peak_RSS_bytes'],
    transport_helper_CPU_seconds=None,transport_scope='separate upload/extract/pack/scp helpers not CPUmetered; knownoneLU total includesallnumeric checks/setup/imports/export',
    old29_ledger_modified=False,old_CMRR_rerun=False,old_predictions_or_scores_modified=False,created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(OUT/'COST.json').write_text(json.dumps(cost,indent=2)+'\n')
print(json.dumps(cost),flush=True)
