"""Offline reporting only, both known-interface predictions already sealed."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import time
CPU,WALL=time.process_time(),time.perf_counter()
from pathlib import Path
import json,datetime
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'RESULTS/vector_extension/interface_mie16';HERE=Path(__file__).resolve().parent
data=json.loads((OUT/'SCORES.json').read_text());parent=json.loads((OUT/'PARENT_RECEIPT.json').read_text());workers=[]
for arm in ('scalar_arithmetic','normal_tensor'):workers.append(json.loads((OUT/arm/'RECEIPT.json').read_text()))
with np.load(ROOT/'RESULTS/vector_extension/validation_refinement16/grid16.npz',allow_pickle=False) as z:theta=z['theta'];ref=z['Mie_far']
with np.load(ROOT/'RESULTS/vector_extension/cell_integrated_mie16/fields.npz',allow_pickle=False) as z:MG=z['far']
curves=[('analytic Mie',ref,'black','--'),('MG + IGT',MG,'tab:blue','-')]
for arm,color in [('scalar_arithmetic','tab:orange'),('normal_tensor','tab:green')]:
    with np.load(OUT/arm/'PREDICTION.npz',allow_pickle=False) as z:curves.append((arm,z['far'].copy(),color,'-'))
etheta=np.c_[np.cos(theta),theta*0,-np.sin(theta)];fig,axes=plt.subplots(1,2,figsize=(12,4.5))
for name,field,color,style in curves:
    a=np.sum(field*etheta,axis=1);axes[0].plot(theta*180/np.pi,a.real,label=name,color=color,ls=style);axes[1].plot(theta*180/np.pi,a.imag,label=name,color=color,ls=style)
for ax,title in zip(axes,['Real F_theta in scattering plane','Imaginary F_theta in scattering plane']):
    ax.set_title(title);ax.set_xlabel('theta (degrees), phi=0');ax.set_ylabel('Far amplitude / lambda0');ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('OFFLINE known-interface diagnostic: same sphere16 / IGT / fractions\nErrors use three-component phi=0 theta-grid norm; total cross section is assessed separately')
fig.tight_layout();fig.savefig(ROOT/'FIGURES/vector_extension/interface_mie16_four_curves.png',dpi=160);plt.close(fig)
rows=data['scores'];lines=['# 已知球界面的两次离线前向诊断','',
    '**OFFLINE_DIAGNOSTIC，非未知材料恢复。** 同16³球、同fraction、同IGT Green/self/cube sinc，仅替换边界材料表示。无新反演，无grid/radius/frequency/material拟合；旧29项、原尺寸9臂和既有负证据均未变。两个worker预测全部封存后才一次离线Mie评分；worker未加载Mie输出。','',
    '## 评分口径','',
    '下述“平面角场误差”为预先指定φ=0、181个均匀θ采样上的三分量Frobenius相对范数：||F(φ=0,θ_l)−F_Mie(φ=0,θ_l)||F / ||F_Mie(φ=0,θ_l)||F。它不是全角球面相对场范数。截面另用24点Gauss-cosθ×32个φ的球面积分计算，两者不可替代。图仅画该散射平面的F_theta实/虚部，三分量误差仍用完整保存场。','',
    '|材料表示|φ=0平面角场误差|全角截面误差|','|---|---:|---:|','|既有MG + IGT|43.3934%|30.1015%|']
for r in rows:lines.append(f'|{r["arm"]}|{100*r["far_complex_relative_L2_error"]:.4f}%|{100*r["scattering_cross_section_relative_error"]:.4f}%|')
lines+=['','Arithmetic大幅改变结果；normal tensor的平面角场略好，但截面明显比arithmetic差。scalar截面误差4.7%不能掩盖其27.8%平面角场误差。两者角场误差仍大，不给连续精度PASS，也不证明某一种界面机制是唯一原因。','',
    '![四曲线](../../../FIGURES/vector_extension/interface_mie16_four_curves.png)','',
    '## 冻结规则与合法性','',
    'Arithmetic：ε=1+8f。Tensor：ε_t=1+8f、ε_n=(1−f+f/9)^−1，εcell=ε_t(I−nnᵀ)+ε_n nnᵀ。仅cutcell使用公开已知球的局部法向n=point/|point|；fullcell精确ε=9。此处局部平面/准静态平均依据切向E连续和法向D连续，不是新理论，不是通过Mie结果选择法向或拟合材料。','',
    '|材料表示|最小/最大ε特征值|转置对称误差Linf|fullcell ε9误差|','|---|---|---:|---:|']
for r in rows:lines.append(f'|{r["arm"]}|{r["eigenvalue_min"]:.12g} / {r["eigenvalue_max"]:.12g}|{r["reciprocal_symmetry_Linf"]:.3g}|{r["fullcell_epsilon9_max_error"]:.3g}|')
lines+=['','两个表示均实SPD、互易，使用相同584个cutcell。改变tensor还改变方向响应及有限材料本构；不能将误差变化唯一归因于体平均、方向或几何。已知法向信息也不能自动迁移给未知界面反演，更未提供给冻结9臂。','',
    '## 先例与方法边界','',
    '切向算术/法向调和平均已有直接先例，见[Farjadpour等2006](https://math.mit.edu/~stevenj/papers/FarjadpourRo06-SPIE.pdf)。该文讨论FDTD/Yee相关的各向异性subpixel处理；本实验只是将静态平面界面effective ε适配到pulse/collocation，不能继承其FDTD二阶收敛结论。',
    '[Piller/Martin的high-permittivity CDA工作](https://www.sciencedirect.com/science/article/abs/pii/S0030401898006452)结合WD与FCD；本实验没有完整复现其WD：IGT self保持冻结，也没有引入principal/secondary场的WD修正。本轮不主张公式原创、完整WD复现或未知界面恢复。','',
    '## 成本与封存','',
    '远端XINAN P-core[8]、一BLAS线程、无GPU；沿用R1 TensorDDA张量E-form I−Gχ与SciPy密集LU，χ为完整3×3块。每臂仅一次LU/一源RHS。两份PREDICTION及receipt先写入PREDICTIONS_SEAL.json，然后才离线读取Mie参考、计算角场/截面和eig合法范围。','']
for r in workers:lines.append(f'- {r["arm"]}: whole-process CPU={r["wholeprocess_CPU_seconds"]:.5f}s，wall={r["script_wall_seconds"]:.5f}s，peakRSS={r["peak_RSS_bytes"]/1e6:.2f}MB，线性残差={r["linear_relative_residual"]:.3g}。')
lines+=['','两worker成本、父进程dispatch/seal/offline assessment、Mac报告费用独立记录于COST.json，均不回写旧Mie/29 canonical。传输小helper CPU未单独计量，null不冒充0。本轮不自主扩展。']
(HERE/'REPORT.md').write_text('\n'.join(lines)+'\n')
post=time.process_time()-CPU;cost=dict(status='COMPLETED_OFFLINE_KNOWN_INTERFACE_DIAGNOSTIC',processes=[dict(arm=r['arm'],platform='XINAN_Windows_i9',CPU_seconds=r['wholeprocess_CPU_seconds'],wall_seconds=r['script_wall_seconds'],source=f'RESULTS/vector_extension/interface_mie16/{r["arm"]}/RECEIPT.json') for r in workers],
    parent_CPU_seconds=parent['parent_CPU_seconds'],parent_scope=parent['CPU_scope'],parent_wall_seconds=parent['parent_wall_seconds'],local_report_CPU_seconds=post,
    known_incremental_CPU_seconds=sum(r['wholeprocess_CPU_seconds'] for r in workers)+parent['parent_CPU_seconds']+post,
    local_report_wall_seconds=time.perf_counter()-WALL,new_LU=2,new_source_RHS=2,new_inverse_runs=0,GPU=False,
    unmetered_transport_helper_CPU_seconds=None,scope='Two disjointwholeworkers +parentonlymetadata/scoring +localreport; no nestedsetup/IGTtimer/oldMGcostreaddition',
    interpretation='OFFLINE_DIAGNOSTIC knownsphere normals, not unknowninterface inversion or FDTD convergence/WD reproduction',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(OUT/'COST.json').write_text(json.dumps(cost,indent=2)+'\n');print(json.dumps(cost),flush=True)
