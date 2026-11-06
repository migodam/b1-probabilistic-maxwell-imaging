"""Offline reporting only; require all requested predictions already sealed."""
import os
for k in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import json,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from research.research_max_struct_r8.experiments.paper_benchmark import ROOT,OUT,CASES,ARMS,score,write
from research.research_max_struct_r8.src.paper_physics import acquisition,make_model

def report():
    for case in CASES:
        for arm in ARMS:
            assert (OUT/case['id']/arm/'PREDICTION_SEAL.json').exists(),(case['id'],arm)
    score(None);rows=json.loads((OUT/'OFFLINE_SCORES.json').read_text());figs=ROOT/'FIGURES/paper_benchmarks';figs.mkdir(parents=True,exist_ok=True)
    with (OUT/'OFFLINE_SCORES.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    physical=[]
    for case in CASES:
        folder=OUT/case['id'];d=np.load(folder/'PUBLIC_DATA.npz');side=int(d['inverse_side']);truth=np.load(ROOT/'data/paper_benchmarks_private'/(case['id']+'.npz'))['m']
        model=make_model(acquisition(side,case['frequency_Hz']));Y=d['Y'];physical.append(dict(case=case['id'],coarse_truth_data_residual=float(np.linalg.norm(model.state(truth)['Y']-Y)/np.linalg.norm(Y))))
        fig,axes=plt.subplots(2,5,figsize=(14,6),constrained_layout=True)
        for ax,title,m in [(axes.ravel()[0],'Reference contrast',truth)]+[(ax,arm,np.load(folder/arm/'PREDICTION.npz')['m']) for ax,arm in zip(axes.ravel()[1:],ARMS)]:
            im=ax.imshow(m.reshape(side,side),origin='lower',extent=(-1,1,-1,1),vmin=0,vmax=case['epsilon_r']-1,cmap='viridis');ax.set_title(title+f' (max={m.max():.2g})',fontsize=9);ax.set_xticks([-1,0,1]);ax.set_yticks([-1,0,1])
        fig.suptitle(case['id']+' | original150 forward /32 inverse;40 receivers,20 plane waves');fig.colorbar(im,ax=list(axes.ravel()),shrink=.7,label='Contrast (above-reference values saturated)',extend='max');fig.savefig(figs/(case['id']+'.png'),dpi=160);plt.close(fig)
    write(OUT/'OFFLINE_DISCRETIZATION.json',physical)
    fig,axes=plt.subplots(1,4,figsize=(16,4),constrained_layout=True)
    colors=plt.cm.tab10(np.linspace(0,1,len(ARMS)))
    for ax,case in zip(axes,CASES):
        for arm,color in zip(ARMS,colors):
            r=next(x for x in rows if x['case']==case['id'] and x['arm']==arm);ax.scatter(r['CPU_seconds'],r['NRMSE'],color=color,label=arm,s=40)
        ax.set_title(case['id'].replace('austria_',''),fontsize=9);ax.set_xlabel('Measured CPU seconds');ax.set_ylabel('Offline relative L2 material error');ax.grid(alpha=.3)
    axes[-1].legend(loc='upper left',bbox_to_anchor=(1,1),fontsize=8);fig.savefig(figs/'material_cost.png',dpi=160);plt.close(fig)
    text=['# R8 原论文二维物理基准：实际执行与边界','',
      '本轮实际生成原2016论文的150×150细网格数据，并完成四个 Austria 物理条件 × 九个求解臂的32×32全像素反演。保持2米成像区域、20个均匀平面波、40个半径3.75米接收点、原频率/介电常数/beta/Fourier阶数与10%复白噪声。反演网格缩小是明确适配；不是把同形状放回历史0.96米/柱面源几何。','',
      '来源：[NIE2016原论文](https://doi.org/10.1109/TAP.2016.2535492)，作者全文§IV pp1794–1797。完整原文合同在 `delegated/paper_baselines/PAPER_CONTRACTS.json`。本轮重新实现的独立paper_nie/paper_lsie以原始状态方程、DPCS和native Fourier APCS/sourcewise PR交替更新；不调用CIE的readout/current loop，也不冒称成熟作者软件。','',
      '## 与原文的实际差异','',
      '- 数据为原150²细网格，主九臂比较为32²反演，附加original64运行单独报告。等面积圆盘半径h/√π近似方形MoM单元，解析disk self-term；原文未给出可逐项比对的单元积分源码。这里不是继承旧.45h占空率。',
      '- 噪声先抽取独立实/虚高斯样本，再整体缩放为noise/cleanY的Frobenius范数比精确0.1；四个公开case使用固定seed8600+i。这是明确的噪声实现合同，不假称拿到原论文的同一随机实现，也不把固定范数噪声套入未核实的IID-CN99成功证书。',
      '- CIEfamily与物理L-BFGS的实材料contrast逐像素独立在[0,20]，有共同smooth+TV prior；原论文通常恢复复材料。没有固定真值support、分区材料或轮廓参数化。',
      '- 独立paper_nie/paper_lsie保留论文无额外材料prior版本，故与带prior六臂的结果不能被解释为只有积分方程不同；其80步physical corrector也保持各自prior。',
      '- 每阶段最多300个auxiliary步骤，使用论文delta阈值：首轮.01、后续.001（代码内部delta²门限对应平方）；每臂相同最多80个L-BFGS-B corrector accepted steps。原论文没有这些统一有限上限。',
      '- sourcewise PR的delta²是逐源相对系数变化平方的平均；global PR采用继承的global Frobenius相对变化平方。二者是明确的算法/停止语义差异，不能仅因使用同一个数值阈值就称为完全相同的paper stop。',
      '- 原文以air/zeroAPIC初始化，本轮亦如此；Dual/global/A/A+的current beta固定首阶段值，material beta按论文下降。CIEsingle current/material都按论文下降。A/A+原global current与新增sourcewiseA+分别明确列出。',
      '- 每臂从阶段端点仅依原physical data+prior选候选；不读GT。预测material、Y和seal写完后才单独offline评分。这里四个条件仍是同一已公开几何家族，不能计作四种新独立拓扑或NN确认集。','',
      '## 全部离线结果','',
      '材料误差的两种分母明确区分：relativeL2=‖mhat−m‖/‖m‖；R7兼容amplitudeNRMSE=‖mhat−m‖/√N/(eps_r−1)。早期JSON的NRMSE键保留为relativeL2别名，不能直接与R7 amplitudeNRMSE合并。IoU同时报告预测自身256-bin Otsu阈值和预定half-contrast阈值；JSON另有绝对contrast>0.3控制。Otsu只是离线评分，不参与优化。','',
      '| 对象 | 方法 | relativeL2 | amplitudeNRMSE | OtsuIoU | half-contrastIoU | 数据残差 | CPU(s) | wall(s) | LU | complex RHS |','|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:text.append('| '+r['case']+' | '+r['arm']+' | '+' | '.join(f'{r[k]:.4g}' for k in ['relative_l2_material_error','amplitude_normalized_nrmse','otsu_iou','IoU_fixed_half_contrast','relative_data_residual','CPU_seconds','wall_seconds','LU','complex_RHS'])+' |')
    text+=['','A+相对原A-material的CPU改善在四个条件约1.19–2.62倍，中位约1.78倍；没有由这些结果建立A+对成熟Dual-CIE的统一材料优势。eps4的global A/A+劣于sourcewiseA+和Dual，current更新选择的影响不能归因于材料读出本身。','', '## 费用与解释','',
      '每臂CPU/wall包括该臂辅助求解、full-wave候选验收、common corrector和prediction写入前计算；不包括共享几何组装、细网格合成数据和训练。共享组装另存在每份receipt，数据生成成本在DATA_RECEIPT。FFT auxiliary与GMRES forward成本分开记录，不能只报LU为零就称没有物理费用。未有NN训练在此分支发生。初步64²forward准备也保留在delegated/forward64_setup_attempt，没有丢掉准备费用。','',
      '下表是封存后才在32²网格代入解析Austria参考材料的同数据残差，用于量化跨网格边界/单元离散影响；它没有参与在线模型或候选选择。','',
      '| 对象 | coarse truth physical residual |','|---|---:|']
    for r in physical:text.append(f"| {r['case']} | {r['coarse_truth_data_residual']:.5g} |")
    text+=['','![材料误差与费用](FIGURES/paper_benchmarks/material_cost.png)','','## 资格与未执行项','',
      '2016的其他有损涂层矩形只引用2011 FFT-TSOM论文并省略结果，原2016全文没有完整材料/边界尺寸合同，未臆造一个同名目标。2019论文给出有损coated cube和四立方体完整数值合同，本轮3D分支另行报告；不能把二维程序当成3D复现。IMSA-CIE的受限全文/作者源码未取得，因此没有用普通网格加密程序冒充IMSA。Fresnel FoamTwinDiel的实测合同与执行见独立Fresnel分支。','',
      '实际代码：`src/paper_physics.py`、`src/paper_nie.py`、`src/paper_cie_fft.py`、`experiments/paper_benchmark.py`、`experiments/paper_original64.py`。paper_cie_fft是保留R7数值算法的R8副本，仅开放相应物理Green对象与逐阶段论文stop；历史R1–R7未修改。','',
      '复现：从worktree根运行指定Python，以`python -m research.research_max_struct_r8.experiments.paper_benchmark build --forward 150`生成；`... paper_benchmark solve`执行九臂；全部seal存在后运行`... paper_report`。original64使用`... paper_original64`。一数值进程、BLAS一线程。','']
    (ROOT/'PAPER_BENCHMARKS.md').write_text('\n'.join(text))

if __name__=='__main__':report()
