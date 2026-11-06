"""Offline original64/long-corrector/NN followups after individual seals."""
import os
for k in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from research.research_max_struct_r8.experiments.paper_benchmark import ROOT,OUT,write
from research.research_max_struct_r8.src.paper_scores import material_scores
from research.research_max_struct_r8.src.paper_physics import acquisition,austria

def run():
    rows=[]
    for group in ['original64','strong_corrector']:
        for seal in (OUT/group).glob('*/*/PREDICTION_SEAL.json'):
            dest=seal.parent;case=dest.parent.name;arm=dest.name;eps=3.5 if '3p5' in case else 4 if 'eps4_' in case else 2
            pred=np.load(dest/'PREDICTION.npz');m=pred['m'];side=int(np.sqrt(len(m)));truth=austria(acquisition(side,400e6,dense=False)['centers'],eps);r=json.loads((dest/'RECEIPT.json').read_text())
            row=dict(group=group,case=case,arm=arm,**material_scores(m,truth,eps-1),cpu_seconds=r['cpu_seconds'],wall_seconds=r['wall_seconds'],relative_data_residual=r['relative_data_residual'],LU=r['physical_cost'].get('LU',0),complex_RHS=r['physical_cost'].get('actual_complex_RHS',0),stop=r.get('stop'),projected_gradient=r.get('terminal_projected_gradient_inf'));rows.append(row);write(dest/'OFFLINE_SCORES.json',row)
    for arm in ['nn_cie_dual','nn_a_plus','aug_nn_cie_dual','aug_nn_a_plus','aug_nn_a_plus_sourcewise']:
        dest=OUT/'austria_eps4_400MHz'/arm
        if not (dest/'PREDICTION_SEAL.json').exists():continue
        p=np.load(dest/'PREDICTION.npz');truth=np.load(ROOT/'data/paper_benchmarks_private/austria_eps4_400MHz.npz')['m'];r=json.loads((dest/'RECEIPT.json').read_text());rows.append(dict(group='NN_DOMAIN_REPAIR_early80' if arm.startswith('aug_') else 'NN_OOD_early80',case='austria_eps4_400MHz',arm=arm,**material_scores(p['m'],truth,3),cpu_seconds=r['cpu_seconds'],wall_seconds=r['wall_seconds'],relative_data_residual=r['relative_data_residual'],LU=r['physical_cost'].get('LU',0),complex_RHS=r['physical_cost'].get('actual_complex_RHS',0),stop='80step physical corrector'))
    write(OUT/'FOLLOWUP_OFFLINE_SCORES.json',rows)
    initial_scores=None;aug_initial_scores=None;init=OUT/'austria_eps4_400MHz/NN_INITIALIZATION/MATERIAL.npz'
    if init.exists() and (OUT/'austria_eps4_400MHz/nn_a_plus/PREDICTION_SEAL.json').exists():
        m0=np.load(init)['m'];truth=np.load(ROOT/'data/paper_benchmarks_private/austria_eps4_400MHz.npz')['m'];rec=json.loads((OUT/'austria_eps4_400MHz/nn_cie_dual/RECEIPT.json').read_text())
        initial_scores=dict(NN=material_scores(m0,truth,3),air=material_scores(np.zeros_like(m0),truth,3),NN_initial_original_Phi=rec['initial_original_Phi'],air_initial_original_Phi=.5,shared_feature_inference=rec['feature_inference_cost_shared'])
        write(init.parent/'OFFLINE_SCORES.json',initial_scores)
    init=OUT/'austria_eps4_400MHz/NN_AUG_INITIALIZATION/MATERIAL.npz'
    if init.exists() and (OUT/'austria_eps4_400MHz/aug_nn_a_plus_sourcewise/PREDICTION_SEAL.json').exists():
        m0=np.load(init)['m'];truth=np.load(ROOT/'data/paper_benchmarks_private/austria_eps4_400MHz.npz')['m'];rec=json.loads((OUT/'austria_eps4_400MHz/aug_nn_cie_dual/RECEIPT.json').read_text())
        aug_initial_scores=dict(NN=material_scores(m0,truth,3),NN_initial_original_Phi=rec['initial_original_Phi'],shared_feature_inference=rec['feature_inference_cost_shared']);write(init.parent/'OFFLINE_SCORES.json',aug_initial_scores)
    target=OUT/'original64/austria_eps4_400MHz';arms=['cie_single','cie_dual','a_plus','a_plus_sourcewise','paper_nie']
    if all((target/a/'PREDICTION_SEAL.json').exists() for a in arms):
        fig,axes=plt.subplots(2,3,figsize=(12,8),constrained_layout=True);truth=austria(acquisition(64,400e6,dense=False)['centers'],4)
        displays=[('Reference',truth)]+[(a,np.load(target/a/'PREDICTION.npz')['m']) for a in arms]
        for ax,(label,m) in zip(axes.ravel(),displays):
            im=ax.imshow(m.reshape(64,64),origin='lower',extent=(-1,1,-1,1),vmin=0,vmax=3,cmap='viridis');ax.set_title(label+f' | max={m.max():.3g}');ax.set_xticks([-1,0,1]);ax.set_yticks([-1,0,1])
        fig.colorbar(im,ax=list(axes.ravel()),shrink=.7,extend='max',label='Contrast; values above3 saturated');fig.suptitle('Original150/64 Austria eps4 |4096 independent material pixels')
        folder=ROOT/'FIGURES/paper_benchmarks';folder.mkdir(parents=True,exist_ok=True);fig.savefig(folder/'original64_eps4.png',dpi=160);plt.close(fig)
        fullmax=max(float(m.max()) for _,m in displays);fig,axes=plt.subplots(2,3,figsize=(12,8),constrained_layout=True)
        for ax,(label,m) in zip(axes.ravel(),displays):
            im=ax.imshow(m.reshape(64,64),origin='lower',extent=(-1,1,-1,1),vmin=0,vmax=fullmax,cmap='viridis');ax.set_title(label+f' | max={m.max():.3g}');ax.set_xticks([-1,0,1]);ax.set_yticks([-1,0,1])
        fig.colorbar(im,ax=list(axes.ravel()),shrink=.7,label='Full contrast range; common scale including reference');fig.suptitle(f'Original150/64 Austria eps4 | full range0–{fullmax:g}; reference max3')
        fig.savefig(folder/'original64_eps4_full_range.png',dpi=160);plt.close(fig)
    budget=[];fig,ax=plt.subplots(figsize=(9,5),constrained_layout=True)
    for arm in ['cie_dual','a_plus','a_plus_sourcewise','nn_cie_dual','nn_a_plus','aug_nn_cie_dual','aug_nn_a_plus','aug_nn_a_plus_sourcewise']:
        d=OUT/'strong_corrector/austria_eps4_400MHz'/arm
        if not (d/'PREDICTION_SEAL.json').exists():continue
        early=OUT/'austria_eps4_400MHz'/arm;r0=json.loads((early/'RECEIPT.json').read_text());r1=json.loads((d/'RECEIPT.json').read_text());prefix=r0['cpu_seconds']+r0.get('feature_inference_cost_shared',{}).get('total_cpu_seconds',0)
        truth=np.load(ROOT/'data/paper_benchmarks_private/austria_eps4_400MHz.npz')['m'];p=np.load(d/'PREDICTION.npz');hist=json.loads((d/'HISTORY.json').read_text())['evaluations'];xy=[]
        m=np.load(early/'PREDICTION.npz')['m'];s=material_scores(m,truth,3);xy.append((prefix,s['amplitude_normalized_nrmse']))
        for ev,m in zip(p['budget_evaluations'],p['budget_material']):
            row=hist[int(ev)-1];s=material_scores(m,truth,3);cpu=prefix+row['cpu_seconds'];xy.append((cpu,s['amplitude_normalized_nrmse']));budget.append(dict(case='austria_eps4_400MHz',arm=arm,continuation_evaluation=int(ev),total_solver_CPU_seconds=cpu,**s))
        s=material_scores(p['m'],truth,3);xy.append((prefix+r1['cpu_seconds'],s['amplitude_normalized_nrmse']));xy=np.array(xy)
        ax.plot(xy[:,0],xy[:,1],label=arm,linestyle=':' if arm.startswith('aug_') else '--' if arm.startswith('nn_') else '-',linewidth=1.5)
    ax.set(xlabel='Cumulative per-query solver CPU seconds (training excluded)',ylabel='Offline amplitude-normalized material RMSE',title='Public eps4 benchmark: paid long-corrector trajectories');ax.grid(alpha=.3);ax.legend(fontsize=8,ncol=2)
    folder=ROOT/'FIGURES/paper_benchmarks';folder.mkdir(parents=True,exist_ok=True);fig.savefig(folder/'long_material_cost_curves.png',dpi=160);plt.close(fig);write(OUT/'LONG_BUDGET_OFFLINE_SCORES.json',budget)
    lines=['# R8 原论文基准后续：64网格、较长物理迭代与NN OOD','',
      '本文件是已封存预测的离线评分；只有存在PREDICTION_SEAL的项目才进入表。没有seal的正在运行/失败项目不是成功或零误差。详细原合同、早期36次比较及分母定义见PAPER_BENCHMARKS。','',
      '| 组别 | 对象 | 方法 | relativeL2 | amplitudeNRMSE | OtsuIoU | halfContrastIoU | 数据残差 | 新增CPU(s) | 新增wall(s) | 新增LU | 新增complexRHS | 停止 |','|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|']
    for r in rows:lines.append('| '+r['group']+' | '+r['case']+' | '+r['arm']+' | '+' | '.join(f'{r[k]:.5g}' for k in ['relative_l2_material_error','amplitude_normalized_nrmse','otsu_iou','IoU_fixed_half_contrast','relative_data_residual','cpu_seconds','wall_seconds','LU','complex_RHS'])+' | '+str(r.get('stop'))+' |')
    if initial_scores:
        lines+=['',f"NN初始材料amplitudeNRMSE={initial_scores['NN']['amplitude_normalized_nrmse']:.6g}，air初值={initial_scores['air']['amplitude_normalized_nrmse']:.6g}；初始原Phi NN={initial_scores['NN_initial_original_Phi']:.6g}而air=0.5。没有因NN初始Phi更高而丢弃它，所有后续物理迭代均实际完成。一次共享feature/SVD+weights装载/推理CPU约{initial_scores['shared_feature_inference']['total_cpu_seconds']:.6g}s（不含Python/PyTorch导入与训练费用）。",'']
    if any(r['arm'].startswith('aug_') for r in rows):
        lines+=['','## 轻量域扩充修复的实际对照','','主线程另行冻结25,025参数12通道网络：在原数据上增加256个随机非Austria训练材料＋64独立验证材料；混合后总1792训练/256验证。新增论文采集为400/750MHz、2米域、20plane/40receiver；没有把Austria GT用作训练或epoch选择。原10feature权重不覆盖，新增两公共通道仅为2k/(2π)电尺寸与plane-wave=1。实际运行同eps4、同300辅助上限、同80步early checkpoint，再以同一fresh L-BFGS实现追加1000步，与uniform/旧NN的80＋1000协议相同。新sourcewise A+组合没有旧10feature NN对应臂；uniform sourcewise A+追加1000作控制。此为公开开发对象上的合理轻量修复，不能称新盲确认。训练CPU约499s及旧训练成本均另加主线程训练/数据生成receipts；不能把这些成本藏在推理之外。','']
        if aug_initial_scores:
            lines += [f"修复后NN材料初值amplitudeNRMSE={aug_initial_scores['NN']['amplitude_normalized_nrmse']:.6g}，初始原Phi={aug_initial_scores['NN_initial_original_Phi']:.6g}；一次共享feature/装载/推理CPU={aug_initial_scores['shared_feature_inference']['total_cpu_seconds']:.6g}s。",'']
    lines+=['','较长corrector的费用是从各自封存early80端点继续的新增费用，必须加回early80和辅助求解成本才能比较完整使用成本。初始direct使用均匀air，不把80步预算未完成的轨迹称为成熟物理baseline终点。新增1000accepted-step ceiling、原R7 scale1e6/ftol1e-12/gtol.01，实际stop和terminal projected gradient保留。','',
      'NN只用于预先指定的eps4/32网格两个材料初始化臂，使用主线程已冻结24,737参数十通道CNN及weights；不在论文对象上训练/选epoch/调参。论文2米/400MHz、20平面波/40接收、air背景与训练.96米域/柱面波等条件不同，属于明确的跨物理OOD。仅把公共坐标特征规范为x/1米、y/1米，之后按原训练均值/方差。初始Phi不用于丢弃NN。DPCS/APCS仍按合法原算法构造，没有把NN材料当成辅助电流。特征SVD、装载/推理一次共享成本在NN_INITIALIZATION/RECEIPT，训练成本另加主训练receipt。','',
      '64网格使用实际完整4096个独立像素、同150²观测、同20plane/40receiver/3.75米以及原beta/mf。物理状态/adjoint用严格残差控制FFT-GMRES，源RHS独立Krylov求解；不声称Krylov复用一个LU。最终五臂统一采用每RHS最多1000 GMRES迭代，数值stagnation时按相同离散算子行块组装完整矩阵并采用SciPy LU，立即共享该匹配材料LU给所有source/adjoint RHS；所有失败Krylov与回退LU/setup/RHS都计费。已发生一次未带LU回退的corrector试点失败、一次短中断，以及一次昂贵通用fallback链的中断修复，独立attempt目录保留，必须计入工程准备/失败成本。最终CIEsingle从已付费auxiliary材料重新开始同80步L-BFGS，属于fresh history，不伪称连续resume。该早期auxiliary的FFT/currentloop细计数未在中断前落盘，已在RECOVERED receipt标明协议缺口；其明确auxiliary_complete进度保留CPU/wall与物理RHS，最终完整使用费用加回一次。停止于80步骤的64结果不是无预算限制物理收敛证书。','']
    if (OUT/'OFFLINE_GRID_ADEQUACY.json').exists():
        j=json.loads((OUT/'OFFLINE_GRID_ADEQUACY.json').read_text());lines+=['## 已封存eps3.5的物理网格诊断','','固定原论文解析材料，不做优化；各网格重新组装同一种equal-area disk算子。','', '| grid | noisy observed residual | clean field difference to150 | clean difference to200 |','|---|---:|---:|---:|']
        for r in j['rows']:lines.append(f"|{r['side']}|{r['observed_noisy_residual']:.6g}|{r['clean_relative_difference_to150']:.6g}|{r['clean_relative_difference_to200']:.6g}|")
        lines+=['','64网格明显减小离散场偏差，但其与150网格的约10%clean场差仍约等于噪声尺度。150→200场差约3.7%，不应把150网格当连续Maxwell真值。网格改善并不自动消除由prior、有限预算或非凸搜索造成的材料误差。上述诊断没有反馈给已封存预测的选取。','']
    lines+=['','## 64网格的幅度与形态必须同时读','','![64grid full common material amplitude range](FIGURES/paper_benchmarks/original64_eps4_full_range.png)','','上图统一从0显示到全部候选的最大contrast20，GT也使用同一尺。真材料contrast3，而反演峰值约19.4–20，幅度偏差不能被较好的轮廓IoU掩盖。下图仅用于形态细节，把0–3以外材料饱和显示，并以箭头和各面板max标记；它不是幅度恢复成功图。','','![64grid morphology detail on clipped0–3 range](FIGURES/paper_benchmarks/original64_eps4.png)','','## 实际性能—费用轨迹','','![Paid long-corrector material-cost trajectories](FIGURES/paper_benchmarks/long_material_cost_curves.png)','','曲线从已付费early80点开始，后续每10个已付费原物理evaluation复制一次在线Phi最优材料，最后加入封存终点；没有为了画曲线增加forward/adjoint。横轴加回各自辅助与early80CPU，NN再加每次输入一次feature/推理费用；训练和训练数据生成费用另计，故此图不是单次总训练摊销收益证明。材料评分只在封存后进行，没有反向选择更好的GT checkpoint。两phase重启与observed不同stop保留，不伪称同一连续L-BFGS轨迹。','','64五臂中sourcewise A+与Dual可比较相同current更新下的材料读出；此64子集没有独立global普通材料读出臂，因此global A+与Dual的差异不能单独归因于A+材料块。32九臂完整矩阵保留该对照。sourcewise/global停止的源平均与global Frobenius语义差异沿用并披露。','']
    (ROOT/'PAPER_FOLLOWUPS.md').write_text('\n'.join(lines))

if __name__=='__main__':run()
