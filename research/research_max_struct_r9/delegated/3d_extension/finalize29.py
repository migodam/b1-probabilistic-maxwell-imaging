"""Seal all29 public predictions before independent offline scoring/figures.
No physics; no outcome-dependent selection or followup.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[key]='1'
import time
WALL,CPU=time.perf_counter(),time.process_time()
from pathlib import Path
import sys,json,resource
ROOT=Path(__file__).resolve().parents[2];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
import numpy as np
from scipy import linalg as la
from research.research_max_struct_r9.experiments import vector_extension as base
CONFIG=ROOT/'configs/vector_extension_all29.json';OUT=base.OUT;FIG=ROOT/'FIGURES/vector_extension'

def score():
    manifest=base.load(CONFIG);tasks=manifest['tasks'];entries=[]
    if (OUT/'SCORES.json').exists():raise FileExistsError('Completed scores preserved')
    evidence=[base.stamp(CONFIG),base.stamp(ROOT/'src/vector_extension.py'),base.stamp(ROOT/'experiments/vector_extension.py'),
        base.stamp(Path(__file__)),base.stamp(ROOT/'delegated/3d_extension/resolved_hollow_control.py'),
        base.stamp(ROOT/'delegated/3d_extension/windows_resolved.py'),base.stamp(ROOT/'delegated/3d_extension/REMOTE_TRANSPORT.json')]
    for config in sorted(set(t['config'] for t in tasks)):evidence.append(base.stamp(ROOT/config))
    for case in dict.fromkeys(t['case_id'] for t in tasks):
        evidence.extend([base.stamp(base.PUBLIC/case/'PUBLIC.json'),base.stamp(base.PUBLIC/case/'observations.npz')])
    for task in tasks:
        folder=OUT/task['case_id']/task['arm'];receipt=base.load(folder/'RECEIPT.json')
        if receipt['case_id']!=task['case_id'] or receipt['arm']!=task['arm']:raise RuntimeError('Task identity mismatch')
        if base.load(folder/'FROZEN_CONFIG.json')!=base.load(ROOT/task['config']):raise RuntimeError('Frozen scientific contract mismatch')
        if receipt['status'] not in ('SEALED_PUBLIC_PREDICTION','FAILED'):raise RuntimeError('Unfinished29task')
        if receipt['status']=='SEALED_PUBLIC_PREDICTION' and not (folder/'prediction.npz').exists():raise RuntimeError('Missing public prediction')
        if receipt['status']=='SEALED_PUBLIC_PREDICTION':
            if receipt['actual_material_dof']!=512 or receipt['vector_dimension']!=1536 or receipt['material_box']!=[0.,8.]:raise RuntimeError('Changed physical unknown space')
            counts=receipt['counts']
            if counts['physical_source_RHS']!=60*counts['physical_LU']:raise RuntimeError('All-source retainedLU accounting mismatch')
            if task['arm']!='direct' and (counts['current_PR_steps']!=240 or receipt['current_metric']!='ordinary_CIE_only'):raise RuntimeError('Changed strict current or PR budget')
            if task['arm'] in ('a_material','a_plus','a_plus_qualified') and receipt['material_metric']!='A':raise RuntimeError('A material readout mismatch')
        entries.append(dict(**task,status=receipt['status'],files=[base.stamp(p) for p in sorted(folder.iterdir()) if p.is_file()]))
    path=OUT/'PUBLIC_PREDICTIONS_SEAL.json'
    if path.exists():raise FileExistsError('Remove/archive explicit root HOLD only after24 completion; never overwrite a true seal')
    base.dump(path,dict(status='ALL29_PUBLIC_PREDICTIONS_SEALED',created_utc=base.utc(),entries=entries,evidence=evidence,truth_access='NONE',hash_checks='NONE'))
    for item in evidence+[v for e in entries for v in e['files']]:
        if base.stamp(REPO/item['path'])!=item:raise RuntimeError('Changed public evidence')
    # First private scoring read occurs only here, after all29 terminal outputs.
    rows=[]
    for case in dict.fromkeys(t['case_id'] for t in tasks):
        with np.load(base.TRUTH/case/'truth.npz',allow_pickle=False) as z:truth=z['chi'].copy()
        denominator=la.norm(truth)
        for task in [t for t in tasks if t['case_id']==case]:
            folder=OUT/case/task['arm'];r=base.load(folder/'RECEIPT.json')
            row=dict(case_id=case,arm=task['arm'],status=r['status'],platform=r.get('platform','Mac_AppleSilicon'),raw_NRMSE=None,final_NRMSE=None,
                original_Phi=r.get('original_Phi'),physical_source_RMS_relative_residual=r.get('physical_source_RMS_relative_residual'),
                CPU_seconds=r['phase2_incremental_CPU_seconds'],wall_seconds=r['phase2_incremental_wall_seconds'],counts=r.get('counts'),
                dpcs=r.get('dpcs'),stop=r.get('optimizer'),truth_access='OFFLINE_AFTER_ALL29_PUBLIC_SEAL')
            if r['status']=='SEALED_PUBLIC_PREDICTION':
                for filename,key in [('initial.npz','initial_NRMSE'),('raw.npz','raw_NRMSE'),('prediction.npz','final_NRMSE'),('checkpoint80_accepted.npz','checkpoint80_NRMSE'),('checkpoint80_evaluation.npz','evaluation80_NRMSE')]:
                    if (folder/filename).exists():
                        with np.load(folder/filename,allow_pickle=False) as z:m=z['m'].copy()
                        row[key]=float(la.norm(m-truth)/denominator)
                        if key=='final_NRMSE':row.update(material_min=float(m.min()),material_max=float(m.max()))
                with np.load(folder/'PAID_BUDGET_MATERIAL.npz',allow_pickle=False) as z:materials=z['m'].copy()
                costs=base.load(folder/'PAID_BUDGET.json')
                if len(materials)!=len(costs):raise RuntimeError('Cost/material trajectory mismatch')
                curve=[dict(index=j,chi_relative_L2_NRMSE=float(la.norm(m-truth)/denominator),**c) for j,(m,c) in enumerate(zip(materials,costs))]
                base.dump(OUT/'offline_curves'/f'{case}__{task["arm"]}.json',curve)
            rows.append(row)
    result=dict(rows=rows,metric='||chi_hat-chi_reference||L2/||chi_reference||L2, explicit denominator per case; NOT paper Eq17 or scalar PhaseI amplitude RMSE',
        selection='Minimum originalPhi among paid states only; no GTiterate/endpoint selection',scoring_after_all29_seal=True,failed_arms_retained=True,score_new_LU=0,score_new_RHS=0)
    base.dump(OUT/'SCORES.json',result);return rows

def report(rows):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    FIG.mkdir(exist_ok=True);cases=list(dict.fromkeys(r['case_id'] for r in rows));allarms=base.load(ROOT/'configs/vector_extension_phase2.json')['arms']
    for view in ('central_slice','maximum_projection'):
        fig,axes=plt.subplots(len(cases),9,figsize=(20,2.4*len(cases)),squeeze=False)
        for i,case in enumerate(cases):
            with np.load(base.TRUTH/case/'truth.npz',allow_pickle=False) as z:truth=z['chi'].reshape((8,)*3)
            for j,arm in enumerate(['truth']+allarms):
                ax=axes[i,j];ax.set_xticks([]);ax.set_yticks([])
                if arm=='truth':m=truth
                elif (OUT/case/arm/'prediction.npz').exists():
                    with np.load(OUT/case/arm/'prediction.npz',allow_pickle=False) as z:m=z['m'].reshape((8,)*3)
                else:ax.text(.5,.5,'NOT_RUN / FAILED',ha='center',fontsize=7);continue
                image=m[:,:,4] if view=='central_slice' else np.max(m,axis=2)
                ax.imshow(image.T,origin='lower',vmin=0,vmax=8,cmap='magma')
                if i==0:ax.set_title(arm.replace('_','\n'),fontsize=8)
                if j==0:ax.set_ylabel(case,fontsize=8)
        fig.suptitle(f'Strict real vector DDA: {view.replace("_"," ")}; chi=epsilon-1, all512unknowns\nScaled finite-model comparison; continuum accuracy NOT certified')
        fig.tight_layout(rect=[0,0,.95,.92]);cax=fig.add_axes([.965,.2,.008,.6]);fig.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(0,8),cmap='magma'),cax=cax)
        fig.savefig(FIG/f'{view}.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(13,9))
    for ax,case in zip(axes.ravel(),cases):
        for row in [r for r in rows if r['case_id']==case]:
            path=OUT/'offline_curves'/f'{case}__{row["arm"]}.json'
            if not path.exists():continue
            curve=base.load(path);ax.plot([p['CPU_seconds'] for p in curve],[p['chi_relative_L2_NRMSE'] for p in curve],label=row['arm'])
        cohort=next(r['platform'] for r in rows if r['case_id']==case)
        ax.set_title(case+'\n'+cohort);ax.set_xlabel('Online CPU s (setup/import/export separately accounted)');ax.set_ylabel('Chi relative L2 NRMSE');ax.grid(alpha=.2);ax.legend(fontsize=6)
    fig.suptitle('Offline-scored paid original-Phi incumbents; no GT endpoint selection');fig.tight_layout();fig.savefig(FIG/'paid_cost_curves.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    paths=[(8,OUT/'validation/grid8.npz'),(10,OUT/'validation/grid10.npz'),(16,OUT/'validation_refinement16/grid16.npz')]
    reference=None
    for n,path in paths:
        if not path.exists():continue
        with np.load(path,allow_pickle=False) as z:theta=z['theta'];far=z['far'];reference=z['Mie_far']
        component=np.sum(far*np.c_[np.cos(theta),theta*0,-np.sin(theta)],axis=1)
        axes[0].plot(theta*180/np.pi,component.real,label=f'grid{n}');axes[1].plot(theta*180/np.pi,component.imag,label=f'grid{n}')
    if reference is not None:
        component=np.sum(reference*np.c_[np.cos(theta),theta*0,-np.sin(theta)],axis=1)
        axes[0].plot(theta*180/np.pi,component.real,'k--',label='analytic Mie');axes[1].plot(theta*180/np.pi,component.imag,'k--',label='analytic Mie')
    for ax,title in zip(axes,['Real far-theta amplitude','Imaginary far-theta amplitude']):ax.set_title(title);ax.set_xlabel('Scattering angle degrees');ax.legend();ax.grid(alpha=.2)
    fig.suptitle('epsilon9 sphere, same physical units; no amplitude/phase fit applied');fig.tight_layout();fig.savefig(FIG/'mie_complex_angles.png',dpi=160);plt.close(fig)
    fmt=lambda x:'NOT_AVAILABLE' if x is None else f'{x:.5f}'
    lines=['# R9 Phase2：严格实材料三维矢量扩展','', '**这是缩小电尺寸的有限DDA比较，不是原2019论文3λ、30³/60³的物理复现。ε9球的连续物理误差很大，不能仅用小线性残差掩盖。** R8与PhaseI保持不变；原24臂加追加5臂全部封存后才离线评分。','',
        '每例8³=512个独立实各向同性χ未知量，1536维完整矢量场，60有限电偶极/60三分量接收，10800复观测。所有方法χ盒[0,8]、初值.25、同离散smooth/l2 prior；无GT/support在线输入，无NN。继承R1/R8 dyadic CM+RR DDA，未增加物理后端。','',
        '单CIE sourcewise current β6→1；Dual sourcewise/GlobalDual β6→3，material β6→1。A/A+仅改变材料读出，current严格与GlobalDual相同。原生projected-IDFT γ在MF1→2按频率/component精确transport；240固定PR、材料内迭代max12/A+两轮八色arc，共同physical corrector300，direct300。两类checkpoint80只复制已付费状态。Eq16 sqrt(mean ratio)仅作诊断；不宣称paper stopping或Eq17 Err。A+只适用于本轮实各向同性材料，不扩称复介电常数。','',
        '|Case|Platform|Arm|Raw χ-NRMSE|Final χ-NRMSE|Accepted80|Residual|Whole-worker CPU s|LU|DPCS rank|','|---|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        lines.append(f'|{r["case_id"]}|{r["platform"]}|{r["arm"]}|{fmt(r.get("raw_NRMSE"))}|{fmt(r.get("final_NRMSE"))}|{fmt(r.get("checkpoint80_NRMSE"))}|{fmt(r.get("physical_source_RMS_relative_residual"))}|{r["CPU_seconds"]:.2f}|{(r.get("counts") or {}).get("physical_LU",0)}|{(r.get("dpcs") or {}).get("rank","—")}|')
    lines+=['','NRMSE=||χhat−χreference||₂/||χreference||₂，与PhaseI幅度归一RMSE及paper Eq17不同；不跨池比较原始Phi。Raw是全部已付费joint物理评估中最小原Phi状态，不一定是最后joint状态；保存的stage文件保留真实联合m/J/Q/γ。','',
        '![最大投影](FIGURES/vector_extension/maximum_projection.png)','', '![中央切片](FIGURES/vector_extension/central_slice.png)','', '![已付费轨迹](FIGURES/vector_extension/paid_cost_curves.png)','', '## 物理、电尺寸和分辨率限制','',
        '原四块在.6λ箱内，side=.16λ（2.13 inverse cells），ε2/3/4/5；原paper side=.8λ。ε5时本轮side=.358内部λ，原paper为1.789内部λ。新ε9双块在.4λ箱内，side=.12λ=.36内部λ、gap=.10λ=2 inverse cells。inverse h=.05λ，6.67 cells/内部λ；fine10为8.33。高contrast不自动比原paper较大对象更难。','',
        '薄壳+双杆沿用R8 literal box几何：外边.30λ、内边.18λ、壳厚.06λ=.8 inverse cell，杆横宽.07λ=.93 cell，确属欠分辨。追加resolved hollow仅立方壳，外边.45λ、中心腔边.15λ、壳厚.15λ=2 inverse cells，ε2.5；是公开声明的R9 resolution control，不是paper复现/新blind win。几何加厚也改变散射，因此不是纯grid-convergence或因果机制证明。','',
        'cube/hollow边界使用继承的几何体积分数乘线性χ；独立球使用固定24³ midpoint的MG静态极化率混合，两者不同。Mie控制不是这些cube/hollow场的误差界。']
    for case in cases:
        r=base.load(OUT/'generation'/f'{case}.json');lines.append(f'- {case}: 独立fine10→coarse8 clean-data relative discrepancy={100*r["crossgrid_clean_relative_difference"]:.3f}%。')
    lines+=['','新三例raw iid CN3%是期望Frobenius噪声，公开σ满足E|n|²=σ²，实际噪声不归一化。legacy四块保留原固定能量10% Y。Primary DPCS历史threshold=.1；qualified在新两例用完整raw噪声球σ²χ²(.99,2nrns)/2、最小支持连续SVD rank，legacy用(.1/.9)||Y||确定界。均不涵盖上述离散模型差异，不做per-r检验/GT调rank。追加resolved hollow只跑五个primary核心臂。','', '## ε9球的独立Mie失败与一次必要细化','']
    q=base.load(OUT/'validation/RECEIPT.json')
    for r in q['rows']:lines.append(f'- grid{r["nside"]}: complex far-vector error={100*r["relative_Mie_far_vector_error"]:.3f}%，linear residual={r["linear_relative_residual"]:.3g}。')
    extra=OUT/'validation_refinement16/RECEIPT.json'
    if extra.exists():
        r=base.load(extra);err=(r.get('grid16_angular_diagnostics') or {}).get('relative_complex_vector_error');lines.append(f'- 唯一追加grid16: status={r["status"]}; relative complex error={100*err:.3f}%（若完成）; CPU={r["CPU_seconds"]:.3f}s; peakRSS={r["process_peak_RSS_bytes"]/1e6:.1f}MB。')
    lines+=['','同一εr9/χ8、r=.145λ、箱.4λ；不拟合材料、半径、归一化或复相位。grid16仅消去严格χ=0背景，使用已有TensorDDA arbitrary-occupied-cell API；保留原8/10全部失败与费用。大Mie误差实质限制ε9对象的连续物理解读；有限DDA算法配对比较仍须与连续准确性分开。','', '![复幅度角分布](FIGURES/vector_extension/mie_complex_angles.png)','',
        '## 成本与执行边界','', '最多两个单BLAS线程数值子进程，有限队列不断流。原24项均Mac；root随后明确资源迁移，厚壳5项全部从相同初值在XINAN i9指定P-core logical[8,10]运行。两条已启动Mac厚壳partial单独保留为资源交接中止，不计算法失败，不混入远端cohort成本。Mie16已在Mac完成，不重复远端。所有CPU比较须在同平台、同对象内配对；跨平台总CPU只作消费账，不能跨池声称速度。远端保持connected SSH及Windows kill-on-close JobObject，只传publicY/config/必要代码，无GT、GPU或新框架。','',
        'Whole-worker CPU包含冷启动imports/setup/SVD/online/导出；表中与已付费online曲线的计时范围不同，后者不包含setup/import/export。Windows使用psutil user+system/GetProcessTimes和peak_wset，未知内存不伪记0。Direct300与joint240+300不是等预算，不能按终点迭代数声称速度优势；应结合真实LU/RHS、CPU轨迹及停止信息。','',
        '每项RECEIPT保留真实optimizer success/message/nit/nfev与允许上限；MAXITER不是收敛。所有material-only/GB/SB缓存动作、SVD、物理LU/全部源RHS/adjoint RHS分别记账；exact-m cache命中不伪称独立求解。终点复制不增加LU/RHS。新增validation/generation与失败费用在独立receipt中，旧四块fine10不重跑。','',
        '原24任务的自动score因root追加控制而由显式SCORE_HOLD_ROOT_AMENDMENT_29预GT阻止；该hold是流程控制记录，不是有效预测seal。之后由all29 manifest统一封存并独立评分。详见delegated/3d_extension及RESULTS/vector_extension；未修改历史R8、PhaseI、原24配置或预测。','',
        '入口与数学范围：[README](delegated/3d_extension/README.md)。配置：[原24](configs/vector_extension_phase2.json)、[resolved hollow](configs/vector_extension_resolved_hollow.json)、[29项合并清单](configs/vector_extension_all29.json)。']
    (ROOT/'VECTOR3D_EXTENSION.md').write_text('\n'.join(lines)+'\n')

if __name__=='__main__':
    rows=score();report(rows);base.dump(OUT/'FINALIZATION_COST.json',dict(CPU_seconds=time.process_time()-CPU,wall_seconds=time.perf_counter()-WALL,
        score_new_LU=0,score_new_RHS=0,scope='cold imports, all29 seal, private offline scoring, figures/report; no online GT selection',peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
