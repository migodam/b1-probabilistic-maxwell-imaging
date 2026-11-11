"""Post-all12-seal scoring of actual true3D high-chi material reconstructions."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
from datetime import datetime,timezone
import sys,json,time,csv
import numpy as np

HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3]
sys.path.insert(0,str(REPO))


def read(path):return json.loads(Path(path).read_text())


def dump(path,value):
    def convert(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value,indent=2,default=convert,allow_nan=False)+'\n')


def metric(x,t):
    x=np.asarray(x).real;t=np.asarray(t).real
    return dict(relative_material_l2=float(np.linalg.norm(x-t)/np.linalg.norm(t)),
        amplitude8_NRMSE=float(np.sqrt(np.mean((x-t)**2))/8),
        material_min=float(np.min(x)),material_max=float(np.max(x)),
        legal_box=bool(np.min(x)>=-1e-12 and np.max(x)<=8+1e-12))


def main():
    cpu,wall=time.process_time(),time.perf_counter();cfg=read(HERE/'FROZEN_PROTOCOL.json')
    seal=read(HERE/'ALL12_PREDICTIONS_SEALED.json');prefixseal=read(HERE/'PREFIX_SELECTION_SEAL.json')
    assert seal['status']=='ALL12_PREDICTIONS_AND_PREFIX_CHOICES_SEALED'
    assert {(r['case_id'],r['method']) for r in seal['rows']}=={(c,m) for c in cfg['cases'] for m in cfg['methods']}
    for row in seal['rows']:
        assert Path(row['prediction']).is_file()
        assert read(row['receipt'])['status']=='SEALED_PUBLIC_PREDICTION'
    # Every prediction and every CPU-prefix choice is sealed before this point.
    private_start=datetime.now(timezone.utc).isoformat()
    truths={};pictures={};witness=[]
    for cid in cfg['cases']:
        with np.load(HERE/'sealed_truth'/f'{cid}.npz',allow_pickle=False) as z:
            t=z['chi'].copy();truths[cid]=t;pictures[(cid,'truth')]=t
            mismatch=float(np.linalg.norm(z['Y_clean']-z['Y_coarse_clean'])/np.linalg.norm(z['Y_clean']))
            witness.append(dict(case_id=cid,generation_side=int(z['generation_side']),inverse_side=6,
                inverse_truth_max=float(np.max(t)),forward_material_max=float(np.max(z['generation_chi'])),
                relative_fine_to_coarse_clean_difference=mismatch,
                truth_kind='actual analytic volumetric mean of same continuous rectilinear geometry',
                paired_with='random_high8' if cid.endswith('crossgrid8') else None))
    assert np.array_equal(truths['random_high8'],truths['random_high8_crossgrid8'])
    # A receiver projection of saved, already-paid physical currents supplies
    # raw global discrepancy alongside the inherited source-normalized score.
    # This is not a new Maxwell solve. Operator assembly and S@J are charged.
    from model_adapter import HighContrastVectorPR3D
    model=HighContrastVectorPR3D(ns=6)
    rows=[];receipts=[];receiver_columns=0;raws=[]
    for cid in cfg['cases']:
        with np.load(HERE/'data_public'/f'{cid}.npz',allow_pickle=False) as z:Y=z['Y'].copy()
        for method in cfg['methods']:
            out=HERE/'RESULTS'/cid/method;rec=read(out/'RECEIPT.json');receipts.append(rec)
            with np.load(out/'prediction.npz',allow_pickle=False) as z:
                x=z['prediction'].real.copy();raw=z['raw'].real.copy();J=z['J'].copy()
            raw_res=float(np.linalg.norm(model.S@J-Y));receiver_columns+=model.P
            row=dict(case_id=cid,method=method,**metric(x,truths[cid]),
                raw_relative_material_l2=metric(raw,truths[cid])['relative_material_l2'],
                physical_phi=rec['physical_phi'],physical_data_residual_source_normalized=rec['physical_data_residual'],
                prior=rec['prior'],raw_global_residual=raw_res,
                raw_global_relative_residual=raw_res/float(np.linalg.norm(Y)),
                raw_global_delta=rec['delta'],raw_CN99_pass=bool(raw_res<=rec['delta']),
                raw_phi=rec['raw_phi'],raw_data_residual=rec['raw_data_residual'],
                CPU_seconds=rec['cpu_seconds'],wall_seconds=rec['wall_seconds'],
                proposal_CPU_seconds=rec['proposal_cpu_seconds'],corrector=rec['corrector'],
                counts=rec['counts'],projection=rec['projection'],setup=rec['setup'])
            rows.append(row);pictures[(cid,method)]=x;raws.append((cid,method,raw))
    prefixes=[]
    with np.load(HERE/'PREFIX_PREDICTIONS.npz',allow_pickle=False) as z:
        for row in prefixseal['rows']:
            prefixes.append({**row,**metric(z[row['array_key']],truths[row['case_id']])})
    primary=cfg['cases'][:2]
    result=dict(status='OFFLINE_SCORED_AFTER_ALL12_PUBLIC_SEAL',
        seal_created_utc=seal['created_utc'],first_private_access_utc=private_start,
        independent_geometries=2,crossgrid_pair_not_independent=True,
        rows=rows,crossgrid_contract=witness,
        matched_two_geometry_mean_relative_l2={method:float(np.mean([r['relative_material_l2'] for r in rows if r['method']==method and r['case_id'] in primary])) for method in cfg['methods']},
        mean_scope='only firsttwo independent matched geometries; crossgrid pair reported separately',
        primary_metric='relative material L2 as earlier R12 vector study; amplitude8NRMSE also retained',
        qualification=cfg['qualification'],cell_phase_max=cfg['cell_phase_max'],
        prior='unchanged inherited quadraticprior; no labels in onlineoptimizer or candidateselection',
        new_scoring_physical_LU=0,new_scoring_physical_RHS=0,
        scoring_receiver_projection_columns=receiver_columns,scoring_operator_setup=model.setup)
    dump(HERE/'OFFLINE_SCORES.json',result);dump(HERE/'OFFLINE_PREFIX_SCORES.json',dict(rows=prefixes,new_physical_LU=0,new_physical_RHS=0))
    flat=[{k:v for k,v in r.items() if not isinstance(v,(dict,list))} for r in rows]
    with (HERE/'OFFLINE_SCORES.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in flat for k in r)));writer.writeheader();writer.writerows(flat)
    ledger=dict(status='COMPLETE',generation=read(HERE/'GENERATION_RECEIPT.json'),
        solver_cpu_seconds=sum(r['cpu_seconds'] for r in receipts),
        solver_wall_seconds_serial_sum=sum(r['wall_seconds'] for r in receipts),
        physical_LU=sum(r['counts']['physical_LU'] for r in receipts),
        physical_source_RHS=sum(r['counts']['source_RHS'] for r in receipts),
        physical_adjoint_RHS=sum(r['counts']['adjoint_RHS'] for r in receipts),
        H_cholesky=sum((r['projection'] or {}).get('counts_including_setup',{}).get('H_cholesky',0) for r in receipts),
        H_setup_RHS=sum((r['projection'] or {}).get('counts_including_setup',{}).get('H_setup_rhs',0) for r in receipts),
        H_iteration_RHS=sum((r['projection'] or {}).get('counts_including_setup',{}).get('H_iteration_rhs',0) for r in receipts),
        rows=[dict(case_id=r['case_id'],method=r['method'],CPU_seconds=r['cpu_seconds'],
                   counts=r['counts'],setup=r['setup'],projection=r['projection'],corrector=r['corrector']) for r in receipts],
        overlapping_scope='phaseCPU overlaps wholeCPU; setup/factors/actions/currentPR/materialcalls/physicalLU all retained; no equal-total-cost claim',
        failed_solver_runs=0,module_imports='before solve clocks; excluded')
    dump(HERE/'COST_LEDGER.json',ledger)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    names=['truth']+cfg['methods'];fig,axes=plt.subplots(9,5,figsize=(15,24),layout='constrained')
    for i,cid in enumerate(cfg['cases']):
        for iz,zindex in enumerate([1,2,3]):
            for j,name in enumerate(names):
                im=axes[i*3+iz,j].imshow(pictures[(cid,name)].reshape(6,6,6)[:,:,zindex],vmin=0,vmax=8,cmap='viridis',origin='lower')
                title=name
                if name!='truth' and iz==0:
                    r=next(r for r in rows if r['case_id']==cid and r['method']==name)
                    title+=f"\nrelL2 {r['relative_material_l2']:.3f}, Phi {r['physical_phi']:.3g}"
                axes[i*3+iz,j].set_title(title,fontsize=10);axes[i*3+iz,j].set_xticks([]);axes[i*3+iz,j].set_yticks([])
                if j==0:axes[i*3+iz,j].set_ylabel(f'{cid}\nz slice {zindex}',fontsize=9)
    fig.colorbar(im,ax=list(axes.ravel()),shrink=.5,pad=.01,label='chi, common [0,8] scale')
    fig.suptitle('High-contrast finite true3D DDA: two independent objects plus paired 8³ to 6³ observation',fontsize=13)
    fig.savefig(HERE/'MATERIAL_SLICES.png',dpi=130);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(16,4.5),layout='constrained')
    for ax,cid in zip(axes,cfg['cases']):
        for method in cfg['methods']:
            h=read(HERE/'RESULTS'/cid/method/'history.json')
            ax.step([r['cpu_seconds'] for r in h],[r['incumbent_phi'] for r in h],where='post',label=method)
        ax.set_yscale('log');ax.set_title(cid);ax.set_xlabel('Paid CPU seconds, incl. proposal');ax.set_ylabel('Best paid original Phi');ax.grid(alpha=.2)
    axes[0].legend(fontsize=8);fig.savefig(HERE/'PHI_CPU.png',dpi=150);plt.close(fig)
    lines=['# R12 高对比度 true3D 补充','',
        '两个独立 seeded 非球形三维对象：chi8，以及混合 chi2/8；再增加第一个对象的8³正演→6³反演配对。配对不是第三个独立对象。seed2026111213，在随机几何和新观察前冻结。',
        '', '几何是带空腔的随机长方体并集/优先级材料分区；6³和8³均通过连续分区的解析体积交计算真实 cell 平均材料，未使用 cell 中心标签近似。12个最终预测及 CPU前缀选择全部封存后才离线读取truth。',
        '', '每个逆问题216个共享实材料、648个矢量电流分量、60个原 R8 电偶极照明、180接收通道。统一初始chi=.25，box[0,8]，保留原source-normalized物理Phi和原quadraticprior。',
        '', 'CIE/Dual保留原R8 current_PR、material、physical与prepare，120cycles、inner12、mf1→2；Dual晚期currentbeta3。AP/DR均3000步、audit50、scale6/omega1、groupedports；projection320CPU秒上限含setup。共同corrector至多1000迭代/60CPU秒。各方案的proposal费用不同，不声称总费用相等。',
        '', 'R8 projected-IFFT basis默认经过SVD range白化；这里完整披露该坐标选择，不能称为逐字2019论文复现或原始gamma坐标。物理值缓存只复用已付费LU得到的场，不改变方程。共享R8[0,4]优化器未改写，本补充独立wrapper向原优化步骤提供[0,8]。',
        '', '|对象|方法|材料相对L2|proposal材料相对L2|原Phi|source-normalized residual|CPU秒|corrector停止|',
        '|---|---|---:|---:|---:|---:|---:|---|']
    for r in rows:lines.append(f"|{r['case_id']}|{r['method']}|{r['relative_material_l2']:.6f}|{r['raw_relative_material_l2']:.6f}|{r['physical_phi']:.8g}|{r['physical_data_residual_source_normalized']:.6g}|{r['CPU_seconds']:.2f}|{r['corrector']['stop']}|")
    lines+=['', '两独立 matched 对象的平均材料相对L2（不混入cross-grid配对）：', '']
    for method,value in result['matched_two_geometry_mean_relative_l2'].items():lines.append(f'- {method}: {value:.6f}')
    lines+=['', '本补充的明确负结果：两个 matched 对象经过共同 corrector 后，四臂最终材料误差与原Phi均非常接近，DR没有最终材料收益且总CPU更高。DR在proposal阶段的材料相对L2为0.061631/0.078315，共同corrector进一步降低原Phi后却变为0.129783/0.227239。cross-grid中DR原Phi略低，材料相对L2反而略高（0.434481 vs 约0.43222）。这三个观察均保留，不能只报告proposal材料优势或更低Phi。',
        '', '该小规模高对比度离散挑战不是连续 Maxwell 精度证明。chi8对应epsilon_r9；最大内部 cell 相位约1.885（6³）/1.414（8³），网格很粗。cross-grid实际 clean数据差异与原始CN99噪声球分别保留；只用sigma定义的球没有包含模型离散误差。',
        '', '|对象|8³/6³ clean相对差异|coarse truth最大chi|forward最大chi|',
        '|---|---:|---:|---:|']
    for w in witness:lines.append(f"|{w['case_id']}|{w['relative_fine_to_coarse_clean_difference']:.6g}|{w['inverse_truth_max']:.6g}|{w['forward_material_max']:.6g}|")
    lines+=['', 'OFFLINE_PREFIX_SCORES 使用先前已按原Phi冻结的10/30/60/120CPU秒候选，没有用GT排名。相对函数容差终止不是局部最优或唯一性证书；低物理残差与材料恢复单独列出。全部负结果、原始完整projection audits、成本和corrector停止原因保留。',
        '', f"Solver CPU合计 {ledger['solver_cpu_seconds']:.3f} 秒，physicalLU {ledger['physical_LU']}，sourceRHS {ledger['physical_source_RHS']}，adjointRHS {ledger['physical_adjoint_RHS']}；projection H分解 {ledger['H_cholesky']:.0f}，setupRHS {ledger['H_setup_RHS']:.0f}，iterationRHS {ledger['H_iteration_RHS']:.0f}。详见COST_LEDGER，嵌套timers/counters不重复相加。", 
        '', '复现：generate_high.py（仅新目录）；run_high.py run；整组seal后score_high.py。无训练、无在线truth、无GPU占用、无push。科学主张由主线程最终评审。']
    (HERE/'RESULTS_ZH.md').write_text('\n'.join(lines)+'\n')
    dump(HERE/'SCORING_RECEIPT.json',dict(status='COMPLETE_POSTSEAL_OFFLINE',
        cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        new_physical_LU=0,new_physical_RHS=0,receiver_projection_columns=receiver_columns,
        public_operator_setup=model.setup))
    print(json.dumps(dict(status=result['status'],matched_means=result['matched_two_geometry_mean_relative_l2'],
                          solverCPU=ledger['solver_cpu_seconds'])),flush=True)


if __name__=='__main__':main()
