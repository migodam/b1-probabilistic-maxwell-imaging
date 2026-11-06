"""Offline exposed-development assessment, after all three new predictions seal.
No physics solve, no GT selection, and no historical file changes.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import time
ENTRY=time.perf_counter();CPU=time.process_time()
import json,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'RESULTS/original_vector3d/physical_refinement';FIG=ROOT/'FIGURES/original_vector3d/physical_refinement'
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def main():
    seal=json.loads((OUT/'SEALED_ALL_THREE.json').read_text());verify=json.loads((OUT/'FETCH_SEAL_VERIFICATION.json').read_text())
    if len(seal['entries'])!=3 or verify['status']!='ALL_REMOTE_BYTES_AND_MTIMES_MATCH_SEAL':raise RuntimeError('Missing allthree verified seals')
    for e in seal['entries']:
        for f in e['files']:
            p=OUT/e['start']/Path(f['path'].replace('\\','/')).name
            if not p.exists() or p.stat().st_size!=f['bytes']:raise RuntimeError('Incomplete sealed download')
    if (OUT/'SCORES.json').exists():raise RuntimeError('Writeonce offline assessment')
    # First private material read is after all-three seal/remote verification guards.
    truth=np.load(ROOT/'RESULTS/original_vector3d/sealed_truth/truth.npz',allow_pickle=False)['inverse_material']
    Y=np.load(ROOT/'RESULTS/original_vector3d/cohorts/matched_grid/public/observations.npz',allow_pickle=False)['Y']
    events=json.loads((OUT/'PAID_EVALUATIONS.json').read_text());whole=json.loads((OUT/'RUN_PROCESS_COST.json').read_text())
    rows=[];allstates=[];FIG.mkdir(parents=True,exist_ok=True)
    fig,axes=plt.subplots(1,2,figsize=(10,3.6))
    figmat,axs=plt.subplots(2,4,figsize=(12,6));maxv=4.;cm='viridis'
    axs[0,0].imshow(truth.reshape(30,30,30).max(2),vmin=0,vmax=maxv,cmap=cm);axs[0,0].set_title('Known truth: max-z')
    axs[1,0].imshow(truth.reshape(30,30,30)[:,:,15],vmin=0,vmax=maxv,cmap=cm);axs[1,0].set_title('Known truth: z=15')
    for col,e in enumerate(seal['entries'],1):
        name=e['start'];d=OUT/name;receipt=json.loads((d/'RECEIPT.json').read_text());states=[]
        for p in [d/'initial.npz']+sorted(d.glob('accepted_*.npz'))+[d/'prediction.npz']:
            with np.load(p,allow_pickle=False) as z:m=z['material'];prediction=z['prediction']
            info=json.loads(p.with_suffix('.json').read_text());phi=float(np.mean(np.sum(abs(prediction-Y)**2,axis=0)/np.sum(abs(Y)**2,axis=0)))
            state=dict(start=name,label=p.stem,state_id=info['state_id'],Phi=phi,stored_Phi_difference=phi-info['Phi'],PG_inf=info['PG_inf'],
                contrast_NRMSE=float(np.linalg.norm(m-truth)/np.linalg.norm(truth)),paper_Err=float(np.sqrt(np.mean((m-truth)**2/(1+truth)**2))),
                contrast_min=float(m.min()),contrast_max=float(m.max()),mean_contrast=float(m.mean()))
            states.append(state);allstates.append(state)
        initial,final=states[0],states[-1];es=receipt['events'];unique={x['state_id'] for x in es if not x['cache_hit']}
        grad=sum(bool(x['gradient_requested']) and (not x['cache_hit'] or x['gradient_upgrade_from_paid_J']) for x in es)
        forward=sum(x['actual_primal_RHS'] for x in es)//60;adj=sum(x['actual_adjoint_RHS'] for x in es)//60
        row=dict(start=name,initial_Phi=initial['Phi'],Phi=final['Phi'],initial_PG_inf=initial['PG_inf'],PG_inf=final['PG_inf'],
            initial_contrast_NRMSE=initial['contrast_NRMSE'],contrast_NRMSE=final['contrast_NRMSE'],paper_Err=final['paper_Err'],
            contrast_min=final['contrast_min'],contrast_max=final['contrast_max'],strict_requests=receipt['requests'],accepted_steps=receipt['accepted_steps'],
            actual_new_state_records=len(unique),actual_new_gradient_evaluations=grad,actual_forward_batches_60=forward,actual_adjoint_batches_60=adj,
            actual_primal_RHS=sum(x['actual_primal_RHS'] for x in es),actual_adjoint_RHS=sum(x['actual_adjoint_RHS'] for x in es),
            cache_hits=sum(x['cache_hit'] and not x['gradient_upgrade_from_paid_J'] for x in es),
            CPU_seconds=receipt['CPU_seconds'],wall_seconds=receipt['wall_seconds'],stop=receipt['stop'],selected_state_id=receipt['selected_state_id'])
        rows.append(row)
        accepted=states[1:-1];axes[0].semilogy(range(len(accepted)+1),[initial['Phi']]+[x['Phi'] for x in accepted],'-o',label=name)
        axes[1].semilogy(range(len(accepted)+1),[initial['PG_inf']]+[x['PG_inf'] for x in accepted],'-o',label=name)
        with np.load(d/'prediction.npz',allow_pickle=False) as z:m=z['material'].reshape(30,30,30)
        axs[0,col].imshow(m.max(2),vmin=0,vmax=maxv,cmap=cm);axs[0,col].set_title(name+' max-z')
        axs[1,col].imshow(m[:,:,15],vmin=0,vmax=maxv,cmap=cm);axs[1,col].set_title('z=15; NRMSE %.3f'%final['contrast_NRMSE'])
    for a,y in zip(axes,['Original physical Phi','Projected gradient infinity norm']):a.set_xlabel('Accepted L-BFGS-B steps');a.set_ylabel(y);a.grid(alpha=.25)
    axes[0].legend(fontsize=7);fig.tight_layout();fig.savefig(FIG/'physical_curves.png',dpi=170);plt.close(fig)
    for a in axs.ravel():a.set_xticks([]);a.set_yticks([])
    figmat.suptitle('Exposed-development CM30 finite-DDA materials; fixed range chi=[0,4]');figmat.tight_layout();figmat.savefig(FIG/'materials.png',dpi=170);plt.close(figmat)
    result=dict(classification='PUBLIC_DEVELOPMENT_EXPOSED_ORIGINAL9_ENDPOINTS_NOT_CONFIRMATION',selection='Paid public Phi only; offline truth never chooses endpoint',
        rows=rows,all_saved_states=allstates,GT_access_after_all_three_seal=True,physical_process_CPU_seconds=whole['CPU_seconds'],physical_process_wall_seconds=whole['wall_seconds'],
        scope='Finite CM+RR DDA refinement, not qualified continuum reconstruction or new CIE run')
    dump(OUT/'SCORES.json',result)
    with (OUT/'TABLE.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    report=['# 原尺寸 matched-grid 短原物理校正（development）','',
        '三个起点是 air、已公开 GlobalDual stage2、已公开 AplusGlobal stage2。原九臂 GT 已暴露，本组不是新的确认实验。单驻留 complex128 packed16 CUDA 算子；GMRES 仍逐源，未使用或宣称共享 LU。',
        '', '共同原目标为 mean_source ||SJ−Y||²/||Y||²，无 1/2、无先验，χ∈[0,4]。L-BFGS-B 固定最多 8 个 accepted 步、24 次请求、12 次 line-search，ftol=1e−12、gtol=1e−8；所有拒绝 trial 都付正/伴随费用。按每臂所有已付状态的最低公开 Phi 选择，不按 GT。', '',
        '|起点|初始→最终 Phi|初始→最终 PG∞|初始→最终 contrast NRMSE|accepted / 请求|新增正/伴随 RHS|CPU s / wall s|','|---|---:|---:|---:|---:|---:|---:|']
    for r in rows:report.append(f"|{r['start']}|{r['initial_Phi']:.6g} → {r['Phi']:.6g}|{r['initial_PG_inf']:.4g} → {r['PG_inf']:.4g}|{r['initial_contrast_NRMSE']:.4g} → {r['contrast_NRMSE']:.4g}|{r['accepted_steps']} / {r['strict_requests']}|{r['actual_primal_RHS']} / {r['actual_adjoint_RHS']}|{r['CPU_seconds']:.3f} / {r['wall_seconds']:.3f}|")
    report+=['','请求数与新求解数不同：initial 在 SciPy fun 前计算，fun 首请求命中相同材料缓存；air 初值已由梯度接线支付，解析 J=0 与伴随无 GMRES RHS。每次非零新 value+gradient 付 60 primal + 60 adjoint；回调读取已付状态，不再计算。上表为每臂此次新增实际归属费用，不能再加到整进程费用。两个 stage2 起点的既有 CIE 运行费用仍在旧九臂账；本表不是包含这些前缀的独立冷启动总成本，也不能把前缀重加为新消费。','',
        '梯度接线只在公开 air 点做预定单向 h=1e−5 检查；无必要 h/2，不复做旧恒等式。三个终点及每个 accepted 材料/预测全部封存后，才统一离线读材料 GT。contrast NRMSE 分母为 ||χ_truth||₂；paper Err 为 sqrt(mean((χ−χ_truth)²/(1+χ_truth)²))。PG∞ 是 ||χ−clip(χ−g,0,4)||∞，使用未重标定的原目标梯度。', '',
        '本组只能检验有限 8 步物理方向是否继续改善这些公开端点，不能证明最终驻点、唯一性、全局最优，或依靠继续迭代拯救原 CIE。小梯度、残差下降与材料恢复分别报告。保留旧 CM+RR 高介电/原尺寸连续误差限制。', '',
        '结果：RESULTS/original_vector3d/physical_refinement/{SCORES.json,TABLE.csv,GRADIENT_WIRING.json,RUN_PROCESS_COST.json,ACCOUNTING_FINAL.json}。图：FIGURES/original_vector3d/physical_refinement/{physical_curves.png,materials.png}。']
    (ROOT/'delegated/3d_extension/physical_refinement/REPORT_ZH.md').write_text('\n'.join(report)+'\n')
    dump(OUT/'LOCAL_REPORT_PROCESS_COST.json',dict(platform='Mac',CPU_seconds=time.process_time()-CPU,wall_seconds=time.perf_counter()-ENTRY,additive=True,included_in=None,operation='OFFLINE_ALLTHREE_SCORE_FIGURE_REPORT',scope='Local imports/sealguard/offline savedarray scoring/two figures/report; zero solves',counts=dict(new_physical_RHS=0)))
if __name__=='__main__':main()
