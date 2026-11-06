"""Historical-only offline scoring, after all21 R8 predictions are sealed.

Only the four old truth entries are opened. Never searches/opens R8 new42 GT.
Old uniform predictions are read-only reused, with their original version/cost.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import json,time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from research.research_max_struct_r8.experiments.run_case import ROOT,REPO,dump,load_public
from research.research_max_struct_r8.experiments.history_warm_prepare import IDS
from research.research_max_struct_r7.experiments.score_primary import material_metrics

def read(p):return json.loads(p.read_text())
def main():
    begin,cpu=time.perf_counter(),time.process_time();out=ROOT/'RESULTS/history_warm';seal=read(out/'ALL_PREDICTIONS_SEALED.json');manifest=read(ROOT/'experiments/history_warm_manifest.json')
    assert seal['status']=='ALL21_HISTORICAL_PREDICTIONS_SEALED' and len(seal['records'])==21
    assert set(t['case_id'] for t in manifest['tasks'])==set(IDS)
    frozen=read(ROOT/'configs/solver.json');audit=[]
    for t in manifest['tasks']:
        rec=read(REPO/t['output_dir']/'RECEIPT.json');assert rec['status']=='SEALED_PUBLIC_PREDICTION'
        assert rec['configuration']==frozen and not t.get('initialize_physical_current',False)
        assert rec['material_dof']==576 and rec['truth_access']=='NONE'
        audit.append(dict(task_id=t['task_id'],same_frozen_R8_configuration=True,material_dof=576,physical_current_projection=False,truth_access=rec['truth_access'],arm=t['arm'],initialization=t['initialization']))
    assert sum('initialization_file' in t for t in manifest['tasks'])==16
    dump(out/'PROTOCOL_AUDIT.json',dict(status='ALL21_FIXED_SOLVER_AND_MATERIAL_ONLY_INITIALIZATION_VERIFIED',records=audit,NN_tasks=16,uniform_sourcewise_controls=4,uniform_global_control_42cb=1,new42_truth_access='NONE'))
    # First historical private opens occur here, after all21 new endpoint seals.
    historical={r['case_id']:r for r in read(REPO/'research/research_max_struct_r6/sealed_truth/scoring/historical_mapping_PRIVATE.json')}
    old_confirmation={r['case_id']:r for r in read(REPO/'research/research_max_struct_r6/sealed_truth/confirmation/mapping_PRIVATE.json')['records']}
    truths={};supports={};truthpaths=[]
    for ident in IDS:
        meta=historical.get(ident,old_confirmation.get(ident));path=REPO/meta['truth_path']
        with np.load(path,allow_pickle=False) as z:truths[ident]=z['m'].ravel().copy();supports[ident]=z['mask'].ravel().astype(bool)
        truthpaths.append(str(path.relative_to(REPO)))
    rows=[];pictures={};artifacts=[]
    def append(ident,label,path,scope):
        receipt=read(path/'RECEIPT.json');assert receipt['status']=='SEALED_PUBLIC_PREDICTION'
        with np.load(path/'prediction.npz',allow_pickle=False) as z:m=z['m'].ravel();yp=z['Y_prediction']
        metrics=material_metrics(m,truths[ident],supports[ident],12.);cost=receipt['physical_model_cost'];task=receipt['task'];Y=load_public(task['public_files'])['Y']
        row=dict(case_id=ident,method=label,cohort='EXPOSED_HISTORICAL_REGRESSION_ONLY',reuse_status=scope,prediction_source=str((path/'prediction.npz').relative_to(REPO)),receipt_source=str((path/'RECEIPT.json').relative_to(REPO)),cpu_seconds=receipt['cpu_seconds'],wall_seconds=receipt['wall_seconds'],LU=cost.get('LU',0),complex_RHS=cost.get('actual_complex_RHS',0),relative_residual=float(np.linalg.norm(yp-Y)/np.linalg.norm(Y)),original_Phi=receipt['physical_objective'],corrector_stop=receipt['corrector']['stop'],corrector_iterations=receipt['corrector']['iterations'],auxiliary_iterations=receipt['auxiliary_iterations'],**metrics)
        rows.append(row);pictures[ident,label]=m.copy();artifacts.append(dict(case_id=ident,method=label,path=str(path.relative_to(REPO)),config=receipt['configuration'],status=receipt['status']))
    for ident in IDS:
        for arm in ['cie_single','cie_dual','global_dual','a_material']:
            p=REPO/'research/research_max_struct_r6/RESULTS/primary'/(ident+'__'+arm)
            append(ident,'R6_uniform_'+arm,p,'REUSED_R6_DENSE_OR_EXPLICIT_BASELINE_SAME3000_AUX3000_CORRECTOR_CAP')
        p=REPO/'research/research_max_struct_r7/delegated/aplus/trajectories'/(ident+'__mobius_arc_candidates')
        if p.exists():append(ident,'R7_uniform_global_Aplus',p,'REUSED_R7_FINAL_ARC_A_PLUS_NATIVEFFT_SAME3000_AUX3000_CORRECTOR_CAP')
        for task in [t for t in manifest['tasks'] if t['case_id']==ident]:
            mode='nn' if 'initialization_file' in task else 'uniform';arm=task['arm']
            append(ident,'R8_'+mode+'_'+arm,REPO/task['output_dir'],'NEW_FIXED_R8_SAME_SOLVER_ORIGINAL_PHYSICS_WEIGHT')
    dump(out/'SCORES.json',dict(status='HISTORICAL_ONLY_SCORED_AFTER_ALL21_SEALED',rows=rows,new_tasks=21,old_predictions_reused=len(rows)-21,objects=4,truth_paths_opened=truthpaths,new42_truth_access='NONE',new_physical_calls_for_scoring=0))
    dump(out/'BASELINE_PROVENANCE.json',artifacts)
    initial=[]
    for ident in IDS:
        for mode,arm in [('nn','cie_dual'),('uniform','a_plus_sourcewise')]:
            task=next(t for t in manifest['tasks'] if t['case_id']==ident and t['arm']==arm and ('initialization_file' in t)==(mode=='nn'));p=REPO/task['output_dir']
            with np.load(p/'initial.npz',allow_pickle=False) as z:m=z['m'].ravel()
            r=read(p/'initial.json');initial.append(dict(case_id=ident,initialization=mode,initial_original_Phi=r['objective'],initial_relative_residual=r['physical_relative_residual'],**material_metrics(m,truths[ident],supports[ident],12.)))
    dump(out/'INITIALIZATION_SCORES.json',dict(rows=initial,scope='OfflineinitialGTdistance only afterALL21seal; originalinitialPhi alreadypaid; no physicalcalls and no NN rejection'))
    # One complete figure keeps allfour objects, onecommon row range inclmax.
    labels=['Reference','R6_uniform_cie_single','R6_uniform_cie_dual','R6_uniform_global_dual','Historical_A_control','R8_uniform_a_plus_sourcewise','R8_nn_cie_dual','R8_nn_a_plus','R8_nn_a_plus_sourcewise','R8_nn_physical_lbfgs']
    fig,axes=plt.subplots(4,len(labels),figsize=(30,13),constrained_layout=True)
    pretty={'R6_uniform_cie_single':'R6 uniform CIE','R6_uniform_cie_dual':'R6 uniform Dual','R6_uniform_global_dual':'R6 uniform GlobalDual','R8_uniform_a_plus_sourcewise':'R8 uniform sourcewise A+','R8_nn_cie_dual':'NN + Dual','R8_nn_a_plus':'NN + global A+','R8_nn_a_plus_sourcewise':'NN + sourcewise A+','R8_nn_physical_lbfgs':'NN + direct physics'}
    for i,ident in enumerate(IDS):
        old='R8_uniform_a_plus' if (ident,'R8_uniform_a_plus') in pictures else 'R7_uniform_global_Aplus' if (ident,'R7_uniform_global_Aplus') in pictures else 'R6_uniform_a_material';mats=[truths[ident]]+[pictures[ident,k] for k in labels[1:4]]+[pictures[ident,old]]+[pictures[ident,k] for k in labels[5:]];vmax=max(float(m.max()) for m in mats)
        for j,(label,m) in enumerate(zip(labels,mats)):
            ax=axes[i,j];im=ax.imshow(m.reshape(24,24),origin='lower',vmin=.02,vmax=vmax,cmap='viridis');title=pretty.get(label,label)
            if j==0:title=ident.replace('obs_','')+' reference'
            else:
                actual=old if label=='Historical_A_control' else label;r=next(r for r in rows if r['case_id']==ident and r['method']==actual);title=('R8 uniform global A+' if old.startswith('R8') else 'R7 uniform global A+' if old.startswith('R7') else 'R6 uniform original A') if label=='Historical_A_control' else title;title+=f"\nN={r['material_nrmse']:.3f} IoU={r['otsu_iou']:.2f}\nr={r['relative_residual']:.3g} CPU={r['cpu_seconds']:.1f}"
            ax.set_title(title,fontsize=8);ax.set_xticks([]);ax.set_yticks([])
        fig.colorbar(im,ax=list(axes[i]),shrink=.7,label='Full common row material scale')
    fig.suptitle('Four exposed historical objects | original10feature NN | all pixels free | old versions explicitly labeled');folder=ROOT/'FIGURES/history_warm';folder.mkdir(parents=True,exist_ok=True);fig.savefig(folder/'all_four_objects.png',dpi=160);plt.close(fig)
    prep=read(ROOT/'data_public/history_warm/PREPROCESS_RECEIPT.json');training=read(ROOT/'RESULTS/training/physics/RECEIPT.json');new=[r for r in rows if r['method'].startswith('R8_')]
    first=read(out/'INITIAL20_BATCH_SEALED.json') if (out/'INITIAL20_BATCH_SEALED.json').exists() else {}
    cost=dict(status='HISTORICAL21_COSTS',new_solver_CPU_seconds=sum(r['cpu_seconds'] for r in new),new_solver_wall_seconds_sum=sum(r['wall_seconds'] for r in new),serial_wrapper_CPU_seconds=first.get('cpu_seconds',0)+seal.get('cpu_seconds_this_invocation',0),serial_wrapper_wall_seconds=first.get('wall_seconds',0)+seal.get('wall_seconds_this_invocation',0),wrapper_scope='First20batch plusserialonecontrolresume; priorreceiptsreused, do notsumalgorithmcosttwice',LU=sum(r['LU'] for r in new),complex_RHS=sum(r['complex_RHS'] for r in new),preprocessing_CPU_seconds=prep['cpu_seconds'],preprocessing_wall_seconds=prep['wall_seconds'],model_load_shared=prep['model_load'],inference_features_records=prep['records'],training_cost=dict(scope='Reused originalphysics training; external to new21 cost and must amortize/add separately; generation cost separate in listed receipts',CPU_seconds=training['cpu_seconds'],wall_seconds=training['wall_seconds'],parameters=training['parameters'],generation_cost_receipts=training['generation_cost_receipts']),old_baseline_costs='Recorded historical costs only; no new charge; differentversion/operatorwiring prevents a purewalltime causalclaim')
    dump(out/'COST_LEDGER.json',cost)
    lines=['# 历史困难对象的Neural Warm Start有限回归','', '本轮仅复查已公开的7b8d、b20d、42cb、40ac四个旧开发对象，21个新任务（含唯一补充42cb uniform global A+配对控制）；不是20个新材料，也不加入新42确认集。第四对象40ac是不同拓扑开缝壳，由父线程在任何NN结果前选定；R6另一历史组合的第四63ab未纳入，R7原四对象还含简单近双岛00f1，集合差异不隐藏。','',
      '原10feature physics/weights.pt冻结，未使用domain_augmented权重；FeatureOperator按旧.96米域坐标和真实公开G/E/S/Y直接计算。没有用初始Phi拒绝NN。所有NN只初始化材料，关闭可选physical-current projection；Dual及sourcewise A+对每一源列分别计算current共轭梯度系数；global A+以所有源的联合current计算一个系数。两种A+均使用同一个Mobius-material arc候选与原Phi选择，sourcewise不是把方法名字改成global，也不是增加可观测信息。单数值进程、每BLAS一个线程串行执行；首批20项与补充1项分两次串行调用，互不并行，未重跑已封存项。NN在前处理一次装载，冷装载与features/推理费用另列。','',
      '原solver.json未改：mf3/6/10各500步＋mf12×1500（共3000辅助上限），共同实材料[.02,20]、原prior、原DPCS和Polak–Ribiere更新；corrector最多3000，direct最多6000。本轮没有重新生成Y、几何或噪声，沿用原离散，不能当成独立网格或真实测量验证。主预测仍按原Phi选取，GT只在21个新端点全部封存后由本评分器打开上述四个历史路径；没有读取新42 GT或新增物理评分调用。','',
      '## 初始材料距离和初始原物理目标','', '| 对象 | 初值 | amplitudeNRMSE | 初始原Phi | 初始residual |','|---|---|---:|---:|---:|']
    for r in initial:lines.append(f"|{r['case_id']}|{r['initialization']}|{r['material_nrmse']:.5g}|{r['initial_original_Phi']:.6g}|{r['initial_relative_residual']:.5g}|")
    lines+=['','NN初始原Phi是否更高不决定其去留；所有固定预算完整执行。初始GT距离只在封存后报告，不能用它在线选择NN或uniform。','', '## 全部端点','', '| 对象 | 版本/方法 | amplitudeNRMSE | OtsuIoU | residual | 原Phi | CPU(s) | LU | complexRHS | corrector步 | stop |','|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|']
    for r in rows:lines.append('| '+r['case_id']+' | '+r['method']+' | '+' | '.join(f'{r[k]:.5g}' for k in ['material_nrmse','otsu_iou','relative_residual','original_Phi','cpu_seconds','LU','complex_RHS','corrector_iterations'])+' | '+r['corrector_stop']+' |')
    lines+=['','R6 uniform CIE/Dual/Global/A为原封存旧版本，R7 final-arc A+在7b8d/b20d/40ac直接复用原封存预测；42cb没有相同预算的R7 final-arc A+文件，因此仅补一个R8 uniform global A+配对控制；原R6 A仍独立列出，不把rational额外校正端点改名作A+。R6 explicit/dense算子wiring与R8 nativeFFT存在实现/长期roundoff路径差异，不能把跨版本差异全部归因为初始化。R7与R8关闭初始current投影后的FFT数值主体一致；原配置与来源在BASELINE_PROVENANCE.json。','',
      '![Allfourhistoricalobjects](FIGURES/history_warm/all_four_objects.png)','','图中每行GT及全部候选使用共同完整幅度范围，没有只显示成功对象或截断高峰材料。NRMSE=RMSE/12；IoU使用预测自身256bin Otsu阈值，原R7材料判据NRMSE≤.15且IoU≥.7并合法，仅作历史可比标签，不是新实验停止gate。','',
      f"新增21任务solverCPU合计{cost['new_solver_CPU_seconds']:.3f}s，LU={cost['LU']}，complexRHS={cost['complex_RHS']}；另加一次公开features/装载/推理CPU={prep['cpu_seconds']:.6g}s。原权重训练另计CPU387.690s、wall283.448s，训练数据生成成本见原训练receipt引用的generation receipts；旧baseline费用不再次计入新花费。完整CPU/wall/LU/RHS和初值费用见RESULTS/history_warm/COST_LEDGER.json。",'']
    (ROOT/'HISTORICAL_WARM_STARTS.md').write_text('\n'.join(lines));dump(out/'SCORING_RECEIPT.json',dict(status='COMPLETE',objects=4,new_tasks=21,reused_old_predictions=len(rows)-21,truth_access='FOUR_OLD_PATHS_AFTER_ALL21_SEALED_ONLY',new42_truth_access='NONE',new_physical_calls=0,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-begin))
    print(json.dumps(dict(status='SCORED',rows=len(rows),new_tasks=21,report=str(ROOT/'HISTORICAL_WARM_STARTS.md'))),flush=True)

if __name__=='__main__':main()
