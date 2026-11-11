"""Offline post-complete-seal material scoring; no physics or optimization."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
from datetime import datetime,timezone
import json
import sys
import csv
import time

HERE=Path(__file__).resolve().parent
STUDY=HERE.parents[1];REPO=STUDY.parents[1]
R8=REPO.parent/'r8_worktree/research/research_max_struct_r8'
sys.path.insert(0,str(REPO))


def read(path):return json.loads(Path(path).read_text())


def main():
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from research.research_max_struct_r7.experiments.score_primary import material_metrics
    wall,cpu=time.perf_counter(),time.process_time()
    protocol=read(HERE/'FROZEN_PROTOCOL.json');seal=read(HERE/'ALL_PREDICTIONS_SEALED.json')
    assert seal['status']=='ALL16_PREDICTIONS_AND4_INITIALIZATIONS_SEALED'
    expected={(i,a) for i in protocol['ids'] for a in protocol['arms']}
    assert {(e['case_id'],e['arm']) for e in seal['entries']}==expected
    for e in seal['entries']:
        receipt=read(e['receipt'])
        assert receipt['status']=='SEALED_PUBLIC_PREDICTION'
        assert Path(e['prediction']).is_file()
    # First private access starts here, after all16 predictions above are sealed.
    first_private_access_utc=datetime.now(timezone.utc).isoformat()
    metadata={r['case_id']:r for r in read(R8/'sealed_truth/warm_test/score_metadata.json')['records']}
    preprocess=read(HERE/'PREPROCESS_RECEIPT.json')
    cost_by_case={r['case_id']:r for r in preprocess['records']}
    rows=[];initials=[];stages=[];prefixes=[];pictures={};pairs=[];dr_audits=[]
    for ident in protocol['ids']:
        with np.load(R8/'sealed_truth/warm_test'/f'{ident}.npz',allow_pickle=False) as z:
            truth=z['m'].ravel().copy();support=z['mask'].ravel().astype(bool)
        amplitude=metadata[ident]['amplitude'];pictures[(ident,'truth')]=truth
        with np.load(HERE/'public'/f'{ident}_nn.npz',allow_pickle=False) as z:mnn=z['m_init'].copy()
        initial=material_metrics(mnn,truth,support,amplitude)
        initials.append(dict(case_id=ident,amplitude=amplitude,**initial));pictures[(ident,'nn')]=mnn
        case_rows={};case_pre={}
        for arm in protocol['arms']:
            out=HERE/'runs'/f'{ident}__{arm}';rec=read(out/'RECEIPT.json')
            with np.load(out/'prediction.npz',allow_pickle=False) as z:m=z['m'].copy()
            with np.load(out/'selected_before_corrector.npz',allow_pickle=False) as z:mpre=z['m'].copy()
            metrics=material_metrics(m,truth,support,amplitude)
            pr=cost_by_case[ident]
            warm_prefix=(pr['features']['cpu_seconds']+pr['inference']['cpu_seconds']) if arm!='uniform_direct' else 0.
            cold_extra=(preprocess['model_load']['cpu_seconds']+pr['feature_setup']['cpu_seconds']) if arm!='uniform_direct' else 0.
            row=dict(case_id=ident,arm=arm,amplitude=amplitude,**metrics,
                initial_material_nrmse=initial['material_nrmse'] if arm!='uniform_direct' else None,
                objective=rec['objective'],initial_objective=rec['initial_objective'],
                precorrector_objective=rec['precorrector_objective'],relative_residual=rec['relative_residual'],
                noise_feasible=rec['noise_feasible'],prior=rec['prior'],
                solver_cpu_seconds=rec['cpu_seconds'],solver_wall_seconds=rec['wall_seconds'],
                warm_deployment_cpu_seconds=rec['cpu_seconds']+warm_prefix,
                cold_perobject_cpu_seconds=rec['cpu_seconds']+warm_prefix+cold_extra,
                proposal_cpu_seconds=rec['proposal_cpu_seconds'],
                corrector_cpu_seconds=rec['corrector']['cpu_seconds'],
                corrector_stop=rec['corrector']['stop'],corrector_iterations=rec['corrector']['accepted_iterations'],
                corrector_evaluations=rec['corrector']['evaluations'],physical_cost=rec['physical_cost'],
                truth_access='OFFLINE_AFTER_COMPLETE16_PUBLIC_SEAL')
            rows.append(row);case_rows[arm]=row;pictures[(ident,arm)]=m
            pre_metrics=material_metrics(mpre,truth,support,amplitude);case_pre[arm]=pre_metrics
            stages.append(dict(case_id=ident,arm=arm,phase='selected_before_corrector',
                               objective=rec['precorrector_objective'],**pre_metrics))
            if arm=='nn_dr':pictures[(ident,'dr_before_corrector')]=mpre
            history=read(out/'history.json')
            with np.load(out/'paid_materials.npz',allow_pickle=False) as z:mats=z['m']
            assert len(history)==len(mats)
            if arm=='nn_dr':
                for j,item in enumerate(history):
                    if item['phase']!='dr_proposal':continue
                    dr_audits.append(dict(case_id=ident,evaluation=j,
                        cpu_seconds=item['cpu_seconds'],objective=item['objective'],
                        incumbent_objective=item['incumbent'],relative_residual=item['relative_residual'],
                        **material_metrics(mats[j],truth,support,amplitude)))
            for limit in (10.,30.,60.,120.):
                eligible=[(r['objective'],j) for j,r in enumerate(history) if r['cpu_seconds']<=limit]
                if not eligible:continue
                _,j=min(eligible)
                prefixes.append(dict(case_id=ident,arm=arm,cpu_limit=limit,
                    paid_snapshot_cpu=history[j]['cpu_seconds'],phase=history[j]['phase'],
                    objective=history[j]['objective'],physical_cost=history[j]['physics_cost'],
                    **material_metrics(mats[j],truth,support,amplitude)))
        direct=case_rows['nn_direct'];dr=case_rows['nn_dr'];cie=case_rows['nn_cie_single']
        pairs.append(dict(case_id=ident,
            DR_minus_NNdirect_material_nrmse=dr['material_nrmse']-direct['material_nrmse'],
            DR_minus_NNdirect_Phi=dr['objective']-direct['objective'],
            DR_to_NNdirect_solver_cpu_ratio=dr['solver_cpu_seconds']/direct['solver_cpu_seconds'],
            DR_minus_initialNN_material_nrmse=dr['material_nrmse']-initial['material_nrmse'],
            DR_pre_minus_initialNN_material_nrmse=case_pre['nn_dr']['material_nrmse']-initial['material_nrmse'],
            CIE_pre_minus_initialNN_material_nrmse=case_pre['nn_cie_single']['material_nrmse']-initial['material_nrmse'],
            CIE_minus_NNdirect_material_nrmse=cie['material_nrmse']-direct['material_nrmse']))
    report=dict(status='SCORED_AFTER_ALL16_SEALED',created_utc=datetime.now(timezone.utc).isoformat(),
        seal_created_utc=seal['created_utc'],first_private_access_utc=first_private_access_utc,
        scope=protocol['scope'],material_metric='R8/R7 amplitude-normalizedRMSE + prediction-onlyOtsuIoU; noGTselection',
        rows=rows,initializations=initials,precorrector=stages,paired= pairs,
        DR_better_than_NNdirect_material_cases=sum(p['DR_minus_NNdirect_material_nrmse']<0 for p in pairs),
        DR_worse_than_NNdirect_material_cases=sum(p['DR_minus_NNdirect_material_nrmse']>0 for p in pairs),
        DR_better_than_initialNN_material_cases=sum(p['DR_minus_initialNN_material_nrmse']<0 for p in pairs),
        mean_material_nrmse={arm:float(np.mean([r['material_nrmse'] for r in rows if r['arm']==arm])) for arm in protocol['arms']},
        mean_initialNN_material_nrmse=float(np.mean([r['material_nrmse'] for r in initials])),
        new_scoring_physical_LU=0,new_scoring_physical_RHS=0,
        cost_scope='solver times include publicload/modelsetup/proposal/projectionarchive/physicalaudits/corrector/output; moduleimports beforeexecute excluded. PREPROCESS total includesTorchimports/load/features/inference sharedonce. cold_perobject addsmodelweightload+featureoperatorsetup, excludes sharedmoduleimports. Nestedcounters/timersoverlap.',
        equal_cost_scope='Samecorrector allowance; DR/CIEproposal extras measured. Prefixes use lowestPhi already-paid state beforeCPUbudget, notGTranking.',
        historical_reference='OriginalR8 NN-CIE uncappedWindows endpoints not treated as equalcost; allfour CIEcontrols rerun here',
        diagnostic_scope='Four exposedcases cannot establish universality or independentNNprior robustness; no deliberatelywrongprior arm was run')
    (HERE/'OFFLINE_SCORES.json').write_text(json.dumps(report,indent=2)+'\n')
    (HERE/'OFFLINE_CPU_PREFIXES.json').write_text(json.dumps(dict(rows=prefixes,new_physical_LU=0,new_physical_RHS=0),indent=2)+'\n')
    (HERE/'OFFLINE_DR_AUDITS.json').write_text(json.dumps(dict(rows=dr_audits,
        scope='Offline candidate material metrics, never used for online selection',
        new_physical_LU=0,new_physical_RHS=0),indent=2)+'\n')
    flat=[{k:v for k,v in row.items() if not isinstance(v,(dict,list))} for row in rows]
    with (HERE/'OFFLINE_SCORES.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in flat for k in r)));writer.writeheader();writer.writerows(flat)
    labels=[('truth','GT (offline)'),('nn','Frozen NN'),('uniform_direct','Uniform → direct'),
            ('nn_direct','NN → direct'),('dr_before_corrector','NN → DR (selected)'),
            ('nn_dr','NN → DR + corrector'),('nn_cie_single','NN → CIE + corrector')]
    fig,axes=plt.subplots(4,7,figsize=(21,12),layout='constrained')
    for i,ident in enumerate(protocol['ids']):
        scale=max(float(np.max(pictures[(ident,key)])) for key,_ in labels)
        for j,(key,label) in enumerate(labels):
            im=axes[i,j].imshow(pictures[(ident,key)].reshape(24,24),vmin=.02,vmax=scale,cmap='viridis')
            title=label
            if key not in ('truth','nn','dr_before_corrector'):
                row=next(r for r in rows if r['case_id']==ident and r['arm']==key)
                title+=f"\nNRMSE {row['material_nrmse']:.3f}, IoU {row['otsu_iou']:.2f}"
            elif key=='nn':
                row=next(r for r in initials if r['case_id']==ident)
                title+=f"\nNRMSE {row['material_nrmse']:.3f}"
            axes[i,j].set_title(title,fontsize=9);axes[i,j].set_xticks([]);axes[i,j].set_yticks([])
            if j==0:axes[i,j].set_ylabel(ident,fontsize=10)
        fig.colorbar(im,ax=list(axes[i,:]),fraction=.012,pad=.006,label='m (shared row scale)')
    fig.suptitle('Four fixed exposed R8 cases: frozen NN information through direct, DR and CIE',fontsize=14)
    fig.savefig(HERE/'materials.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for ax,ident in zip(axes.ravel(),protocol['ids']):
        for arm in protocol['arms']:
            hist=read(HERE/'runs'/f'{ident}__{arm}'/'history.json')
            ax.step([r['cpu_seconds'] for r in hist],[r['incumbent'] for r in hist],where='post',label=arm)
        ax.set_yscale('log')
        ax.set_title(ident);ax.set_xlabel('CPU seconds, including proposal');ax.set_ylabel('Best paid original Phi');ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    fig.savefig(HERE/'objective_cpu.png',dpi=160);plt.close(fig)
    lines=['# R12 冻结 NN 接口四对象比较','',
           '四个 ID 在新推断/新反演前冻结，全部 16 个预测及四个 NN 初值封存后才打开 truth。无训练、无在线 GT、无物理 backend 替换。',
           '','实际 R8 ID 为 0001、0004、0022、0027；0001 为 6 源，其余 12 源。0027 是既有暴露的 heldout-family 反例。它们不是 R11 同编号对象，也不是独立确认集。',
           '', '原 physics NN（epoch86）及原十通道 features 的四个输出与历史初值完全一致。所有比较使用原 [.02,20] box、原 Maxwell Phi 和原 smooth/TV prior。',
           '', 'DR 使用完整物理 NN 双场开始，scale6/omega1，3000 步，每50步物理审计，global CN99 数据球及 grouped_ports=True。CIE 为原 R8 cie_single/原DPCS、Fourier和3000步调度；不初始化物理 current，保留原接口。所有四臂使用同一个最多3000迭代、120CPU秒的原 corrector；proposal 费用另外实测，不声称总费用相同。',
           '', '|对象|方法|幅度NRMSE|IoU|原Phi|solverCPU秒|corrector停止|', '|---|---|---:|---:|---:|---:|---|']
    for row in rows:
        lines.append(f"|{row['case_id']}|{row['arm']}|{row['material_nrmse']:.6f}|{row['otsu_iou']:.4f}|{row['objective']:.7g}|{row['solver_cpu_seconds']:.2f}|{row['corrector_stop']}|")
    lines+=['',f"DR 相对 NN→direct 的材料幅度NRMSE：改善 {report['DR_better_than_NNdirect_material_cases']}/4，变差 {report['DR_worse_than_NNdirect_material_cases']}/4。四臂均保留原物理 Phi 的 incumbent，材料标签从未参与选择。",'',
            '|对象|DR−NNdirect NRMSE|DR−NNdirect Phi|DR/NNdirect CPU|DR−初始NN NRMSE|',
            '|---|---:|---:|---:|---:|']
    for p in pairs:lines.append(f"|{p['case_id']}|{p['DR_minus_NNdirect_material_nrmse']:.6f}|{p['DR_minus_NNdirect_Phi']:.7g}|{p['DR_to_NNdirect_solver_cpu_ratio']:.3f}|{p['DR_minus_initialNN_material_nrmse']:.6f}|")
    lines+=['','投影本身的 NN 保留诊断（共同 corrector 之前）：',
            '', '|对象|初始 NN NRMSE|DR最低Phi候选 NRMSE|CIE最低Phi候选 NRMSE|',
            '|---|---:|---:|---:|']
    for ident in protocol['ids']:
        nr=next(r['material_nrmse'] for r in initials if r['case_id']==ident)
        drpre=next(r['material_nrmse'] for r in stages if r['case_id']==ident and r['arm']=='nn_dr')
        ciepre=next(r['material_nrmse'] for r in stages if r['case_id']==ident and r['arm']=='nn_cie_single')
        lines.append(f'|{ident}|{nr:.6f}|{drpre:.6f}|{ciepre:.6f}|')
    lines+=['','成本与解释边界：PREPROCESS_RECEIPT 给出实际 feature、NN加载和推断费用；各 RECEIPT 的 physical_cost 与 projection/auxiliary_cost 分开保留。OFFLINE_CPU_PREFIXES 按已付费原Phi候选给10/30/60/120CPU秒前缀，不按GT挑选。材料误差下降与物理残差下降分别报告。',
            '', '该诊断只回答四个暴露对象上 NN 信息通过不同接口后的保留/损伤。未运行人为错误 prior 的新臂，故不能推出任意错误 prior 的可靠性。科学主张与 novelty 由主线程独立评审。',
            '', '复现：`run_nn_interface.py prepare`（会重写协议/初值，仅新目录复现），随后 `run_nn_interface.py run`，整组 seal 后运行 `score_nn_interface.py`。原历史来源只读；无 push。']
    (HERE/'RESULTS_ZH.md').write_text('\n'.join(lines)+'\n')
    (HERE/'SCORING_RECEIPT.json').write_text(json.dumps(dict(status='COMPLETE_POSTSEAL_OFFLINE_SCORING',
        cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,new_physical_LU=0,new_physical_RHS=0),indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],mean_material_nrmse=report['mean_material_nrmse'],paired=pairs)),flush=True)


if __name__=='__main__':main()
