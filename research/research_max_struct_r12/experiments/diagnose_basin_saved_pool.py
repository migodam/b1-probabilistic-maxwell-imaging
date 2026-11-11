"""Authorized post-seal material oracle of saved basin pools; zero FW/optimizer.

The oracle is an OFFLINE diagnostic and never a prediction or stopping rule.
"""
from pathlib import Path
from datetime import datetime,timezone
import sys,json,time

ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
BASE=ROOT/'RESULTS/basin_remote'
sys.path.insert(0,str(REPO))


def read(path):return json.loads(Path(path).read_text())


def dump(path,value):
    def convert(x):
        import numpy as np
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value,indent=2,default=convert,allow_nan=False)+'\n')


def main():
    cpu,wall=time.process_time(),time.perf_counter()
    import numpy as np
    from research.research_max_struct_r12.experiments.score_basin import metrics
    cases=['r10_twod_0001','r10_twod_0004','r10_twod_0022']
    methods=['direct','multistart','alternating_projection','douglas_rachford']
    seal=read(BASE/'PREFIX_SEAL.json')
    assert len({r['run'] for r in seal['records']})==12
    assert all((BASE/f'{cid}__{method}'/'RECEIPT.json').is_file() for cid in cases for method in methods)
    config=dict(status='ROOT_AUTHORIZED_POSTSEAL_EXPOSED_POOL_DIAGNOSTIC',
        created_utc=datetime.now(timezone.utc).isoformat(),cases=cases,methods=methods,
        maximum_body_CPU_seconds=240.,CPU_prefixes=[30.,60.,120.,240.],
        K_prefixes=[1,10,50,100,500,1000],
        online_rule='minimum saved originalPhi, first recorded candidate resolves ties',
        oracle_rule='minimum offline materialNRMSE within already-paid eligible pool; tieoriginalPhi thenindex',
        oracle_scope='diagnostic only, never online prediction/earlystop/newalgorithm',
        metric_source='unchanged experiments/score_basin.py::metrics and exactR10 truth mapping',
        new_forward=0,new_optimization=0)
    dump(BASE/'SAVED_POOL_DIAGNOSTIC_CONFIG.json',config)
    official=read(BASE/'OFFLINE_SCORES.json')['rows']
    truthbase=REPO.parent/'r10_worktree/research/research_max_struct_r10/sealed_truth/twod'
    rows=[];prefixrows=[];krows=[];coverage=[];oracle_arrays={};truth_reads=0
    for cid in cases:
        with np.load(truthbase/f'{cid}.npz',allow_pickle=False) as z:truth=z['m'].copy()
        truth_reads+=1
        for method in methods:
            folder=BASE/f'{cid}__{method}'
            if not (folder/'prediction.npz').is_file() or not (folder/'history.json').is_file():
                coverage.append(dict(case=cid,method=method,status='MISSING_ARRAY_OR_HISTORY'));continue
            h=read(folder/'history.json')
            with np.load(folder/'prediction.npz',allow_pickle=False) as z:
                if 'candidates' not in z.files:
                    coverage.append(dict(case=cid,method=method,status='MISSING_CANDIDATE_ARRAY'));continue
                mats=z['candidates'].copy()
            ids=[r.get('candidate') for r in h]
            aligned=(mats.ndim==2 and mats.shape[1]==576 and len(h)==len(mats)
                     and ids==list(range(len(h))))
            if not aligned:
                coverage.append(dict(case=cid,method=method,status='METADATA_ARRAY_MISMATCH',
                    history_count=len(h),array_shape=list(mats.shape)));continue
            if any(not np.isfinite(mats[i]).all() or not np.isfinite(r['Phi']) or not np.isfinite(r['cpu']) for i,r in enumerate(h)):
                coverage.append(dict(case=cid,method=method,status='NONFINITE_SAVED_PAYLOAD'));continue
            metric_rows=[metrics(m,truth) for m in mats]
            eligible=[i for i,r in enumerate(h) if r['cpu']<=240.]
            if not eligible:
                coverage.append(dict(case=cid,method=method,status='NO_PAID_CANDIDATE_UNDER240'));continue
            def chosen(indices):
                online=min(indices,key=lambda i:(h[i]['Phi'],i))
                oracle=min(indices,key=lambda i:(metric_rows[i]['NRMSE'],h[i]['Phi'],i))
                def item(i):
                    return dict(candidate=i,body_CPU_seconds=h[i]['cpu'],Phi=h[i]['Phi'],
                        saved_relative_data_residual=h[i].get('relative_data_residual'),
                        saved_prior=h[i].get('prior'),counts=h[i].get('cost'),**metric_rows[i])
                return online,oracle,item(online),item(oracle)
            online,oracle,on,ora=chosen(eligible)
            with np.load(folder/'paid_prefix.npz',allow_pickle=False) as z:
                match=np.array_equal(z['m_240'],mats[online]) if 'm_240' in z.files else None
            assert match is not False,'Online pool selection differs from sealed240prefix'
            score=next(r for r in official if r['case']==cid and r['method']==method and r['view']=='paid_CPU_240')
            assert abs(score['NRMSE']-on['NRMSE'])<1e-13
            row=dict(case=cid,method=method,saved_K=len(h),eligible_K=len(eligible),
                excluded_over240=len(h)-len(eligible),online=on,offline_oracle=ora,
                online_matches_sealed_240prefix=match,
                oracle_is_prediction=False,K_includes_repeated_paid_states=True)
            rows.append(row);coverage.append(dict(case=cid,method=method,status='FULL_INDEX_ALIGNMENT',
                saved_K=len(h),eligible_K=len(eligible),online_matches_sealed_prefix=match))
            key=f'{cid}__{method}';oracle_arrays[key+'__OFFLINE_ORACLE']=mats[oracle]
            oracle_arrays[key+'__ONLINE_MIN_PHI']=mats[online]
            for limit in config['CPU_prefixes']:
                ids0=[i for i in eligible if h[i]['cpu']<=limit]
                if not ids0:continue
                _,_,a,b=chosen(ids0);prefixrows.append(dict(case=cid,method=method,
                    maximum_body_CPU_seconds=limit,K=len(ids0),online=a,offline_oracle=b,oracle_is_prediction=False))
            for K in sorted(set([k for k in config['K_prefixes'] if k<=len(eligible)]+[len(eligible)])):
                _,_,a,b=chosen(eligible[:K]);krows.append(dict(case=cid,method=method,
                    K=K,maximum_body_CPU_seconds=240.,online=a,offline_oracle=b,
                    last_considered_candidate_CPU=h[eligible[K-1]]['cpu'],oracle_is_prediction=False))
    comparisons=[]
    for cid in cases:
        subset=[r for r in rows if r['case']==cid]
        dr=next((r for r in subset if r['method']=='douglas_rachford'),None)
        if dr is None:continue
        for r in subset:
            if r['method']=='douglas_rachford':continue
            comparisons.append(dict(case=cid,control=r['method'],
                control_online_NRMSE=r['online']['NRMSE'],control_oracle_NRMSE=r['offline_oracle']['NRMSE'],
                DR_online_NRMSE=dr['online']['NRMSE'],DR_oracle_NRMSE=dr['offline_oracle']['NRMSE'],
                control_pool_contains_better_material_than_DR_online=bool(r['offline_oracle']['NRMSE']<dr['online']['NRMSE']),
                oracle_is_prediction=False))
    result=dict(status='POSTSEAL_OFFLINE_POOL_ORACLE_COMPLETE',rows=rows,
        comparisons=comparisons,coverage=coverage,new_forward=0,new_optimization=0,
        GT_reads='authorized already-exposed R10 material truth after all12seal, for diagnostics only',
        limits='Only savedphysicalaudits/callbacks. Missing inneriterates/unauditedshadows are not reconstructed. Repeats are not independentdraws. No unboundedrun claims.')
    dump(BASE/'SAVED_POOL_ORACLE_DIAGNOSTIC.json',result)
    dump(BASE/'SAVED_POOL_CPU_PREFIXES.json',dict(rows=prefixrows,oracle_is_prediction=False,new_forward=0))
    dump(BASE/'SAVED_POOL_BEST_AT_K.json',dict(rows=krows,oracle_is_prediction=False,new_forward=0,
        K_definition='number of recorded cap-eligible paid states in chronological order, including repeats; notequalcomputebudget'))
    np.savez_compressed(BASE/'OFFLINE_ORACLE_DIAGNOSTIC_MATERIALS.npz',**oracle_arrays)
    lines=['# 已封存basin候选池：online与离线oracle','',
        '本项为root明确授权的封存后GT诊断，0新forward、0优化。oracle只描述已付费、已保存、bodyCPU≤240的有限池，不是预测、early-stop或新方法。材料/归一化/IoU完全复用原score_basin.py，online仍取最小原Phi。',
        '', '|对象|方法|savedK / eligibleK|online NRMSE|offline best@K NRMSE|oracle候选CPU|oracle数据残差|',
        '|---|---|---:|---:|---:|---:|---:|']
    for r in rows:
        o=r['offline_oracle'];lines.append(f"|{r['case']}|{r['method']}|{r['saved_K']}/{r['eligible_K']}|{r['online']['NRMSE']:.6f}|{o['NRMSE']:.6f}|{o['body_CPU_seconds']:.3f}|{o['saved_relative_data_residual']:.6g}|")
    lines+=['', '|对象|控制|控制pool oracle|DR online|控制pool是否含比DR online材料误差更小的已保存点|',
            '|---|---|---:|---:|---|']
    for r in comparisons:lines.append(f"|{r['case']}|{r['control']}|{r['control_oracle_NRMSE']:.6f}|{r['DR_online_NRMSE']:.6f}|{r['control_pool_contains_better_material_than_DR_online']}|")
    lines+=['', '所有12臂都核查候选行号与history一一对应，online材料与已封存240前缀逐元素一致；任何缺失/错位会列于coverage而不补跑。K含重复已付callback，不是独立起点数。raw未审计inner/shadow没有恢复，不能据有限池断言方法从未到过某个区域。',
        '', '30/60/120/240 bodyCPU前缀和firstK视图分别保存在SAVED_POOL_CPU_PREFIXES及SAVED_POOL_BEST_AT_K。两者均按公开Phi产生online结果，oracle只在离线另列，不与在线胜率混合。',
        '', '结果应区分：控制的Phi-selected估计较差；控制已保存pool材料oracle是否也较差；DR自己的pool是否含被Phi选择遗漏的更好材料。即使有限池oracle仍差，也不能推出全部未保存迭代或其他restart分布都失败。完整运行超过240的候选均排除，完整费用仍保留。']
    (BASE/'SAVED_POOL_ORACLE_DIAGNOSTIC_ZH.md').write_text('\n'.join(lines)+'\n')
    dump(BASE/'SAVED_POOL_DIAGNOSTIC_RECEIPT.json',dict(status='COMPLETE_0FW_POSTSEAL',
        cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        truth_files_read=truth_reads,new_forward=0,new_Maxwell_LU=0,new_source_RHS=0,new_adjoint_RHS=0,
        new_optimization=0,coverage=len(rows),scope='newdiagnostic only; originalpredictions/seal/scores unchanged'))
    print(json.dumps(dict(rows=[dict(case=r['case'],method=r['method'],K=r['eligible_K'],
        online=r['online']['NRMSE'],oracle=r['offline_oracle']['NRMSE']) for r in rows],comparisons=comparisons)),flush=True)


if __name__=='__main__':main()
