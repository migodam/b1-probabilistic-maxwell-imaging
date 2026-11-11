"""Decompose already sealed high3D selections; no Maxwell forward/adjoint.

Uses receipts, paid-state materials, existing offline metrics and saved physical
current audits. Missing baseline raw physical currents remain explicitly missing.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
from types import SimpleNamespace
from datetime import datetime,timezone
import sys,json,time
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


def quadratic_prior(x):
    image=np.asarray(x).real.reshape(6,6,6)
    d=np.concatenate([np.diff(image,axis=a).ravel() for a in range(3)])
    return float(2e-4*np.mean(d*d)+1e-6*np.mean(image*image))


def receiver_operator_only():
    """Reuse inherited external electric operator, without constructing G/LU."""
    from research.research_max_struct_r1.src.vector_backend import TensorDDA
    side=6;h=.6/side;ijk=np.indices((side,)*3).reshape(3,-1).T
    points=(ijk-(side-1)/2)*h;a=2*np.pi*np.arange(20)/20;circles=[]
    for fixed,center in [(2,[.04,0,-.02]),(0,[.02,0,-.03]),(1,[-.01,.02,0])]:
        q=np.zeros((20,3));axes=[j for j in range(3) if j!=fixed]
        q[:,axes[0]]=.6*np.cos(a);q[:,axes[1]]=.6*np.sin(a);q+=center;circles.append(q)
    # The inherited method uses only these public acquisition/grid attributes.
    spec=SimpleNamespace(points=points,volume=h**3,k=2*np.pi,n=3*len(points))
    return TensorDDA.external_electric_operator(spec,np.concatenate(circles))


def main():
    cpu,wall=time.process_time(),time.perf_counter()
    cfg=read(HERE/'FROZEN_PROTOCOL.json');seal=read(HERE/'ALL12_PREDICTIONS_SEALED.json')
    assert seal['status']=='ALL12_PREDICTIONS_AND_PREFIX_CHOICES_SEALED'
    scored=read(HERE/'OFFLINE_SCORES.json')
    metrics={(r['case_id'],r['method']):r for r in scored['rows']}
    old_prefix=read(HERE/'OFFLINE_PREFIX_SCORES.json')['rows']
    S=receiver_operator_only();receiver_columns=0
    rows=[];prefixes=[]
    for cid in cfg['cases']:
        with np.load(HERE/'data_public'/f'{cid}.npz',allow_pickle=False) as z:Y=z['Y'].copy()
        energy=np.sum(abs(Y)**2,axis=0)
        for method in cfg['methods']:
            out=HERE/'RESULTS'/cid/method;r=read(out/'RECEIPT.json')
            h=read(out/'history.json');met=metrics[(cid,method)]
            with np.load(out/'prediction.npz',allow_pickle=False) as z:
                raw=z['raw'].real.copy();final=z['prediction'].real.copy()
            with np.load(out/'paid_materials.npz',allow_pickle=False) as z:mats=z['m'].copy()
            proposal_indices=[i for i,item in enumerate(h) if item['phase']!='corrector']
            raw_index=min(proposal_indices,key=lambda i:(h[i]['phi'],i))
            final_index=min(range(len(h)),key=lambda i:(h[i]['phi'],i))
            assert np.array_equal(mats[raw_index],raw)
            assert np.array_equal(mats[final_index],final)
            raw_data=r['raw_data_residual']**2;final_data=r['physical_data_residual']**2
            raw_prior=quadratic_prior(raw);final_prior=quadratic_prior(final)
            closure=max(abs(r['raw_phi']-raw_data-raw_prior),abs(r['physical_phi']-final_data-final_prior))
            if closure>1e-13:raise ArithmeticError('Saved data/prior/Phi decomposition does not close')
            raw_iteration=None;raw_cycle=None;raw_audit=None;fixed_cost=None
            raw_noise=dict(exact_global_residual=None,raw_CN99_pass=None,
                           availability='MISSING: baseline raw physical current/prediction was not saved',
                           data_term_source_normalized=raw_data)
            if method in ('DR','AP'):
                ph=read(out/'projection/history.json');raw_iteration=ph['incumbent_iteration']
                raw_audit=next(i for i,record in enumerate(ph['history']) if record['iteration']==raw_iteration)
                fixed_cost=ph['history'][raw_audit]['costs']
                with np.load(out/'projection/states.npz',allow_pickle=False) as states:
                    saved_m=states[f'audit_{raw_audit:05d}_incumbent_material']
                    assert np.array_equal(saved_m,raw)
                    current=states[f'audit_{raw_audit:05d}_incumbent_J']
                rn=float(np.linalg.norm(S@current-Y));receiver_columns+=Y.shape[1]
                raw_noise=dict(exact_global_residual=rn,raw_CN99_pass=bool(rn<=r['delta']),
                    raw_global_delta=r['delta'],residual_over_CN99_delta=rn/r['delta'],
                    availability='S@stored physical incumbent J; no Maxwell solve',
                    data_term_source_normalized=raw_data)
            elif h[raw_index]['phase']=='auxiliary':
                raw_cycle=sum(item['phase']=='auxiliary' for item in h[:raw_index+1])
            # Bounds follow P*f0=sum ||r_p||²/||Y_p||². They do not fill in
            # an unsaved baseline exact raw residual or manufacture feasibility.
            bounds=np.sqrt(Y.shape[1]*raw_data*np.array([energy.min(),energy.max()]))
            raw_noise.update(global_residual_lower_bound=float(bounds[0]),
                global_residual_upper_bound=float(bounds[1]),
                bounds_status='guaranteed_infeasible' if bounds[0]>r['delta'] else
                              'guaranteed_feasible' if bounds[1]<=r['delta'] else 'undetermined')
            if raw_noise['raw_CN99_pass'] is not None:
                raw_noise.update(CN99_feasibility_from_available_data=raw_noise['raw_CN99_pass'],
                                 CN99_feasibility_evidence='exact S@stored physical current residual')
            else:
                certified=(False if bounds[0]>r['delta'] else True if bounds[1]<=r['delta'] else None)
                raw_noise.update(CN99_feasibility_from_available_data=certified,
                    CN99_feasibility_evidence='weighted norm lower/upper bound; exact baseline raw residual remains missing',
                    raw_global_delta=r['delta'])
            corrector_eval=sum(item['phase']=='corrector' for item in h[:final_index+1])
            ddata=final_data-raw_data;dprior=final_prior-raw_prior
            kind=('both_data_and_prior_improve' if ddata<0 and dprior<0 else
                  'data_worsens_prior_improves' if ddata>0 and dprior<0 else
                  'data_improves_prior_worsens' if ddata<0 and dprior>0 else
                  'other_or_tie')
            row=dict(case_id=cid,method=method,
                raw=dict(Phi=r['raw_phi'],data_term=raw_data,prior=raw_prior,
                    source_normalized_residual=r['raw_data_residual'],
                    relative_material_l2=met['raw_relative_material_l2'],noise=raw_noise,
                    selected_history_index=raw_index,selected_projection_iteration=raw_iteration,
                    selected_auxiliary_cycle_1based=raw_cycle,selected_projection_audit_index=raw_audit,
                    selected_paid_CPU_seconds=h[raw_index]['cpu_seconds'],
                    selected_physical_and_auxiliary_counts=h[raw_index]['counts'],
                    selected_fixed_background_counts=fixed_cost,
                    full_proposal_CPU_seconds_spent=r['proposal_cpu_seconds'],
                    full_proposal_fixed_background_counts=(r['projection'] or {}).get('counts_including_setup')),
                corrected=dict(Phi=r['physical_phi'],data_term=final_data,prior=final_prior,
                    source_normalized_residual=r['physical_data_residual'],
                    relative_material_l2=met['relative_material_l2'],
                    raw_CN99_pass=met['raw_CN99_pass'],exact_raw_global_residual=met['raw_global_residual'],
                    raw_global_delta=r['delta'],residual_over_CN99_delta=met['raw_global_residual']/r['delta'],
                    selected_history_index=final_index,selected_corrector_evaluation_1based=corrector_eval,
                    selected_accepted_iteration=None,
                    selected_accepted_iteration_availability='NOT_RECORDED; evaluationindex is exact, totalacceptednit below',
                    total_accepted_iterations=r['corrector']['nit'],total_evaluations=r['corrector']['evaluations'],
                    selected_paid_CPU_seconds=h[final_index]['cpu_seconds'],
                    selected_counts=h[final_index]['counts'],full_CPU_seconds_spent=r['cpu_seconds'],
                    corrector_CPU_seconds_spent=r['corrector']['cpu_seconds'],stop=r['corrector']['stop']),
                delta=dict(Phi=r['physical_phi']-r['raw_phi'],data_term=ddata,prior=dprior,
                    relative_material_l2=met['relative_material_l2']-met['raw_relative_material_l2'],
                    classification=kind),decomposition_closure_absolute_error=closure)
            rows.append(row)
            for limit in cfg['CPU_prefix_seconds']:
                eligible=[i for i in proposal_indices if h[i]['cpu_seconds']<=limit]
                if not eligible:continue
                index=min(eligible,key=lambda i:(h[i]['phi'],i));error=None;metric_source=None
                if np.array_equal(mats[index],raw):
                    error=met['raw_relative_material_l2'];metric_source='existing offline rawselected metric'
                else:
                    matching=[q for q in old_prefix if q['case_id']==cid and q['method']==method
                              and q['CPU_limit_seconds']==limit and q['index']==index]
                    if matching:
                        error=matching[0]['relative_material_l2'];metric_source='existing offline sealed CPU-prefix metric'
                prefixes.append(dict(case_id=cid,method=method,CPU_limit_seconds=limit,
                    phase=h[index]['phase'],selected_history_index=index,Phi=h[index]['phi'],
                    data_term=h[index]['physical_data_residual']**2,
                    prior=quadratic_prior(mats[index]),relative_material_l2=error,
                    material_metric_source=metric_source or 'MISSING: not previously scored; no fresh private read',
                    selected_paid_CPU_seconds=h[index]['cpu_seconds'],counts=h[index]['counts'],
                    full_fixed_proposal_CPU_seconds_spent=r['proposal_cpu_seconds'],
                    selection='min saved originalPhi among proposal-only states available before CPU limit, neverGT',
                    status='supplemental post-score diagnostic using predeclared Phi rule, not a new primary freeze'))
    result=dict(status='SAVED_STATE_DECOMPOSITION_NO_NEW_MAXWELL',created_utc=datetime.now(timezone.utc).isoformat(),
        data_term_formula='f0=(1/P)sum_p ||S J_p-Y_p||²/||Y_p||²; no factor1/2',
        prior_formula='2e-4 mean_neighbour_edges(chi_i-chi_j)² + 1e-6 mean_cells chi²',
        rows=rows,new_physical_LU=0,new_physical_source_RHS=0,new_physical_adjoint_RHS=0,
        stored_current_receiver_projection_columns=receiver_columns,receiver_operator_only_assemblies=1,
        new_material_truth_reads=0,scope='receipts/savedstates/existing offlinemetrics only; no optimizer or algorithm change')
    dump(HERE/'RAW_CORRECTOR_DECOMPOSITION.json',result)
    dump(HERE/'RAW_CPU_PREFIX_COMPARISON.json',dict(rows=prefixes,
        scope='proposal-only availability at common CPU ceilings, not equal expenditure; fixed proposal may finish early; compare also full spentcost',
        new_physical_LU=0,new_material_truth_reads=0))
    lines=['# 高对比度 true3D：已保存状态的 corrector 分解','',
        '本补充不改变算法、不跑新反演、不做新 Maxwell 正演/伴随，也不重新读取材料truth。只使用已保存材料、Phi、source-normalized residual、原quadraticR和已有离线材料指标。AP/DR原始物理current保存在audit里，因此可用原接收算子S读取raw残差；CIE/Dual的raw物理current未保存，精确raw全局CN99状态如实标缺。',
        '', '原目标为 Phi=D+R，D=(1/P)Σ ||r_p||²/||Y_p||²（无1/2），R=2e-4 mean_edge(chi_i-chi_j)²+1e-6 mean_cell chi²。',
        '', '|对象|方法|raw D|corrected D|raw R|corrected R|raw→corrected relL2|分解|',
        '|---|---|---:|---:|---:|---:|---|---|']
    for r in rows:
        a,b=r['raw'],r['corrected'];lines.append(f"|{r['case_id']}|{r['method']}|{a['data_term']:.8g}|{b['data_term']:.8g}|{a['prior']:.8g}|{b['prior']:.8g}|{a['relative_material_l2']:.6f}→{b['relative_material_l2']:.6f}|{r['delta']['classification']}|")
    lines+=['','|对象|方法|raw Phi|corrected Phi|raw CN99|corrected CN99|',
            '|---|---|---:|---:|---|---|']
    for r in rows:
        a,b=r['raw'],r['corrected'];feas=a['noise']['raw_CN99_pass']
        text=(f"{a['noise']['CN99_feasibility_from_available_data']} (norm界认证；exact残差MISSING)"
              if feas is None else f"{feas}, ratio={a['noise']['residual_over_CN99_delta']:.5f}")
        lines.append(f"|{r['case_id']}|{r['method']}|{a['Phi']:.9g}|{b['Phi']:.9g}|{text}|{b['raw_CN99_pass']}, ratio={b['residual_over_CN99_delta']:.5f}|")
    lines+=['','|对象|方法|raw选择迭代/周期|raw候选CPU|完整proposal实际CPU|corrected选择eval / 总acceptednit|完整CPU|',
            '|---|---|---|---:|---:|---|---:|']
    for r in rows:
        a,b=r['raw'],r['corrected'];where=(f"projection{a['selected_projection_iteration']}" if a['selected_projection_iteration'] is not None else f"auxcycle{a['selected_auxiliary_cycle_1based']}")
        lines.append(f"|{r['case_id']}|{r['method']}|{where}|{a['selected_paid_CPU_seconds']:.3f}|{a['full_proposal_CPU_seconds_spent']:.3f}|{b['selected_corrector_evaluation_1based']} / {b['total_accepted_iterations']}|{b['full_CPU_seconds_spent']:.3f}|")
    lines+=['', 'CIE/Dual的exact raw global residual缺失，但其CN99不可行可由已知source-normalized D严格推出：P D min_p||Y_p||² ≤ ||r||² ≤ P D max_p||Y_p||²。六个baseline raw下界均已大于delta，故无需新forward即可认证球外；精确残差仍不补造。',
        '', 'selected candidate出现时CPU与跑完整个固定proposal实际CPU是不同数字，均保留。corrector所选paid evaluation序号准确；该evaluation对应accepted iteration未保存，不能臆造。每项已付费物理/aux计数及AP/DR fixed-H计数见JSON。',
        '', '同CPU前缀的raw-only比较（只考虑proposal，按原Phi选，不按材料误差选）：',
        '', '|对象|CPU ceiling|方法|已付费候选CPU|Phi|prior|材料relL2|',
        '|---|---:|---|---:|---:|---:|---:|']
    for p in prefixes:
        error='MISSING' if p['relative_material_l2'] is None else f"{p['relative_material_l2']:.6f}"
        lines.append(f"|{p['case_id']}|{p['CPU_limit_seconds']}|{p['method']}|{p['selected_paid_CPU_seconds']:.3f}|{p['Phi']:.8g}|{p['prior']:.8g}|{error}|")
    lines+=['', '上述是相同CPU上限下已可获得的固定proposal状态，未花完的时间不虚计为开销；已完成CIE/Dual后其raw状态仍可取得。完整proposal费用与此前包含corrector的预封存前缀另列。该后续诊断不产生新的独立确认或新GT选择。',
        '', f"新增成本：Maxwell LU/sourceRHS/adjointRHS全为0；仅原外部接收算子装配1次和对已保存current做{receiver_columns}个接收列动作。原始数据、预测与seal未改写。"]
    (HERE/'RAW_CORRECTOR_DECOMPOSITION_ZH.md').write_text('\n'.join(lines)+'\n')
    dump(HERE/'DECOMPOSITION_RECEIPT.json',dict(status='COMPLETE_SAVED_STATE_ONLY',
        cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        new_Maxwell_LU=0,new_source_RHS=0,new_adjoint_RHS=0,new_material_truth_reads=0,
        receiver_projection_columns=receiver_columns,external_receiver_operator_assemblies=1))
    print(json.dumps(dict(status=result['status'],DR_rows=[r for r in rows if r['method']=='DR'])),flush=True)


if __name__=='__main__':main()
