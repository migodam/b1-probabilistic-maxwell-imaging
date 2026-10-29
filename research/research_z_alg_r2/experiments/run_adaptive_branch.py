"""Frozen data-only anchored slice study with one exact global Schur kernel."""
import os
for _variable in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_variable]='1'
import sys,time,json,itertools
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from scipy.optimize import least_squares
P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P/'src'));sys.path.insert(0,str(P/'experiments'))
from physics import Model,ReducedKernel
from validate import CountingBook,load_case
from branch_enumerator import enumerate_bilinear_dual,filter_bilinear_dual_roots
from run_branch_study import write,starts,common_validate,qualified

CASES=['E1_LOW6_00','E1_LOW6_01','E1_LOW6_02']
OFFSETS=[-.6,-.3,0.,.3,.6]

def read_anchor(case):
    path=P/'results'/f'blind_{case}_frozen_J_TRF.json'
    row=json.loads(path.read_text())
    if not row['validation']['response_fit_pass']:raise RuntimeError('Public data-only anchor invalid')
    return dict(path=str(path.relative_to(P)),m=np.asarray(row['m']),validation=row['validation'],
                original_cost={k:row['accounting'][k] for k in ['total_complex_rhs','wall_seconds','CPU_seconds']})

def grid(anchor):
    unique=[];nominal=[]
    for a,b in itertools.product(OFFSETS,repeat=2):
        free=np.clip(anchor['m'][4:]+[a,b],.5,4.)
        match=next((i for i,old in enumerate(unique) if np.array_equal(free,old)),None)
        if match is None:match=len(unique);unique.append(free)
        nominal.append(dict(offset=[a,b],free=free,unique_slice=match))
    return dict(unique_free_pairs=unique,nominal_grid=nominal)

def baseline_starts(anchor):
    result=[]
    for candidate in [anchor['m'][:4],*starts()]:
        if all(np.linalg.norm(candidate-old)>1e-12 for old in result):result.append(candidate.copy())
        if len(result)==8:break
    return result

def setup_contract():
    anchors={case:read_anchor(case) for case in CASES}
    return dict(schema='Z-ALG-R2-ADAPTIVE-BRANCH-v1',frozen_utc=datetime.now(timezone.utc).isoformat(),
        cases=CASES,anchors=anchors,grids={case:grid(anchor) for case,anchor in anchors.items()},
        fixed_free_indices=[4,5],active=[0,1,2,3],offsets=OFFSETS,box=[.5,4.],
        grid_policy='75 nominal clipped points;identical clipped free pairs run once;70 distinct slices',
        anchor_policy='Existing public frozen-J predictions only;validated observation fits;no scores/private',
        setup=dict(reference_material=[2.25]*6,global_setups=1,actual_complex_RHS=8,
                   reuse='one6-material kernel shared across3 observations and every slice;F11 Schur restriction'),
        six=dict(seed=2026102904,tolerance=2e-11,wall_seconds=60),
        baseline=dict(starts={case:baseline_starts(a) for case,a in anchors.items()},
                      start_policy='Data-only anchor active coords,midpoint,and first6 other distinct frozen Sobol starts',
                      max_nfev=120,ftol=1e-12,xtol=1e-12,gtol=1e-12,x_scale=1.,loss='linear'),
        validation=dict(response_relative=1e-8,state_relative=1e-10,sigma_min=1e-6,merge_RMS=1e-7,
                        independent='fresh full nonlinear numpy solve at every unique legal candidate'),
        accounting=dict(anchor_cost_charged_to_both=True,global_setup_charged_to_six_once=True,
                        reduced_schur_solve_separate=True,single_thread=True),
        scope='Material coverage on common data-adapted slices;no connected-component,global-fiber,posterior,or rigorous completeness claim')

def six_case(model,global_kernel,case,y,slice_id,fixed):
    book=CountingBook();model.book=book;before=model.counts.copy()
    output=dict(case=case,slice=slice_id,method='six_shared_Schur',fixed_material=fixed.copy())
    try:
        kernel=global_kernel.restrict([0,1,2,3],fixed)
        result=enumerate_bilinear_dual(kernel,y,wall_seconds=60)
        filtered=filter_bilinear_dual_roots(result,validate=False)
        system=result['system']
        output.update(status='COMPLETE',homotopy={k:v for k,v in result.items() if k!='system'},
            physical_filter=filtered,independently_validated=common_validate(model,y,filtered['roots']),
            kernel={k:getattr(kernel,k) for k in ['m0','active','K','L','q','y0','restriction']},
            dual_system={k:getattr(system,k) for k in ['w0','N','W','V','coefficients','scale']})
    except BaseException as exc:output.update(status='FAILED_PAID',error=repr(exc))
    output['accounting']=book.summary(model.ns,model.nr)
    output['model_counts']=dict(model.counts-before)
    return output

def baseline_case(case,g,y,slice_id,fixed,initials):
    book=CountingBook();model=Model(book,**g);target=model.pack(y);scale=max(np.linalg.norm(target),1e-300)
    output=dict(case=case,slice=slice_id,method='FullJ8_anchor',fixed_material=fixed.copy(),starts=[])
    candidates=[]
    try:
        for number,x0 in enumerate(initials):
            cache={}
            def state(x):
                if 'x' not in cache or not np.array_equal(x,cache['x']):
                    m=fixed.copy();m[:4]=x;cache.update(x=x.copy(),s=model.forward(m))
                return cache['s']
            def fun(x):return (model.pack(state(x)['y'])-target)/scale
            def jac(x):return model.jacobian(state(x))[:,:4]/scale
            fit=least_squares(fun,x0,jac=jac,method='trf',bounds=(.5,4.),max_nfev=120,
                              ftol=1e-12,xtol=1e-12,gtol=1e-12,x_scale=1.,loss='linear')
            m=fixed.copy();m[:4]=fit.x;residual=float(np.linalg.norm(fit.fun))
            output['starts'].append(dict(start=number,x0=x0,m=m,response_relative=residual,
                scipy_success=bool(fit.success),scipy_status=int(fit.status),message=fit.message,
                nfev=fit.nfev,njev=fit.njev,root_candidate=residual<=1e-8))
            if residual<=1e-8:candidates.append(dict(status='CANDIDATE',m=m,path=number))
        output.update(status='COMPLETE',independently_validated=common_validate(model,y,candidates))
    except BaseException as exc:output.update(status='FAILED_PAID',error=repr(exc))
    output.update(accounting=book.summary(model.ns,model.nr),model_counts=dict(model.counts))
    return output

def summarize(outputs,contract,global_setup):
    anchor_rhs=sum(a['original_cost']['total_complex_rhs'] for a in contract['anchors'].values())
    anchor_wall=sum(a['original_cost']['wall_seconds'] for a in contract['anchors'].values())
    anchor_cpu=sum(a['original_cost']['CPU_seconds'] for a in contract['anchors'].values())
    totals={}
    for method in ['six_shared_Schur','FullJ8_anchor']:
        rows=[r for r in outputs if r['method']==method]
        rhs=sum(r['accounting']['total_complex_rhs'] for r in rows)
        wall=sum(r['accounting']['wall_seconds'] for r in rows)
        cpu=sum(r['accounting']['CPU_seconds'] for r in rows)
        setup=global_setup if method=='six_shared_Schur' else dict(total_complex_rhs=0,wall_seconds=0,CPU_seconds=0)
        totals[method]=dict(slices=len(rows),complete=sum(r['status']=='COMPLETE' for r in rows),
            covered_slices=sum(bool(qualified(r)) for r in rows),qualified_materials=sum(len(qualified(r)) for r in rows),
            anchor_reused_cost_RHS=anchor_rhs,anchor_reused_wall_seconds=anchor_wall,
            incremental_RHS=rhs+setup['total_complex_rhs'],
            total_RHS=rhs+setup['total_complex_rhs']+anchor_rhs,
            total_wall_seconds=wall+setup['wall_seconds']+anchor_wall,
            total_CPU_seconds=cpu+setup['CPU_seconds']+anchor_cpu,
            resolved_paths=sum(r.get('homotopy',{}).get('resolved_paths',0) for r in rows),
            unresolved_paths=sum(r.get('homotopy',{}).get('unresolved_paths',0) for r in rows),
            collision_slices=sum(not r['homotopy']['all_paths_numerically_resolved'] for r in rows if 'homotopy' in r))
    matches=[];coverage=[]
    for case in CASES:
        sixrows=[r for r in outputs if r['case']==case and r['method']=='six_shared_Schur']
        materials=[(r['slice'],m) for r in sixrows for m in qualified(r)]
        distances=[dict(left_slice=a,right_slice=b,material_RMS=float(np.linalg.norm(m-n)/np.sqrt(len(m))))
                   for (a,m),(b,n) in itertools.combinations(materials,2) if a!=b]
        coverage.append(dict(case=case,nominal_slices=25,distinct_slices=len(sixrows),qualified_materials=len(materials),
            across_slice_distances=distances,max_across_slice_material_RMS=max([d['material_RMS'] for d in distances],default=None),
            connected_components='UNKNOWN'))
        for six in sixrows:
            baseline=next(r for r in outputs if r['case']==case and r['slice']==six['slice'] and r['method']=='FullJ8_anchor')
            left,right=qualified(six),qualified(baseline)
            matched=sum(any(np.linalg.norm(a-b)/np.sqrt(len(a))<=1e-7 for b in right) for a in left)
            matches.append(dict(case=case,slice=six['slice'],six_roots=len(left),FullJ8_roots=len(right),matched=matched,
                                same_root_set=len(left)==len(right)==matched))
    return dict(schema='Z-ALG-R2-ADAPTIVE-BRANCH-SUMMARY-v1',method_totals=totals,case_coverage=coverage,
        matches=matches,all_root_sets_match=all(r['same_root_set'] for r in matches),
        global_setup=global_setup,anchor_cost_reused_RHS=anchor_rhs,
        new_batch_actual_RHS=sum(r['accounting']['total_complex_rhs'] for r in outputs)+global_setup['total_complex_rhs'],
        scope=contract['scope'])

def plot(summary):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    methods=['six_shared_Schur','FullJ8_anchor'];labels=['Shared Schur + six','Anchor Full-J8']
    for ax,key,title in zip(axes,['qualified_materials','total_wall_seconds','total_RHS'],
            ['Qualified materials on common slices','Total wall including anchor + setup','Full RHS including anchor + setup']):
        vals=[summary['method_totals'][m][key] for m in methods]
        ax.bar(labels,vals,color=['#007f86','#6054a3'])
        for i,v in enumerate(vals):ax.text(i,v,str(round(v,3)),ha='center',va='bottom')
        ax.set_title(title);ax.tick_params(axis='x',rotation=12)
    fig.suptitle('Data-only anchored slice coverage; connected components remain unknown')
    fig.tight_layout();fig.savefig(P/'figures/adaptive_branch_comparison.png',dpi=170);plt.close(fig)

def main():
    path=P/'experiments/ADAPTIVE_BRANCH_CONTRACT.json'
    if '--freeze' in sys.argv:
        if path.exists():raise RuntimeError('Adaptive contract already frozen')
        write(path,setup_contract());print('ADAPTIVE_CONTRACT_FROZEN',path,flush=True);return
    if not path.exists():raise RuntimeError('Run --freeze and announce contract before execution')
    if list((P/'results').glob('adaptive_branch_*.json')):
        raise RuntimeError('Adaptive batch exists;preserve frozen initial run')
    contract=json.loads(path.read_text())
    g=load_case(CASES[0])[1];setupbook=CountingBook();model=Model(setupbook,**g)
    global_kernel=ReducedKernel(model,np.full(6,2.25),active=np.arange(6))
    global_setup=setupbook.summary(model.ns,model.nr)
    write(P/'results/adaptive_branch_global_setup.json',dict(accounting=global_setup,
          kernel={k:getattr(global_kernel,k) for k in ['m0','active','K','L','q','y0']},model_counts=dict(model.counts)))
    outputs=[]
    for case in CASES:
        _,g,y=load_case(case)
        for number,free in enumerate(contract['grids'][case]['unique_free_pairs']):
            fixed=contract['anchors'][case]['m'].copy();fixed[4:]=free;slice_id=f'slice_{number:02d}'
            for method in ['six_shared_Schur','FullJ8_anchor']:
                row=six_case(model,global_kernel,case,y,slice_id,fixed) if method=='six_shared_Schur' else baseline_case(case,g,y,slice_id,fixed,contract['baseline']['starts'][case])
                outputs.append(row);write(P/'results'/f'adaptive_branch_{case}_{slice_id}_{method}.json',row)
                print(case,slice_id,method,row['status'],'roots',len(qualified(row)),
                      'RHS',row['accounting']['total_complex_rhs'],flush=True)
    summary=summarize(outputs,contract,global_setup);write(P/'results/adaptive_branch_summary.json',summary)
    plot(summary);print('ADAPTIVE_COMPLETE',json.dumps(summary['method_totals']),flush=True)

if __name__=='__main__':main()
