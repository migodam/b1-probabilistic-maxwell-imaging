"""One frozen observation-only batch: six dual paths,24 paths,and8-start TRF."""
import os
for _variable in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_variable]='1'
import sys
from pathlib import Path
import json,time,itertools
from datetime import datetime,timezone
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import qmc
P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P/'src'))
from physics import Model,ReducedKernel,model_for
from algebraic_solver import fixed_active_polynomial,validate_candidate
from branch_enumerator import (enumerate_multiaffine,filter_real_box_roots,
    enumerate_bilinear_dual,filter_bilinear_dual_roots,classify_determinantal_endpoints)
from validate import CountingBook,load_case

def encode(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,complex):return {'real':value.real,'imag':value.imag}
    raise TypeError(type(value).__name__)

def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,default=encode,indent=2)+'\n');temp.replace(path)

def starts():
    points=qmc.Sobol(4,scramble=False).random_base2(4)
    chosen=[np.full(4,2.25)]
    for point in points:
        candidate=.5+3.5*point
        if all(np.linalg.norm(candidate-old)>1e-12 for old in chosen):chosen.append(candidate)
        if len(chosen)==8:break
    return np.asarray(chosen)

def contract():
    return dict(schema='Z-ALG-R2-BRANCH-v1',frozen_utc=datetime.now(timezone.utc).isoformat(),
        batch='ONE_INITIAL_BATCH',case_policy='E2 public reference;three existing blind E1 observations;never private/scores',
        E2=dict(stage='E2',public_reference_material=[1.,3.,2.,1.4],active=[0,1,2,3],slices=[[]]),
        E1=dict(cases=['E1_LOW6_00','E1_LOW6_01','E1_LOW6_02'],active=[0,1,2,3],free=[4,5],slice_values=[.75,2.25,3.75],
                choice='Fixed geometric cell order;no weak-coordinate selection or GT'),
        material_box=[.5,4.],primary='BILINEAR_DUAL_SIX_PATH_HOMOTOPY',standard='MULTIAFFINE_24_PATH_HOMOTOPY',
        homotopy=dict(six_seed=2026102904,standard_seed=2026102903,tolerance=2e-11,max_steps=4000,
                      initial_step=.01,min_step=1e-9,max_step=.04,max_newton=12,wall_seconds=60),
        baseline=dict(method='BOUNDED_FULLJ_TRF_8_STARTS',starts=starts(),max_nfev=120,
                      ftol=1e-12,xtol=1e-12,gtol=1e-12,loss='linear',x_scale=1.,
                      residual='all labeled observations divided by fixed target norm',
                      Jacobian='full material exact receiver transpose solves then slice columns',
                      start_policy='midpoint plus first7 distinct unscrambled Sobol4 points;midpoint duplicate skipped'),
        validation=dict(independent='fresh numpy full nonlinear solve at every unique legal candidate',
                        response_relative=1e-8,state_relative=1e-10,resolvent_sigma_min=1e-6,
                        merge_material_RMS=1e-7),
        denominator_audit=dict(relative_singular_gap_floor=1e-7,
                               conjugate='conjugate coefficients,keep same complex x'),
        scope=dict(root_coverage='unique independently qualified legal roots per fixed slice',
                   path_counts='numerical finite endpoint resolution only',
                   completeness='NOT_RIGOROUSLY_CERTIFIED',connected_components='UNKNOWN',
                   posterior='NOT_COMPUTED',global_fiber='NOT_COVERED_BY_NINE_SLICES'),
        accounting=dict(single_BLAS_thread=True,prior_unsaved_development_actual_complex_RHS=28,
                        development_detail='first24 run10;second24 collision inspection8;first6 run10;no private input',
                        full_RHS='source or receiver adjoint RHS',reduced='algebraic counts separately',
                        per_method_slice_wall_seconds=120,per_method_slice_RHS=15000))

def geometry(case):
    if case=='E2_PUBLIC':
        model=model_for(stage='E2')
        return dict(centers=model.centers,sources=model.sources,receivers=model.receivers)
    return load_case(case)[1]

def common_validate(model,y,rows):
    unique,output=[],[]
    for row in rows:
        if row['status']!='CANDIDATE':continue
        m=np.asarray(row['m'],dtype=float)
        if any(np.linalg.norm(m-old)/np.sqrt(len(m))<=1e-7 for old in unique):continue
        unique.append(m)
        validation=validate_candidate(model,m,y,response_rtol=1e-8)
        output.append(dict(m=m,validation=validation,source_path=row.get('path'),
                           status='QUALIFIED' if validation['qualified'] else 'REJECTED_INDEPENDENT_FORWARD'))
    return output

def run_algebraic(case,g,y,slice_id,fixed,method):
    book=CountingBook();model=Model(book,**g);begin=time.perf_counter();cpu=time.process_time()
    out=dict(case=case,slice=slice_id,method=method,fixed_material=fixed.copy(),active=[0,1,2,3])
    try:
        if method=='six':
            m0=fixed.copy();m0[:4]=.5
            kernel=ReducedKernel(model,m0,active=[0,1,2,3])
            result=enumerate_bilinear_dual(kernel,y,wall_seconds=60)
            filtered=filter_bilinear_dual_roots(result,validate=False)
            raw={k:v for k,v in result.items() if k!='system'}
        else:
            polynomial=fixed_active_polynomial(model,y,[0,1,2,3],fixed_material=fixed)
            result=enumerate_multiaffine(polynomial['coefficients'],model=model,wall_seconds=60)
            filtered=filter_real_box_roots(result,polynomial,validate=False)
            raw=result
            out['denominator_classification']=classify_determinantal_endpoints(result,polynomial)
            out['polynomial']={k:v for k,v in polynomial.items() if k!='kernel'}
        out.update(status='COMPLETE',homotopy=raw,physical_filter=filtered,
                   independently_validated=common_validate(model,y,filtered['roots']))
    except BaseException as exc:out.update(status='FAILED_PAID',error=repr(exc))
    out.update(method_wall_seconds=time.perf_counter()-begin,method_CPU_seconds=time.process_time()-cpu,
               accounting=book.summary(model.ns,model.nr),model_counts=dict(model.counts),model_timings=dict(model.timings))
    return out

def run_baseline(case,g,y,slice_id,fixed):
    book=CountingBook();model=Model(book,**g);begin=time.perf_counter();cpu=time.process_time()
    output=dict(case=case,slice=slice_id,method='FullJ8',fixed_material=fixed.copy(),active=[0,1,2,3],starts=[])
    target=model.pack(y);scale=max(np.linalg.norm(target),1e-300);candidates=[]
    try:
        for number,x0 in enumerate(starts()):
            cache={}
            def state(x):
                if 'x' not in cache or not np.array_equal(x,cache['x']):
                    m=fixed.copy();m[:4]=x
                    cache.update(x=x.copy(),s=model.forward(m))
                return cache['s']
            def fun(x):return (model.pack(state(x)['y'])-target)/scale
            def jac(x):return model.jacobian(state(x))[:,:4]/scale
            fit=least_squares(fun,x0,jac=jac,method='trf',bounds=(.5,4.),max_nfev=120,
                              ftol=1e-12,xtol=1e-12,gtol=1e-12,x_scale=1.,loss='linear')
            m=fixed.copy();m[:4]=fit.x
            residual=float(np.linalg.norm(fit.fun));is_root=residual<=1e-8
            output['starts'].append(dict(start=number,x0=x0,m=m,scipy_success=bool(fit.success),
                scipy_status=int(fit.status),message=fit.message,nfev=fit.nfev,njev=fit.njev,
                response_relative=residual,root_candidate=is_root))
            if is_root:candidates.append(dict(status='CANDIDATE',m=m,path=number))
        output.update(status='COMPLETE',independently_validated=common_validate(model,y,candidates))
    except BaseException as exc:output.update(status='FAILED_PAID',error=repr(exc))
    output.update(method_wall_seconds=time.perf_counter()-begin,method_CPU_seconds=time.process_time()-cpu,
                  accounting=book.summary(model.ns,model.nr),model_counts=dict(model.counts),model_timings=dict(model.timings))
    return output

def qualified(output):return [np.asarray(row['m']) for row in output.get('independently_validated',[]) if row['status']=='QUALIFIED']

def summarize(outputs,observation_cost):
    methods=('six','standard24','FullJ8');rows=[]
    for case,slice_id in dict.fromkeys((r['case'],r['slice']) for r in outputs):
        selected={r['method']:r for r in outputs if r['case']==case and r['slice']==slice_id}
        roots={method:qualified(selected[method]) for method in methods}
        comparisons={}
        for left,right in itertools.combinations(methods,2):
            dist=np.array([[np.linalg.norm(a-b)/np.sqrt(len(a)) for b in roots[right]] for a in roots[left]])
            comparisons[left+'_vs_'+right]=dict(left_roots=len(roots[left]),right_roots=len(roots[right]),
                 matches=sum(any(np.linalg.norm(a-b)/np.sqrt(len(a))<=1e-7 for b in roots[right]) for a in roots[left]),
                 distance_matrix=dist)
        rows.append(dict(case=case,slice=slice_id,qualified_counts={m:len(roots[m]) for m in methods},
                         matched=comparisons))
    by_method={}
    for method in methods:
        subset=[r for r in outputs if r['method']==method]
        by_method[method]=dict(slices=len(subset),complete=sum(r['status']=='COMPLETE' for r in subset),
            qualified_roots=sum(len(qualified(r)) for r in subset),
            covered_slices=sum(bool(qualified(r)) for r in subset),
            total_wall_seconds=sum(r['method_wall_seconds'] for r in subset),
            total_CPU_seconds=sum(r['method_CPU_seconds'] for r in subset),
            actual_complex_RHS=sum(r['accounting']['total_complex_rhs'] for r in subset),
            full_forward_calls=sum(r['accounting']['forward']['calls'] for r in subset),
            full_J_builds=sum(r['accounting']['adjoint_J']['calls'] for r in subset),
            algebraic_linear_solves=sum(r['model_counts'].get('algebraic_linear_solves',0) for r in subset),
            resolved_paths=sum(r.get('homotopy',{}).get('resolved_paths',0) for r in subset),
            unresolved_paths=sum(r.get('homotopy',{}).get('unresolved_paths',0) for r in subset))
    cross=[]
    for case in dict.fromkeys(r['case'] for r in outputs):
        material_rows=[(r['slice'],m) for r in outputs if r['case']==case and r['method']=='six' for m in qualified(r)]
        pairs=[dict(left_slice=s1,right_slice=s2,material_RMS=float(np.linalg.norm(m1-m2)/np.sqrt(len(m1))))
               for (s1,m1),(s2,m2) in itertools.combinations(material_rows,2) if s1!=s2]
        cross.append(dict(case=case,qualified_six_roots=len(material_rows),across_slice_distances=pairs,
                          connected_components='UNKNOWN;distance does not establish connectivity'))
    locus={key:sum(r.get('denominator_classification',{}).get('counts',{}).get(key,0) for r in outputs)
           for key in ('REGULAR_BOTH','ONE_BASE_LOCUS','BOTH_BASE_LOCUS')}
    return dict(schema='Z-ALG-R2-BRANCH-SUMMARY-v1',method_totals=by_method,slices=rows,
        across_slices=cross,standard24_denominator_counts=locus,
        development_prior_actual_complex_RHS=28,public_observation_generation=observation_cost,
        combined_actual_complex_RHS=28+observation_cost['total_complex_rhs']+sum(r['accounting']['total_complex_rhs'] for r in outputs),
        scope='Numerical fixed-slice root coverage;no connected-component/global-fiber/posterior/completeness claim')

def plot(summary):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    methods=['six','standard24','FullJ8'];labels=['Dual6 paths','Cleared24 paths','Full-J8 starts']
    totals=summary['method_totals'];fig,axes=plt.subplots(1,3,figsize=(13,4))
    for ax,key,title in zip(axes,['qualified_roots','total_wall_seconds','actual_complex_RHS'],
                            ['Unique qualified roots over28 slices','Total measured wall seconds','Full source/adjoint RHS']):
        vals=[totals[m][key] for m in methods];ax.bar(labels,vals,color=['#007f86','#d68d2b','#6054a3'])
        for i,val in enumerate(vals):ax.text(i,val,str(round(val,3)),ha='center',va='bottom')
        ax.set_title(title);ax.tick_params(axis='x',rotation=15)
    fig.suptitle('Fixed slices;blind E1 observations and public E2;finite numerical coverage')
    fig.tight_layout();(P/'figures').mkdir(exist_ok=True)
    fig.savefig(P/'figures/branch_comparison.png',dpi=170);plt.close(fig)

def main():
    contract_path=P/'experiments/BRANCH_CONTRACT.json'
    if contract_path.exists():raise RuntimeError('Frozen initial branch batch exists;no automatic rerun')
    if list((P/'results').glob('branch_*.json')):raise RuntimeError('Existing branch outputs;preserve paid run')
    write(contract_path,contract())
    book=CountingBook();reference=model_for(book,stage='E2')
    public_y=reference.forward([1.,3.,2.,1.4],independent=True)['y']
    observation_cost=book.summary(reference.ns,reference.nr)
    write(P/'results/branch_public_observation.json',dict(case='E2_PUBLIC',y=public_y,accounting=observation_cost))
    batches=[('E2_PUBLIC',geometry('E2_PUBLIC'),public_y,'all4',np.full(4,2.25))]
    for case in ['E1_LOW6_00','E1_LOW6_01','E1_LOW6_02']:
        _,g,y=load_case(case)
        for a,b in itertools.product([.75,2.25,3.75],repeat=2):
            fixed=np.full(6,2.25);fixed[4:]=[a,b]
            batches.append((case,g,y,f'free_{a:.2f}_{b:.2f}',fixed))
    outputs=[]
    for case,g,y,slice_id,fixed in batches:
        for method in ('six','standard24','FullJ8'):
            out=run_baseline(case,g,y,slice_id,fixed) if method=='FullJ8' else run_algebraic(case,g,y,slice_id,fixed,method)
            outputs.append(out);write(P/'results'/f'branch_{case}_{slice_id}_{method}.json',out)
            print(case,slice_id,method,out['status'],'roots',len(qualified(out)),
                  'rhs',out['accounting']['total_complex_rhs'],'wall',round(out['method_wall_seconds'],3),flush=True)
    summary=summarize(outputs,observation_cost);write(P/'results/branch_summary.json',summary)
    plot(summary)
    print('BRANCH_BATCH_COMPLETE',json.dumps(summary['method_totals']),flush=True)

if __name__=='__main__':main()
