"""Exposed-data quotient coverage and a frozen six-task five-repeat timing."""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import sys,json,time,itertools
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from scipy.optimize import least_squares
P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P/'src'));sys.path.insert(0,str(P/'experiments'))
from physics import Model,ReducedKernel,model_for
from quotient_solver import enumerate_quotient
from branch_enumerator import enumerate_bilinear_dual,filter_bilinear_dual_roots
from algebraic_solver import validate_candidate
from validate import CountingBook,load_case
from run_branch_study import write,starts

def decode(value):
    if isinstance(value,dict):
        if set(value)=={'real','imag'}:return complex(value['real'],value['imag'])
        return {key:decode(v) for key,v in value.items()}
    if isinstance(value,list):return [decode(v) for v in value]
    return value

def qualified(row):return [np.asarray(r['m']) for r in row.get('independently_validated',[]) if r['status']=='QUALIFIED']

def data_rows():
    rows=[];adaptive=json.loads((P/'experiments/ADAPTIVE_BRANCH_CONTRACT.json').read_text())
    high=json.loads((P/'experiments/HIGH4_CONTRACT.json').read_text())
    public_y=np.asarray(decode(json.loads((P/'results/branch_public_observation.json').read_text()))['y'])
    low4=model_for(stage='E2');g4={k:getattr(low4,k) for k in ['centers','sources','receivers']}
    for path in sorted((P/'results').glob('branch_*_six.json')):
        saved=json.loads(path.read_text());case=saved['case'];slice_id=saved['slice']
        if case=='E2_PUBLIC':g,y,pool=g4,public_y,'LOW4'
        else:_,g,y=load_case(case);pool='LOW6'
        rows.append(dict(id=f'fixed_{case}_{slice_id}',study='fixed',case=case,slice=slice_id,pool=pool,
            geometry=g,y=y,box=[.5,4.],fixed=np.asarray(saved['fixed_material']),starts=starts(),
            saved_six=str(path.relative_to(P)),saved_FullJ8=str(path.relative_to(P)).replace('_six.json','_FullJ8.json')))
    for path in sorted((P/'results').glob('adaptive_branch_*_six_shared_Schur.json')):
        saved=json.loads(path.read_text());case=saved['case'];slice_id=saved['slice'];_,g,y=load_case(case)
        rows.append(dict(id=f'adaptive_{case}_{slice_id}',study='adaptive',case=case,slice=slice_id,pool='LOW6',
            geometry=g,y=y,box=[.5,4.],fixed=np.asarray(saved['fixed_material']),starts=np.asarray(adaptive['baseline']['starts'][case]),
            saved_six=str(path.relative_to(P)),saved_FullJ8=str(path.relative_to(P)).replace('_six_shared_Schur.json','_FullJ8_anchor.json')))
    for row in json.loads((P/'experiments/high4_public/manifest.json').read_text())['cases']:
        case=row['case'];y=np.load(P/'experiments/high4_public'/row['observation_file'])['y']
        rows.append(dict(id=f'high4_{case}',study='high4',case=case,slice='all4',pool='HIGH4',geometry=high['geometry'],
            y=y,box=high['material_box'],fixed=np.full(4,15.25),starts=np.asarray(high['initialization']['FullJ8_starts']),
            saved_six=f'results/high4_{case}_six.json',saved_FullJ8=f'results/high4_{case}_FullJ8.json'))
    if len(rows)!=101:raise RuntimeError(f'Expected101 frozen tasks,found{len(rows)}')
    return rows,adaptive

def freeze():
    path=P/'experiments/QUOTIENT_BENCHMARK_CONTRACT.json'
    if path.exists():raise RuntimeError('Quotient benchmark already frozen')
    rows,adaptive=data_rows()
    selected=['fixed_E2_PUBLIC_all4','adaptive_E1_LOW6_02_slice_12','adaptive_E1_LOW6_01_slice_00',
              'high4_E2_HIGH4_00','high4_E2_HIGH4_01','high4_E2_HIGH4_02']
    contract=dict(schema='Z-ALG-R2-QUOTIENT-BENCH-v1',frozen_utc=datetime.now(timezone.utc).isoformat(),
        exposure='EXPOSED_DATA_DIAGNOSTIC_NOT_FRESH_BLIND',coverage_tasks=101,
        coverage='One quotient pass;compare qualified legal material sets against saved six andFullJ8 receipts',
        timing=dict(task_ids=selected,repeats=5,methods=['quotient','six','FullJ8','ReducedFullJ8'],
                    order='cyclic method rotation byrepeat',aggregation='all5 samples,median;no favorable repeat selection',
                    total='whole sharedpool includes geometry,reference setup,restriction,allroots,independent validation,JSONlogging'),
        quotient=dict(separating=.371+.193j,tolerance=2e-10,max_polish=5,leading_condition_limit=1e12,
                      commutator_relative_limit=1e-7,eigenvector_condition_limit=1e12,
                      separating_relative_eigengap_min=1e-10,backsub_rank=2,backsub_condition_limit=1e12,
                      simple_endpoint_J_condition_limit=1e12,certificate_flag='NOT_RIGOROUSLY_CERTIFIED',
                      fallback='Anyfailure,gates,unresolved orcolliding endpoints invokes counted sixpath;retain failedattempt'),
        shared_setup=dict(pools=['LOW6','LOW4','HIGH4'],reference='box midpoint eachmaterial',
                          full_setup_RHS={'LOW6':8,'LOW4':6,'HIGH4':6},
                          all_slices='F11 Schur restriction;one geometry/kernel setup perpool perrun'),
        baseline=dict(starts='Exact original8 starts for eachfrozen task;no new support orGT use',
                      max_nfev=120,ftol=1e-12,xtol=1e-12,gtol=1e-12,loss='linear',x_scale=1.),
        reduced_baseline='Same8starts andTRF with exact samecached Schur kernel andanalytic reducedJacobian;full validation eachunique candidate',
        validation=dict(response_relative=1e-8,state_relative=1e-10,sigma_min=1e-6,root_merge_RMS=1e-7,
                        fresh_full_solve='everyunique legal candidate on everybenchmark run'),
        anchor_costs=adaptive['anchors'],
        scope='Finite roots on declared slices;no rigorous completeness,component,posterior,or freshblindclaim',
        task_inputs=[{k:v for k,v in row.items() if k not in ('y','geometry')} for row in rows])
    write(path,contract);print('QUOTIENT_CONTRACT_FROZEN',flush=True)

def qualify(model,y,rows,box):
    out=[];seen=[]
    for row in rows:
        if row['status']!='CANDIDATE':continue
        m=np.asarray(row['m'])
        if any(np.linalg.norm(m-old)/np.sqrt(len(m))<=1e-7 for old in seen):continue
        seen.append(m);validation=validate_candidate(model,m,y,box=box,response_rtol=1e-8)
        out.append(dict(m=m,validation=validation,source_path=row.get('path'),
                        status='QUALIFIED' if validation['qualified'] else 'REJECTED_INDEPENDENT_FORWARD'))
    return out

def algebraic(model,kernel,row,method):
    output=dict(id=row['id'],method=method,status='RUNNING',fallback_used=False)
    try:
        if method=='quotient':
            try:
                result=enumerate_quotient(kernel,row['y'])
                if not result['all_paths_numerically_resolved']:
                    output['failed_quotient_attempt']={k:v for k,v in result.items() if k!='system'}
                    raise ValueError('UNRESOLVED_OR_COLLIDING_QUOTIENT_ENDPOINTS')
            except (ValueError,np.linalg.LinAlgError,FloatingPointError) as exc:
                output.update(fallback_used=True,fallback_reason=repr(exc),
                              failed_quotient_diagnostic=getattr(exc,'diagnostic',{}))
                result=enumerate_bilinear_dual(kernel,row['y'],wall_seconds=60)
        else:result=enumerate_bilinear_dual(kernel,row['y'],wall_seconds=60)
        filtered=filter_bilinear_dual_roots(result,box=row['box'],validate=False)
        output.update(status='COMPLETE',result={k:v for k,v in result.items() if k!='system'},
                      physical_filter=filtered,independently_validated=qualify(model,row['y'],filtered['roots'],row['box']))
    except BaseException as exc:output.update(status='FAILED_PAID',error=repr(exc))
    return output

def fullj(model,row,kernel=None):
    target=model.pack(row['y']);scale=max(np.linalg.norm(target),1e-300);initials=row['starts'];candidates=[]
    output=dict(id=row['id'],method='ReducedFullJ8' if kernel is not None else 'FullJ8',starts=[],status='RUNNING')
    try:
        for index,x0 in enumerate(initials):
            cache={}
            def state(x):
                if 'x' not in cache or not np.array_equal(cache['x'],x):
                    if kernel is None:
                        m=row['fixed'].copy();m[:4]=x;cache.clear();cache.update(x=x.copy(),s=model.forward(m))
                    else:
                        cache.clear();cache.update(x=x.copy(),s=dict(y=kernel.evaluate(x-kernel.m0[:4])))
                return cache['s']
            def fun(x):return (model.pack(state(x)['y'])-target)/scale
            def jac(x):
                current=state(x)
                if kernel is None:return model.jacobian(current)[:,:4]/scale
                if 'J' not in current:current['J']=kernel.jacobian(x-kernel.m0[:4])
                return current['J']/scale
            fit=least_squares(fun,x0,jac=jac,method='trf',bounds=tuple(row['box']),max_nfev=120,
                              ftol=1e-12,xtol=1e-12,gtol=1e-12,x_scale=1.,loss='linear')
            m=row['fixed'].copy();m[:4]=fit.x;relative=float(np.linalg.norm(fit.fun))
            output['starts'].append(dict(start=index,x0=x0,m=m,nfev=fit.nfev,njev=fit.njev,
                scipy_status=int(fit.status),scipy_success=bool(fit.success),message=fit.message,response_relative=relative))
            if relative<=1e-8:candidates.append(dict(status='CANDIDATE',m=m,path=index))
        output.update(status='COMPLETE',independently_validated=qualify(model,row['y'],candidates,row['box']))
    except BaseException as exc:output.update(status='FAILED_PAID',error=repr(exc))
    return output

def compare(left,right):
    a,b=qualified(left),qualified(right)
    matches=sum(any(np.linalg.norm(x-y)/np.sqrt(len(x))<=1e-7 for y in b) for x in a)
    return dict(left_roots=len(a),right_roots=len(b),matches=matches,same_root_set=len(a)==len(b)==matches)

def run_pool(rows,method,prefix,anchors):
    begin=time.perf_counter();cpu=time.process_time();outputs=[];pool_receipts=[]
    for pool in ['LOW6','LOW4','HIGH4']:
        selected=[row for row in rows if row['pool']==pool]
        if not selected:continue
        book=CountingBook();model=Model(book,**selected[0]['geometry'])
        global_kernel=None
        if method!='FullJ8':global_kernel=ReducedKernel(model,np.full(model.n,np.mean(selected[0]['box'])),active=np.arange(model.n))
        for row in selected:
            counts=model.counts.copy();event_start=len(book.events);rhs=book.rhs;task_begin=time.perf_counter()
            if method=='FullJ8':out=fullj(model,row)
            else:
                kernel=global_kernel.restrict([0,1,2,3],row['fixed']) if model.n==6 else global_kernel
                out=fullj(model,row,kernel) if method=='ReducedFullJ8' else algebraic(model,kernel,row,method)
            out.update(study=row['study'],case=row['case'],slice=row['slice'],box=row['box'],
                       counts=dict(model.counts-counts),incremental_RHS=book.rhs-rhs,
                       events=book.events[event_start:],task_wall_before_logging=time.perf_counter()-task_begin)
            outputs.append(out);write(P/'results'/f'{prefix}_{row["id"]}_{method}.json',out)
        pool_receipts.append(dict(pool=pool,accounting=book.summary(model.ns,model.nr),model_counts=dict(model.counts)))
    actual_elapsed=time.perf_counter()-begin
    anchor_cases=sorted({r['case'] for r in rows if r['study']=='adaptive'})
    anchor_rhs=sum(anchors[case]['original_cost']['total_complex_rhs'] for case in anchor_cases)
    anchor_wall=sum(anchors[case]['original_cost']['wall_seconds'] for case in anchor_cases)
    anchor_cpu=sum(anchors[case]['original_cost']['CPU_seconds'] for case in anchor_cases)
    receipt=dict(method=method,tasks=len(rows),new_actual_RHS=sum(r['accounting']['total_complex_rhs'] for r in pool_receipts),
        anchor_reused_RHS=anchor_rhs,total_RHS=sum(r['accounting']['total_complex_rhs'] for r in pool_receipts)+anchor_rhs,
        incremental_wall_including_logging=actual_elapsed,total_wall_including_logging_and_anchor=actual_elapsed+anchor_wall,
        total_CPU_seconds=time.process_time()-cpu+anchor_cpu,pools=pool_receipts,
        qualified_roots=sum(len(qualified(r)) for r in outputs),fallbacks=sum(r.get('fallback_used',False) for r in outputs),
        failures=sum(r['status']!='COMPLETE' for r in outputs))
    write(P/'results'/f'{prefix}_{method}_pool_receipt.json',receipt)
    return outputs,receipt

def plot(summary):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    methods=['quotient','six','FullJ8','ReducedFullJ8'];labels=['Quotient','Six-path','Full-J8','Reduced Full-J8']
    fig,axes=plt.subplots(1,3,figsize=(13,4.2))
    for ax,key,title in zip(axes,['median_wall_seconds','qualified_roots_per_repeat','total_RHS_per_repeat'],
                            ['Median total wall over five repeats','Qualified roots in six tasks','Full RHS including anchor + setup']):
        values=[summary['timing'][m][key] for m in methods];ax.bar(labels,values,color=['#007f86','#d68d2b','#6054a3','#4774a5'])
        for i,v in enumerate(values):ax.text(i,v,str(round(v,4)),ha='center',va='bottom')
        ax.set_title(title);ax.tick_params(axis='x',rotation=18)
    fig.suptitle('Previously exposed observations; fixed six-task comparison; no completeness certificate')
    fig.tight_layout();fig.savefig(P/'figures/quotient_benchmark_comparison.png',dpi=170);plt.close(fig)

def main():
    if '--freeze' in sys.argv:freeze();return
    contract=json.loads((P/'experiments/QUOTIENT_BENCHMARK_CONTRACT.json').read_text())
    if list((P/'results').glob('quotient_bench_*.json')):raise RuntimeError('Preserve first quotient benchmark;no automatic rerun')
    rows,adaptive=data_rows();anchors=adaptive['anchors']
    coverage,coverage_receipt=run_pool(rows,'quotient','quotient_bench_coverage',anchors)
    comparison=[]
    for row,out in zip([r for pool in ['LOW6','LOW4','HIGH4'] for r in rows if r['pool']==pool],coverage):
        oldsix=json.loads((P/row['saved_six']).read_text());oldfull=json.loads((P/row['saved_FullJ8']).read_text())
        comparison.append(dict(id=row['id'],vs_saved_six=compare(out,oldsix),vs_saved_FullJ8=compare(out,oldfull)))
    selected=[r for r in rows if r['id'] in contract['timing']['task_ids']]
    repeats=[]
    for repeat in range(5):
        methods=['quotient','six','FullJ8','ReducedFullJ8'];order=methods[repeat%4:]+methods[:repeat%4]
        for method in order:
            output,receipt=run_pool(selected,method,f'quotient_bench_timing_r{repeat}',anchors)
            receipt['repeat']=repeat;repeats.append(receipt)
            print('TIMING',repeat,method,'wall',round(receipt['total_wall_including_logging_and_anchor'],5),
                  'roots',receipt['qualified_roots'],'RHS',receipt['total_RHS'],flush=True)
    timing={}
    for method in ['quotient','six','FullJ8','ReducedFullJ8']:
        samples=[r for r in repeats if r['method']==method];walls=[r['total_wall_including_logging_and_anchor'] for r in samples]
        timing[method]=dict(wall_samples=walls,median_wall_seconds=float(np.median(walls)),
            min_wall_seconds=min(walls),max_wall_seconds=max(walls),
            qualified_roots_per_repeat=samples[0]['qualified_roots'],root_counts=[r['qualified_roots'] for r in samples],
            total_RHS_per_repeat=samples[0]['total_RHS'],RHS_samples=[r['total_RHS'] for r in samples],
            fallbacks=sum(r['fallbacks'] for r in samples),failures=sum(r['failures'] for r in samples))
    summary=dict(schema='Z-ALG-R2-QUOTIENT-BENCH-SUMMARY-v1',exposure=contract['exposure'],coverage=coverage_receipt,
        coverage_comparisons=comparison,all101_match_saved_six=all(r['vs_saved_six']['same_root_set'] for r in comparison),
        timing=timing,repeats=repeats,new_total_actual_RHS=coverage_receipt['new_actual_RHS']+sum(r['new_actual_RHS'] for r in repeats),
        certificate_flag='NOT_RIGOROUSLY_CERTIFIED',scope=contract['scope'])
    write(P/'results/quotient_bench_summary.json',summary);plot(summary)
    print('QUOTIENT_BENCH_COMPLETE',json.dumps(timing),flush=True)

if __name__=='__main__':main()
