"""Frozen six-material/three-source quotient versus64-start reduced Full-J."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import sys,json,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from scipy.linalg import lu_factor,lu_solve
from scipy.optimize import least_squares
P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P/'src'))
from physics import Model,ReducedKernel
from algebraic_solver import validate_candidate
from validate import CountingBook
from general_quotient_solver import enumerate_general_quotient,filter_general_roots

SOURCES=['src/general_quotient_solver.py','experiments/run_six_material_holdout.py',
         'src/physics.py','src/algebraic_solver.py','src/validate.py']

def encode(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,complex):return dict(real=value.real,imag=value.imag)
    raise TypeError(type(value).__name__)

def write(path,value):
    path=Path(path);temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,default=encode)+'\n');temporary.replace(path)

def snapshot(paths):
    return [dict(path=p,bytes=(P/p).stat().st_size,mtime_ns=(P/p).stat().st_mtime_ns) for p in paths]

def freeze():
    path=P/'results/general6_SOLVER_FREEZE_RECEIPT.json'
    if path.exists():raise RuntimeError('General6 solver already frozen')
    write(path,dict(status='CORE_SOLVER_AND_PARAMETERS_FROZEN_BEFORE_GENERATION',
        frozen_utc=datetime.now(timezone.utc).isoformat(),source_snapshot=snapshot(SOURCES),
        solver='n6,s3,one receiver;15quartic maximal minors;15x15leadingblock;20dimensional quotient',
        parameters=dict(reference_material=[.5]*6,separating=[1.,.371+.193j,-.227+.419j],
            residual_tolerance=2e-10,max_Newton_polish=5,leading_condition_limit=1e12,
            commutator_relative_limit=1e-7,eigenvector_condition_limit=1e12,separating_relative_eigengap_min=1e-10,
            backsub_rank=3,backsub_condition_limit=1e12,endpoint_J_condition_limit=1e12),
        baseline='Strong ReducedFullJ64 exact samekernel;cached reducedstate LU;3sourceforward RHS and1receiveradjoint RHS pernewstate/J',
        optimizer=dict(max_nfev=120,ftol=1e-12,xtol=1e-12,gtol=1e-12,x_scale=1.,loss='linear'),
        starts='exact64declared inGENERAL6_CONTRACT.json;noGTinitialization',
        validation=dict(relative_response=1e-8,state_relative=1e-10,sigma_min=1e-6,root_merge_RMS=1e-7),
        resource=dict(single_BLAS_thread=True,shared_full_setup_RHS_per_arm=9,targets=3,repeats=1),
        stop_policy='Stop onfailedleadingblock,noncommuting ornonseparating quotient,orunresolved/colliding endpoints;savepaidattempts;noGTtuning',
        certificate_flag='NOT_RIGOROUSLY_CERTIFIED',GT_access='NONE'))
    print('GENERAL6_CORE_FROZEN',path,flush=True)

def qualify(model,y,rows,box):
    output=[];seen=[]
    for row in rows:
        if row['status']!='CANDIDATE':continue
        m=np.asarray(row['m'])
        if any(np.linalg.norm(m-old)/np.sqrt(len(m))<=1e-7 for old in seen):continue
        seen.append(m);validation=validate_candidate(model,m,y,box=box,response_rtol=1e-8)
        output.append(dict(m=m,validation=validation,source_path=row.get('path'),
                           status='QUALIFIED' if validation['qualified'] else 'REJECTED_INDEPENDENT_FORWARD'))
    return output

def quotient(model,kernel,case,y,box):
    output=dict(case=case,method='general_quotient20',status='RUNNING')
    try:
        result=enumerate_general_quotient(kernel,y)
        filtered=filter_general_roots(result,box=box)
        system=result['system']
        output.update(result={k:v for k,v in result.items() if k!='system'},physical_filter=filtered,
            dual_system={key:getattr(system,key) for key in ['w0','N','W','V','coefficients','scale']},
            independently_validated=qualify(model,y,filtered['roots'],box),
            status='COMPLETE' if result['all_paths_numerically_resolved'] else 'NUMERICAL_GATE_FAILURE')
    except BaseException as exc:
        output.update(status='MATH_GATE_FAILURE_PAID',error=repr(exc),diagnostic=getattr(exc,'diagnostic',{}))
    return output

def reduced64(model,kernel,case,y,box,initials):
    target=model.pack(y);scale=max(np.linalg.norm(target),1e-300)
    output=dict(case=case,method='ReducedFullJ64',status='RUNNING',starts=[]);candidates=[]
    try:
        for number,x0 in enumerate(initials):
            cache={}
            def state(x):
                if 'x' not in cache or not np.array_equal(cache['x'],x):
                    delta=x-kernel.m0
                    A,b=kernel.coefficients(delta)
                    with model.operation('GENERAL6_REDUCED_FORWARD',reduced_LU=1,reduced_complex_RHS=model.ns):
                        factor=lu_factor(A,check_finite=False);w=lu_solve(factor,b,check_finite=False)
                        current=dict(y=kernel.y0+kernel.L@w,w=w,factor=factor)
                    cache.clear();cache.update(x=x.copy(),state=current)
                return cache['state']
            def fun(x):return (model.pack(state(x)['y'])-target)/scale
            def jac(x):
                current=state(x)
                if 'J' not in current:
                    with model.operation('GENERAL6_REDUCED_J',reduced_complex_RHS=model.nr,reduced_adjoint_RHS=model.nr,reduced_J_builds=1):
                        transfer=lu_solve(current['factor'],kernel.L.T,trans=1,check_finite=False).T
                        field=kernel.q+kernel.K@current['w']
                        J=np.einsum('ri,is->sri',transfer,field).reshape(model.ns*model.nr,kernel.r)
                        current['J']=np.stack([J.real,J.imag],axis=1).reshape(2*model.ns*model.nr,kernel.r)
                return current['J']/scale
            fit=least_squares(fun,x0,jac=jac,bounds=tuple(box),method='trf',max_nfev=120,
                              ftol=1e-12,xtol=1e-12,gtol=1e-12,x_scale=1.,loss='linear')
            residual=float(np.linalg.norm(fit.fun))
            output['starts'].append(dict(start=number,x0=x0,m=fit.x,response_relative=residual,
                scipy_success=bool(fit.success),scipy_status=int(fit.status),message=fit.message,nfev=fit.nfev,njev=fit.njev))
            if residual<=1e-8:candidates.append(dict(status='CANDIDATE',m=fit.x,path=number))
        output.update(status='COMPLETE',independently_validated=qualify(model,y,candidates,box))
    except BaseException as exc:output.update(status='FAILED_PAID',error=repr(exc))
    return output

def roots(output):return [np.asarray(r['m']) for r in output.get('independently_validated',[]) if r['status']=='QUALIFIED']

def main():
    if '--freeze' in sys.argv:freeze();return
    freeze_path=P/'results/general6_SOLVER_FREEZE_RECEIPT.json';frozen=json.loads(freeze_path.read_text())
    if snapshot(SOURCES)!=frozen['source_snapshot']:raise RuntimeError('Frozen source changed;stop')
    if list((P/'results').glob('general6_GENERAL6_*.json')):raise RuntimeError('General6 predictions alreadyexist;no repeats')
    contract=json.loads((P/'experiments/GENERAL6_CONTRACT.json').read_text())
    manifest=json.loads((P/'experiments/general6_public/manifest.json').read_text())
    outputs=[];receipts=[];stopped=False;begin=time.perf_counter();prediction_files=[]
    for method in ['general_quotient20','ReducedFullJ64']:
        if stopped:break
        tick=time.perf_counter();cpu=time.process_time();book=CountingBook();model=Model(book,**contract['geometry'])
        kernel=ReducedKernel(model,np.full(6,.5),active=np.arange(6))
        for row in manifest['cases']:
            case=row['case'];y=np.load(P/'experiments/general6_public'/row['observation_file'])['y']
            before=model.counts.copy();rhs=book.rhs;ev=len(book.events);task=time.perf_counter()
            output=quotient(model,kernel,case,y,contract['material_box']) if method=='general_quotient20' else reduced64(model,kernel,case,y,contract['material_box'],contract['initialization']['ReducedFullJ64_starts'])
            output.update(box=contract['material_box'],incremental_RHS=book.rhs-rhs,events=book.events[ev:],
                model_counts=dict(model.counts-before),task_wall_before_logging=time.perf_counter()-task,
                exposure='FRESH_GENERAL6_TRUTH_WITHHELD_FROM_SOLVER',certificate_flag='NOT_RIGOROUSLY_CERTIFIED')
            filename=f'results/general6_{case}_{method}.json';write(P/filename,output)
            outputs.append(output);prediction_files.append(filename)
            print('GENERAL6',case,method,output['status'],'qualified',len(roots(output)),flush=True)
            if output['status']!='COMPLETE':stopped=True;break
        receipt=dict(method=method,accounting=book.summary(model.ns,model.nr),model_counts=dict(model.counts),
            total_wall_including_logging=time.perf_counter()-tick,total_CPU_seconds=time.process_time()-cpu,
            setup_RHS=9,qualified_roots=sum(len(roots(o)) for o in outputs if o['method']==method))
        receipts.append(receipt);write(P/'results'/f'general6_{method}_accounting.json',receipt)
    comparison=[]
    for row in manifest['cases']:
        selected={o['method']:o for o in outputs if o['case']==row['case']}
        if len(selected)==2:
            left,right=roots(selected['general_quotient20']),roots(selected['ReducedFullJ64'])
            matches=sum(any(np.linalg.norm(a-b)/np.sqrt(6)<=1e-7 for b in right) for a in left)
            comparison.append(dict(case=row['case'],quotient_roots=len(left),ReducedFullJ64_roots=len(right),
                                   matches=matches,same_root_set=len(left)==len(right)==matches))
    summary=dict(status='STOPPED_NO_GO_DIAGNOSTIC' if stopped else 'ALL_PREDICTIONS_COMPLETE_UNSCORED',
        methods=receipts,cases=comparison,new_solver_actual_RHS=sum(r['accounting']['total_complex_rhs'] for r in receipts),
        total_solver_wall_including_logging=time.perf_counter()-begin,generation_RHS_declared=9,
        source_snapshot=snapshot(SOURCES),truth_access='NONE',certificate_flag='NOT_RIGOROUSLY_CERTIFIED')
    write(P/'results/general6_summary.json',summary)
    if snapshot(SOURCES)!=frozen['source_snapshot']:raise RuntimeError('Frozen sources changed duringrun')
    write(P/'results/general6_PREDICTION_FREEZE_RECEIPT.json',dict(status='PARTIAL_PREDICTIONS_FROZEN_NO_GO' if stopped else 'ALL_PREDICTIONS_FROZEN_BEFORE_SCORING',
        frozen_utc=datetime.now(timezone.utc).isoformat(),predictions=snapshot(prediction_files),source_snapshot=snapshot(SOURCES),
        truth_access='NONE',summary='results/general6_summary.json'))
    print('GENERAL6_PREDICTIONS_FROZEN',json.dumps({r['method']:dict(roots=r['qualified_roots'],RHS=r['accounting']['total_complex_rhs'],wall=r['total_wall_including_logging']) for r in receipts}),flush=True)

if __name__=='__main__':main()
