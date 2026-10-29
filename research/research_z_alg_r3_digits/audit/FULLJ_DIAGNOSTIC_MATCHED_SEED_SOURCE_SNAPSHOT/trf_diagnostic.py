#!/usr/bin/env python3
"""POST-MAIN-GT development diagnostic: mature bounded SciPy TRF, no truth input."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse,json,time,shutil,traceback
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from threadpoolctl import threadpool_limits
from src.solver import Run,load_public,json_safe
ROOT=Path(__file__).resolve().parent
CONFIG=ROOT/'audit/FULLJ_DIAGNOSTIC_CONFIG.json'
FREEZE=ROOT/'audit/FULLJ_DIAGNOSTIC_SOURCE_SNAPSHOT'


class TRFRun(Run):
    def __init__(self,public,config):super().__init__(public,'fullj_native',config)
    def gn(self,start,B,seconds_end,reciprocal=False,origin='TRF',max_iterations=None):
        if reciprocal:raise ValueError('Diagnostic is native bounded TRF')
        self.phase_deadline=seconds_end
        identity=B.shape[1]==self.model.n
        DB=self.reg.D.toarray() if identity else self.reg.D@B
        cache={};calls=jac_calls=0;begun=self.elapsed();status='INTERRUPTED_BY_STAGE_OR_WALL_OR_ERROR'
        cs=np.sqrt(2*self.reg.smooth/self.reg.edges)
        ct=np.sqrt(2*self.reg.tv/self.reg.edges)
        def state(x):
            if 'x' not in cache or not np.array_equal(cache['x'],x) or cache['state'] is not self.model.state:
                current=self.model.forward(x if identity else B@x,label='mature_scipy_trf')
                self.score(current,'mature_scipy_trf',x)
                d=self.reg.D@current['m'];q=np.sqrt(d*d+self.reg.epsilon**2)
                cache.update(x=x.copy(),state=current,d=d,q=q)
            return cache['state']
        def fun(x):
            nonlocal calls
            self.guard();calls+=1;s=state(x);d,q=cache['d'],cache['q']
            # .5||rTV||² = lambdaTV mean(sqrt(d²+eps²)-eps), using signed stable residual.
            return np.r_[self.model.pack(s['y']-self.target)/np.sqrt(self.scale2),cs*d,ct*d/np.sqrt(q+self.reg.epsilon)]
        def jac(x):
            nonlocal jac_calls
            self.guard();jac_calls+=1;s=state(x);q=cache['q']
            J=self.model.jacobian(s)
            self.model.counts['TRF_augmented_J_builds']+=1
            J=J if identity else J@B
            # d[d/sqrt(q+eps)]/dd = sqrt(q+eps)/(2q), exact also at d=0.
            tv_derivative=ct*np.sqrt(q+self.reg.epsilon)/(2*q)
            return np.vstack((J/np.sqrt(self.scale2),cs*DB,tv_derivative[:,None]*DB))
        try:
            fit=least_squares(fun,np.clip(start,self.cfg['lower'],self.cfg['upper']),jac=jac,
                bounds=(self.cfg['lower'],self.cfg['upper']),method='trf',tr_solver='exact',x_scale='jac',
                loss='linear',max_nfev=1000000000,ftol=1e-10,xtol=1e-10,gtol=1e-10)
            status='SCIPY_TERMINATED'
            self.events.append(dict(kind='SCIPY_TRF_TERMINATION',stage=self.active_stage,status=int(fit.status),
                 message=fit.message,nfev=int(fit.nfev),njev=int(fit.njev),cost=float(fit.cost),
                 optimality=float(fit.optimality),active_bound_count=int(np.count_nonzero(fit.active_mask))))
            return fit.x,dict(iterations=None,accepted=None,optimizer='SciPy least_squares TRF',nfev=int(fit.nfev),njev=int(fit.njev))
        finally:
            self.events.append(dict(kind='SCIPY_TRF_PHASE',stage=self.active_stage,phase_status=status,
                residual_calls=calls,jacobian_calls=jac_calls,wall_seconds=self.elapsed()-begun,
                internal_dense_svd='SciPy exact trust-region subproblem; optimizer internal factorizations not separately instrumented',
                all_source_LU_and_FullJ_RHS='counted by shared physical ledger'))
    def finish(self):
        result=super().finish();result['method']='scipy_trf_diagnostic'
        result['development_status']='POST_MAIN_GT_DEVELOPMENT; cases selected by already-visible public high residual; not fresh blind evidence'
        result['truth_access']='NONE_IN_THIS_SCRIPT; public operators and anonymous Y only'
        result['optimizer']='scipy.optimize.least_squares(method=trf,tr_solver=exact,x_scale=jac)'
        result['TV_residual']='sqrt(2lambdaTV/nedge)*d/sqrt(sqrt(d²+epsilon²)+epsilon); exact same objective'
        return result


def write(path,data):
    path=Path(path)
    if path.exists():raise FileExistsError('Preserving '+str(path))
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(json_safe(data),indent=2,allow_nan=False)+'\n')


def freeze():
    if FREEZE.exists():raise FileExistsError('Diagnostic source snapshot exists')
    FREEZE.mkdir(parents=True)
    names=['trf_diagnostic.py','audit/FULLJ_DIAGNOSTIC_CONFIG.json','src/solver.py','src/physics.py','src/multiscale.py','src/__init__.py']
    rows=[]
    for name in names:
        out=FREEZE/name;out.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,out)
        rows.append(dict(path=name,bytes=out.stat().st_size))
    write(FREEZE/'FREEZE.json',dict(frozen_utc=datetime.now(timezone.utc).isoformat(),files=rows,
         status='POST_MAIN_GT_DEVELOPMENT_NOT_FRESH_BLIND',truth_input_access='NONE',
         source_freeze_before_diagnostic_predictions=True,formula_experiments='NONE_ADDED',hashes='NOT_COMPUTED'))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--freeze',action='store_true');p.add_argument('--public');p.add_argument('--operators');p.add_argument('--output')
    a=p.parse_args()
    if a.freeze:freeze();return
    if not all((a.public,a.operators,a.output)):p.error('Supply public/operators/output')
    manifest=json.loads((FREEZE/'FREEZE.json').read_text())
    for row in manifest['files']:
        if (ROOT/row['path']).read_bytes()!=(FREEZE/row['path']).read_bytes():raise RuntimeError('Diagnostic source changed afterfreeze')
    config=json.loads(CONFIG.read_text())['config'];run=None;begun,cpu=time.perf_counter(),time.process_time()
    try:
        with threadpool_limits(limits=1):
            run=TRFRun(load_public(a.public,a.operators),config);run.started,run.cpu_started=begun,cpu
            try:result=run.run()
            except BaseException as exc:
                run.status='FAILED_PAID_UNEXPECTED';run.stop_reason=type(exc).__name__;run.failures[type(exc).__name__]+=1
                run.events.append(dict(kind='UNEXPECTED_FAILURE_PAID',error=repr(exc),seconds=run.elapsed()));result=run.finish()
    except BaseException as exc:
        result=dict(method='scipy_trf_diagnostic',status='FAILED_INITIALIZATION',error=repr(exc),traceback=traceback.format_exc(),
                    wall_seconds=time.perf_counter()-begun,cpu_seconds=time.process_time()-cpu,accounting={'counts':None},top1=None)
    result['diagnostic_source_freeze_utc']=manifest['frozen_utc'];result['public_input_path']=str(Path(a.public).resolve())
    write(a.output,result)
    print(json.dumps(dict(output=a.output,status=result['status'],wall_seconds=result['wall_seconds'],
        residual=(result.get('top1') or {}).get('relative_residual'),accounting=result.get('accounting'))),flush=True)


if __name__=='__main__':main()
