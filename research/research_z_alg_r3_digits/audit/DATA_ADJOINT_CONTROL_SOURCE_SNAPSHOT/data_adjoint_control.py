#!/usr/bin/env python3
"""Post-generation, pre-GT ordinary data-adjoint L-BFGS + Full-J attribution control.

Only the virtual proposal objective differs from GlobalDual. Original solver and
budget-fill sources are immutable. No truth, class, support or template is read.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import argparse,json,shutil,time,traceback
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.linalg import lu_solve
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from src.solver import Run,load_public,json_safe
ROOT=Path(__file__).resolve().parent
CONTRACT=ROOT/'configs/DATA_ADJOINT_CONTROL.json'
SNAPSHOT=ROOT/'audit/DATA_ADJOINT_CONTROL_SOURCE_SNAPSHOT'


def utc():return datetime.now(timezone.utc).isoformat()


def write_json(path,data):
    path=Path(path)
    if path.exists():raise FileExistsError('Preserving existing output '+str(path))
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(json_safe(data),indent=2,allow_nan=False)+'\n')


class DataAdjointRun(Run):
    def __init__(self,public,config):
        # Keep the inherited GlobalDual stage dispatcher and replace only its virtual dual() hook.
        super().__init__(public,'globaldual',config)

    def dual(self,start,B,seconds_end):
        self.phase_deadline=seconds_end
        identity_basis=B.shape[1]==self.model.n
        accepted_iterations=0
        def accepted_callback(c):
            nonlocal accepted_iterations
            accepted_iterations+=1
            self.model.counts['data_adjoint_optimizer_accepted_iterations']+=1
            self.events.append(dict(kind='DATA_ADJOINT_ACCEPTED_ITERATION',stage=self.active_stage,
                                    accepted_iterations=accepted_iterations,seconds=self.elapsed(),additional_RHS=0))
        def objective(c):
            self.guard()
            state=self.model.forward(c if identity_basis else B@c,label='ordinary_data_adjoint_lbfgs')
            row=self.score(state,'ordinary_data_adjoint_lbfgs',c)
            residual=state['y']-self.target
            # f_data=.5||R||²/||Y||², so its real material gradient has no leading2.
            # dF=S M^-1 diag(dm)Q; M^H Lambda=S^H R gives g_i=Re(sum conj(Lambda_is)Q_is)/||Y||².
            with self.model.operation('ORDINARY_DATA_ADJOINT',rhs=self.model.ns,
                                      adjoint_RHS=self.model.ns,data_adjoint_gradients=1):
                adj=lu_solve(self.model.factor,self.model.SH@residual,trans=2,check_finite=False)
                g=np.real(np.sum(adj.conj()*state['E'],axis=1))/self.scale2
                value=.5*float(np.vdot(residual,residual).real)/self.scale2
            reg,greg,_=self.reg.evaluate(state['m'])
            gradient=g+greg
            if not np.isfinite(value+reg) or not np.all(np.isfinite(gradient)):
                raise FloatingPointError('Nonfinite ordinary-data adjoint objective/gradient')
            self.events.append(dict(kind='DATA_ADJOINT_EVALUATION',stage=self.active_stage,seconds=self.elapsed(),
                                    proposal_objective=value+reg,ranking_score=row['objective'],Gram_calls=0))
            return value+reg,(gradient if identity_basis else B.T@gradient)
        fit=minimize(objective,np.asarray(start),jac=True,method='L-BFGS-B',callback=accepted_callback,
                     bounds=[(self.cfg['lower'],self.cfg['upper'])]*B.shape[1],
                     options=dict(maxiter=int(self.cfg['dual_maxiter']),maxfun=int(self.cfg['dual_maxfun']),
                                  ftol=float(self.cfg['dual_ftol']),gtol=float(self.cfg['dual_gtol']),maxls=25))
        self.events.append(dict(kind='DATA_ADJOINT_OPTIMIZER',stage=self.active_stage,status=int(fit.status),
                                message=str(fit.message),nfev=int(fit.nfev),iterations=int(fit.nit),
                                proposal_objective=float(fit.fun)))
        return fit.x

    def finish(self):
        result=super().finish()
        result['method']='data_adjoint'
        result['supplemental_control']='ORDINARY_DATA_ADJOINT_LBFGS_PLUS_FULLJ'
        result['source_freeze_after_public_Y_generation']=True
        result['source_freeze_fully_pre_generation']=False
        result['truth_access']='NONE_PUBLIC_INPUT_ONLY_BEFORE_MAIN_GT_SCORING'
        result['proposal_objective']='same normalized original-data objective and generic material regularizer as final online ranking'
        result['gradient_formula']='M^H Lambda=S^H(F-Y); g_i=Re(sum_s conj(Lambda_is)Q_is)/||Y||² + reg_gradient'
        result['matching']='same common coarse320, .45 proposal fraction, LBFGS options, NativeFullJ tail,90s; GN cap1e9'
        result['inherited_origin_label_note']='dual_fullj_polish is the inherited generic virtual-hook tail name; this control performs zero Gram objectives'
        for group in ('stages','events'):
            for row in result[group]:
                if row.get('kind')=='MULTISCALE_STAGE' and row.get('stage')==24:
                    row['optimizer']='data_adjoint_lbfgs_plus_native_fullj'
        return result


def freeze():
    if SNAPSHOT.exists():raise FileExistsError('Data-adjoint snapshot already exists')
    SNAPSHOT.mkdir(parents=True)
    names=['data_adjoint_control.py','configs/DATA_ADJOINT_CONTROL.json','src/solver.py','src/physics.py','src/multiscale.py','src/__init__.py']
    files=[]
    for name in names:
        target=SNAPSHOT/name;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,target);files.append(dict(path=name,bytes=target.stat().st_size))
    write_json(SNAPSHOT/'FREEZE.json',dict(frozen_utc=utc(),files=files,
        scope='Independent ordinary-data optimizer attribution control before any main GT scoring',
        after_public_Y_generation=True,fully_pre_generation_source_freeze=False,
        original_primary_and_budget_fill_sources_unchanged=True,truth_access='NONE',hashes='NOT_COMPUTED',
        validation='Derived standard exact adjoint analytically; syntax/static checks only, no additional formula experiments'))


def check_freeze():
    manifest=json.loads((SNAPSHOT/'FREEZE.json').read_text())
    for row in manifest['files']:
        if (ROOT/row['path']).read_bytes()!=(SNAPSHOT/row['path']).read_bytes():
            raise RuntimeError('Source changed after data-adjoint control freeze: '+row['path'])
    return manifest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--freeze',action='store_true');p.add_argument('--public');p.add_argument('--operators');p.add_argument('--output')
    a=p.parse_args()
    if a.freeze:freeze();return
    if not all((a.public,a.operators,a.output)):p.error('Supply --public --operators --output')
    manifest=check_freeze();contract=json.loads(CONTRACT.read_text())
    live=None;begun,cpu=time.perf_counter(),time.process_time()
    try:
        with threadpool_limits(limits=1):
            live=DataAdjointRun(load_public(a.public,a.operators),contract['config'])
            #90s includes public loading and construction as well as all physical/numerical work.
            live.started,live.cpu_started=begun,cpu
            try:result=live.run()
            except BaseException as exc:
                live.status='FAILED_PAID_UNEXPECTED';live.stop_reason=type(exc).__name__
                live.failures[type(exc).__name__]+=1
                live.events.append(dict(kind='UNEXPECTED_FAILURE_PAID',error=repr(exc),seconds=live.elapsed()))
                result=live.finish()
    except BaseException as exc:
        result=dict(method='data_adjoint',status='FAILED_BEFORE_RUN_INITIALIZATION',error=repr(exc),
                    traceback=traceback.format_exc(),wall_seconds=time.perf_counter()-begun,cpu_seconds=time.process_time()-cpu,
                    top1=None,m=None,candidates=[],archive=[],all_unknowns_released=False,
                    accounting=dict(counts=None,reason='Run initialization failed; counters unavailable'),truth_access='NONE')
    result['control_contract']=contract;result['control_source_freeze_utc']=manifest['frozen_utc']
    result['public_input_path']=str(Path(a.public).resolve())
    write_json(a.output,result)
    print(json.dumps(dict(output=a.output,method='data_adjoint',status=result['status'],wall_seconds=result['wall_seconds'],
                          cpu_seconds=result['cpu_seconds'],relative_residual=(result.get('top1') or {}).get('relative_residual'),
                          all_unknowns_released=result.get('all_unknowns_released'),accounting=result.get('accounting'))))


if __name__=='__main__':main()
