"""Full FWI command line; observations contain only the public complex y."""
from pathlib import Path
import argparse
import json
import os
import sys
import time
sys.dont_write_bytecode=True
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS']:
    os.environ[key]='1'
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from gsda.fem_model import FEMModel
from gsda.fwi import invert,FWIOptions,regularization,dump,complete_outputs
from gsda.guard import OnlineGuard
from gsda.accounting import BudgetBook


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--config',type=Path,required=True)
    ap.add_argument('--measurement','--observation',dest='measurement',type=Path)
    ap.add_argument('--output','--out',dest='output',type=Path,required=True)
    ap.add_argument('--options',type=Path)
    ap.add_argument('--job',default='fwi')
    ap.add_argument('--audit-gradient',action='store_true')
    ap.add_argument('--audit-assembly',action='store_true')
    ap.add_argument('--assembly-pair',nargs=2,default=['legacy_two_pass','cached_combined'])
    ap.add_argument('--heldout-endpoint',type=Path,help='predict declared heldout frequencies from an already sealed endpoint')
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    if sum([a.audit_gradient,a.audit_assembly,a.heldout_endpoint is not None])>1:raise ValueError('choose one audit/prediction action')
    guard=None if a.audit_gradient else OnlineGuard(ROOT).install()
    config=json.loads(a.config.read_text())
    freeze=json.loads(a.options.read_text()) if a.options else {}
    options=freeze.get('fwi_options',freeze)
    book=BudgetBook(root=ROOT,job=a.job);status='FAILED';result=None;model=None
    try:
        book.check(new_job=True)
        if a.heldout_endpoint is None:
            if a.measurement is None:raise ValueError('--measurement is required for inversion/audit')
            with np.load(a.measurement,allow_pickle=False) as obs:
                if set(obs.files)!={'y'}:raise ValueError('online NPZ must contain only y')
                y=obs['y'].copy()
        else:
            seal=a.heldout_endpoint.parent/'ONLINE_SEALED.json'
            if not seal.exists():raise RuntimeError('heldout prediction requires sealed endpoint')
            if a.heldout_endpoint.stat().st_mtime>seal.stat().st_mtime:raise RuntimeError('endpoint modified after seal')
            config['frequencies']=config['heldout_frequencies']
        model=FEMModel(config,book,phase='R4_FULL_FWI')
        if a.heldout_endpoint is not None:
            with np.load(a.heldout_endpoint,allow_pickle=False) as endpoint:z=endpoint['z'].copy()
            pred=model.predict(z)
            np.savez_compressed(a.output/'heldout_prediction.npz',data=pred.data,z=z,frequencies=config['frequencies'])
            result=dict(status='HELDOUT_PREDICTION_SAVED',sealed_endpoint=str(a.heldout_endpoint),
                truth_access=False,frequencies=config['frequencies'],diagnostics=pred.diagnostics)
            dump(a.output/'HELDOUT_PREDICTION.json',result)
        elif a.audit_assembly:
            z=model.initial(0);o=FWIOptions(**options);prior,pg,_=regularization(z,model.nwall,model.target_shape,model.wall_shape,o)
            comparisons={}
            for strategy in a.assembly_pair:
                model.config['assembly_strategy']=strategy;cpu=time.process_time();wall=time.perf_counter()
                value,g,pred=model.gradient(z,y)
                comparisons[strategy]=dict(value=value+prior,g=g+pg,data=pred.data,
                    cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
            old,new=[comparisons[strategy] for strategy in a.assembly_pair]
            data_absolute=float(np.max(np.abs(old['data']-new['data'])))
            gradient_relative=float(np.linalg.norm(old['g']-new['g'])/max(np.linalg.norm(old['g']),1e-15))
            objective_relative=abs(old['value']-new['value'])/max(abs(old['value']),1.)
            passed=data_absolute<1e-10 and gradient_relative<1e-9 and objective_relative<1e-10
            np.savez_compressed(a.output/'ASSEMBLY_EQUIVALENCE.npz',z=z,legacy_data=old['data'],cached_data=new['data'],legacy_gradient=old['g'],cached_gradient=new['g'])
            result=dict(status='PASS' if passed else 'FAIL',assembly_pair=a.assembly_pair,data_max_absolute_error=data_absolute,
                gradient_relative_error=gradient_relative,objective_relative_error=objective_relative,
                legacy_objective=old['value'],cached_objective=new['value'],
                legacy_cpu_seconds=old['cpu_seconds'],cached_cpu_seconds=new['cpu_seconds'],
                legacy_wall_seconds=old['wall_seconds'],cached_wall_seconds=new['wall_seconds'],
                scope='one same public uniform material state; identical physical quadrature, operator and objective; not reconstruction evidence')
            dump(a.output/'ASSEMBLY_AUDIT.json',result)
            if not passed:raise RuntimeError('cached assembly equivalence failed')
        elif a.audit_gradient:
            z=model.initial(0);_,_,scale=model.bounds_and_scale()
            rng=np.random.default_rng(4103001);v=rng.normal(size=model.nparam)*scale
            v/=np.linalg.norm(v/scale);step=1e-5
            o=FWIOptions(**options)
            value,g,pred=model.gradient(z,y)
            prior,pg,_=regularization(z,model.nwall,model.target_shape,model.wall_shape,o)
            def score(zz):
                dv,_,_=model._run(zz,y,False)
                return dv+regularization(zz,model.nwall,model.target_shape,model.wall_shape,o)[0]
            plus=score(z+step*v);minus=score(z-step*v)
            fd=(plus-minus)/(2*step);analytic=float((g+pg)@v)
            error=abs(fd-analytic)/max(abs(fd),abs(analytic),1e-12)
            result=dict(status='PASS' if error<1e-3 else 'FAIL',relative_error=error,
                finite_difference=fd,weighted_adjoint=analytic,objective=value+prior,
                scope='one necessary mixed wall/material directional derivative',
                cost=book.snapshot(pid=os.getpid()),diagnostics=pred.diagnostics)
            dump(a.output/'GRADIENT_AUDIT.json',result)
            if error>=1e-3:raise RuntimeError('weighted adjoint audit failed')
        else:
            result=invert(model,y,a.output,options)
        status=result['status']
    except BaseException as error:
        book.failure('R4_FULL_FWI','run_fwi',error)
        dump(a.output/'FAILURE.json',dict(type=type(error).__name__,message=str(error)))
        raise
    finally:
        artifact_error=None
        if model is not None and a.heldout_endpoint is None and not a.audit_gradient and (a.output/'endpoint.npz').exists() and (a.output/'ONLINE_SEALED.json').exists():
            try:complete_outputs(model,a.output)
            except Exception as error:
                artifact_error=str(error);book.failure('R4_FULL_FWI','output aliases',error)
        receipt=book.receipt(status=status,output=str(a.output),options=str(a.options),artifact_error=artifact_error)
        dump(a.output/('HELDOUT_RECEIPT.json' if a.heldout_endpoint is not None else 'RECEIPT.json'),receipt)
        if a.heldout_endpoint is None:dump(a.output/'runtime_cost.json',receipt)
        if guard is not None:guard.save(a.output/('HELDOUT_FILE_ACCESS_SUMMARY.json' if a.heldout_endpoint is not None else 'FILE_ACCESS_SUMMARY.json'))
    print(json.dumps(dict(status=status,output=str(a.output),cost=book.snapshot(pid=os.getpid()))),flush=True)

if __name__=='__main__':main()
