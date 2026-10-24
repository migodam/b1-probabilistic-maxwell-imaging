"""CLOSED-only snapshot + independent algebra checks, zero physical actions."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import csv,json,time,zipfile,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import torch
from p3r4.core import Book,OUT,pack
from p3r4.frozen_geometry_pcn import FrozenMatrixPCN
from p3r4.probability_physics import stable_importance

P=OUT/'AUDIT/independent_probability_review'
MODULES={'p3r4.train_density','p3r4.train_anchor_density','p3r4.train_risk','p3r4.train_risk_bin','p3r4.train_risk_anchorbin','p3r4.train_field_moments','p3r4.train_risk_fusion'}

def load(p):return json.loads(p.read_text())
def opt(argv,key,default=None):
    try:return argv[argv.index(key)+1]
    except (ValueError,IndexError):return default

def input_count(path):
    if path is None:return None
    path=Path(path);path=path if path.is_absolute() else ROOT/path
    try:
        with zipfile.ZipFile(path) as archive:
            for name in ('parent.npy','parents.npy','y.npy','observation.npy'):
                if name in archive.namelist():
                    with archive.open(name) as f:
                        version=np.lib.format.read_magic(f)
                        shape,_,_=np.lib.format.read_array_header_1_0(f) if version==(1,0) else np.lib.format.read_array_header_2_0(f)
                        return int(shape[0])
    except (OSError,ValueError):return None
    return None


def main():
    job='independent_probability_review_data_'+str(time.time_ns())
    b=Book('I',job,analysis_only=True,max_wall=600)
    status='COMPLETE';error=None
    try:
        torch.set_num_threads(1);started=time.time_ns()
        closed=[];open_jobs=[]
        for command_path in sorted(OUT.glob('*/jobs/*/COMMAND.json')):
            command=load(command_path);argv=command.get('argv',[])
            module=opt(argv,'-m','')
            if module not in MODULES:continue
            folder=command_path.parent;receipt=folder/'CLOSED.json'
            if not receipt.exists():
                open_jobs.append(dict(job=folder.name,module=module,status='OPEN_EXCLUDED_NO_METRICS_READ',output=opt(argv,'--out')));continue
            receipt_before=receipt.stat().st_mtime_ns;closure=load(receipt)
            out=Path(opt(argv,'--out',opt(argv,'--output','')))
            out=out if out.is_absolute() else ROOT/out
            result_path=out/'RESULT.json';contract_path=out/'RUN_CONTRACT.json'
            # Never read live/partial results. Failed CLOSED attempts are recorded
            # with no selected metric, even if residual artifacts are present.
            valid=closure['status']=='COMPLETE' and result_path.exists()
            result=load(result_path) if valid else {}
            contract=load(contract_path) if contract_path.exists() else {}
            if receipt.stat().st_mtime_ns!=receipt_before:raise RuntimeError('CLOSED receipt changed during snapshot')
            mc=result.get('model_config',contract.get('model_config',{}))
            ac=result.get('anchor_contract',contract.get('anchor_contract',{}))
            kind=result.get('kind',contract.get('kind',opt(argv,'--kind','')))
            metric_name=result.get('criterion',contract.get('criterion'))
            value=result.get('selected_validation_metric')
            if value is None:
                if 'selected_validation_full_complex_mse' in result:metric_name='validation_full576_complex_MSE';value=result['selected_validation_full_complex_mse']
                elif 'validation_full_complex_mse' in result:metric_name='validation_full576_complex_MSE';value=result['validation_full_complex_mse']
                elif 'validation_nll' in result:metric_name='validation_nll_raw_x66';value=result['validation_nll']
            train=result.get('train_size',contract.get('train_size'))
            val=result.get('validation_size',contract.get('validation_size'))
            train_source='RESULT/RUN_CONTRACT' if train is not None else 'requested/header only; failed attempt not a completed fit'
            requested_n=opt(argv,'--train-size')
            if train is None:train=int(requested_n) if requested_n else input_count(opt(argv,'--train'))
            if val is None:val=input_count(opt(argv,'--validation'))
            if module=='p3r4.train_density':scope='raw_y_only_joint'
            elif module=='p3r4.train_anchor_density':scope='explicit_anchor_conditional_joint'
            elif module=='p3r4.train_field_moments':scope='raw_y_only_mean_or_pointwise_moments'
            else:scope='anchor_risk_context' if opt(argv,'--train-context',opt(argv,'--anchor-train')) else result.get('information_set','risk/control; see command')
            if train==32768 and val==2048:panel='joint32768/2048'
            elif train==1024 and val==256:panel='strong_anchor1024/256'
            elif train==512 and val==128:panel='historical_anchor512/128'
            else:panel=f'{train}/{val}'
            history=result.get('history',[]);initial=history[0] if history else {}
            initial_metric=initial.get('validation_metric',initial.get('validation_full_complex_mse',initial.get('validation_component_mse')))
            initial_name='validation_metric (same criterion)' if 'validation_metric' in initial else 'validation_full_complex_mse' if 'validation_full_complex_mse' in initial else 'validation_component_mse (ROI/component units)' if 'validation_component_mse' in initial else ''
            row=dict(job=folder.name,stage=folder.parent.parent.name,module=module,CLOSED_status=closure['status'],final_result_eligible=valid,
                output=str(out),result_path=str(result_path) if result_path.exists() else '',receipt=str(receipt),closure_error=closure.get('error'),
                train_count=train,validation_count=val,count_source=train_source,data_panel=panel,input_scope=scope,
                kind=kind,seed=result.get('seed',contract.get('seed',opt(argv,'--seed'))),parameter_count=result.get('parameter_count',contract.get('parameter_count')),
                selected_update=result.get('selected_update',result.get('best_update')),
                selected_initialization=result.get('selected_initialization'),
                updates_requested=result.get('updates_requested',result.get('updates',opt(argv,'--updates'))),
                completed_updates=result.get('completed_updates'),stop_reason=result.get('stop_reason'),
                selected_metric_name=metric_name,selected_metric_value=value,initial_metric_name=initial_name,initial_metric_value=initial_metric,
                condition_dim=mc.get('context_dim'),CFM_velocity_hidden=mc.get('velocity_hidden'),encoder_width=mc.get('width',result.get('encoder_width',opt(argv,'--width'))),
                coordinate_preconditioner=ac.get('coordinate_preconditioner','none' if module=='p3r4.train_anchor_density' else ''),
                fixed_frame=result.get('frame',contract.get('frame',opt(argv,'--frame'))),
                triangular_first_block_dim=mc.get('weak_dim') if kind=='triangular' else None,
                CFM_metric=ac.get('metric') if kind=='cfm' else None,learning_rate=result.get('learning_rate',opt(argv,'--lr',opt(argv,'--learning-rate'))),
                controlled_wall_seconds=closure.get('controlled_wall_seconds'),CPU_seconds=closure.get('profile',{}).get('CPU_seconds'),
                RSS_peak_bytes=closure.get('profile',{}).get('RSS_peak_bytes'),
                interpretation='frozen historical reference, no new optimization' if kind=='r2_cnn' else 'classical label-fit/control; zero neural parameter count does not mean no fitted coefficients' if kind in ('affine','sharpen','spatial') else '')
            closed.append(row)
        fields=list(closed[0]) if closed else []
        with (P/'CLOSED_FITS.csv').open('w',newline='') as file:
            writer=csv.DictWriter(file,fieldnames=fields);writer.writeheader();writer.writerows(closed)
        (P/'CLOSED_FITS.json').write_text(json.dumps(closed,indent=2,allow_nan=False)+'\n')
        (P/'OPEN_JOBS_EXCLUDED.json').write_text(json.dumps(open_jobs,indent=2)+'\n')
        # No Target is constructed: check frozen arithmetic/kernels only.
        rng=np.random.default_rng(105);noise=rng.normal(size=(32,8))+1j*rng.normal(size=(32,8));sigma=.00119211798514519
        phi=float(np.vdot(noise,noise).real/sigma**2);white=pack(noise)*np.sqrt(2)/sigma
        whiten_relative=abs(.5*float(white@white)-phi)/phi
        q,_=np.linalg.qr(rng.normal(size=(66,66)));Q=(q*np.linspace(.001,.049,66))@q.T;Q=(Q+Q.T)/2
        kernel=FrozenMatrixPCN(Q);a=rng.normal(size=66);z=kernel.propose(a,rng.normal(size=66))
        cancellation=-.5*(z@z-a@a)+kernel.log_proposal(a,z)-kernel.log_proposal(z,a)
        points=rng.normal(size=(9,66));logq=-.5*(points*points).sum(1)-33*np.log(2*np.pi);phis=np.full(9,7.)
        imp=stable_importance(phis,points,logq)
        assert whiten_relative<1e-14 and abs(cancellation)<1e-10
        assert np.max(abs(imp['normalized_weights']-1/9))<1e-12
        assert abs(imp['log_mean_unnormalized_weight']-(-7+33*np.log(2*np.pi)))<1e-12
        (P/'INDEPENDENT_ALGEBRA.json').write_text(json.dumps(dict(status='PASS',whiten_Phi_relative_error=whiten_relative,
            pCN_prior_cancellation_error=float(abs(cancellation)),flat_likelihood_IS_ESS=imp['ESS'],
            log_mean_unnormalized_weight=imp['log_mean_unnormalized_weight'],
            normalized_prior_log_evidence_for_fixture=-7.,omitted_Gaussian_normalizing_constant=33*np.log(2*np.pi),
            interpretation='ISnormalized weights correct; logmean unnormalized is not normalized Bayesian evidence without Gaussian constant',
            physics_calls=0,optimizer_updates=0),indent=2)+'\n')
        files=['src/p3r4/core.py','src/p3r4/density.py','src/p3r4/anchor_density.py','src/p3r4/conditional_coordinates.py','src/p3r4/probability_physics.py','src/p3r4/frozen_geometry_pcn.py','src/p3r4/probability_eval.py','src/p3r4/field_moments.py','tools/run_r4_probability_audit.py','tools/r4_geometry_probability.py']
        snapshots=[]
        for relative in files:
            path=ROOT/relative;stat=path.stat();snapshots.append(dict(path=str(path),size=stat.st_size,mtime_ns=stat.st_mtime_ns))
        summary=dict(snapshot_time_ns=started,closed_fit_attempts=len(closed),eligible_completed_fits=sum(r['final_result_eligible'] for r in closed),
            closed_failed_attempts=sum(r['CLOSED_status']!='COMPLETE' for r in closed),open_fits_excluded=len(open_jobs),
            complete_data_panels={panel:sum(r['final_result_eligible'] and r['data_panel']==panel for r in closed) for panel in sorted(set(r['data_panel'] for r in closed))},
            source_snapshot=snapshots,physics_calls=0,optimizer_updates=0,Book_receipt=str(b.jobdir/'CLOSED.json'),
            scaffold_receipt=load(P/'BOOK_SCOPE.json').get('receipt'))
        (P/'SNAPSHOT_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
        print(json.dumps({k:v for k,v in summary.items() if k!='source_snapshot'},indent=2))
    except BaseException as exc:status='FAILED';error=repr(exc);raise
    finally:b.close(status,error)

if __name__=='__main__':main()
