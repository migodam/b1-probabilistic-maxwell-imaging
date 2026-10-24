"""Create inert post-fit protocol and verify guards/ratio algebra; no inference."""
from pathlib import Path
import json,sys,time
from unittest.mock import patch
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'tools'))
from p3r4.core import Book,OUT,write
from p3r4.corrected_transport import log_acceptance,component_for,pcn_proposal,multichain_diagnostics
import r4_postfit_dev as evaluator


def main():
    job='postfit_dev_readiness_'+str(time.time_ns())
    book=Book('I',job,analysis_only=True,max_wall=120)
    folder=OUT/'AUDIT'/job;folder.mkdir(parents=True)
    status,error='COMPLETE',None
    try:
        p=evaluator.build_protocol()
        if evaluator.PROTOCOL.exists():
            if evaluator.load(evaluator.PROTOCOL)!=p:raise FileExistsError('existing postfit protocol differs; no overwrite')
        else:write(evaluator.PROTOCOL,p)
        evaluator.validate_protocol(p)
        identities=[]
        original_identity=evaluator.identity
        def tracked(path):
            identities.append(str(path));return original_identity(path)
        with patch.object(evaluator,'identity',side_effect=tracked):
            closed,pending=evaluator.closed_methods(p,require_all=False)
        pending_ids={row['id'] for row in pending}
        assert not any(m['checkpoint'] in identities for m in p['methods'] if m['id'] in pending_ids)
        inert=folder/'UNRELEASED.json';write(inert,dict(released=False))
        with patch.object(evaluator.subprocess,'run',side_effect=AssertionError('no child launch allowed')) as launch:
            for phase in ('generator','paid_probability','corrected_MH'):
                try:evaluator.phase_release(inert,p,phase)
                except ValueError:pass
                else:raise AssertionError('inert phase accepted')
            assert launch.call_count==0
        try:evaluator.validate_binding({},closed)
        except ValueError:pass
        else:raise AssertionError('missing validation-only preinference MH binding accepted')
        rng=np.random.default_rng(171);x=rng.normal(size=66);z=rng.normal(size=66)
        lp=lambda a:float(-.5*(a@a)-33*np.log(2*np.pi))
        prior_cancel=abs(log_acceptance('q',x,7.,z,7.,lp(x),lp(z)))
        assert prior_cancel<1e-12
        assert log_acceptance('pcn',x,7.,z,9.)==-2.
        assert component_for('mixture_half',.499)=='q' and component_for('mixture_half',.5)=='pcn'
        beta=.002
        e=rng.normal(size=66);proposal=pcn_proposal(x,e)
        ratio_forward=-.5*np.sum((proposal-np.sqrt(1-beta*beta)*x)**2)/beta**2
        ratio_backward=-.5*np.sum((x-np.sqrt(1-beta*beta)*proposal)**2)/beta**2
        pcn_cancel=abs(-.5*(proposal@proposal-x@x)+ratio_backward-ratio_forward)
        assert pcn_cancel<1e-10
        diagnostic=multichain_diagnostics(np.ones((4,128,2)))
        assert diagnostic['constant_or_stuck_chain_function_count']==8
        assert all(np.isnan(v) for row in diagnostic['ESS_per_chain'] for v in row)
        assert p['paid_probability']['planned_RHS']==9*(128+8*63)*8==45504
        assert p['mh']['planned_RHS']==4*3*4*129*8==49536
        write(folder/'FROZEN_INERT_PROTOCOL.json',p)
        write(folder/'READINESS.json',dict(status='PASS_INERT_PENDING_CLOSED_FITS',protocol=original_identity(evaluator.PROTOCOL),
            closed_methods=[m['id'] for m in closed],pending=pending,no_live_checkpoint_identity_reads=True,
            no_unreleased_stage_execution=True,no_missing_MH_binding=True,
            independence_prior_cancellation_error=prior_cancel,pcn_prior_cancellation_error=pcn_cancel,
            mixture='constant half selector; each component individually Metropolized',
            stuck_chain_not_assigned_ESS_N=True,paid_probability_RHS=45504,MH_RHS=49536,
            actual_Target_instances=0,physical_RHS=0,optimizer_updates=0,inference_samples=0,GT_reads=0,child_launches=0,
            no_frozen_fit_source_modified=True,receipt=str(book.jobdir/'CLOSED.json')))
        print(folder)
    except BaseException as exc:status,error='FAILED',repr(exc);raise
    finally:book.close(status,error)


if __name__=='__main__':main()
