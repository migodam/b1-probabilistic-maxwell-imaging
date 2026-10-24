"""Create the inert campaign registry and audit its gate without any fit."""
from pathlib import Path
import copy,json,sys,time
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'tools'))
from p3r4.core import Book,OUT,write
import queue_r4_physics_transport as controller


def fixture(composite=94.,phi=79.,nll=109.,selected=300):
    history=[]
    for update in range(0,3001,300):
        row=dict(update=update,validation_composite=110.,validation_Phi=110.,validation_full_NLL=110.)
        if update==0:row.update(validation_composite=100.,validation_Phi=100.,validation_full_NLL=100.)
        if update==selected:row.update(validation_composite=composite,validation_Phi=phi,validation_full_NLL=nll)
        history.append(row)
    return dict(status='COMPLETE',completed_updates=3000,mode='physics',selected_update=selected,
                selected_validation_metric=composite,history=history)


def main():
    job='physics_transport_campaign_readiness_'+str(time.time_ns())
    book=Book('I',job,analysis_only=True,max_wall=120)
    out=OUT/'AUDIT'/job;out.mkdir(parents=True)
    status,error='COMPLETE',None
    try:
        registry=controller.build_registry()
        if controller.REGISTRY.exists():
            if controller.load(controller.REGISTRY)!=registry:raise FileExistsError('existing registry differs; never overwrite a frozen campaign')
        else:write(controller.REGISTRY,registry)
        plans=controller.validate_registry(registry,execution=False)
        assert registry['released'] is False and len(plans)==18
        for seed in controller.SEEDS:
            assert sum(p['seed']==seed and p['mode']=='physics' for p in plans)==3
            assert sum(p['seed']==seed and p['mode']=='supervised' for p in plans)==3
        supervised=dict(status='COMPLETE',completed_updates=3000,mode='supervised')
        cases=[]
        for name,result,expected in [
            ('all_pass',fixture(),True),
            ('composite_fails',fixture(composite=96.),False),
            ('Phi_fails',fixture(phi=81.),False),
            ('NLL_fails',fixture(nll=111.),False),
            ('exact_boundaries_pass',fixture(composite=95.,phi=80.,nll=110.),True)]:
            gate=controller.replication_gate(result,supervised)
            assert gate['eligible']==expected
            cases.append(dict(name=name,expected=expected,actual=gate['eligible'],conditions=gate['conditions']))
        no_reselection=fixture(composite=96.)
        no_reselection['history'][2].update(validation_composite=50.,validation_Phi=50.,validation_full_NLL=50.)
        assert not controller.replication_gate(no_reselection,supervised)['eligible']
        # Neither execution nor subprocess creation is allowed by an inert registry.
        with patch.object(controller.subprocess,'run',side_effect=AssertionError('unexpected process launch')) as launch:
            try:controller.run_campaign(registry,controller.REGISTRY)
            except ValueError as exc:assert 'inert' in str(exc)
            else:raise AssertionError('inert registry executed')
            assert launch.call_count==0
        fail=fixture();fail['status']='FAILED'
        try:controller.replication_gate(fail,supervised)
        except ValueError:pass
        else:raise AssertionError('failed fit admitted to gate')
        write(out/'READINESS.json',dict(status='PASS_READY_INERT',registry=controller.identity(controller.REGISTRY),
            released=False,conditional_plans=plans,mock_gate_cases=cases,alternate_checkpoint_reselection_rejected=True,
            failed_fit_gate_rejected=True,unreleased_execution_rejected=True,subprocess_launches=0,
            actual_Target_instances=0,actual_physics_RHS=0,optimizer_updates=0,trainer_math_changes=False,
            max_physics_fits=9,max_supervised_fits=9,max_wall_per_fit=2400,
            maximum_planned_RHS=1930752,maximum_hard_RHS=2484000,
            preflight_scope='saved actual64RHS preflight FAIL retained; corrected SAME_INPUT audit PASS;0newphysics',
            posterior_qualification='NOT_ESTABLISHED',receipt=str(book.jobdir/'CLOSED.json')))
        print(out)
    except BaseException as exc:status,error='FAILED',repr(exc);raise
    finally:book.close(status,error)


if __name__=='__main__':main()
