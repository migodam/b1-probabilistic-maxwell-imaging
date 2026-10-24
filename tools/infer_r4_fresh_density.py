"""Authorized observation-only fresh inference, after physical/context closure.

Never invokes score/truth/Full/training. Waits60s between marker checks; each
cohort runs once under an analysis-only StageIV Book with2GiB RSS cap.
"""
from pathlib import Path
import json,os,sys,subprocess,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,write
from p3r4.train_physics_transport import identity,verify_identity

GLOBAL=OUT/'configs/FROZEN_DENSITY_BOTH_COHORTS.json'
PPC=OUT/'configs/INERT_FRESH_DENSITY_FULL_PPC_REQUEST.json'
START=OUT/'results/FRESH_DENSITY_INFERENCE_STARTED.json'
COMPLETE=OUT/'results/FRESH_DENSITY_ONLINE_BOTH_COMPLETE.json'
PHYSICAL=OUT/'results/FRESH_PHYSICS_AND_CONTEXTS_COMPLETE.json'


def load(path):return json.loads(Path(path).read_text())


def validate_frozen():
    g=load(GLOBAL)
    if g['status']!='FROZEN_BEFORE_EITHER_NEW_COHORT' or len(g['methods'])!=26:raise ValueError('original26-model preobservations freeze required')
    for source in g['source_identities']:verify_identity(source)
    for method in g['methods']:
        verify_identity(method['checkpoint_identity']);verify_identity(method['source_receipt'])
        if load(method['source_receipt']['path'])['status']!='COMPLETE':raise ValueError('all original model fits must remain CLOSED COMPLETE')
    freezes={}
    for cohort,n,prefix in [('locked',512,'R4L1'),('confirmation',256,'R4L2')]:
        path=OUT/f'results/fresh/{cohort}/FROZEN_DENSITY.json';f=load(path)
        if f['status']!='FROZEN_BEFORE_TRUTH' or f['truth_read'] is not False or len(f['methods'])!=26 or f['sample_count']!=128:
            raise ValueError('sealed fresh inference freeze required')
        if f['parents']!=[f'{prefix}_{i:05d}' for i in range(n)]:raise ValueError('frozen independent parent roster changed')
        for m in f['methods']:
            verify_identity(m['checkpoint_identity'])
            if m['kind']=='cfm' and m.get('sample_nfes')!=[32]:raise ValueError('fresh CFM NFE32 frozen only')
        freezes[cohort]=(path,f)
    return g,freezes


def inert_ppc(g,freezes):
    expected=['raw_spline_s2026102501','anchor_cfm_none_s2026102501__nfe32',
              'anchor_cfm_geometry_s2026102501__nfe32','anchor_cfm_random_s2026102501__nfe32']
    if g['full_checks']['groups']!=expected or g['full_checks']['max_RHS']!=56832:raise ValueError('preobservations PPC scope mismatch')
    generated=[];evaluations=0
    for cohort,(path,f) in freezes.items():
        for group in expected:
            extended=f['parents'][:16] if cohort=='locked' else []
            evaluations+=len(f['parents'])+63*len(extended)
            generated.append(dict(id=cohort+'__'+group,directory=str(OUT/f'results/fresh/{cohort}/density_online'/group),
                online_gate=str(OUT/f'results/fresh/{cohort}/density_online/ONLINE_COMPLETE.json'),
                observations=str(OUT/'inputs/online'),parents=f['parents'],extended_parents=extended,
                density='explicit' if group=='raw_spline_s2026102501' else 'none'))
    if evaluations*8!=56832:raise AssertionError('PPC fixed index cost changed')
    request=dict(status='INERT_REQUEST_NOT_A_PHYSICS_RELEASE',released=False,preobservations_plan=identity(GLOBAL),
        full_checks=g['full_checks'],generated=generated,likelihood_samples=64,importance_weights=True,
        planned_Full_evaluations=evaluations,planned_max_RHS=56832,required_Maxwell_RHS=56832,
        fixed_sample_index=0,no_residual_selection=True,duplicates='keep every frozen index',
        CFM_density='none; numerical density not treated as exact',rawSpline_density='explicit algebraic full66 rawdx66',
        first_locked_GT_gate='main must close BOTH risk+density inference for BOTH cohorts before any first heldout truth score',
        no_execution_authorized=True,posterior_qualification='NOT_ESTABLISHED')
    if PPC.exists():
        if load(PPC)!=request:raise FileExistsError('existing inert PPC request differs; preserve history')
    else:write(PPC,request)


def main():
    if START.exists() or COMPLETE.exists():raise FileExistsError('fresh inference controller already registered; no replay')
    g,freezes=validate_frozen();inert_ppc(g,freezes)
    bound=dict(global_freeze=identity(GLOBAL),cohorts={k:identity(v[0]) for k,v in freezes.items()})
    write(START,dict(status='WAITING_FOR_PHYSICS_AND_CONTEXTS_ONLY',pid=os.getpid(),poll_seconds=60,
        bindings=bound,truth_read=False,physics_RHS=0,optimizer_updates=0))
    while not PHYSICAL.exists():time.sleep(60)
    if load(PHYSICAL).get('status')!='COMPLETE':raise ValueError('physical/context marker not COMPLETE; no inference')
    verify_identity(bound['global_freeze']);validate_frozen()
    gates={};profiles={};status,error='COMPLETE',None
    try:
        for cohort,(frozen,f) in freezes.items():
            verify_identity(bound['cohorts'][cohort])
            out=OUT/f'results/fresh/{cohort}/density_online'
            if out.exists():raise FileExistsError('no cohort inference replay: '+cohort)
            export=OUT/f'inputs/anchor_conditions/{cohort}/EXPORT_COMPLETE.json'
            if load(export).get('status')!='COMPLETE':raise ValueError('explicit observation-derived context exports not complete')
            job=f'fresh_density_{cohort}_infer_NFE32'
            args=[sys.executable,'tools/run_r4_job.py','--stage','IV','--job',job,'--max-wall','7200',
                '--rss-limit-mib','2048','--',sys.executable,'tools/evaluate_r4_density_frozen_nfe.py','infer',
                '--freeze',str(frozen),'--contexts',str(export.parent),'--out',str(out),'--device','cpu','--batch-size','16']
            env=os.environ.copy();env.update(PYTHONPATH=str(ROOT/'src'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',R4_ARRAY_CODEC='stored')
            child=subprocess.run(args,cwd=ROOT,env=env)
            if child.returncode:raise RuntimeError('cohort inference failed; no rerun: '+cohort)
            receipt=OUT/'IV/jobs'/job/'CLOSED.json';closed=load(receipt);gate=load(out/'ONLINE_COMPLETE.json')
            if closed['status']!='COMPLETE' or gate['status']!='ONLINE_COMPLETE' or gate['truth_read'] is not False:
                raise ValueError('closed observation-only inference gate required')
            if len(gate['records'])!=26*len(f['parents']):raise ValueError('missing frozen method/parent slots')
            gates[cohort]=dict(gate=identity(out/'ONLINE_COMPLETE.json'),receipt=identity(receipt),
                contexts_complete=identity(export),failures=gate['failures'],records=len(gate['records']),truth_read=False)
            profiles[cohort]=dict(wall_seconds=closed['controlled_wall_seconds'],**closed['profile'])
            write(OUT/'results/FRESH_DENSITY_ONLINE_PROGRESS.json',dict(status='RUNNING',gates=gates,profiles=profiles,truth_read=False))
        validate_frozen()
    except BaseException as exc:status,error='FAILED_NO_AUTORERUN',repr(exc);raise
    finally:
        write(COMPLETE,dict(status=status,error=error,bindings=bound,physical_context_gate=identity(PHYSICAL),
            gates=gates,profiles=profiles,all_failures_preserved=True,truth_read=False,GT_scoring='NOT_RUN_NOT_AUTHORIZED',
            physics_RHS=0,optimizer_updates=0,Maxwell_PPC='INERT_REQUEST_ONLY',PPC_request=identity(PPC)))
    print(COMPLETE,flush=True)


if __name__=='__main__':main()
