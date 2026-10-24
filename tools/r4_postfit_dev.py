"""Inert post-fit DEV plan; separate released generator/Full/MH stages.

No fit or physical action is launched by check/build. Every model must have a
CLOSED COMPLETE receipt before freeze; live checkpoints are never inspected.
"""
from __future__ import annotations
import argparse,copy,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'tools'))
from p3r4.core import OUT,write,save
from p3r4.train_physics_transport import identity,verify_identity
from queue_r4_physics_transport import fit_paths

PROTOCOL=OUT/'configs/POSTFIT_DEV_V2_INERT.json'
ROOTOUT=OUT/'results/postfit_dev_v2'
PARENTS=[f'R4DEV_{i:05d}' for i in range(128)]
MH_PARENTS=PARENTS[64:68]
SOURCES=['tools/r4_postfit_dev.py','tools/evaluate_r4_density.py','src/p3r4/probability_eval.py',
 'src/p3r4/evaluation.py','tools/run_r4_probability_audit.py','src/p3r4/probability_physics.py',
 'tools/r4_corrected_transport.py','src/p3r4/corrected_transport.py','tools/run_r4_job.py',
 'src/p3r4/core.py','src/p3r4/density.py','src/p3r4/anchor_density.py','src/p3r4/conditional_coordinates.py']


def load(path):return json.loads(Path(path).read_text())


def build_protocol():
    registry_path=OUT/'configs/PHYSICS_TRANSPORT_PAIRED_V2_RUNTIME_RECOVERY.json'
    registry=load(registry_path)
    if registry.get('released') is not True:raise ValueError('parent-released V2 campaign required')
    methods=[]
    for family in registry['families']:
        source=family['sources'][0]
        for arm in ('initial','physics','supervised'):
            if arm=='initial':
                checkpoint=source['checkpoint']['path'];receipt=source['original_receipt']['path'];result=source['original_result']['path']
            else:
                _,out,jobdir=fit_paths(registry,family,source,arm)
                checkpoint=str(out/'BEST.pt');receipt=str(jobdir/'CLOSED.json');result=str(out/'RESULT.json')
            methods.append(dict(id=family['id']+'__'+arm,family=family['id'],arm=arm,seed=source['seed'],
                checkpoint=checkpoint,receipt=receipt,result=result,train_size=family['train_size'],validation_size=family['validation_size']))
    conditions=OUT/'inputs/anchor_conditions/dev'
    exported=load(conditions/'EXPORT_COMPLETE.json')
    if exported.get('status')!='COMPLETE':raise ValueError('explicit DEV export not COMPLETE')
    old_job='transport_raw_spline_physics_s2026102501'
    old_receipt=OUT/'I/jobs'/old_job/'CLOSED.json'
    old_result=OUT/'models/physics_transport/paired_v1'/old_job/'RESULT.json'
    return dict(format='p3r4_postfit_dev_v2_inert_v1',released=False,
        parent_release='separate explicit phase release, bound to this immutable protocol identity',
        paired_campaign=identity(registry_path),methods=methods,parents=PARENTS,
        contexts=str(conditions),contexts_complete=identity(conditions/'EXPORT_COMPLETE.json'),
        context_identities={p:identity(conditions/(p+'.npz')) for p in PARENTS},
        original_observations={p:identity(OUT/'inputs/online'/(p+'.npz')) for p in PARENTS},
        target_config=identity(OUT/'configs/TARGET.json'),
        source_identities=[identity(ROOT/p) for p in SOURCES],
        generator=dict(K=128,fixed_sample_index=0,sample_seed=20261025031,device='cpu',batch_size=16,
            full66_proper_NLL=True,material_energy=True,raw_sample_mean_risk=True,
            projections='fixed16 real material+4 spatial',coverage=[.5,.8,.9,.95],phase_events=[.25,.5,.75],
            bootstrap=dict(draws=4000,seed=20261025099),
            primary_comparisons='physics vs supervised material energy, separately in each3 family; other scores descriptive',
            no_GT_until_all_inference_complete=True,separate_PID_score=True),
        paid_probability=dict(sample0_parents=PARENTS,first64_parents=PARENTS[:8],sample_indices=list(range(64)),
            evaluations_per_model=632,models=9,planned_Full_evaluations=5688,planned_RHS=45504,
            source='sealed original K128 draws; explicit exact algebraic full66 logq only',
            keep_current_field=True,no_residual_selection=True,duplicate_policy='retain every fixed draw',
            max_wall=2400,stage='I',job='postfit_dev_v2_paid_probability',posterior_qualification='NOT_ESTABLISHED'),
        mh=dict(parents=MH_PARENTS,chains=4,starts='one online anchor+3 independent q draws common across all kernels',
            warmup=0,production=128,kernels=['pcn','independence_q','mixture_half'],beta=.002,mixture_q_probability=.5,
            seed=20261025141,planned_Full_evaluations=6192,planned_RHS=49536,max_wall=2400,stage='III',
            job='postfit_dev_v2_corrected_MH',checkpoint_selection='parent-frozen validation-only binding BEFORE DEV inference',
            cross_panel_caution='raw properNLL VAL2048 vs anchorVAL256; no automatic cross-panel winner or DEV winner',
            exact_ratio_scope='formal rawdx66 algebraic q and prior pCN; floating-point precision limits not a machine-exact certificate',
            repeats='retain all rejection states/endpoints; no thinning, endpoint replacement or residual selection',
            posterior_qualification='NOT_ESTABLISHED'),
        root=str(ROOTOUT),truth=str(OUT/'inputs/offline_only'),
        V1_failure_retained=dict(job=old_job,receipt=identity(old_receipt),result=identity(old_result),
            controlled_wall_seconds=load(old_receipt)['controlled_wall_seconds'],
            costs=load(old_result)['costs'],completed_updates=load(old_result)['completed_updates'],
            scope='pre-physics MPS float64 reduction startup failure; retained, never merged with V2'),
        no_CFM=True,no_locked_GT=True,no_training=True,no_automatic_followups=True)


def validate_protocol(p):
    if p.get('format')!='p3r4_postfit_dev_v2_inert_v1' or p.get('released') is not False:
        raise ValueError('immutable inert protocol required; release in a separate phase file')
    if p['parents']!=PARENTS or len(p['methods'])!=9 or p['mh']['parents']!=MH_PARENTS:
        raise ValueError('fixed9 models/DEV128/MH64..67 required')
    verify_identity(p['paired_campaign']);verify_identity(p['contexts_complete']);verify_identity(p['target_config'])
    for r in p['source_identities']:verify_identity(r)
    for pnt in PARENTS:
        verify_identity(p['context_identities'][pnt]);verify_identity(p['original_observations'][pnt])
    return p


def closed_methods(p,*,require_all):
    """Receipt before checkpoint/result; no live/partial data access."""
    resolved=[];pending=[]
    for row in p['methods']:
        path=Path(row['receipt'])
        if not path.exists():pending.append(dict(id=row['id'],status='PENDING_NO_LIVE_CHECKPOINT_READ'));continue
        receipt=load(path)
        if receipt['status']!='COMPLETE':
            pending.append(dict(id=row['id'],status='HOLD_NONCOMPLETE_RECEIPT',receipt_status=receipt['status']));continue
        r=load(row['result'])
        if row['arm']!='initial' and (r.get('status')!='COMPLETE' or r.get('completed_updates')!=3000):
            raise ValueError('refinement not complete3000 updates: '+row['id'])
        if (r['train_size'],r['validation_size'])!=(row['train_size'],row['validation_size']) or r['kind'] not in ('spline','diag_gaussian'):
            raise ValueError('closed source data/family mismatch')
        item=copy.deepcopy(row)
        item.update(checkpoint_identity=identity(row['checkpoint']),receipt_identity=identity(path),result_identity=identity(row['result']),
            kind=r['kind'],selected_update=r['selected_update'],validation_criterion=r.get('criterion'),
            validation_selected_metric=r.get('selected_validation_metric',r.get('validation_nll')))
        resolved.append(item)
    if require_all and pending:raise ValueError('all9 CLOSED COMPLETE sources required; pending='+str(pending))
    return resolved,pending


def phase_release(path,p,phase):
    r=load(path)
    if r.get('released') is not True or r.get('phase')!=phase or r.get('protocol')!=identity(PROTOCOL):
        raise ValueError('explicit parent phase release bound to immutable protocol required')
    validate_protocol(p)
    return r


def call(command):
    env=os.environ.copy();env.update(PYTHONPATH=str(ROOT/'src'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',R4_ARRAY_CODEC='stored')
    child=subprocess.run([sys.executable,*command],cwd=ROOT,env=env)
    if child.returncode:raise RuntimeError('stage failed; no automatic replay: '+str(command))


def analysis(job,wall,command):
    call(['tools/run_r4_job.py','--stage','I','--job',job,'--max-wall',str(wall),'--rss-limit-mib','2048','--',sys.executable,*command])


def validate_binding(release,methods):
    binding=release.get('MH_validation_binding',{})
    if binding.get('validation_only') is not True or binding.get('before_DEV_inference') is not True or not binding.get('reason'):
        raise ValueError('parent must freeze validation-only MH choice and rationale BEFORE inference')
    match=[m for m in methods if m['id']==binding.get('method_id')]
    if len(match)!=1:raise ValueError('MH choice must be one of9 CLOSED seed1 models')
    return match[0]


def run_generator(p,release_path):
    release=phase_release(release_path,p,'generator')
    methods,_=closed_methods(p,require_all=True)
    chosen=validate_binding(release,methods)
    root=Path(p['root'])
    if root.exists():raise FileExistsError('fresh postfit generator output only; no replay')
    root.mkdir(parents=True)
    final=dict(format='p3r4_postfit_final_preinference_v1',status='FROZEN_BEFORE_INFERENCE_AND_NEW_GT_SCORE',
        protocol=identity(PROTOCOL),generator_release=identity(release_path),resolved_methods=methods,
        MH_validation_binding=release['MH_validation_binding'],MH_selected_checkpoint=chosen,
        freeze_time_ns=time.time_ns(),truth_read=False,physics_RHS=0,optimizer_updates=0)
    write(root/'FINAL_PREINFERENCE_SPEC.json',final)
    spec=dict(methods=[dict(id=m['id'],family=m['id'],seed=m['seed'],checkpoint=m['checkpoint']) for m in methods],
        parents=PARENTS,sample_seed=p['generator']['sample_seed'],comparisons=[
            dict(candidate_family=f+'__physics',reference_family=f+'__supervised',metric='material_energy_score')
            for f in ['raw_spline','raw_diag_gaussian','anchor_spline_none']])
    write(root/'GENERATOR_SPEC.json',spec)
    frozen=root/'GENERATOR_FROZEN.json'
    analysis('postfit_dev_v2_freeze',180,['tools/evaluate_r4_density.py','freeze','--methods',str(root/'GENERATOR_SPEC.json'),'--out',str(frozen)])
    analysis('postfit_dev_v2_infer',1800,['tools/evaluate_r4_density.py','infer','--freeze',str(frozen),
        '--contexts',p['contexts'],'--out',str(root/'online'),'--device','cpu','--batch-size','16'])
    gate=load(root/'online/ONLINE_COMPLETE.json')
    if gate['status']!='ONLINE_COMPLETE' or gate['failures']:raise RuntimeError('HOLD failed generator slots; no source/model replacement')
    # Separate subprocess, after the COMPLETE inference Book, performs first
    # new GT access. Controller never reads truth arrays itself.
    analysis('postfit_dev_v2_score',1800,['tools/evaluate_r4_density.py','score','--freeze',str(frozen),
        '--online',str(root/'online'),'--truth',p['truth'],'--out',str(root/'scores'),'--device','cpu'])
    write(root/'GENERATOR_COMPLETE.json',dict(status='COMPLETE',truth_read_by_controller=False,
        final_spec=identity(root/'FINAL_PREINFERENCE_SPEC.json'),online_gate=identity(root/'online/ONLINE_COMPLETE.json')))


def check_final(p):
    root=Path(p['root']);final=load(root/'FINAL_PREINFERENCE_SPEC.json')
    if final['protocol']!=identity(PROTOCOL):raise ValueError('postfit protocol changed')
    validate_protocol(p)
    for m in final['resolved_methods']:
        for key in ('checkpoint_identity','receipt_identity','result_identity'):verify_identity(m[key])
        if load(m['receipt'])['status']!='COMPLETE':raise ValueError('source fit not complete')
    return root,final


def run_probability(p,release_path):
    release=phase_release(release_path,p,'paid_probability')
    root,final=check_final(p)
    if release.get('final_spec')!=identity(root/'FINAL_PREINFERENCE_SPEC.json'):
        raise ValueError('paid release must bind final preinference spec')
    if load(root/'GENERATOR_COMPLETE.json')['status']!='COMPLETE':raise ValueError('generator evaluation incomplete')
    folder=root/'paid_probability'
    if folder.exists():raise FileExistsError('no paid replay')
    folder.mkdir()
    spec=dict(generated=[dict(id=m['id'],directory=str(root/'online'/m['id']),online_gate=str(root/'online/ONLINE_COMPLETE.json'),
            observations=str(OUT/'inputs/online'),parents=PARENTS,extended_parents=PARENTS[:8],density='explicit') for m in final['resolved_methods']],
        likelihood_samples=64,importance_weights=True,cohort='FIXED_EXPOSED_DEV128_POSTFIT_V2',
        no_residual_selection=True,required_Maxwell_RHS=45504,posterior_qualification='NOT_ESTABLISHED')
    write(folder/'SPEC.json',spec)
    frozen=folder/'FROZEN.json'
    analysis('postfit_dev_v2_paid_freeze',180,['tools/run_r4_probability_audit.py','freeze','--spec',str(folder/'SPEC.json'),'--out',str(frozen)])
    if load(frozen)['planned_max_RHS']!=45504:raise ValueError('fixed paid roster cost changed')
    call(['tools/run_r4_probability_audit.py','run','--freeze',str(frozen),'--out',str(folder/'physics'),
        '--stage','I','--job',p['paid_probability']['job'],'--max-rhs','45504','--max-wall','2400'])
    analysis('postfit_dev_v2_paid_score',600,['tools/run_r4_probability_audit.py','score','--freeze',str(frozen),
        '--online',str(folder/'physics'),'--truth',p['truth'],'--out',str(folder/'scores')])


def prepare_mh(p,final,folder):
    """Pure freeze using the same explicit allowlisted DEV condition exports."""
    import numpy as np
    from r4_corrected_transport import get_model
    method=final['MH_selected_checkpoint']
    loader='anchor' if method['family']=='anchor_spline_none' else 'raw'
    model=get_model(method['checkpoint'],loader)
    cfg=load(OUT/'configs/TARGET.json');rows=[]
    folder.mkdir(parents=True)
    for parent in MH_PARENTS:
        context=Path(p['contexts'])/(parent+'.npz')
        with np.load(context,allow_pickle=False) as z:
            y=z['observation'];x0=z['context_x']
        original=p['original_observations'][parent]
        with np.load(original['path'],allow_pickle=False) as z:
            if set(z.files)!={'observation','target_id','sigma'} or not np.array_equal(y,z['observation']):raise ValueError('original exact observation mismatch')
            if y.dtype!=np.complex128 or str(z['target_id'])!=cfg['target_id'] or float(z['sigma'])!=cfg['noise']['sigma']:raise ValueError('strict original observation target')
        if x0.shape!=(66,) or not np.isfinite(x0).all():raise ValueError('full66 explicit anchor needed')
        destination=folder/'conditions'/(parent+'.npz')
        save(destination,parent=np.array(parent),observation=y,context_x=x0)
        rows.append(dict(parent=parent,export=identity(destination),observation=original,
            online_context=p['context_identities'][parent],online_sidecar=p['contexts_complete'],
            q_condition='raw512+explicit frozen context_x66' if loader=='anchor' else 'raw512 only'))
    frozen=dict(status='FROZEN_CORRECTED_TRANSPORT',experiment='corrected_transport_full66_v1',checkpoint=method['checkpoint_identity'],
        loader=loader,kind=model.payload['kind'],parents=MH_PARENTS,rows=rows,chains=4,
        kernels=['pcn','independence_q','mixture_half'],warmup=0,production=128,seed=20261025141,beta=.002,mixture_q_probability=.5,
        planned_Full_evaluations=6192,planned_RHS=49536,target_config=cfg,
        initialization='same online anchor+3q starts common across kernels; starts Full paid separately per kernel',
        selection=final['MH_validation_binding'],GT_read=False,new_J_RHS=0,no_autoexpand=True,
        q_density='full66 explicit rawdx66, formal MH ratios; float32 arithmetic is not machine-exact certification',
        posterior_qualification='NOT_ESTABLISHED')
    write(folder/'FROZEN.json',frozen)


def run_mh(p,release_path):
    release=phase_release(release_path,p,'corrected_MH')
    root,final=check_final(p)
    if release.get('final_spec')!=identity(root/'FINAL_PREINFERENCE_SPEC.json'):raise ValueError('MH release must bind final preinference selection')
    if release.get('experiment')!='corrected_transport_full66_v1' or release.get('campaign_id')!='P3LR4-20261025':raise ValueError('existing exact-scope MH runner release fields required')
    prepared=root/'mh_prepared'
    if prepared.exists():raise FileExistsError('no MH pilot replay')
    analysis('postfit_dev_v2_mh_freeze',180,['tools/r4_postfit_dev.py','freeze_MH','--release',str(release_path)])
    call(['tools/r4_corrected_transport.py','run','--prepared',str(prepared),'--release',str(release_path),
        '--out',str(root/'mh_physics'),'--stage','III','--job',p['mh']['job'],'--max-rhs','49536','--max-wall','2400'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['check','generator','paid_probability','corrected_MH','freeze_MH'])
    parser.add_argument('--release',type=Path)
    args=parser.parse_args();p=validate_protocol(load(PROTOCOL))
    if args.command=='check':
        methods,pending=closed_methods(p,require_all=False)
        print(json.dumps(dict(status='INERT',CLOSED_COMPLETE=len(methods),pending=pending,subprocesses=0,physics_RHS=0)))
    elif args.command=='generator':run_generator(p,args.release)
    elif args.command=='paid_probability':run_probability(p,args.release)
    elif args.command=='freeze_MH':
        release=phase_release(args.release,p,'corrected_MH');root,final=check_final(p)
        if release.get('final_spec')!=identity(root/'FINAL_PREINFERENCE_SPEC.json'):raise ValueError('MH frozen selection mismatch')
        if (root/'mh_prepared').exists():raise FileExistsError('fresh frozen MH pilot only')
        prepare_mh(p,final,root/'mh_prepared')
    else:run_mh(p,args.release)


if __name__=='__main__':main()
