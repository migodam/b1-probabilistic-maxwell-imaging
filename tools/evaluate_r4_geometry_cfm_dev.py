"""Authorized zero-physics exposed-DEV comparison, after six CLOSED fits.

Freeze -> CPU1 inference -> different-PID offline truth scoring. No new fits,
Full calls, CFM density audits, candidate filtering or checkpoint mutation.
"""
from pathlib import Path
import csv,json,os,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,write
from p3r4.train_physics_transport import identity,verify_identity
SEEDS=[2026102501,2026102502,2026102503]
DEST=OUT/'results/geometry_weighted_cfm_dev_v1'


def load(path):return json.loads(Path(path).read_text())


def call(job,seconds,args):
    cmd=[sys.executable,'tools/run_r4_job.py','--stage','I','--job',job,'--max-wall',str(seconds),
         '--rss-limit-mib','2048','--',sys.executable,*args]
    env=os.environ.copy();env.update(PYTHONPATH=str(ROOT/'src'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',R4_ARRAY_CODEC='stored')
    result=subprocess.run(cmd,cwd=ROOT,env=env)
    if result.returncode:raise RuntimeError('evaluation failed, no replay: '+job)
    return identity(OUT/'I/jobs'/job/'CLOSED.json')


def closed_models():
    methods=[];sources=[];validation=[]
    for arm in ('none','geometry','random'):
        for seed in SEEDS:
            job=f'anchor_cfm_none_s{seed}' if arm=='none' else f'geometry_cfm_{arm}_s{seed}'
            directory=OUT/'models'/('anchor_density' if arm=='none' else 'geometry_weighted_cfm_v1')/job
            receipt=OUT/'I/jobs'/job/'CLOSED.json'
            if not receipt.exists() or load(receipt)['status']!='COMPLETE':raise ValueError('CLOSED COMPLETE required before any checkpoint read: '+job)
            r=load(directory/'RESULT.json')
            if r.get('completed_updates')!=3000 or r['kind']!='cfm' or (r['train_size'],r['validation_size'])!=(1024,256):raise ValueError('matched3000/full66 anchorCFM source required')
            methods.append(dict(id=f'{arm}_s{seed}',family='cfm_'+arm,seed=seed,checkpoint=str(directory/'BEST.pt')))
            sources.append(dict(id=methods[-1]['id'],checkpoint=identity(directory/'BEST.pt'),result=identity(directory/'RESULT.json'),receipt=identity(receipt)))
            validation.append(dict(arm=arm,seed=seed,parameter_count=r['parameter_count'],selected_update=r['selected_update'],
                selected_raw_latent_energy=r['selected_validation_metric'],selected_initialization=r['selected_initialization']))
    return methods,sources,validation


def summarize(validation,receipts):
    import numpy as np
    from p3r4.probability_eval import bootstrap_mean
    results=load(DEST/'scores/RESULTS.json')
    rows=[];phase=[]
    for family,data in results['families'].items():
        arm=family.split('__')[0].removeprefix('cfm_');nfe=int(family.split('__nfe')[1])
        if sorted(data['training_seeds'])!=SEEDS:raise ValueError('family must average exactly3 seed losses per parent')
        row=dict(arm=arm,NFE=nfe,TRAIN=1024,VAL=256,DEV=128,seeds=';'.join(map(str,SEEDS)),
                 validation_raw_latent_energy_seed_mean=float(np.mean([v['selected_raw_latent_energy'] for v in validation if v['arm']==arm])))
        for key,summary in data['metrics'].items():
            row[key]=summary.get('mean_all_parents')
        rows.append(row)
        variant_arrays=[]
        for seed in SEEDS:
            with np.load(DEST/f'scores/{arm}_s{seed}__nfe{nfe}_per_parent.npz',allow_pickle=False) as z:
                variant_arrays.append(z['phase_brier_per_threshold'])
        values=np.stack(variant_arrays).mean(0)
        for j,t in enumerate([.25,.5,.75]):phase.append(dict(arm=arm,NFE=nfe,threshold=t,**bootstrap_mean(values[:,j])))
    rows.sort(key=lambda r:(r['NFE'],r['arm']))
    columns=sorted(set().union(*(r.keys() for r in rows)))
    with (DEST/'SUMMARY.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader();writer.writerows(rows)
    write(DEST/'SUMMARY.json',dict(status='COMPLETE_EXPOSED_DEV_DIAGNOSTICS',validation=validation,families=rows,
        phase_Brier_by_threshold=phase,comparisons=results['comparisons'],all_online_failures=results['all_online_failures'],
        successful_parents_per_variant={k:v['successful_parents'] for k,v in results['methods'].items()},
        aggregation='average per-parent seed losses; no independent-seed sample-size inflation or image ensemble',
        CFM_exact_density='NOT_EVALUATED_NOT_AVAILABLE',Full_Phi='NOT_RUN',Gamma_SBC='NOT_RUN',
        posterior_qualification='NOT_ESTABLISHED',physics_RHS=0,optimizer_updates=0,receipts=receipts))
    selected=[r for r in rows if r['NFE']==32]
    lines=['# Bounded geometry CFM: actual exposed-DEV diagnostics','',
        'All9 checkpoints were CLOSED COMPLETE before freeze. DEV128 is exposed; no new Maxwell likelihood, Gamma/SBC, exact CFM density, or posterior qualification was evaluated.',
        '', 'Three seed losses are averaged within each parent;128 parents are the inference unit. All failures remain explicit. NFE16/32/64, K128 and seed20261025031 were fixed before inference.','',
        '| NFE32 arm | VAL raw66 energy (seed mean) | DEV material energy | DEV sample-mean MSE | ROI90 empirical coverage | ROI90 randomized rank coverage |',
        '|---|---:|---:|---:|---:|---:|']
    def f(value):return 'NA' if value is None else f'{value:.8g}'
    for r in selected:
        keys=['validation_raw_latent_energy_seed_mean','material_energy_score','raw_sample_mean_full_mse',
              'material_real_roi__empirical_interval_coverage_90','material_real_roi__randomized_rank_coverage_90']
        lines.append('| '+r['arm']+' | '+' | '.join(f(r.get(k)) for k in keys)+' |')
    lines+=['','Full NFE/phase/coverage metrics, intervals, fixed projections and paired comparisons are in SUMMARY.json/CSV and scores/RESULTS.json. Pointwise diagnostic intervals do not provide simultaneous certification across these metrics.','',
        'Inference and offline truth scoring used separate PIDs after the sealed ONLINE_COMPLETE gate. Each stage has an analysis-only Book; zero new physical RHS and no checkpoint modification.']
    (DEST/'REPORT.md').write_text('\n'.join(lines)+'\n')


def main():
    started=time.monotonic()
    if DEST.exists():raise FileExistsError('no automatic evaluation replay')
    DEST.mkdir(parents=True)
    marker=OUT/'models/geometry_weighted_cfm_v1/COMPLETE.json'
    write(DEST/'START.json',dict(status='WAITING_ONLY_FOR_CLOSED_MARKER',authorization='parent explicit zero-physics DEV evaluation',
        wait_check_seconds=60,maximum_wait_seconds=1800,physics_RHS=0,optimizer_updates=0))
    while not marker.exists():
        if time.monotonic()-started>1800:raise TimeoutError('no six-fit COMPLETE marker; no inference')
        time.sleep(60)
    if load(marker).get('status')!='COMPLETE':raise ValueError('six-fit controller not complete')
    methods,sources,validation=closed_models()
    parents=[f'R4DEV_{i:05d}' for i in range(128)]
    contexts=OUT/'inputs/anchor_conditions/dev'
    compared=[dict(candidate_family=f'cfm_geometry__nfe{nfe}',reference_family=f'cfm_{reference}__nfe{nfe}',metric=metric)
        for nfe in (16,32,64) for reference in ('none','random') for metric in ('material_energy_score','raw_sample_mean_full_mse')]
    spec=dict(methods=methods,parents=parents,sample_seed=20261025031,comparisons=compared,
        cohort='EXPOSED_DEV128',no_CFM_density_audit=True,zero_physics=True,
        aggregation='3 seed losses within parent; no seed as independent inference unit')
    write(DEST/'METHODS_SPEC.json',spec)
    sources += [identity(ROOT/p) for p in ['tools/evaluate_r4_density.py','src/p3r4/probability_eval.py','src/p3r4/density.py','src/p3r4/anchor_density.py']]
    write(DEST/'PREINFERENCE_FREEZE.json',dict(status='FROZEN_BEFORE_INFERENCE',methods=sources,validation=validation,
        contexts={p:identity(contexts/(p+'.npz')) for p in parents},spec=identity(DEST/'METHODS_SPEC.json'),
        truth_read=False,physics_RHS=0,CFM_numeric_density_not_frozen=True))
    receipts=[]
    receipts.append(call('geometry_cfm_DEV_freeze_v1',180,['tools/evaluate_r4_density.py','freeze','--methods',str(DEST/'METHODS_SPEC.json'),'--out',str(DEST/'FROZEN.json')]))
    receipts.append(call('geometry_cfm_DEV_infer_v1',3000,['tools/evaluate_r4_density.py','infer','--freeze',str(DEST/'FROZEN.json'),
        '--contexts',str(contexts),'--out',str(DEST/'online'),'--device','cpu','--batch-size','16']))
    receipts.append(call('geometry_cfm_DEV_score_v1',1800,['tools/evaluate_r4_density.py','score','--freeze',str(DEST/'FROZEN.json'),
        '--online',str(DEST/'online'),'--truth',str(OUT/'inputs/offline_only'),'--out',str(DEST/'scores'),'--device','cpu']))
    for row in sources:
        if 'checkpoint' in row:verify_identity(row['checkpoint'])
    summarize(validation,receipts)
    write(DEST/'COMPLETE.json',dict(status='COMPLETE',elapsed_controller_wall_seconds=time.monotonic()-started,
        physics_RHS=0,optimizer_updates=0,receipts=receipts,posterior_qualification='NOT_ESTABLISHED'))
    print(DEST,flush=True)


if __name__=='__main__':main()
