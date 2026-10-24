"""Inert six-fit bounded-geometry CFM plan; no waits or implicit release."""
from pathlib import Path
import json,os,sys,subprocess
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,write
from p3r4.train_physics_transport import identity,verify_identity
REGISTRY=OUT/'configs/GEOMETRY_WEIGHTED_CFM_V1.json'
SEEDS=[2026102501,2026102502,2026102503]
SOURCES=['src/p3r4/train_geometry_cfm.py','tools/queue_r4_geometry_cfm.py',
 'src/p3r4/train_anchor_density.py','src/p3r4/train_density.py','src/p3r4/anchor_density.py',
 'src/p3r4/density.py','src/p3r4/conditional_coordinates.py','tools/run_r4_job.py']


def build_registry(cache):
    baselines=[]
    for seed in SEEDS:
        name=f'anchor_cfm_none_s{seed}'
        receipt=OUT/'I/jobs'/name/'CLOSED.json';checkpoint=OUT/'models/anchor_density'/name/'BEST.pt'
        if json.loads(receipt.read_text())['status']!='COMPLETE':raise ValueError('original baseline not complete')
        baselines.append(dict(seed=seed,checkpoint=identity(checkpoint),receipt=identity(receipt)))
    return dict(format='p3r4_geometry_weighted_CFM_release_v1',released=False,
        paired_VI_complete=str(OUT/'results/physics_transport_paired_v2/CONTROLLER_RESULT.json'),
        train=identity(OUT/'inputs/strong_context_1024_256/density_train.npz'),
        validation=identity(OUT/'inputs/strong_context_1024_256/density_validation.npz'),
        metrics=str(cache),metric_cache=identity(cache/'CACHE_COMPLETE.json'),baselines=baselines,
        source_identities=[identity(ROOT/p) for p in SOURCES],
        fits=[dict(seed=s,kind=k) for s in SEEDS for k in ('geometry','random')],
        updates=3000,learning_rate=5e-5,batch_size=128,parameter_count=691022,
        validation_selection='same baseline rawlatent energy; K4,NFE16,every100,update0eligible; weighted FM diagnostic not selection',
        loss='networkvelocity ERROR.T W_r(y) ERROR /66; Jr=Jx*T; SPD eigenvalues1..4',
        raw_metric='T^{-T}W_rT^{-1}, bounded relative to baselineT^{-T}T^{-1}; no raw-Euclidean1..4 claim',
        random_seed=20261025191,initialization='fresh, identical seed/model/base/normalizer/noise streams as original CFM_none',
        maximum_new_Maxwell_RHS=0,maximum_fits=6,wall_cap_per_fit=600,inner_max_seconds=540,
        inference='existing full66 anchor CFM checkpoint format; numeric density remains unqualified; no exactMH',
        output_root=str(OUT/'models/geometry_weighted_cfm_v1'),no_locked_GT=True,no_autoexpansion=True)


def validate(r,execute=False):
    if r['format']!='p3r4_geometry_weighted_CFM_release_v1' or r['fits']!=[dict(seed=s,kind=k) for s in SEEDS for k in ('geometry','random')]:raise ValueError('fixed six-fit protocol mismatch')
    for key in ('train','validation','metric_cache'):verify_identity(r[key])
    for row in r['source_identities']:verify_identity(row)
    for row in r['baselines']:
        verify_identity(row['checkpoint']);verify_identity(row['receipt'])
    if execute:
        if r.get('released') is not True:raise ValueError('inert until explicit parent release')
        gate=Path(r['paired_VI_complete'])
        if not gate.exists() or json.loads(gate.read_text()).get('status')!='COMPLETE':raise ValueError('paired VI not closed COMPLETE; no blocking wait')
    return r


def main():
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--execute',action='store_true');a=p.parse_args()
    r=validate(json.loads(REGISTRY.read_text()),a.execute)
    if not a.execute:
        print(json.dumps(dict(status='INERT' if not r['released'] else 'RELEASED_CHECK_ONLY',fits=6,physical_RHS=0,optimizer_updates=0)));return
    out=Path(r['output_root']);out.mkdir(parents=True,exist_ok=True)
    marker=out/'CONTROLLER_START.json'
    if marker.exists():raise FileExistsError('no automatic weighted-CFM replay')
    write(marker,dict(registry=identity(REGISTRY)))
    frozen=identity(REGISTRY)
    for fit in r['fits']:
        verify_identity(frozen);validate(r,True)
        job=f"geometry_cfm_{fit['kind']}_s{fit['seed']}"
        cmd=[sys.executable,'tools/run_r4_job.py','--stage','I','--job',job,'--training','--gpu','--max-wall','600',
            '--rss-limit-mib','2048','--',sys.executable,'-m','p3r4.train_geometry_cfm','--release',str(REGISTRY),
            '--train',r['train']['path'],'--validation',r['validation']['path'],'--metrics',r['metrics'],'--out',str(out/job),
            '--kind',fit['kind'],'--seed',str(fit['seed']),'--updates','3000','--lr','5e-5','--batch-size','128','--device','mps','--max-seconds','540']
        env=os.environ.copy();env['R4_ARRAY_CODEC']='stored'
        child=subprocess.run(cmd,cwd=ROOT,env=env)
        if child.returncode:raise RuntimeError('fit failed; no retry')
        result=json.loads((out/job/'RESULT.json').read_text())
        if result['completed_updates']!=3000:raise RuntimeError('partial fit; no expansion/retry')
    write(out/'COMPLETE.json',dict(status='COMPLETE',fits=6,new_Maxwell_RHS=0,registry=frozen))


if __name__=='__main__':main()
