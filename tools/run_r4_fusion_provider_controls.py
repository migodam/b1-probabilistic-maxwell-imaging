"""Test whether expensive Flow moments add anything over ordinary global features."""
from pathlib import Path
import json,sys,subprocess,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r4_conditional_transport'
def call(args):
    p=subprocess.run([sys.executable,*args],cwd=ROOT)
    if p.returncode:raise SystemExit(p.returncode)
def main():
    cfg=OUT/'configs/FUSION_PROVIDER_CONTROLS_01.json'
    if cfg.exists():raise FileExistsError('provider control release already registered')
    providers={
        'diag_gaussian':OUT/'models/density/extended_diag_gaussian_s2026102501/BEST.pt',
        'field_moments':OUT/'models/field/extended_mlp_small_moments_s2026102501/BEST.pt'}
    plan=[dict(provider=p,seed=s) for p in providers for s in [2026102501,2026102502,2026102503]]
    cfg.write_text(json.dumps(dict(released=True,plan=plan,provider_checkpoints={k:str(v) for k,v in providers.items()},
        initial_observation='Spline fusion validation .007126 versus strongreoptimizedCNN .007447; this may be ordinary global context rather than probability benefit',
        exact_K=dict(diag_gaussian=32,field_moments=0),same_risk_head_parameters=6418,train=1024,validation=256,
        updates=4000,material_target_unchanged=True,GT_not_used_for_features=True),indent=2)+'\n')
    while not (OUT/'results/RISK_FUSION_COMPLETE.json').exists():time.sleep(10)
    packed=OUT/'inputs/strong_context_1024_256'
    for provider,checkpoint in providers.items():
        dest=OUT/'inputs'/('fusion_'+provider);dest.mkdir(parents=True,exist_ok=False)
        for group,n in [('train',1024),('validation',256)]:
            parents=dest/f'{group}_parents.json';parents.write_text(json.dumps(dict(parents=[f'R2_{group}_{k:05d}' for k in range(n)]),indent=2)+'\n')
            call(['tools/run_r4_job.py','--stage','I','--job',f'fusion_{provider}_features_{group}','--max-wall','600','--rss-limit-mib','2048','--',
                sys.executable,'tools/r4_fusion_moments.py','--provider',provider,'--checkpoint',str(checkpoint),
                '--observations',str(OUT/'inputs/online'),'--parents',str(parents),'--out',str(dest/group),'--device','cpu'])
        for row in [r for r in plan if r['provider']==provider]:
            job=f"riskfusion_{provider}_s{row['seed']}"
            command=[sys.executable,'-m','p3r4.train_risk_fusion','--train',str(OUT/'inputs/joint/train32768.npz'),
                '--validation',str(OUT/'inputs/joint/validation.npz'),'--anchor-train',str(packed/'anchor_train.npz'),
                '--anchor-validation',str(packed/'anchor_validation.npz'),'--moments-train',str(dest/'train'),
                '--moments-validation',str(dest/'validation'),'--out',str(OUT/'models/risk'/job),'--control','moments',
                '--seed',str(row['seed']),'--updates','4000','--refit-feature-normalization','--device','cpu']
            call(['tools/run_r4_job.py','--stage','I','--job',job,'--training','--cpu-parallel','--max-wall','1800','--rss-limit-mib','2048','--',*command])
    (OUT/'results/FUSION_PROVIDER_CONTROLS_COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',fits=len(plan)),indent=2)+'\n')
if __name__=='__main__':main()
