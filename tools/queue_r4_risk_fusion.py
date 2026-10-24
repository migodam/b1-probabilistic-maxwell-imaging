"""One bounded material-moment fusion and information-matched spatial control."""
from pathlib import Path
import json,sys,time,subprocess
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r4_conditional_transport'
def call(args):
    p=subprocess.run([sys.executable,*args],cwd=ROOT)
    if p.returncode:raise SystemExit(p.returncode)
def main():
    cfg=OUT/'configs/RISK_FUSION_01.json'
    if cfg.exists():raise FileExistsError('fusion sequence already registered')
    plan=[dict(control=c,seed=s) for c in ['moments','permuted'] for s in [2026102501,2026102502,2026102503]]
    cfg.write_text(json.dumps(dict(plan=plan,train=1024,validation=256,updates=4000,parameters=6418,
        hypothesis='coarse full66 model material moments provide useful context to local physics-risk estimator despite unqualified model posterior',
        base_model='raw-y full66 Spline seed2026102501, frozen before this feature generation',K=32,
        target='supervised material risk only, not calibrated posterior',
        controls='fixed perparent spatial permutation of same four generated moment channels',
        feature_normalization='same TRAIN-only online-anchor refit as strongrisk baseline; four extra channels TRAIN-only',
        no_locked_GT=True,extra_teacher_full_Maxwell_RHS=0),indent=2)+'\n')
    for marker in ['STRONG_CONTEXT_PACK_READY.json','STRONG_RISK_COMPLETE.json']:
        while not (OUT/'results'/marker).exists():time.sleep(10)
    packed=OUT/'inputs/strong_context_1024_256';dest=OUT/'inputs/fusion_moments'
    dest.mkdir(parents=True,exist_ok=False)
    checkpoint=OUT/'models/density/density_spline_32768_s2026102501/BEST.pt'
    for group,n in [('train',1024),('validation',256)]:
        parents=dest/f'{group}_parents.json';parents.write_text(json.dumps(dict(parents=[f'R2_{group}_{k:05d}' for k in range(n)]),indent=2)+'\n')
        call(['tools/run_r4_job.py','--stage','I','--job','fusion_moments_'+group,'--max-wall','600','--rss-limit-mib','2048','--',
            sys.executable,'tools/r4_fusion_moments.py','--checkpoint',str(checkpoint),'--observations',str(OUT/'inputs/online'),
            '--parents',str(parents),'--out',str(dest/group),'--device','cpu'])
    for row in plan:
        job=f"riskfusion_{row['control']}_s{row['seed']}"
        command=[sys.executable,'-m','p3r4.train_risk_fusion','--train',str(OUT/'inputs/joint/train32768.npz'),
            '--validation',str(OUT/'inputs/joint/validation.npz'),'--anchor-train',str(packed/'anchor_train.npz'),
            '--anchor-validation',str(packed/'anchor_validation.npz'),'--moments-train',str(dest/'train'),
            '--moments-validation',str(dest/'validation'),'--out',str(OUT/'models/risk'/job),
            '--control',row['control'],'--seed',str(row['seed']),'--updates','4000',
            '--refit-feature-normalization','--device','cpu']
        call(['tools/run_r4_job.py','--stage','I','--job',job,'--training','--cpu-parallel','--max-wall','1800','--rss-limit-mib','2048','--',*command])
    (OUT/'results/RISK_FUSION_COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',fits=len(plan)),indent=2)+'\n')
if __name__=='__main__':main()
