"""Matched1024/256 strong-anchor risk fits with controlled geometry ablations."""
from pathlib import Path
import sys,json,subprocess,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r4_conditional_transport'
def main():
    cfg=OUT/'configs/STRONG_RISK_01.json'
    if cfg.exists():raise FileExistsError('risk sequence already registered')
    plan=[dict(mode=mode,seed=seed) for mode in ['spectral','decoder_only','random'] for seed in [2026102501,2026102502,2026102503]]
    cfg.write_text(json.dumps(dict(plan=plan,train=1024,validation=256,updates=4000,parameters=5842,
        feature_scaling='refit on identical TRAIN contexts for every frame',
        same_initialization_optimizer_and_information_budget=True,no_locked_GT=True),indent=2)+'\n')
    for marker in ['STRONG_CONTEXT_PACK_READY.json','ANCHORBIN_REVISION2_COMPLETE.json']:
        while not (OUT/'results'/marker).exists():time.sleep(10)
    while not (OUT/'results/II/interval_chain_tight_pair/COMPLETE.json').exists():time.sleep(10)
    packed=OUT/'inputs/strong_context_1024_256'
    for row in plan:
        job=f"strongrisk_{row['mode']}_s{row['seed']}"
        command=[sys.executable,'-m','p3r4.train_risk','--kind','local_cnn','--context-mode',row['mode'],
            '--refit-feature-normalization','--train',str(OUT/'inputs/joint/train32768.npz'),
            '--validation',str(OUT/'inputs/joint/validation.npz'),'--train-context',str(packed/'anchor_train.npz'),
            '--validation-context',str(packed/'anchor_validation.npz'),'--out',str(OUT/'models/risk'/job),
            '--seed',str(row['seed']),'--updates','4000','--device','cpu']
        p=subprocess.run([sys.executable,'tools/run_r4_job.py','--stage','I','--job',job,'--training','--cpu-parallel',
            '--max-wall','1800','--rss-limit-mib','2048','--',*command],cwd=ROOT)
        print(json.dumps(dict(job=job,returncode=p.returncode)),flush=True)
        if p.returncode:raise SystemExit(p.returncode)
    (OUT/'results/STRONG_RISK_COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',fits=len(plan)),indent=2)+'\n')
if __name__=='__main__':main()
