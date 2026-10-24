"""Build matched strong-context supervision and fit full66 anchor models."""
from pathlib import Path
import json,sys,subprocess,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r4_conditional_transport'
SEEDS=[2026102501,2026102502,2026102503]

def call(command):
    p=subprocess.run([sys.executable,*command],cwd=ROOT)
    if p.returncode:raise SystemExit(p.returncode)

def main():
    path=OUT/'configs/ANCHOR_NEURAL_01.json'
    if path.exists():raise FileExistsError('anchor sequence already registered')
    plan=[]
    for kind in ['spline','diag_gaussian','cfm']:
        for seed in SEEDS:plan.append(dict(kind=kind,coordinate='none',seed=seed))
        for coordinate in ['gn','random']:plan.append(dict(kind=kind,coordinate=coordinate,seed=SEEDS[0]))
    path.write_text(json.dumps(dict(plan=plan,train=1024,validation=256,updates=3000,learning_rate=5e-5,
        anchor_budget='native96+probit64 for every TRAIN/VAL object; paid StageI',
        all66=True,j_summary='none',conditional_map='raw512+anchor66; optional exact full-rank GN preconditioner',
        geometry_control='same GN singular scale, fixed rotated latent alignment',
        selection='independent validation before new locked generation; update0 eligible',no_locked_GT=True),indent=2)+'\n')
    while not (OUT/'results/CONTEXT_EXPANSION_COMPLETE.json').exists():time.sleep(10)
    packed=OUT/'inputs/strong_context_1024_256'
    call(['tools/run_r4_job.py','--stage','I','--job','pack_strong_context_1024_256','--max-wall','300','--rss-limit-mib','2048','--',
        sys.executable,'tools/pack_r4_context_training.py','--train-count','1024','--validation-count','256','--out',str(packed),'--base-chol-scale','1'])
    (OUT/'results/STRONG_CONTEXT_PACK_READY.json').write_text(json.dumps(dict(status='READY',path=str(packed)),indent=2)+'\n')
    while not (OUT/'results/OPERATOR_TRAIN_COMPLETE.json').exists():time.sleep(10)
    for row in plan:
        job=f"anchor_{row['kind']}_{row['coordinate']}_s{row['seed']}"
        call(['tools/run_r4_job.py','--stage','I','--job',job,'--training','--gpu','--max-wall','2100','--',
            sys.executable,'-m','p3r4.train_anchor_density',
            '--train',str(packed/'density_train.npz'),'--validation',str(packed/'density_validation.npz'),
            '--out',str(OUT/'models/anchor_density'/job),'--kind',row['kind'],'--seed',str(row['seed']),
            '--coordinate-preconditioner',row['coordinate'],'--updates','3000','--lr','5e-5',
            '--eval-every','100','--batch-size','128','--validation-batch-size','128','--max-seconds','1950','--device','mps'])
    (OUT/'results/ANCHOR_NEURAL_COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',fits=len(plan)),indent=2)+'\n')
if __name__=='__main__':main()
