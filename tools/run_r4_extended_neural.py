"""Prelocked additional capacity, measurement-risk and TRAIN-frame controls.

Serial training slot; no online solver calls. Does not open locked GT or change
the target. Every fit has its own append-only resource receipt.
"""
from pathlib import Path
import sys,json,subprocess,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r4_conditional_transport'
SEEDS=[2026102501,2026102502,2026102503]

def main():
    plan=[]
    for kind in ['diag_gaussian','cfm_matched','triangular_maxwell','triangular_random']:
        for seed in SEEDS:plan.append(dict(group='density',kind=kind,seed=seed,n=32768,updates=6000))
    for kind in ['mlp_small','cnn_small','mlp_small_moments','cnn_small_moments']:
        for seed in SEEDS:plan.append(dict(group='field',kind=kind,seed=seed,n=32768,updates=4000))
    cfg=OUT/'configs/EXTENDED_NEURAL_01.json'
    if cfg.exists():raise FileExistsError('never reset/replay a registered sequence')
    cfg.write_text(json.dumps(dict(status='FROZEN_PRELOCKED_DEVELOPMENT',plan=plan,selection='existing independent validation only',
        reason='capacity-match full66 baselines; true TRAIN Maxwell versus random triangular frame; genuine full512 measurement-to-material mean/moment heads',
        field_loss='proper squared material loss; second moments supervised separately, not FULL posterior',
        target_unchanged=True),indent=2)+'\n')
    while not (OUT/'results/DENSITY_FIRST_COMPLETE.json').exists():time.sleep(10)
    p=subprocess.run([sys.executable,'tools/build_r4_train_frames.py'],cwd=ROOT)
    if p.returncode:raise SystemExit(p.returncode)
    for row in plan:
        job=f"extended_{row['kind']}_s{row['seed']}"
        command=[sys.executable,'-m','p3r4.train_density' if row['group']=='density' else 'p3r4.train_field_moments',
            '--train',str(OUT/'inputs/joint/train32768.npz'),'--validation',str(OUT/'inputs/joint/validation.npz'),
            '--out',str(OUT/'models'/('density' if row['group']=='density' else 'field')/job),
            '--seed',str(row['seed']),'--updates',str(row['updates']),'--train-size',str(row['n']),
            '--batch-size','256','--validation-batch-size','256','--eval-every','300',
            '--max-seconds','1950','--device','mps']
        kind=row['kind']
        if kind=='cfm_matched':command+=['--kind','cfm','--velocity-hidden','547']
        elif kind.startswith('triangular_'):
            frame='maxwell_frame.npz' if kind.endswith('maxwell') else 'random_frame.npz'
            command+=['--kind','triangular','--frame',str(OUT/'inputs/train_frames_v1'/frame)]
        else:command+=['--kind',kind]
        p=subprocess.run([sys.executable,'tools/run_r4_job.py','--stage','I','--job',job,'--training','--gpu',
            '--max-wall','2100','--rss-limit-mib','6144','--',*command],cwd=ROOT)
        print(json.dumps(dict(job=job,returncode=p.returncode)),flush=True)
        if p.returncode:raise SystemExit(p.returncode)
    (OUT/'results/EXTENDED_NEURAL_COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',fits=len(plan)),indent=2)+'\n')
if __name__=='__main__':main()
