"""Final phase-head revision: exact bins and matched anchor initialization."""
from pathlib import Path
import json,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r4_conditional_transport'
OLD=ROOT.parent/'b1-p3l-r2-risk-posterior/research/p3l_r2_risk_posterior_exploration'
def main():
    path=OUT/'configs/ANCHORBIN_REVISION2_ACTUAL_RUNS.json'
    if path.exists():raise FileExistsError('no phase-family replay')
    plan=[dict(kind=kind,seed=seed,updates=3000) for kind in ['phase8','phase16'] for seed in [2026102501,2026102502,2026102503]]
    path.write_text(json.dumps(dict(plan=plan,revision='same exact bin labels, anchor-based full-mass softmax initialization',
        original_DATA_and_model_hyperparameters_preserved=True,phase_family_final_mechanism_revision=True,
        no_locked_GT=True,CPU_slot='same bounded6k-parameter2GiB risk slot; no extra concurrency slot'),indent=2)+'\n')
    while not (OUT/'results/BIN_REVISION1_COMPLETE.json').exists():time.sleep(10)
    for row in plan:
        job=f"risk_anchorbin_{row['kind']}_s{row['seed']}"
        command=[sys.executable,'-m','p3r4.train_risk_anchorbin',
            '--train',str(OUT/'inputs/joint/train32768.npz'),'--validation',str(OUT/'inputs/joint/validation.npz'),
            '--anchor-train',str(OLD/'A/anchor_data/train.npz'),'--anchor-validation',str(OLD/'A/anchor_data/validation.npz'),
            '--out',str(OUT/'models/risk'/job),'--kind',row['kind'],'--seed',str(row['seed']),
            '--updates',str(row['updates']),'--device','cpu']
        p=subprocess.run([sys.executable,'tools/run_r4_job.py','--stage','I','--job',job,'--training','--cpu-parallel',
            '--max-wall','1200','--rss-limit-mib','2048','--',*command],cwd=ROOT)
        print(json.dumps(dict(job=job,returncode=p.returncode)),flush=True)
        if p.returncode:raise SystemExit(p.returncode)
    (OUT/'results/ANCHORBIN_REVISION2_COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',fits=6),indent=2)+'\n')
if __name__=='__main__':main()
