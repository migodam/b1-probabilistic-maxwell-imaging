"""Exact-bin hypothesis revision, isolated from the original point-rule fits."""
from pathlib import Path
import subprocess,sys,json,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r4_conditional_transport'
OLD=ROOT.parent/'b1-p3l-r2-risk-posterior/research/p3l_r2_risk_posterior_exploration'
def main():
 cfg=OUT/'configs/BIN_REVISION1_ACTUAL_RUNS.json'
 plan=[dict(kind=k,seed=s,updates=3000) for k in ['phase8','phase16'] for s in [2026102501,2026102502,2026102503]]
 if cfg.exists():raise FileExistsError('never replay an existing sequence')
 cfg.write_text(json.dumps(dict(plan=plan,revision='point-threshold quadrature -> exact equal-h-mass-bin conditional moments',target_unchanged=True,
   old_results_retained=True,decision_before_locked=True,initialization_same_as_original_point_phase=True),indent=2)+'\n')
 # Keep the first two independent physical timing arms under unchanged load.
 while not (OUT/'results/II/interval_chain_pair/PAIR_COMPLETE.json').exists():time.sleep(10)
 for p in plan:
  job=f"risk_bin_{p['kind']}_s{p['seed']}"
  cmd=[sys.executable,'tools/run_r4_job.py','--stage','I','--job',job,'--training','--cpu-parallel','--rss-limit-mib','2048','--max-wall','1200',
       '--',sys.executable,'-m','p3r4.train_risk_bin','--train',str(OUT/'inputs/joint/train32768.npz'),
       '--validation',str(OUT/'inputs/joint/validation.npz'),'--anchor-train',str(OLD/'A/anchor_data/train.npz'),
       '--anchor-validation',str(OLD/'A/anchor_data/validation.npz'),'--out',str(OUT/'models/risk'/job),
       '--kind',p['kind'],'--seed',str(p['seed']),'--updates',str(p['updates']),'--device','cpu']
  result=subprocess.run(cmd,cwd=ROOT);print(json.dumps(dict(job=job,returncode=result.returncode)),flush=True)
  if result.returncode:raise SystemExit(result.returncode)
 (OUT/'results/BIN_REVISION1_COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',fits=len(plan)),indent=2)+'\n')
if __name__=='__main__':main()
