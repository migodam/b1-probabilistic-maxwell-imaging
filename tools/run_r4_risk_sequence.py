"""Predeclared first matched512/128 risk competition; all jobs billed by runner."""
from pathlib import Path
import subprocess,sys,json
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'research/p3l_r4_conditional_transport'
OLD=ROOT.parent/'b1-p3l-r2-risk-posterior/research/p3l_r2_risk_posterior_exploration'
def main():
 plan=[]
 for kind in ['r2_cnn','affine','sharpen','spatial','local_cnn','phase8','phase16','no_geometry_cnn']:
  for seed in ([2026102501,2026102502,2026102503] if kind in ['local_cnn','phase8','phase16','no_geometry_cnn'] else [2026102501]):
   plan.append(dict(kind=kind,seed=seed,updates=3000 if kind in ['local_cnn','phase8','phase16','no_geometry_cnn'] else 0))
 manifest=OUT/'configs/RISK_FIRST_COMPETITION.json'
 if manifest.exists():raise FileExistsError('sequence already registered; no implicit replay')
 manifest.write_text(json.dumps(dict(phase='MATCHED_512_128',plan=plan,GT_scope='train/validation labels only',physics_RHS=0),indent=2)+'\n')
 for p in plan:
  job=f"risk_{p['kind']}_s{p['seed']}"
  cmd=[sys.executable,'tools/run_r4_job.py','--stage','I','--job',job,'--training','--max-wall','1200','--',sys.executable,'-m','p3r4.train_risk',
    '--train',str(OUT/'inputs/joint/train32768.npz'),'--validation',str(OUT/'inputs/joint/validation.npz'),
    '--train-context',str(OLD/'A/anchor_data/train.npz'),'--validation-context',str(OLD/'A/anchor_data/validation.npz'),
    '--out',str(OUT/'models/risk'/job),'--kind',p['kind'],'--seed',str(p['seed']),'--updates',str(p['updates']),'--device','cpu']
  done=subprocess.run(cmd,cwd=ROOT)
  print(json.dumps(dict(job=job,returncode=done.returncode)),flush=True)
  if done.returncode:raise SystemExit(done.returncode)
if __name__=='__main__':main()
