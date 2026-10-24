"""Explicit recovery of the MPS watermark launch error; no budget reset."""
from pathlib import Path
import subprocess,sys,json
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r4_conditional_transport'
def main():
 source=OUT/'configs/DENSITY_FIRST_COMPETITION.json';plan=json.loads(source.read_text())['plan']
 failed=OUT/'I/jobs/density_LR_5e-05/CLOSED.json'
 if json.loads(failed.read_text())['status']!='FAILED':raise ValueError('expected exact failed attempt')
 dest=OUT/'configs/DENSITY_RUNTIME_RECOVERY_01.json'
 if dest.exists():raise FileExistsError('recovery already registered')
 dest.write_text(json.dumps(dict(reason='MPS defaultlow1.4 exceeded explicit high0.5; set explicitlow0.4',
    old_receipt=str(failed),old_budget_counts_retained=True,old_logs_retained=True,optimizer_updates_before_failure=0,
    changed_scientific_target=False,changed_model=False,changed_data=False,changed_hyperparameters=False),indent=2)+'\n')
 for i,p in enumerate(plan):
  job='density_'+p['label']+('_runtimefix1' if i==0 else '')
  if (OUT/'I/jobs'/job/'CLOSED.json').exists():raise FileExistsError('no silent replay/skip of scientific fits')
  cmd=[sys.executable,'tools/run_r4_job.py','--stage','I','--job',job,'--training','--gpu','--max-wall','2100','--',sys.executable,'-m','p3r4.train_density',
    '--train',str(OUT/'inputs/joint/train32768.npz'),'--validation',str(OUT/'inputs/joint/validation.npz'),
    '--out',str(OUT/'models/density'/job),'--kind',p['kind'],'--seed',str(p['seed']),'--train-size',str(p['n']),
    '--updates',str(p['updates']),'--lr',str(p['lr']),'--width',str(p['width']),'--batch-size','256',
    '--eval-every','300','--validation-batch-size','256','--max-seconds','1950','--device','mps']
  done=subprocess.run(cmd,cwd=ROOT);print(json.dumps(dict(job=job,returncode=done.returncode)),flush=True)
  if done.returncode:raise SystemExit(done.returncode)
 (OUT/'results/DENSITY_FIRST_COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',jobs=len(plan),runtime_recovery=str(dest)),indent=2)+'\n')
if __name__=='__main__':main()
