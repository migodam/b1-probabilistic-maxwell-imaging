"""Bounded full-dimensional density competition; CPU risk job finishes first."""
from pathlib import Path
import subprocess,sys,json,time,fcntl
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r4_conditional_transport'
def main():
 plan=[]
 # The LR range is a real validation experiment, not an assertion of stability.
 for lr in [5e-5,2e-4,8e-4]:
  plan.append(dict(label=f'LR_{lr}',kind='gaussian',seed=2026102501,n=4096,updates=1000,lr=lr,width=128))
 # Fixed-family first pass and independent training seeds. One training job at
 # a time keeps MPS and CPU memory from accumulating across child lifetimes.
 for kind in ['gaussian','spline','triangular','cfm']:
  for seed in [2026102501,2026102502,2026102503]:
   plan.append(dict(label=f'{kind}_32768_s{seed}',kind=kind,seed=seed,n=32768,updates=6000,lr=2e-4,width=128))
 for n in [4096,16384]:
  for kind in ['spline','cfm']:
   plan.append(dict(label=f'{kind}_{n}_s2026102501',kind=kind,seed=2026102501,n=n,updates=6000,lr=2e-4,width=128))
 manifest=OUT/'configs/DENSITY_FIRST_COMPETITION.json'
 if manifest.exists():raise FileExistsError('no implicit sequence replay')
 manifest.write_text(json.dumps(dict(plan=plan,selection='validation only',device='mps',sampler='full66; no strong coordinates replaced by prior'),indent=2)+'\n')
 expected=[f'risk_{k}_s{s}' for k in ['r2_cnn','affine','sharpen','spatial','local_cnn','phase8','phase16','no_geometry_cnn'] for s in ([2026102501,2026102502,2026102503] if k in ['local_cnn','phase8','phase16','no_geometry_cnn'] else [2026102501])]
 start=time.monotonic()
 while True:
  receipt=[OUT/'I/jobs'/j/'CLOSED.json' for j in expected]
  failed=[str(p) for p in receipt if p.exists() and json.loads(p.read_text())['status']!='COMPLETE']
  if failed:raise RuntimeError('risk execution error needs review: '+str(failed))
  if all(p.exists() for p in receipt):break
  if time.monotonic()-start>3600:raise TimeoutError('risk prequeue wait; no model silently skipped')
  time.sleep(10)
 for p in plan:
  job='density_'+p['label']
  cmd=[sys.executable,'tools/run_r4_job.py','--stage','I','--job',job,'--training','--gpu','--max-wall','2100','--',sys.executable,'-m','p3r4.train_density',
    '--train',str(OUT/'inputs/joint/train32768.npz'),'--validation',str(OUT/'inputs/joint/validation.npz'),
    '--out',str(OUT/'models/density'/job),'--kind',p['kind'],'--seed',str(p['seed']),'--train-size',str(p['n']),
    '--updates',str(p['updates']),'--lr',str(p['lr']),'--width',str(p['width']),'--batch-size','256',
    '--eval-every','300','--validation-batch-size','256','--max-seconds','1950','--device','mps']
  done=subprocess.run(cmd,cwd=ROOT)
  print(json.dumps(dict(job=job,returncode=done.returncode)),flush=True)
  if done.returncode:raise SystemExit(done.returncode)
 (OUT/'results/DENSITY_FIRST_COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',jobs=len(plan)),indent=2)+'\n')
if __name__=='__main__':main()
