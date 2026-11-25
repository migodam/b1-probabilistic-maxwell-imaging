"""Rebuild all six supplied mechanism scripts in an isolated working copy."""
from pathlib import Path
import shutil,subprocess,sys,time,os,json
p=Path(__file__).resolve().parents[1]; src=p/'inputs/Z_ALG_Research_Package/checks'; dest=p/'results/mechanisms'
dest.mkdir(exist_ok=True)
for path in src.glob('*.py'):shutil.copy2(path,dest/path.name)
# active_probe generates the 16-cell operators; no precomputed results serve as inputs.
receipts=[]
for name in ('check_fixed_space','active_probe','validate_active_tm','validate_clean_response','background_compensation_exact','active_ambiguity_closed_loop'):
 start=time.perf_counter()
 with (dest/(name+'.log')).open('w') as log:
  run=subprocess.run([sys.executable,'-B',str(dest/(name+'.py'))],stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,VECLIB_MAXIMUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
 receipts.append(dict(script=name,returncode=run.returncode,seconds=time.perf_counter()-start))
 (dest/'runs.json').write_text(json.dumps(receipts,indent=2))
 if run.returncode:raise SystemExit(name)
