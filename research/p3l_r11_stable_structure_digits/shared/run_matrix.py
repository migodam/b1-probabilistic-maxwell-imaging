from core import *
import subprocess,sys
jobs=[('N','AGENT_B'),('N','AGENT_C'),('P','FULL256'),('P','AGENT_A'),('P','AGENT_B'),('P','AGENT_C')]
for ds,m in jobs:
 sums=0
 for path in (P/'results').glob('*/*/*/PREDICTION_SEALED.json'):sums+=json.loads(path.read_text())['physics_RHS']
 # Each next9job worstcase declared240forward*40=96k plus80adjoint? cap reservation conservative.
 if sums>285000:
  write(P/'results'/f'{ds}_{m}_NOT_RUN.json',dict(status='NOT_RUN_BUDGET',previous_main_RHS=sums));continue
 print('START',ds,m,'previousRHS',sums,flush=True)
 log=(P/'results'/f'{ds}_{m}.log').open('w')
 pr=subprocess.run([sys.executable,str(P/'shared/invert_digits.py'),'--dataset',ds,'--method',m],stdout=log,stderr=subprocess.STDOUT,cwd=REPO)
 log.close()
 print('END',ds,m,'code',pr.returncode,flush=True)
 if pr.returncode:break
