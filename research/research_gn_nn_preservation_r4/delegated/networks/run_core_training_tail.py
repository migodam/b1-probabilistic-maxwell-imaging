"""Finite predeclared confirmation and single-loss tasks; one MPS lane."""
from pathlib import Path
import json,os,subprocess,sys,time
study=Path(__file__).resolve().parents[2]
receipt=study/'audit/CORE_TRAINING_TAIL.json'
state={'status':'WAIT_EXISTING_JFULL_CONFIRMATION','pid':os.getpid(),'test_truth_access':'NONE','constrained_arm':'g_weak','tasks':[]}
receipt.write_text(json.dumps(state,indent=2)+'\n')
for _ in range(180):
    ready=all((study/'models'/f'j_full_seed{seed}_cascade2/summary.json').exists()
              and len(list((study/'predictions'/f'j_full_seed{seed}_cascade2').glob('obs_*.json')))==156
              for seed in (17,27))
    if ready:break
    time.sleep(5)
else:raise RuntimeError('Current confirmation did not close; preserve its evidence')
tasks=[(['g_weak'],[17,27],'material'),(['v_physics'],[27],'material'),
       (['j_full','g_weak'],[7],'const'),(['j_full','g_weak'],[7],'internal')]
for arms,seeds,loss in tasks:
    command=[sys.executable,'-u',str(study/'experiments/train_models.py'),'--data-ready',
             '--arms',*arms,'--seeds',*map(str,seeds),'--loss',loss,'--device','mps','--predict-val']
    state['status']='TRAINING';state['active']={'arms':arms,'seeds':seeds,'loss':loss,'started_unix':time.time()}
    receipt.write_text(json.dumps(state,indent=2)+'\n')
    process=subprocess.run(command,cwd=study)
    state['tasks'].append(dict(state['active'],returncode=process.returncode,finished_unix=time.time()))
    receipt.write_text(json.dumps(state,indent=2)+'\n')
    if process.returncode:raise SystemExit(process.returncode)
state.update(status='COMPLETED_CORE_CONFIRMATION_AND_E2',active=None)
receipt.write_text(json.dumps(state,indent=2)+'\n')
