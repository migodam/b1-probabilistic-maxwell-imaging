"""Finite initial five-arm task: wait for full release, then train and public VAL."""
from pathlib import Path
from datetime import datetime, timezone
import json, os, sys, time
study=Path(__file__).resolve().parents[2]
marker=study/'audit/DATA_READY.json'
protocol=json.loads((study/'audit/PROTOCOL.json').read_text())
deadline=datetime.fromisoformat(protocol['deadline_utc'].replace('Z','+00:00'))
receipt=study/'audit/training_initial5_process.json'
state={'status':'WAITING_FULL_STRICT_TRAIN_VAL_RELEASE','pid':os.getpid(),'marker':str(marker),
       'python':sys.executable,'arms':['v_physics','j_full','j_null','g_weak','g_dark'],
       'seeds':[7],'device':'mps','test_truth_access':'NONE','soft_alpha_stage':'WAIT_PARENT',
       'loss_ablation_stage':'WAIT_PARENT','full_vie_finetune_stage':'WAIT_PARENT'}
receipt.write_text(json.dumps(state,indent=2)+'\n')
print(json.dumps(state),flush=True)
sys.path.insert(0,str(study/"src"))
from train import strict_training_release

while not strict_training_release(study):
    progress_path=study/'audit/data_progress.json'
    if progress_path.exists():
        try:
            progress=json.loads(progress_path.read_text())
            failed=(progress.get('wrong_primary') or any(
                any(word in name.upper() for word in ('FAILED','ERROR')) and count
                for name,count in progress.get('primary_statuses',{}).items())
                or any(word in str(progress.get('status','')).upper() for word in ('FAILED','ERROR')))
            if failed:
                state.update(status='NOT_RUN_ANCHOR_FAILURE_REPORTED',data_progress=progress,
                             missing_by_split={k:v.get('expected',0)-v.get('ready',0)
                                               for k,v in progress.get('splits',{}).items()})
                receipt.write_text(json.dumps(state,indent=2)+'\n')
                print(json.dumps(state),flush=True);sys.exit(3)
        except (OSError,ValueError): pass
    if datetime.now(timezone.utc)>=deadline:
        state['status']='NOT_RUN_RELEASE_NOT_READY_BY_DEADLINE'
        receipt.write_text(json.dumps(state,indent=2)+'\n')
        print(json.dumps(state),flush=True);sys.exit(2)
    time.sleep(30)
state['status']='RELEASE_OBSERVED_LAUNCHING_INITIAL5'
receipt.write_text(json.dumps(state,indent=2)+'\n')
print(json.dumps(state),flush=True)
args=[sys.executable,'-u',str(study/'experiments/train_models.py'),'--data-ready','--device','mps',
      '--arms','v_physics','j_full','j_null','g_weak','g_dark','--seeds','7','--predict-val']
os.execv(sys.executable,args)
