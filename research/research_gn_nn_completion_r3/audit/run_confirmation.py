"""Finite, predeclared final training and public TEST prediction queue.

No TEST truth is opened. Physics workers seal all predictions before scoring.
"""
from pathlib import Path
import json, os, subprocess, sys, time

S=Path(__file__).resolve().parents[1]
R=S.parents[1]
sys.path.insert(0,str(S/'src'))
ENV=os.environ.copy()
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    ENV[k]='1';os.environ[k]='1'
ARMS=['v_image','v_physics','j_full','j_null','j_null_current']
PLAN=[dict(arm=a,seed=seed,passes=1,run=f'confirm1024_{a}_seed{seed}') for seed in (7,17,27) for a in ARMS]
PLAN += [dict(arm=a,seed=7,passes=2,run=f'confirm1024_{a}_seed7_cascade2') for a in ('v_physics','j_null')]
ONLINE=json.loads((S/'data/manifests/online.json').read_text())['records']
TRAIN=[r['obs_id']for r in ONLINE if r['split']=='train']
TEST=[r['obs_id']for r in ONLINE if r['split']=='test']
FROZEN=dict(frozen_unix=time.time(),train_parents=1024,train_observations=len(TRAIN),test_observations=len(TEST),plans=PLAN,controls=['continue_trf','weakreg45_anchor'],test_truth_access='NONE',no_test_based_model_choice=True,main_result='All five seed7 arms; seed17/27 confirmation; cascade2 separately labeled',optional_softfilter='NOT_RUN until core result completion',wall_deadline_utc='2026-11-07 03:26:37 UTC')
if len(sys.argv)==1:
    (S/'audit/CONFIRMATION_PLAN_FROZEN.json').write_text(json.dumps(FROZEN,indent=2))

def worker(index):
    cli=S/'experiments/run_evaluation.py'
    commands=[['--controls','continue_trf','weakreg45_anchor']]+[['--runs',p['run']]for p in PLAN]
    for command in commands:
        if command[0]=='--runs':
            run=command[1]
            while not (S/'predictions'/run/'PREDICTION_COMPLETE.json').exists():
                if (S/'audit/CONFIRMATION_TRAINING_DONE.json').exists():
                    break
                time.sleep(15)
        with (S/f'evaluation/CONFIRMATION_P{index}.log').open('a')as f:
            p=subprocess.run([sys.executable,str(cli),'evaluate',*command,'--splits','test','--partition-count','2','--partition-index',str(index),'--partition-block-size','8'],cwd=R,env=ENV,stdout=f,stderr=subprocess.STDOUT)
        if p.returncode:
            (S/f'audit/CONFIRMATION_P{index}_ERROR.json').write_text(json.dumps(dict(command=command,returncode=p.returncode)))
            return
    (S/f'audit/CONFIRMATION_P{index}_DONE.json').write_text(json.dumps(dict(unix=time.time())))

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='worker':
        worker(int(sys.argv[2]));sys.exit(0)
    while not all((S/'anchors'/f'{i}.npz').exists()and(S/'anchors'/f'{i}.json').exists()for i in TRAIN):
        time.sleep(30)
    from train import train_model
    from predict import predict_anchors
    physics=[subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'worker',str(i)],cwd=S,env=ENV)for i in (0,1)]
    rows=[]
    start=time.time()
    for plan in PLAN:
        tick=time.time()
        suffix='_cascade2'if plan['passes']==2 else''
        output=S/'models/confirm1024'/f"{plan['arm']}_seed{plan['seed']}{suffix}"
        train_model(S/'data',S/'data/manifests/evaluation.json',S/'data/operators.npz',output,plan['arm'],plan['seed'],100,64,15,'mps',plan['passes'],S/'anchors',0,False)
        trained=time.time()
        dest=S/'predictions'/plan['run']
        predict_anchors(output/'best.pt',[S/'anchors'/f'{i}.npz'for i in TEST],S/'data/operators.npz',dest,'mps')
        row=dict(**plan,training_wall_seconds=trained-tick,prediction_wall_seconds=time.time()-trained,checkpoint=str(output/'best.pt'),test_truth_access='NONE')
        (dest/'PREDICTION_COMPLETE.json').write_text(json.dumps(row,indent=2))
        rows.append(row)
        (S/'audit/CONFIRMATION_TRAINING_PROGRESS.json').write_text(json.dumps(rows,indent=2))
    (S/'audit/CONFIRMATION_TRAINING_DONE.json').write_text(json.dumps(dict(rows=rows,wall_seconds=time.time()-start)))
    codes=[p.wait()for p in physics]
    (S/'audit/CONFIRMATION_PHYSICS_DONE.json').write_text(json.dumps(dict(returncodes=codes,unix=time.time())))
