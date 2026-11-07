"""Finite local fetch monitor; no GT loading and no physics work."""
from pathlib import Path
import json,time,subprocess,sys,collections
STUDY=Path(__file__).resolve().parents[2]
PYTHON=sys.executable
online=json.loads((STUDY/'data/manifests/online.json').read_text())['records']
ids={s:{r['obs_id'] for r in online if r['split']==s} for s in ('train','val','test')}

def dump(path,value):
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)

def fetch(split):
    with (STUDY/f'delegated/data/fetch_{split}_progress.log').open('a') as log:
        result=subprocess.run([PYTHON,str(STUDY/'src/remote_jobs.py'),'fetch','--split',split],stdout=log,stderr=subprocess.STDOUT)
    return result.returncode

def summarize():
    all_rows={};costs=collections.Counter();stops=collections.Counter();solve=0.;wrong=[]
    for split,expected in ids.items():
        records=[]
        for obsid in expected:
            path=STUDY/'anchors'/f'{obsid}.json'
            if not path.exists() or not path.with_suffix('.npz').exists():continue
            try:r=json.loads(path.read_text())
            except (ValueError,OSError):continue
            if r.get('config',{}).get('discrepancy_stop',False):wrong.append(obsid);continue
            records.append((obsid,r));costs.update(r['counts']);stops[r['status']]+=1;solve+=r['solve_wall_seconds']
        all_rows[split]=dict(ready=len(records),expected=len(expected),full24=sum(r['final_stage']==24 and r['all_unknowns_released'] for _,r in records),nonfull=[o for o,r in records if r['final_stage']!=24 or not r['all_unknowns_released']])
    return dict(status='PUBLIC_ANCHOR_PREPARATION',time_unix=time.time(),splits=all_rows,primary_counts=dict(costs),primary_statuses=dict(stops),primary_sum_solve_seconds=solve,wrong_primary=wrong,test_truth_status='SEALED_NOT_SCORED',remote_affinity_mask=20560,remote_cores=[4,6,12,14])

if __name__=='__main__':
    began=time.time();loops=0;baselines=None;published=False
    while time.time()-began<6*3600:
        for split in ('train','val','test'):
            if not (STUDY/'audit/data_progress.json').exists() or split!='test' or published:
                result=fetch(split)
                if result:print(json.dumps(dict(fetch_error=split,code=result)),flush=True)
        status=summarize();dump(STUDY/'audit/data_progress.json',status)
        print(json.dumps(status),flush=True)
        trainval=all(status['splits'][s]['ready']==len(ids[s]) and status['splits'][s]['full24']==len(ids[s]) for s in ('train','val')) and not status['wrong_primary']
        if trainval and not published:
            dump(STUDY/'audit/DATA_READY.json',dict(status='TRAIN_VAL_STRICT_ANCHORS_READY',train=1008,val=156,anchors=str(STUDY/'anchors'),online=str(STUDY/'data/manifests/online.json'),evaluator=str(STUDY/'data/manifests/evaluation.json'),test_truth_status='SEALED_NOT_SCORED',all_primary_discrepancy_stop=False,all_primary_full24=True))
            published=True;print('DATA_READY',flush=True)
        all_base=all(status['splits'][s]['ready']==len(ids[s]) for s in ids)
        if all_base and baselines is None:
            # The existing base queue has finished its actual tasks before this
            # independent four-worker baseline queue starts.
            log=(STUDY/'delegated/data/remote_baselines_session.log').open('w')
            baselines=subprocess.Popen([PYTHON,str(STUDY/'src/remote_jobs.py'),'run','--kind','baselines'],stdout=log,stderr=subprocess.STDOUT)
            print(json.dumps(dict(baseline_local_process=baselines.pid)),flush=True)
        if baselines is not None:
            fetch('val');fetch('test')
            complete=sum((STUDY/'anchors'/f'{obsid}{suffix}.json').exists() and (STUDY/'anchors'/f'{obsid}{suffix}.npz').exists() for split in ('val','test') for obsid in ids[split] for suffix in ('_discrepancy','_weakreg1e-5'))
            dump(STUDY/'audit/data_baseline_progress.json',dict(complete=complete,expected=744,process_returncode=baselines.poll(),test_truth_status='SEALED_NOT_SCORED'))
            if complete==744:
                print('BASELINES_COMPLETE',flush=True);break
            if baselines.poll() is not None:
                print('BASELINE_PROCESS_STOPPED_INCOMPLETE',flush=True);break
        loops+=1
        for _ in range(2):time.sleep(60)
