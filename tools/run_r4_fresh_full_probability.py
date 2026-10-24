"""Parent's fixed Full-query release, after BOTH online cohorts have closed.

No heldout truth scoring, re-sampling, candidate ranking or checkpoint changes.
"""
from pathlib import Path
import json,subprocess,sys,os
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,Book,write
from p3r4.evaluation import checkpoint_identity

def main():
    guards=[OUT/'results/FRESH_DENSITY_ONLINE_BOTH_COMPLETE.json',OUT/'results/fresh/risk_online_queue_v1/COMPLETE.json',
            OUT/'results/fresh/risk_physical_queue_v1/PHYSICAL_COMPLETE.json']
    expected=['COMPLETE','BOTH_COHORT_RISK_ONLINE_COMPLETE','ALL_THREE_PHYSICAL_PROCESSES_CLOSED']
    for p,status in zip(guards,expected):
        if json.loads(p.read_text())['status']!=status:raise ValueError('matching closed producer status required')
    for name in ['fresh_raw_locked_v1','fresh_raw_confirmation_v1','fresh_locked64_representatives_v1']:
        if json.loads((OUT/'IV/jobs'/name/'CLOSED.json').read_text())['status']!='COMPLETE':raise ValueError('physical process failed')
    if json.loads(guards[0].read_text())['truth_read'] or json.loads(guards[1].read_text())['heldout_GT_opened'] or json.loads(guards[2].read_text())['GT_read']:
        raise ValueError('first-heldout-truth gate violated')
    root=OUT/'results/fresh/full_probability'
    if root.exists():raise FileExistsError('no paid replay')
    root.mkdir()
    inert=OUT/'configs/INERT_FRESH_DENSITY_FULL_PPC_REQUEST.json'
    spec=json.loads(inert.read_text());spec.update(status='PARENT_RELEASED_FIXED_FULL_QUERIES',released=True,no_execution_authorized=False,
        phase='IV',online_guards=[checkpoint_identity(p) for p in guards],new_truth_read=False,
        source_inert_request=checkpoint_identity(inert))
    write(root/'SPEC.json',spec)
    command=[sys.executable,'tools/run_r4_probability_audit.py','freeze','--spec',str(root/'SPEC.json'),'--out',str(root/'FROZEN.json')]
    b=Book('IV','fresh_full_probability_roster_freeze',analysis_only=True,max_wall=180)
    try:
        p=subprocess.run(command,cwd=ROOT)
        if p.returncode:raise RuntimeError('roster freeze rejected beforephysics')
        if json.loads((root/'FROZEN.json').read_text())['planned_max_RHS']!=56832:raise ValueError('frozen56832RHS roster changed')
    except BaseException as exc:b.close('FAILED',repr(exc));raise
    else:b.close()
    environment=dict(os.environ,R4_ARRAY_CODEC='stored')
    p=subprocess.run([sys.executable,'-u','tools/run_r4_probability_audit.py','run','--freeze',str(root/'FROZEN.json'),
        '--out',str(root/'physics'),'--stage','IV','--job','fresh_full_probability_fixed_draws_v1',
        '--max-rhs','56832','--max-wall','2400'],cwd=ROOT,env=environment)
    if p.returncode:raise SystemExit(p.returncode)
    write(root/'FIXED_FULL_QUERIES_COMPLETE.json',dict(status='COMPLETE',new_truth_read=False,planned_RHS=56832))

if __name__=='__main__':main()
