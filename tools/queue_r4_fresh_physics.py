"""One serial fresh-data/anchor/quiet-kernel campaign; no truth scoring.

Both cohort decisions must already be frozen. All child producers own their
append-only Books. Failed or existing stages are not silently restarted.
"""
from pathlib import Path
import json,subprocess,sys,os,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,write
from p3r4.evaluation import checkpoint_identity

def main():
    if json.loads((OUT/'results/FRESH_METHODS_FROZEN_BEFORE_OBSERVATIONS.json').read_text())['status']!='COMPLETE':raise ValueError('freeze before data')
    ledger=json.loads((OUT/'RESOURCE_LEDGER.json').read_text())
    if ledger['active_meta']:raise ValueError('start fresh physical queue only after other jobs close')
    root=OUT/'results/fresh_physics_queue_v1'
    if root.exists():raise FileExistsError('no replay; explicit partial recovery must retain all prior debits')
    root.mkdir();write(root/'START.json',dict(pid=os.getpid(),time_ns=time.time_ns(),new_GT_scoring=False,
        freezes={g:checkpoint_identity(OUT/f'configs/FROZEN_{g.upper()}.json') for g in ('locked','confirmation')}))
    environment=dict(os.environ,R4_ARRAY_CODEC='stored',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    calls=[];status,error='COMPLETE',None
    def call(name,command):
        path=root/name;path.mkdir();tick=time.perf_counter()
        write(path/'COMMAND.json',dict(argv=[sys.executable,*command],cwd=str(ROOT),truth_scoring=False,array_codec='stored'))
        with (path/'stdout.log').open('x') as output,(path/'stderr.log').open('x') as errors:
            p=subprocess.run([sys.executable,'-u',*command],cwd=ROOT,env=environment,stdout=output,stderr=errors)
        calls.append(dict(name=name,returncode=p.returncode,external_wall=time.perf_counter()-tick))
        write(root/'PROGRESS.json',dict(status='RUNNING',steps=calls))
        if p.returncode:raise RuntimeError(name+' failed; keep partial arrays/ledger, no automatic retry')
    try:
        call('generate_locked',['tools/generate_r4_data.py','locked','--count','512'])
        call('generate_confirmation',['tools/generate_r4_data.py','confirmation','--count','256'])
        call('contexts_locked_0',['tools/run_r4_contexts.py','locked','--start','0','--count','64'])
        release=OUT/'configs/INTERVAL_CHAIN_LOCKED_RELEASE.json'
        if release.exists():raise FileExistsError(release)
        r=json.loads((OUT/'configs/INTERVAL_CHAIN_LOCKED_RELEASE_TEMPLATE.json').read_text())
        r.update(released=True,context_receipt=str(OUT/'IV/jobs/contexts_locked_0_64/CLOSED.json'),
            scope='pre-frozen independent eight-parent fixed-kernel validation, no other numerical jobs scheduled; not a posterior reference')
        write(release,r)
        call('quiet_locked_preflight',['tools/run_r4_interval_chain_locked.py','preflight','--release',str(release)])
        call('quiet_locked_four_arm',['tools/run_r4_interval_chain_locked.py','pair','--release',str(release)])
        for group,n in [('locked',512),('confirmation',256)]:
            for start in range(64 if group=='locked' else 0,n,64):
                call(f'contexts_{group}_{start}',['tools/run_r4_contexts.py',group,'--start',str(start),'--count','64'])
            call('export_'+group,['tools/export_r4_anchor_conditions.py',group,'--count',str(n)])
    except BaseException as exc:status,error='FAILED',repr(exc);raise
    finally:write(root/'RESULT.json',dict(status=status,error=error,steps=calls,GT_scoring=False,all_child_costs_ledgered=True))
    write(OUT/'results/FRESH_PHYSICS_AND_CONTEXTS_COMPLETE.json',dict(status='COMPLETE',new_parents=768,
        independent_cohorts=[512,256],GT_scoring=False,source_queue=str(root)))

if __name__=='__main__':main()
