"""Finite post-scoring sequence: independent SLICE refinement and evidence audit.

No algorithm is reselected after fresh truth access. These tasks read fixed
outputs or refine an historically exposed conditional posterior.
"""
from pathlib import Path
import json, subprocess, sys, time

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import OUT,write


def execute(label,args):
    path=OUT/('post_'+label+'.log')
    write(OUT/'POST_SCORING_CHECKPOINT.json',dict(status='RUNNING',task=label))
    with path.open('x') as log:
        result=subprocess.run([sys.executable,*args],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
    if result.returncode:
        write(OUT/'POST_SCORING_CHECKPOINT.json',dict(status='STOPPED_ERROR',task=label,
            exit_code=result.returncode,log=str(path),automatic_replay=False))
        raise RuntimeError('Post-scoring task stopped: '+label)
    print(label+' COMPLETE',flush=True)


def analysis(label,script,extra=()):
    execute(label,['tools/run_p3r2_accounted.py','--stage','D','--job',label,
        '--max-wall','3600','--',sys.executable,script,*extra])


def main():
    while True:
        path=OUT/'FRESH_CHECKPOINT.json'
        state=json.loads(path.read_text()) if path.exists() else {}
        if state.get('status')=='STOPPED_ERROR':raise RuntimeError('Fresh sequence needs diagnosis; no paid automatic replay')
        if state.get('status')=='COMPLETE':break
        time.sleep(15)
    if (OUT/'POST_SCORING_CHECKPOINT.json').exists():raise FileExistsError('Review existing post-scoring execution')
    execute('field_cdf',['tools/run_p3r2_field_cdf_reference.py','--job','B_field2_direct_CDF_DEV030_r0'])
    analysis('independent_locked_statistics','tools/audit_p3r2_locked_statistics.py')
    analysis('final_figures','tools/plot_p3r2_results.py',['--output-directory',str(OUT/'figures/final')])
    write(OUT/'POST_SCORING_CHECKPOINT.json',dict(status='COMPLETE',no_locked_reselection=True))


if __name__=='__main__':main()
