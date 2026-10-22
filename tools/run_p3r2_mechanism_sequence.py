"""One finite, separately frozen confirmation; no automatic paid replay."""
from pathlib import Path
import json,subprocess,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import OUT,write
from p3r2_mechanism.contract import validate_freeze


def execute(name,args):
    log=OUT/('confirmation_'+name+'.log')
    write(OUT/'MECHANISM_CHECKPOINT.json',dict(status='RUNNING',task=name))
    with log.open('x') as handle:
        result=subprocess.run([sys.executable,*args],cwd=ROOT,stdout=handle,stderr=subprocess.STDOUT)
    if result.returncode:
        write(OUT/'MECHANISM_CHECKPOINT.json',dict(status='STOPPED_ERROR',task=name,exit_code=result.returncode,
            log=str(log),automatic_replay=False,paid_outputs_preserved=True))
        raise RuntimeError('Confirmation stopped: '+name)
    print(name+' COMPLETE',flush=True)


def analysis(name,script):
    execute(name,['tools/run_p3r2_accounted.py','--stage','L','--job',name,'--max-wall','1800','--',sys.executable,script])


def main():
    validate_freeze(OUT,OUT/'MECHANISM_CONFIRMATION_FROZEN.json')
    if (OUT/'MECHANISM_CHECKPOINT.json').exists():raise FileExistsError('Review any previous execution before continuation')
    execute('data',['tools/generate_p3r2_data.py','locked','--start','256','--count','256'])
    execute('anchors',['tools/run_p3r2_anchors.py','locked','--start','256','--count','256','--native','96','--probit','64'])
    analysis('mechanism_confirmation_inference','tools/infer_p3r2_mechanism_confirmation.py')
    analysis('mechanism_confirmation_evaluation','tools/evaluate_p3r2_mechanism_confirmation.py')
    write(OUT/'MECHANISM_CHECKPOINT.json',dict(status='COMPLETE',count=256,models=21,
        full_online_completion_before_truth=True,first_cohort_untouched=True,no_test_reselection=True))


if __name__=='__main__':main()
