"""First heldout truth gate and isolated scoring, only on explicit main release.

No training, inference, physical query, re-ranking or new method is permitted.
All cold/warm costs stay in their original owning receipts.
"""
from pathlib import Path
import argparse,json,subprocess,sys,time,os
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,write
from p3r4.evaluation import checkpoint_identity

def main():
    p=argparse.ArgumentParser();p.add_argument('--release',type=Path,required=True);a=p.parse_args()
    r=json.loads(a.release.read_text())
    if r.get('released') is not True or r.get('experiment')!='fresh_r4_isolated_scores_v1':raise ValueError('main explicit firstGT release required')
    for identity in r['closed_gates']:
        if checkpoint_identity(identity['path'])!=identity:raise ValueError('release gate changed')
    if not r.get('all_online_both_cohorts_closed_before_truth') or not r.get('pretruth_optimizer_diagnostics_closed'):raise ValueError('chronology gate missing')
    root=OUT/'results/fresh/scoring_queue_v1'
    if root.exists():raise FileExistsError('no score replay; preservepartialfailure')
    root.mkdir()
    write(root/'FIRST_HELDOUT_TRUTH_GATE.json',dict(release=checkpoint_identity(a.release),time_ns=time.time_ns(),pid=os.getpid(),
        no_further_model_or_algorithm_selection=True,all_model_predictions_frozen=True))
    component=OUT/'AUDIT/risk/fresh_freeze_v3/RISK_COMPONENT_FROZEN.json'
    records=[];status,error='COMPLETE',None
    def run(name,command,cap=1200):
        folder=root/name;folder.mkdir();cmd=[sys.executable,'tools/run_r4_job.py','--stage','IV','--job',name,'--max-wall',str(cap),'--rss-limit-mib','3072','--',sys.executable,*command]
        write(folder/'COMMAND.json',dict(argv=cmd,scope='isolatedofflineevaluationonly'))
        with (folder/'stdout.log').open('x') as stdout,(folder/'stderr.log').open('x') as stderr:
            q=subprocess.run(cmd,cwd=ROOT,stdout=stdout,stderr=stderr)
        records.append(dict(job=name,returncode=q.returncode));write(root/'PROGRESS.json',dict(status='RUNNING',records=records))
        if q.returncode:raise RuntimeError(name+' failed; no automaticretry')
    try:
        for group in ['locked','confirmation']:
            folder=OUT/'results/fresh'/group
            run(f'fresh_risk_{group}_score_v1',['tools/r4_fresh_risk.py','score','--cohort',group,'--component',str(component),
                '--top-freeze',str(OUT/f'configs/FROZEN_{group.upper()}.json'),'--freeze',str(folder/'BOUND_RISK.json'),
                '--online',str(folder/'risk45/online'),'--truth',str(OUT/'inputs/offline_only'),'--out',str(folder/'risk45/scores')])
        for group in ['locked','confirmation']:
            folder=OUT/'results/fresh'/group
            run(f'fresh_density_{group}_score_NFE32',['tools/evaluate_r4_density_frozen_nfe.py','score','--freeze',str(folder/'FROZEN_DENSITY.json'),
                '--online',str(folder/'density_online'),'--truth',str(OUT/'inputs/offline_only'),'--out',str(folder/'density_scores'),
                '--device','cpu','--batch-size','16'],cap=7200)
        full=OUT/'results/fresh/full_probability'
        run('fresh_full_probability_offline_score',['tools/run_r4_probability_audit.py','score','--freeze',str(full/'FROZEN.json'),
            '--online',str(full/'physics'),'--truth',str(OUT/'inputs/offline_only'),'--out',str(full/'scores')])
        rep=OUT/'results/fresh/locked'
        run('fresh_locked64_representatives_offline_score',['tools/r4_fresh_representatives.py','score','--stage','IV',
            '--prepared-request',str(rep/'representatives_prepared/FROZEN_REQUEST.json'),
            '--physical',str(rep/'representatives'),'--receipt',str(OUT/'IV/jobs/fresh_locked64_representatives_v1/CLOSED.json'),
            '--truth',str(OUT/'inputs/offline_only'),'--out',str(rep/'representative_scores')])
    except BaseException as exc:status,error='FAILED',repr(exc);raise
    finally:write(root/'RESULT.json',dict(status=status,error=error,records=records,physical_queries=0,training=0,checkpoints_unchanged=True))
    write(OUT/'results/FRESH_ISOLATED_SCORES_COMPLETE.json',dict(status='COMPLETE',cohorts=[512,256],first_truth_gate=str(root/'FIRST_HELDOUT_TRUTH_GATE.json')))

if __name__=='__main__':main()
