"""Continue only the CLI-rejected third seed and remaining JVP arms.

The pre-frozen parent plan already requested three seeds; an older helper enum
only listed two. No data, model, objective, optimizer or paid result is changed.
"""
from pathlib import Path
import json,sys,subprocess
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,Book,write
def main():
    release=OUT/'configs/DINO_SEED_REGISTRY_RECOVERY.json'
    if release.exists():raise FileExistsError('no repeated recovery')
    prior=OUT/'II/jobs/forward_value_66_s2026102502_evaluation/CLOSED.json';log=OUT/'logs/fast_and_dino.log'
    if json.loads(prior.read_text())['status']!='COMPLETE' or 'invalid choice' not in log.read_text():raise ValueError('expected exact pre-Book parser failure')
    envelope=(log.stat().st_mtime_ns-prior.stat().st_mtime_ns)/1e9
    if not 0<=envelope<10:raise ValueError('cannot safely bound parser overhead; reconcile first')
    b=Book('II','dino_seed_parser_overhead_reconciliation',analysis_only=True,max_wall=60)
    b._update('FAILED_COMMAND_OVERHEAD_UPPER_BOUND',controlled=envelope,physical_actions=0,optimizer_updates=0,
        scope='filetime envelope after last successful child CLOSED until parser stderr, conservative wall upper bound not exact CPU time',
        source_log=str(log),no_previous_receipts_or_caps_reset=True)
    b.close()
    plan=[dict(mode='value',seed=2026102503)]+[dict(mode='jvp',seed=s) for s in [2026102501,2026102502,2026102503]]
    write(release,dict(reason='helper enum omitted the third seed already present in pre-frozen6-fit parent plan',
        unchanged_data_model_objective_hyperparameters=True,existing_fastpair_and_two_valuefits_retained=True,
        parser_before_Book=True,actual_new_physics=0,parser_wall_upper_bound=envelope,remaining_plan=plan))
    for row in plan:
        job=f"forward_{row['mode']}_66_s{row['seed']}"
        commands=[['tools/train_forward_jvp.py','--train','--mode',row['mode'],'--coordinates','66','--seed',str(row['seed']),
            '--updates','1000','--job',job,'--max-wall','900'],
            ['tools/train_forward_jvp.py','--evaluate','--checkpoint',str(OUT/'results/II/forward_jvp_training'/job/'BEST.pt'),
             '--job',job+'_evaluation','--max-wall','300']]
        for cmd in commands:
            p=subprocess.run([sys.executable,*cmd],cwd=ROOT)
            if p.returncode:raise SystemExit(p.returncode)
    write(OUT/'results/FAST_AND_DINO_COMPLETE.json',dict(status='COMPLETE',fits=6,recovery=str(release),fast_pair_NOT_replayed=True))
if __name__=='__main__':main()
