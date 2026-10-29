"""Explicit fresh deployment budget; never reopens the expired R3 campaign.

Usage: python tools/deploy_r3.py --deployment-budget-seconds 1800
  --measurement observation.npz --config config.json --output NEW_DIRECTORY
  [--method FULL_GN]

This separate entry point uses frozen R3 scientific source, keeps the global
R3 heavy-FEM/ledger locks for concurrency, and writes all deployment charges
and its fresh RUN_CLOCK under the requested output. No campaign reservation,
clock, authorization or scientific method is extended or modified.
"""
from pathlib import Path
import argparse,datetime,importlib.util,json,os,sys,time
ROOT=Path(__file__).resolve().parents[1]
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--deployment-budget-seconds',required=True,type=int);ap.add_argument('--measurement',required=True);ap.add_argument('--config',required=True);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--method',choices=['FULL_GN','FIXED_CURRENT_CHART','ADAPTIVE_CURRENT_CHART'],default='FULL_GN');a=ap.parse_args()
    if not 1<=a.deployment_budget_seconds<=10800:raise ValueError('Explicit deployment budget must be1..10800seconds')
    out=a.output.resolve()
    if out.exists() and any(out.iterdir()):raise FileExistsError('Fresh deployment requires a new or empty output directory')
    out.mkdir(parents=True,exist_ok=True);accounting=out/'deployment_accounting';accounting.mkdir(exist_ok=False)
    sys.path.insert(0,str(ROOT/'tools'))
    if (ROOT/'vendor').is_dir():sys.path.insert(0,str(ROOT/'vendor'))
    os.environ['R3_SOURCE_ROOT']=str(ROOT/'frozen/src');os.environ['PYTHONDONTWRITEBYTECODE']='1';sys.dont_write_bytecode=True
    entry=ROOT/'tools/run_r3.py';spec=importlib.util.spec_from_file_location('r3_deployment_entry',entry);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    BaseBook=module.BudgetBook
    class DeploymentBook(BaseBook):
        def __init__(self,root=None):
            # root is global only for the shared physical and ledger mutexes.
            # The path and deadline belong exclusively to this deployment.
            self.root=Path(os.environ.get('R3_SHARED_LOCK_ROOT',str(ROOT))).resolve();self.path=accounting/'PHYSICS_LEDGER.jsonl';self.cpu0=time.process_time();now=time.time();end=now+a.deployment_budget_seconds
            utc=lambda stamp:datetime.datetime.fromtimestamp(stamp,datetime.timezone.utc).isoformat()
            self.clock=dict(mode='EXPLICIT_USER_DEPLOYMENT_SEPARATE_FROM_CLOSED_RESEARCH_CAMPAIGN',started_utc=utc(now),deadline_utc=utc(end),no_new_heavy_after_utc=utc(end),wall_seconds=a.deployment_budget_seconds,single_job_rss_bytes=4294967296,heavy_FEM_concurrency=1,global_heavy_lock=str(self.root/'.fem.lock'),campaign_budget_extended=False)
            self.deadline=end;self.heavy_after=end;(accounting/'RUN_CLOCK.json').write_text(json.dumps(self.clock,indent=2)+'\n')
            self.event('DEPLOYMENT','explicit_new_deployment_budget',budget_seconds=a.deployment_budget_seconds,campaign_budget_extended=False,accounting_root=str(accounting))
    module.BudgetBook=DeploymentBook
    sys.argv=[str(entry),'run','--measurement',str(Path(a.measurement).resolve()),'--config',str(Path(a.config).resolve()),'--output',str(out),'--method',a.method]
    try:module.main()
    finally:
        (out/'DEPLOYMENT_RECEIPT.json').write_text(json.dumps(dict(mode='explicit_deployment',budget_seconds=a.deployment_budget_seconds,campaign_budget_extended=False,accounting_root=str(accounting),source_root=str(ROOT/'frozen/src')),indent=2)+'\n')
if __name__=='__main__':main()
