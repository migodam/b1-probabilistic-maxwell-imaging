"""Execute only a pre-frozen method contract; never reads private target files."""
import argparse
import datetime
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--contract',required=True,type=Path)
    a=p.parse_args()
    contract=json.loads(a.contract.read_text())
    if not contract.get('frozen_utc') or not contract.get('method_configs'):
        raise ValueError('Frozen method_configs contract required')
    source_freeze=json.loads((ROOT/contract['source_manifest']).read_text())
    def assert_source_frozen():
        for item in source_freeze['files']:
            if (ROOT/item['source']).read_bytes()!=(ROOT/item['snapshot']).read_bytes():
                raise RuntimeError('Source changed after frozen method contract: '+item['source'])
    assert_source_frozen()
    public=sorted((ROOT/'public_data'/'holdout').glob('*.npz'))
    if len(public)!=contract['holdout_count']:
        raise ValueError('Expected exactly frozen target count')
    receipt=ROOT/'audit'/'HOLDOUT_EXECUTION_RECEIPT.json'
    if receipt.exists():raise FileExistsError(receipt)
    config_dir=ROOT/'experiments'/'frozen_configs'
    config_dir.mkdir(exist_ok=True)
    env=os.environ.copy()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
        env[key]='1'
    starts=datetime.datetime.now(datetime.timezone.utc).isoformat()
    rows=[]
    for target in public:
        for method in contract['method_configs']:
            assert_source_frozen()
            config=config_dir/(method['name']+'.json')
            content=json.dumps(method['config'],indent=2)+'\n'
            if config.exists() and config.read_text()!=content:raise RuntimeError('Frozen config changed')
            if not config.exists():config.write_text(content);config.chmod(0o444)
            output=ROOT/'results'/('holdout_'+target.stem+'_'+method['name']+'.json')
            if output.exists():raise FileExistsError(output)
            started=time.perf_counter()
            run=subprocess.run([sys.executable,str(ROOT/'run_solver.py'),'--public',str(target),
                                '--method',method['method'],'--config',str(config),'--output',str(output)],
                               cwd=ROOT,env=env,capture_output=True,text=True)
            row=dict(case_id=target.stem,name=method['name'],method=method['method'],
                     output=str(output.relative_to(ROOT)),process_wall_seconds=time.perf_counter()-started,
                     returncode=run.returncode,stdout=run.stdout,stderr=run.stderr)
            rows.append(row)
            receipt.write_text(json.dumps(dict(started_utc=starts,contract=str(a.contract),
                                               rows=rows,complete=False,truth_read=False),indent=2)+'\n')
            print(json.dumps({k:row[k] for k in ('case_id','name','returncode','process_wall_seconds')}),flush=True)
            if run.returncode:
                raise RuntimeError('Failed solver; retained execution receipt, no silent retries')
    receipt.write_text(json.dumps(dict(started_utc=starts,
                                      finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                      contract=str(a.contract),rows=rows,complete=True,truth_read=False),indent=2)+'\n')


if __name__=='__main__':main()
