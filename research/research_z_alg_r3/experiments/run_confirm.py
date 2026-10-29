"""Run only a separately frozen V7 confirm contract; no private target input."""
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
    p=argparse.ArgumentParser();p.add_argument('--contract',type=Path,required=True);a=p.parse_args()
    contract=json.loads(a.contract.read_text())
    source_manifest=json.loads((ROOT/contract['source_manifest']).read_text())
    def check_sources():
        for item in source_manifest['files']:
            if (ROOT/item['source']).read_bytes()!=(ROOT/item['snapshot']).read_bytes():
                raise RuntimeError('Source changed after confirm freeze: '+item['source'])
    check_sources()
    targets=sorted((ROOT/'public_data'/'confirm').glob('*.npz'))
    if len(targets)!=contract['confirm_count']:raise ValueError('Wrong confirm target count')
    receipt=ROOT/'audit'/'CONFIRM_V7_EXECUTION_RECEIPT.json'
    if receipt.exists():raise FileExistsError(receipt)
    env=os.environ.copy()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):env[key]='1'
    configs=ROOT/'experiments'/'confirm_v7_configs';configs.mkdir(exist_ok=True)
    rows=[];started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat();batch=time.perf_counter()
    for target in targets:
        for method in contract['method_configs']:
            check_sources()
            config=configs/(method['name']+'.json');content=json.dumps(method['config'],indent=2)+'\n'
            if config.exists() and config.read_text()!=content:raise RuntimeError('Confirm config changed')
            if not config.exists():config.write_text(content);config.chmod(0o444)
            output=ROOT/'results'/('confirm_'+target.stem+'_'+method['name']+'.json')
            if output.exists():raise FileExistsError(output)
            started=time.perf_counter()
            run=subprocess.run([sys.executable,str(ROOT/contract['runner']),'--public',str(target),
                                '--method',method['method'],'--config',str(config),'--output',str(output)],
                               cwd=ROOT,env=env,capture_output=True,text=True)
            row=dict(case_id=target.stem,name=method['name'],method=method['method'],output=str(output.relative_to(ROOT)),
                     process_wall_seconds=time.perf_counter()-started,returncode=run.returncode,stdout=run.stdout,stderr=run.stderr)
            rows.append(row)
            receipt.write_text(json.dumps(dict(started_utc=started_utc,contract=str(a.contract),rows=rows,complete=False,truth_read=False),indent=2)+'\n')
            print(json.dumps({k:row[k] for k in ('case_id','name','returncode','process_wall_seconds')}),flush=True)
            if run.returncode:raise RuntimeError('Confirm process failed; preserve paid receipt; no silent retry')
    receipt.write_text(json.dumps(dict(started_utc=started_utc,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                      batch_wall_seconds=time.perf_counter()-batch,contract=str(a.contract),rows=rows,complete=True,truth_read=False),indent=2)+'\n')


if __name__=='__main__':main()
