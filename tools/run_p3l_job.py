from pathlib import Path
import sys,json,argparse,time,traceback
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,Target,write,save,plain

def main():
    p=argparse.ArgumentParser();p.add_argument('operation',choices=['anchor','reference']);p.add_argument('parent');p.add_argument('--method',default='full_TRF');p.add_argument('--chain',type=int,default=0);p.add_argument('--warmup',type=int,default=1024);p.add_argument('--production',type=int,default=4096);p.add_argument('--max-evaluations',type=int,default=96);p.add_argument('--polish',action='store_true');p.add_argument('--job-suffix',default='');p.add_argument('--stage',default='A');a=p.parse_args()
    if not (OUT/'AUDIT/development_DATA_GENERATOR_EXIT_REQUIRED.json').exists():raise RuntimeError('dedicated observation job missing')
    contract=json.loads((OUT/'AUDIT/PHYSICS_PROBABILITY_CONTRACT.json').read_text())
    if not contract['PASS']:raise RuntimeError('probability/physics interface incomplete')
    job=f'{a.operation}_{a.parent}_{a.method}'+(f'_c{a.chain}' if a.operation=='reference' else '')+a.job_suffix
    b=Book(a.stage,job,max_wall=7200 if a.operation=='reference' else 1200);status='COMPLETE';error=None
    try:
        T=Target(b,keep_fields=a.operation=='anchor')
        if a.operation=='anchor':
            from p3l.inversion import reconstruct
            best,report=reconstruct(T,T.observation(a.parent),a.parent,method=a.method,max_evaluations=a.max_evaluations,polish=a.polish)
            dest=OUT/f'{a.stage}/anchor_runs/{a.parent}_{a.method}{a.job_suffix}.npz';T.persist_cache(best,dest)
            write(dest.with_suffix('.json'),report)
            print(json.dumps(plain(report)),flush=True)
        else:
            from p3l.reference import run_reference
            run_reference(T,a.parent,a.chain,method=a.method,warmup=a.warmup,production=a.production)
    except BaseException as e:
        status='ERROR';error=repr(e);traceback.print_exc();raise
    finally:b.close(status,error)

if __name__=='__main__':main()
