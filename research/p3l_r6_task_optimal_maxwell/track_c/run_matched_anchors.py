"""Root-only central paid generation of historical TRAIN512 low48/32 contexts."""
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
import argparse,gc,json,time,sys
from pathlib import Path
import numpy as np
OUT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(OUT))
from shared.physics import Book,Target,write
from track_c.anchors import run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=OUT/'track_c/matched_inputs_TRAIN512')
    parser.add_argument('--output',type=Path,default=OUT/'track_c/matched_anchors_TRAIN512')
    args=parser.parse_args()
    manifest=json.loads((args.input/'PARENTS.json').read_text())
    parents=[str(row['parent']) for row in manifest['parents']]
    if manifest['group']!='HISTORICAL_TRAIN512' or len(parents)!=512:raise ValueError('frozen TRAIN512 export required')
    if not args.output.resolve().is_relative_to((OUT/'track_c').resolve()):raise ValueError('Track C outputs only')
    pending=[p for p in parents if not (args.output/p/'native48_probit32/REPORT.json').exists()]
    auth=json.loads((OUT/'configs/AUTHORIZATION.json').read_text());ledger=Book.totals()
    available=min(auth['stage_caps']['C']-ledger['stages'].get('C',{}).get('RHS',0),auth['total_complex_RHS']-ledger['RHS'])
    reserve=len(pending)*3280
    write(args.output/'PRE_RUN_BUDGET.json',dict(pending_parents=len(pending),worst_RHS_per_parent=3280,
        worst_remaining_RHS=reserve,available_RHS=available,decision='RUN' if reserve<=available else 'HOLD',
        full_training_contexts_reused=True,labels_read=False))
    if reserve>available:raise RuntimeError('HOLD: insufficient remaining authorized C/total RHS for full pending TRAIN512')
    for index,parent in enumerate(parents):
        destination=args.output/parent/'native48_probit32'
        if (destination/'REPORT.json').exists():continue
        with np.load(args.input/'online'/f'{parent}.npz',allow_pickle=False) as z:
            if set(z.files)!={'observation','sigma','parent','target_id'} or str(z['parent'])!=parent:raise ValueError('TRAIN online allowlist')
            y,sigma=z['observation'].copy(),float(z['sigma'])
        book=Book('C',f'matched_anchor_{parent}_native48_probit32',parent,'HISTORICAL_TRAIN_ANCHOR',max_wall=120)
        status,error='COMPLETE',None;tick=time.perf_counter()
        try:
            target=Target(book,sigma=sigma)
            context,report=run(target,y,parent,destination,'native48_probit32')
            write(destination/'REPORT.json',dict(report,GT_read=False,sigma=sigma,
                                                end_to_end_anchor_wall_seconds=time.perf_counter()-tick))
            target.backend()._transfers.clear();del target,context,report;gc.collect()
        except BaseException as exc:status,error='FAILED',repr(exc);raise
        finally:book.close(status,error)
        if (index+1)%16==0:print(json.dumps(dict(event='MATCHED_ANCHORS_PROGRESS',completed=index+1,RHS=Book.totals()['RHS'])),flush=True)
    write(args.output/'ANCHORS_COMPLETE.json',dict(status='COMPLETE',n=512,parents=parents,config='native48_probit32',GT_read=False))


if __name__=='__main__':main()
