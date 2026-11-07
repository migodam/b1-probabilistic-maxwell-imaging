"""Infer every frozen NN on anonymous TEST inputs, with one operator upload."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
import argparse,json,sys,time
from pathlib import Path
S=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(S/'src'))
from networks import OperatorCache,choose_device
from predict import predict_anchors,_sync
from predict_recurrent import predict_recurrent_anchors

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--device',default='mps');a=parser.parse_args()
    freeze=json.loads((S/'audit/METHOD_FREEZE.json').read_text())
    final=json.loads((S/'audit/FINAL_RUNS.json').read_text())
    records=json.loads((S/'data/manifests/online.json').read_text())['records']
    paths=[S/'anchors'/f"{r['obs_id']}.npz"for r in records if r['split']=='test']
    device=choose_device(a.device);started=time.perf_counter()
    operators=OperatorCache(S/'data/operators.npz').to(device);_sync(device)
    setup_seconds=time.perf_counter()-started;rows=[]
    receipt=S/f'predictions/prediction_runs_{time.time_ns()}.json';receipt.parent.mkdir(exist_ok=True)
    for p in final['plans']:
        if p['kind']!='NN':continue
        directory=S/'predictions'/p['run']
        ready=[v for v in paths if (directory/v.name).exists()and(directory/v.with_suffix('.json').name).exists()]
        pending=[v for v in paths if v not in ready]
        predictor=predict_recurrent_anchors if p['checkpoint_family']=='recurrent'else predict_anchors
        if pending:
            predictor(S/p['checkpoint'],pending,S/'data/operators.npz',directory,a.device,operators=operators)
        rows.append(dict(run=p['run'],public_TEST_observations=len(paths),new_predictions=len(pending),reused_predictions=len(ready)))
        receipt.write_text(json.dumps(dict(status='RUNNING',global_operator_setup_seconds=setup_seconds,operator_uploads=1,runs=rows,truth_access='NONE',method_freeze='audit/METHOD_FREEZE.json'),indent=2)+'\n')
        print(json.dumps(rows[-1]),flush=True)
    receipt.write_text(json.dumps(dict(status='COMPLETE_FROZEN_TEST_NN_PREDICTIONS',global_operator_setup_seconds=setup_seconds,operator_uploads=1,runs=rows,truth_access='NONE',method_freeze='audit/METHOD_FREEZE.json',test_observations=freeze['test_observations']),indent=2)+'\n')
    print(str(receipt),flush=True)

if __name__=='__main__':main()
