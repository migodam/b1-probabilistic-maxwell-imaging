"""Explicit closure and offline scoring; TEST requires a method freeze record."""
import argparse,json,sys,os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
S=Path(__file__).resolve().parents[1];sys.path.insert(0,str(S.parents[1]))
from research.research_gn_nn_preservation_r4.src.score import seal_predictions,score
p=argparse.ArgumentParser();p.add_argument('command',choices=['seal','score']);p.add_argument('--runs',nargs='+');p.add_argument('--split',choices=['val','test'],required=True);p.add_argument('--seal',type=Path);p.add_argument('--output',type=Path);a=p.parse_args()
seal=a.seal or S/'evaluation'/f'{a.split.upper()}_SEALED.json'
actual_splits=json.loads(seal.read_text())['splits']if a.command=='score'else[a.split]
if set(actual_splits)!={a.split}:
    raise ValueError('Requested split must match the actual sealed cohort')
if 'test' in actual_splits and not (S/'audit/METHOD_FREEZE.json').exists():
    raise RuntimeError('TEST scoring starts after final method freeze, not during validation selection')
if a.command=='seal':
    print(seal_predictions(S,S/'data/manifests/evaluation.json',a.runs,splits=(a.split,),output=seal))
else:print(score(seal,a.output or S/'results'/f'{a.split.upper()}_SCORES.json'))
