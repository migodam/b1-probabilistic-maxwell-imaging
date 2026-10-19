"""Separate dataset process. Closed P2 still refuses every new paid call."""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from material_fiber.core import OUT,Book,Target,write,save

def main():
    if (OUT/'A/OBSERVATIONS_GENERATED.json').exists():raise RuntimeError('observations already generated; no repeats')
    book=Book('A');target=Target(book);error=None
    try:
        parents=json.loads((OUT/'configs/PARENTS.json').read_text())
        for r in parents:
            parent=r['parent']
            if (OUT/f'inputs/online/{parent}.npz').exists():raise RuntimeError('existing observation; reconcile rather than repeat')
            with np.load(OUT/f'inputs/offline_only/{parent}.npz') as a:q=a['q'].copy();noise=a['noise'].copy()
            truth=target.evaluate(q,'DATASET_TRUTH_'+parent)
            save(OUT/f'inputs/online/{parent}.npz',observation=truth.state.data+noise,target_id=np.array(target.target_id),sigma=target.sigma)
            write(OUT/f'offline_only/{parent}_truth.json',dict(key=truth.key,split=r['split']))
        write(OUT/'A/OBSERVATIONS_GENERATED.json',dict(parent_count=len(parents),process='dataset-only; inversion must start in another process',costs_in_same_A_ledger=True))
    except BaseException as e:error=repr(e);raise
    finally:book.close('OBSERVATIONS_COMPLETE' if error is None else 'ERROR',error,final_phase=error is not None)
if __name__=='__main__':main()
