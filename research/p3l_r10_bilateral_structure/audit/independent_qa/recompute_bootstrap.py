"""Persisted producer for the independent one-off paired bootstrap calculation.

Defaults to a new replay artifact so the original canonical receipt is preserved.
Pure statistics on already independently rescored parent rows; no physics/GT IO.
"""
from pathlib import Path
import argparse,json,time
import numpy as np

A=Path(__file__).resolve().parent
P=A.parents[1]

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=A/'INDEPENDENT_PAIRED_BOOTSTRAP_REPLAY.json')
    args=parser.parse_args()
    rows=json.loads((A/'FRESH10_SCORE_RECOMPUTATION.json').read_text())['rows']
    config=json.loads((P/'config/candidate_fresh10_statistics.json').read_text())
    rng=np.random.default_rng(np.random.SeedSequence(config['bootstrap_seed']))
    indices=rng.integers(0,10,size=(20000,10))
    output=[]
    for pair in config['primary_family']:
        metric=pair['metric']
        first=sorted([r for r in rows if r['arm']==pair['candidate'] and r['K']==8],key=lambda r:r['parent'])
        second=sorted([r for r in rows if r['arm']==pair['reference'] and r['K']==8],key=lambda r:r['parent'])
        assert len(first)==len(second)==10
        assert [r['parent'] for r in first]==[r['parent'] for r in second]
        delta=np.array([float(x[metric])-float(y[metric]) for x,y in zip(first,second)])
        means=delta[indices].mean(axis=1)
        output.append(dict(**pair,parents=[r['parent'] for r in first],delta=delta.tolist(),mean_delta=float(delta.mean()),CI9875=np.quantile(means,[.00625,.99375]).tolist(),wins=int(np.sum(delta>0)) if metric.startswith('coverage') else int(np.sum(delta<0)),losses=int(np.sum(delta<0)) if metric.startswith('coverage') else int(np.sum(delta>0))))
    result=dict(new_RHS=0,utc=time.time(),repeats=20000,seed=config['bootstrap_seed'],paired_draws_shared_across_family=True,direction='candidate minus reference; positive coverage or negative MSE favors candidate',results=output)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(args.output)

if __name__=='__main__':main()
