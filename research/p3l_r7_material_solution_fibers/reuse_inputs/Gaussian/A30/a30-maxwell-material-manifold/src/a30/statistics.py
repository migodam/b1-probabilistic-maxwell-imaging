"""Exactly two conditions averaged before exploratory parent paired statistics."""
from collections import defaultdict
import numpy as np


COMPLEX={'shell','asymmetric_piecewise','mixed_scale'}


def paired(rows, candidate, *, baseline='N3_nrmse', families=COMPLEX):
    group=defaultdict(list)
    for r in rows:
        if r['family'] in families:group[int(r['parent_id'])].append(r)
    expected=[p for p in range(250065,250081) if (
        ('shell' if p<250069 else 'asymmetric_piecewise' if p<250073 else 'smooth' if p<250077 else 'mixed_scale') in families)]
    out=[]
    for p in expected:
        a=group.get(p,[])
        if len(a)!=2 or {x['background_id'] for x in a}!={'weak','ring'}:raise ValueError('PAIRED_PARENTS_INCOMPLETE')
        x=np.mean([float(r[candidate]) for r in a]);b=np.mean([float(r[baseline]) for r in a])
        if not np.isfinite(x+b) or b<=0:raise ValueError('UNDEFINED_PARENT_METRIC')
        out.append(dict(parent_id=p,baseline=float(b),candidate=float(x),gain=float((b-x)/b)))
    v=np.asarray([x['gain'] for x in out]);rng=np.random.default_rng(291011)
    bs=v[rng.integers(0,len(v),(2000,len(v)))].mean(axis=1)
    return dict(parent_count=len(v),point=float(v.mean()),median=float(np.median(v)),ci95=np.quantile(bs,[.025,.975]).tolist(),parents=out,
                statistical_unit='EXPOSED_PARENT_EXPLORATORY',bootstrap_samples=2000,seed=291011)


def paired_difference(rows,candidate,baseline, *, families=COMPLEX):
    out=[]
    for pid in range(250065,250081):
        pair=[r for r in rows if int(r['parent_id'])==pid and r['family'] in families]
        if not pair:continue
        if len(pair)!=2:raise ValueError('MISSING_BACKGROUND_FOR_DIFFERENCE')
        out.append(float(np.mean([float(x[candidate])-float(x[baseline]) for x in pair])))
    v=np.asarray(out);rng=np.random.default_rng(291011)
    bs=v[rng.integers(0,len(v),(2000,len(v)))].mean(1)
    return dict(point=float(v.mean()),ci95=np.quantile(bs,[.025,.975]).tolist(),parent_count=len(v))
