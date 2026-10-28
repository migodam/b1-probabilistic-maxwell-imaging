"""Whitelist-only TRAIN freeze / exposed E0 array access. No locked loader."""
from pathlib import Path
import csv,json
import numpy as np


def safe(root, member):
    root=Path(root).resolve();p=root/member
    if any('locked' in x.lower() or x.lower() in {'test','id_test','ood_test'} for x in p.parts):
        raise ValueError('LOCKED_PATH_FORBIDDEN')
    if not p.is_file() or p.is_symlink() or not p.resolve().is_relative_to(root):
        raise ValueError('INPUT_PATH_MISSING_OR_UNSAFE:'+str(member))
    return p


def source_root(root):return safe(root,'inputs/a29_release/research/A29_PRE_RUN_PROTOCOL.md').parents[1]


def release_root(root):
    p=Path(root).resolve()/'inputs/a29_release'
    if p.is_symlink() or not p.is_dir():raise ValueError('AUTHENTIC_A29_RELEASE_REQUIRED')
    return p


def case_metadata(root):
    r=release_root(root)
    with safe(r,'results/a29/COVERAGE_CASES.csv').open() as f:rows=list(csv.DictReader(f))
    with np.load(safe(r,'results/a29/OFFLINE_TARGETS.npz'),allow_pickle=False) as z:
        a={k:z[k].copy() for k in ('case_id','parent_id','background_id','family')}
    if len(rows)!=32 or set(map(int,a['parent_id']))!=set(range(250065,250081)) or len(set(a['case_id']))!=32:
        raise ValueError('FULL32_EXPOSED_CASES_REQUIRED')
    cases={r['case_id']:r for r in rows}
    if set(cases)!=set(a['case_id']):raise ValueError('CASE_METADATA_MISMATCH')
    for pid in range(250065,250081):
        ix=np.flatnonzero(a['parent_id']==pid)
        if len(ix)!=2 or set(a['background_id'][ix])!={'weak','ring'}:raise ValueError('PAIRED_BACKGROUND_REQUIRED')
    return [dict(case_id=str(c),parent_id=int(p),background_id=str(b),family=str(f))
            for c,p,b,f in zip(a['case_id'],a['parent_id'],a['background_id'],a['family'])]


def online_e0(root):
    """Read only baseline + known backgrounds; not OFFLINE_TARGETS truth."""
    r=release_root(root);rows=case_metadata(root);backgrounds={}
    for bg in ('weak','ring'):
        with np.load(safe(r,'inputs/a25/background/'+bg+'.npz'),allow_pickle=False) as z:
            backgrounds[bg]=z['chi_b'].copy().reshape(576)
    baseline=[]
    for row in rows:
        with np.load(safe(r,'results/a29/pools/'+row['case_id']+'.npz'),allow_pickle=False) as z:
            baseline.append(z['chi'][0].copy())
    return rows,np.asarray(baseline),np.asarray([backgrounds[x['background_id']] for x in rows])


def train_deltas(root):
    r=release_root(root);path=safe(r,'inputs/a25/train.npz')
    with np.load(path,allow_pickle=False) as z:
        p=z['parent_id'].copy();split=z['split'].copy()
        if set(map(int,p))!=set(range(250001,250065)) or set(np.asarray(split).astype(str).ravel())!={'train'}:
            raise ValueError('TRAIN64_IDENTITY_REQUIRED_BEFORE_LABEL_READ')
        b=z['background_id'].astype(str);amp=z['amplitude'].copy();cases=z['case_id'].astype(str)
        selected=[]
        for pid in range(250001,250065):
            ix=np.flatnonzero((p==pid)&(b=='weak')&(amp==1.2))
            if len(ix)!=2:raise ValueError('TRAIN_CASE_LAYOUT')
            selected.append(sorted(ix,key=lambda i:cases[i])[0])
        chi=z['chi'][selected].copy()
    with np.load(safe(r,'inputs/a25/background/weak.npz'),allow_pickle=False) as z:bg=z['chi_b'].copy()
    return (chi-bg[None])/1.2,dict(parents=list(range(250001,250065)),selected_case_ids=cases[selected].tolist(),
                independent_parents=64,background='weak',amplitude=1.2,normalization='(chi-chi_b)/saved amplitude')


def offline_targets(root,case_ids):
    r=release_root(root)
    with np.load(safe(r,'results/a29/OFFLINE_TARGETS.npz'),allow_pickle=False) as z:
        ids=list(map(str,z['case_id']))
        if ids!=list(case_ids):raise ValueError('OFFLINE_CASE_ORDER_MISMATCH')
        truth=z['truth'].copy();bg=z['background'].copy();baseline=z['N3'].copy()
    if truth.shape!=(32,576) or not np.isfinite(truth).all():raise ValueError('REAL32_SAVED_TARGETS_REQUIRED')
    return truth,bg,baseline
