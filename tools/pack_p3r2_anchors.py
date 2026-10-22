"""Permitted supervised train/validation labels; never accepts locked parents."""
from pathlib import Path
import sys,argparse
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import OUT,save,write
def main():
    p=argparse.ArgumentParser();p.add_argument('group',choices=['train','validation']);p.add_argument('--count',type=int,required=True);a=p.parse_args();rows=[]
    for k in range(a.count):
        parent=f'R2_{a.group}_{k:05d}'
        with np.load(OUT/f'A/contexts/{parent}.npz') as d:r={key:d[key].copy() for key in d.files}
        with np.load(OUT/f'inputs/offline_only/{parent}.npz') as d:r.update(chi_GT=d['chi'].copy(),q_GT=d['q'].copy())
        r['parent']=np.array(parent);rows.append(r)
    save(OUT/f'A/anchor_data/{a.group}.npz',**{key:np.stack([r[key] for r in rows]) for key in rows[0]})
    write(OUT/f'A/anchor_data/{a.group}.json',dict(count=a.count,GT_scope='explicitly authorized train/validation material risk supervision',locked_access=False,physics_calls=0))
if __name__=='__main__':main()
