"""Prove disconnection from an empty projection interval and feasible material witnesses."""
import json
import math
from pathlib import Path
import numpy as np
from data14 import ROOT, load_cases
from material_controls import GRID,LO,HI,K

P=ROOT/'results'
def run():
    pool=json.loads((P/'pivot_diameter_witnesses.json').read_text())['cases']
    records={}
    for case in load_cases():
        cid=case.case_id
        points=[]
        for row in pool[cid]['actions'].values():
            if row['witness']:
                points.extend(row['witness']['materials'])
        pts=np.unique(points,axis=0)
        boxes=np.load(P/'material_covers_seeded'/(cid+'.npz'))['boxes']
        dirs=list(np.eye(6))
        if len(pts)>1:
            dirs.extend(np.linalg.svd(pts-pts.mean(0),full_matrices=False)[2])
            dirs.extend(pts[i]-pts[j] for i,j in zip(*np.triu_indices(len(pts),1)))
        best=None
        for d in dirs:
            if max(abs(d))==0:continue
            v=np.rint(512*d/max(abs(d))).astype(np.int64)
            vp,vm=np.maximum(v,0),np.minimum(v,0)
            za=boxes[:,0]@vp+boxes[:,1]@vm
            zb=boxes[:,1]@vp+boxes[:,0]@vm
            order=np.argsort(za)
            a,b=za[order],np.maximum.accumulate(zb[order])
            gaps=np.flatnonzero(a[1:]>b[:-1])
            if not len(gaps):continue
            # Exact binary64 witness inputs give exact rational projections here.
            from fractions import Fraction as Q
            truth=[sum(int(vi)*(Q(float(x))-Q(1,2))*GRID/Q(59,2) for x,vi in zip(pt,v)) for pt in pts]
            for j in gaps:
                left=[i for i,f in enumerate(truth) if f<=int(b[j])]
                right=[i for i,f in enumerate(truth) if f>=int(a[j+1])]
                if not left or not right:continue
                width=(int(a[j+1])-int(b[j]))*(HI-LO)/GRID/np.linalg.norm(v)
                if best is None or width>best['empty_slab_material_L2_display']:
                    best=dict(v=v.tolist(),gap_integer=[int(b[j]),int(a[j+1])],
                              old_materials=[pts[left[0]].tolist(),pts[right[0]].tolist()],
                              empty_slab_material_L2_display=float(width))
        records[cid]=best
    (P/'disconnected_old.json').write_text(json.dumps(records,indent=2)+'\n')
    print([(k, round(v['empty_slab_material_L2_display'],3)) for k,v in records.items() if v])
if __name__=='__main__':run()
