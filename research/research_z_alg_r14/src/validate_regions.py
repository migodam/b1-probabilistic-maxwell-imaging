"""Independent materials at two matched noise levels; use only the public acquisition file."""
import argparse
import json
from time import perf_counter

import numpy as np

from material_controls import GRID, K, RESULTS, cover
from pivot_experiment import run_arm
from posterior import run


def validate():
    from measure_validation import load_validation_cases, load_validation_reads
    cases=load_validation_cases();reads=load_validation_reads()
    folder=RESULTS/"validation_regions";folder.mkdir(exist_ok=True)
    for case in cases:
        old,boxes=run_arm(case,"r10_scalar",.70,4)
        seed,shared=cover(case,.12,16000,start_boxes=boxes)
        info,arrays=run(case,reads[case.case_id+"/split_0.02"],4500,numerator=True,seed=seed)
        np.savez_compressed(folder/(case.case_id+".npz"),old=seed,**arrays)
        (folder/(case.case_id+".json")).write_text(json.dumps(dict(r10=old,old=shared,posterior=info),indent=2)+"\n")
        print(case.case_id,[(x,round(k["D_upper"],4)) for x,k in info.items()],flush=True)


def start_ablations():
    from data14 import load_cases
    first_by_geom=[f"spectral-{v:02d}" for v in range(0,36,6)]
    names=set(first_by_geom+["fresh6-00","fresh6-06","fresh6-12","confirm6-00","confirm6-05"])
    folder=RESULTS/"posterior_starts";folder.mkdir(exist_ok=True)
    with np.load(RESULTS/"measure_new_reads.npz") as reads, np.load(RESULTS/"pivot_old_regions.npz") as initial:
        for case in load_cases():
            if case.case_id not in names: continue
            lo=np.zeros(K,dtype=np.int64);hi=np.full(K,GRID,dtype=np.int64)
            whole=np.array([[lo,hi]])
            row={}
            for label,seed in (("whole_box",whole),("r10_only",initial[case.case_id+"/r10_scalar"])):
                if label=="r10_only":
                    from material_controls import HI,LO
                    x=(seed-LO)*(GRID/(HI-LO))
                    seed=np.stack((np.maximum(0,np.floor(x[:,0])-1),np.minimum(GRID,np.ceil(x[:,1])+1)),axis=1).astype(np.int64)
                info,_=run(case,reads[case.case_id+"/split_0.02"],4500,numerator=True,seed=seed)
                row[label]=info
            (folder/(case.case_id+".json")).write_text(json.dumps(row,indent=2)+"\n")
            print(case.case_id,[(x,round(k["first_half"]["D_upper"],4)) for x,k in row.items()],flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--starts",action="store_true");a=p.parse_args()
    if a.starts:start_ablations()
    else:validate()
