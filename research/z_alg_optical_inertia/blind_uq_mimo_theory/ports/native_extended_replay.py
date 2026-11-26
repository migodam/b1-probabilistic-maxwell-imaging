"""Re-evaluate stable native blind pairs on acquired array configurations."""
import argparse
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.ndimage import label

from model import Green, physical_solve
from sweep import tag

HERE=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--band",type=int,choices=(400,460),required=True)
    p.add_argument("--counts",type=int,nargs="+",choices=(6,8,12,16),required=True)
    p.add_argument("--cases",nargs="+",default=("n96_c01","n96_c02"))
    args=p.parse_args()
    started=perf_counter(),process_time()
    rows,predictions=[],{}
    for case_id in args.cases:
        ext=HERE.parent/"results"/"extended"/f"{case_id}.json"
        receipt=json.loads(ext.read_text())
        zero=(HERE.parent/receipt["zero_file"]).with_suffix(".npz")
        with np.load(ext.with_suffix(".npz")) as a,np.load(zero) as b:
            present,absent=a["material"].ravel(),b["material"].ravel()
            roi=a["roi"].ravel()
            patch=a["patch"]
        region=np.zeros_like(roi,dtype=bool)
        region[patch]=True
        assert receipt["upper_residual_over_r"] < 1 and np.all(absent[roi]==0)
        assert np.min(present[patch])>=2 and label(region.reshape(96,96))[1]==1
        with np.load(HERE/"public"/"cases"/f"{case_id}_mimo.npz") as d:
            sweep=dict(d)
        with np.load(HERE/"public"/"cases"/f"{case_id}.npz") as d:
            primary=dict(d)
        for n in args.counts:
            base=(n,args.band)==(8,400)
            op_id=(str(primary["operator_id"]) if base else f"n96_f{args.band}_a{n:02}_c048")
            with np.load(HERE/"public"/"operators"/f"{op_id}.npz") as d:
                op=dict(d)
            green=Green(op)
            key=tag(args.band,n,48,0.)
            data=primary if base else {k:sweep[key+"_"+k] for k in ("Y","Y_hold","radius","radius_hold")}
            N=[]; L=[]; charges=[]
            for m in (present,absent):
                j,c=physical_solve(m,green.apply,op["E0"],green.diagonal)
                N.append(op["S"]@j)
                L.append(op["S_hold"]@j)
                charges.append(c)
            train=[float(np.linalg.norm(y-data["Y"])/data["radius"]) for y in N]
            hold=[float(np.linalg.norm(y-data["Y_hold"])/data["radius_hold"]) for y in L]
            r=float(data["radius"])
            sep=float(np.linalg.norm(N[0]-N[1])/(2*r))
            result=dict(case=case_id,count=n,band_MHz=args.band,source_total_available_W=1.0,
                        operator_id=op_id,measurement="original blind" if base else "same target paired acquisition",
                        connected_patch_cells=len(patch),patch_average_difference=float(present[patch].mean()),
                        roi_average_difference=float(present[roi].mean()),noise_radius_sqrtW=r,
                        new_train_endpoint_residual_over_r=train,new_hold_endpoint_residual_over_r=hold,
                        pair_separation_over_two_r=sep,costs=charges)
            rows.append(result)
            k=f"{case_id}_f{args.band}_a{n:02}"
            predictions[k+"_fit_present"],predictions[k+"_fit_absent"]=N
            predictions[k+"_hold_present"],predictions[k+"_hold_absent"]=L
            print(case_id,args.band,n,train,hold,sep,flush=True)
    filename=f"native_extended_replay_f{args.band}"
    np.savez_compressed(HERE/"results"/f"{filename}.npz",**predictions)
    (HERE/"results"/f"{filename}.json").write_text(json.dumps(dict(
        rows=rows,wall=perf_counter()-started[0],cpu=process_time()-started[1],
        scope="No optimization or truth access; physically recomputed existing derived material maps",
    ),indent=2))


if __name__=="__main__":
    main()
