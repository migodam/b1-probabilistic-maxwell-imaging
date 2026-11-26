"""Independent MIMO replay of complete image-generated connected-patch pairs."""
import argparse
import csv
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.ndimage import label

from model import Green, physical_solve
from sweep import tag
from witness import pair_grid

HERE=Path(__file__).resolve().parent
COHORT=[f"n32_c{i:02}" for i in range(1,9)]


def fixed_audit(extended,zero,observation,op):
    values=[]
    green=Green(op)
    for m in (extended,zero):
        j,cost=physical_solve(m,green.apply,op["E0"],green.diagonal,op["G"])
        vals=dict(train_r=float(np.linalg.norm(op["S"]@j-observation["Y"])/observation["radius"]),
                  hold_r=float(np.linalg.norm(op["S_hold"]@j-observation["Y_hold"])/observation["radius_hold"]),
                  residual_max=max(cost["state_residual"]),RHS=cost["physical_rhs"],
                  LU=cost["LU"],seconds_wall=cost["seconds_wall"])
        values.append(vals)
    return values


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cases", nargs="+", default=COHORT)
    args=p.parse_args()
    t,tc=perf_counter(),process_time()
    with np.load(HERE/"public/operators/n32_f400_a08.npz") as op:
        reference=dict(op)
    all_rows,cases=[],[]
    for oid in args.cases:
        extended_path=HERE.parent/"results"/"extended"/f"{oid}.json"
        if not extended_path.exists():
            cases.append(dict(case=oid,status="no complete extension receipt"))
            continue
        pinfo=json.loads(extended_path.read_text())
        if pinfo["upper_residual_over_r"] is None or pinfo["upper_residual_over_r"]>=1:
            cases.append(dict(case=oid,status="extension not certified",receipt=str(extended_path.relative_to(HERE.parent))))
            continue
        zero_file=(HERE.parent/pinfo["zero_file"]).with_suffix(".npz")
        with np.load(extended_path.with_suffix(".npz")) as e,np.load(zero_file) as z:
            x,y=e["material"].ravel(),z["material"].ravel()
            patch=e["patch"]
            roi=e["roi"].ravel()
        patchmask=np.zeros_like(roi,dtype=bool)
        patchmask[patch]=True
        component=label(patchmask.reshape(32,32))[1]
        assert np.all(roi[patch]) and component==1 and np.min(x[patch])>=2 and np.all(y[roi]==0)
        assert min(x.min(),y.min())>=0 and max(x.max(),y.max())<=4.5
        with np.load(HERE/"public/cases"/f"{oid}.npz") as obs:
            fixed=fixed_audit(x,y,obs,reference)
            sigma=float(obs["sigma"])
        if any(item["train_r"]>=1 for item in fixed):
            cases.append(dict(case=oid,status="reference pair does not share actual fit ball",first=fixed[0],second=fixed[1]))
            continue
        s0=perf_counter()
        rows,arrays=pair_grid(x,y,sigma)
        with np.load(HERE/"public"/"cases"/f"{oid}_mimo.npz") as observed:
            for row in rows:
                key=tag(row["frequency_MHz"],row["q"],row["radiation_scale_ohms"],row["calibration_eta"])
                columns={}
                for state in ("with","without"):
                    for split in ("","_hold"):
                        ports=arrays[key+f"_F_{state}_region"+split]
                        residual=np.linalg.norm(ports-observed[key+"_Y"+split])/float(observed[key+"_radius"+split])
                        columns[f"observed_{state}_patch{'_hold' if split else '_train'}_over_r"]=residual
                all_rows.append(dict(case_id=oid,roi_index=pinfo["roi_index"],patch_cells=len(patch),
                                     patch_area_cm2=len(patch)*39.0625,patch_minimum_difference=float(np.min(x[patch])),
                                     patch_average_difference=float(x[patch].mean()),
                                     roi_average_difference=float(x[roi].mean()),
                                     primary_ext_train_r=fixed[0]["train_r"],primary_zero_train_r=fixed[1]["train_r"],
                                     **columns,**row))
        cases.append(dict(case=oid,status="paired arrays replayed",patch_cells=len(patch),roi_cells=int(roi.sum()),
                          mean_patch_gap=float(x[patch].mean()),mean_roi_gap=float(x[roi].mean()),
                          proof_ref_ext=pinfo["upper_residual_over_r"],first=fixed[0],second=fixed[1],
                          first_material_file=str(extended_path.with_suffix(".npz").relative_to(HERE.parent)),
                          zero_material_file=str(zero_file.relative_to(HERE.parent)),array_replay_wall=perf_counter()-s0,
                          paid_bare_LU=4,paid_bare_rhs=288))
        print(oid,len(patch),float(x[patch].mean()),float(x[roi].mean()),
              fixed[0]["train_r"],fixed[1]["train_r"],fixed[0]["hold_r"],fixed[1]["hold_r"],flush=True)
    if not all_rows:
        raise RuntimeError("No available completed alternative")
    with (HERE/"results"/"extended_array_replay.csv").open("w") as f:
        writer=csv.DictWriter(f,fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    (HERE/"results"/"extended_array_replay_cost.json").write_text(json.dumps(dict(
        cohort=cases,evaluated_rows=len(all_rows),total_wall=perf_counter()-t,total_CPU=process_time()-tc,
        acquisition_replay="All arrays reuse the same source union for each frequency, no inverse refits",
        semantics="Each paired chosen material and its original measured reference case; multiple measurements per material remain paired",
    ),indent=2))


if __name__=="__main__":
    main()
