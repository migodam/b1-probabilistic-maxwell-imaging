"""Read exposed frozen truth after scoring and check physical alternate witnesses."""
from collections import Counter
import json

from common import PRIVATE, PUBLIC, REGION, RESULTS, load, save, write_json
from bounds import Acquisition, ctx
import numpy as np

ctx.prec = 192


def main():
    assert (RESULTS/"DECISION_SEAL.json").exists() and (RESULTS/"summary.json").exists()
    source, truth = load(PUBLIC),load(PRIVATE)
    records = []
    exemplar = None
    whole_exemplar = None
    for pair in range(24):
        prefix=f"pair_{pair+1:02d}"
        raw=json.loads((RESULTS/"raw"/f"{prefix}.json").read_text())
        frames=load(RESULTS/"raw"/f"{prefix}_arrays.npz")
        for row in raw["rows"]:
            target,bg,w,slot=row["box_positions"]
            if slot==truth["generating_branch"][pair,target] or row["status"]["C"]!="WITNESS_FOUND":
                continue
            true_mat=truth["material"][pair,target,bg]
            ell,upp=source["boxes"][pair,bg,w,slot]
            if np.all(ell<=true_mat) and np.all(true_mat<=upp):
                continue
            mat=frames[row["witness"]["material"]]
            acq=Acquisition(source["G"],source["S"],source["B"],source["Y"][pair,target,bg],
                            float(source["radii"][pair,target,bg]))
            witnesses=[acq.strict_witness(mat,ell,upp),
                       acq.strict_witness(true_mat,np.zeros(acq.n),np.full(acq.n,32.))]
            if not all(item["proved"] for item in witnesses):
                raise ArithmeticError("Stored witness did not independently replay")
            tpeak,mpeak=map(float,(np.max(true_mat[REGION]),np.max(mat[REGION])))
            same_local=np.array_equal(true_mat[REGION]>0,mat[REGION]>0)
            record=dict(pair=pair+1,target_id=row["target_id"],target=target,background_cap=row["background_cap"],
                        relative_width=row["width_ratio"],branch_candidate=slot,
                        true_peak=tpeak,witness_peak=mpeak,peak_reduction=tpeak-mpeak,
                        same_actual_local_support=bool(same_local),
                        same_actual_whole_support=bool(np.array_equal(true_mat>0,mat>0)),
                        family=str(source["families"][pair]),contract_radius=acq.radius,
                        truth_residual_upper=witnesses[1]["residual_upper"],
                        original_residual_lower=witnesses[0]["residual_lower"],
                        original_residual_upper=witnesses[0]["residual_upper"],
                        radius_minus_residual_upper=acq.radius-witnesses[0]["residual_upper"],
                        upper_over_r=witnesses[0]["residual_upper"]/acq.radius,
                        witness_array=row["witness"]["material"])
            records.append(record)
            if same_local and tpeak>mpeak and (exemplar is None or tpeak-mpeak>exemplar[0]["peak_reduction"]):
                exemplar=record,mat,true_mat,ell,upp,acq.Y
            if record["same_actual_whole_support"] and tpeak>mpeak and (
                    whole_exemplar is None or tpeak-mpeak>whole_exemplar[0]["peak_reduction"]):
                whole_exemplar=record,mat,true_mat,ell,upp,acq.Y
    assert len(records)==32
    counts=Counter((r["family"],"lower" if r["peak_reduction"]>0 else "higher_or_equal",
                    "same_local" if r["same_actual_local_support"] else "different_local") for r in records)
    report=dict(scope="Existing exposed post-score witnesses, no new acquisitions or proposals",
                totals={":".join(k):v for k,v in counts.items()},
                independent_pairs=sorted({r["pair"] for r in records}),
                paired_examples=records,largest_same_local_peak_reduction=exemplar[0] if exemplar else None,
                largest_same_whole_support_peak_reduction=whole_exemplar[0] if whole_exemplar else None)
    for example,name in ((exemplar,"alternate_witness_example"),
                         (whole_exemplar,"alternate_whole_support_example")):
        if example is None:
            continue
        record,mat,true_mat,ell,upp,Y=example
        save(RESULTS/f"{name}.npz",material_witness=mat,material_generator=true_mat,
             radius=np.array(record["contract_radius"]),
             ell=ell,upp=upp,Y=Y,G=source["G"],S=source["S"],B=source["B"])
    write_json(RESULTS/"alternate_witness_audit.json",report)
    print(json.dumps({k:v for k,v in report.items() if k!="paired_examples"},indent=2))


if __name__=="__main__":
    main()
