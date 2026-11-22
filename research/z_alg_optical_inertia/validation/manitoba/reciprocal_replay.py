"""Apply the ordinary reciprocal projection and replay the original strict profiles."""
import argparse
import json
from pathlib import Path

import numpy as np
from flint import arb, ctx
from threadpoolctl import threadpool_limits

from physical import FREQUENCIES, HERE, MASK, PREPARED, TARGETS, USED, elapsed, timer
from scan import hi, lo, verify

PROFILE_FIELDS=("mu","nu","numeric_lower","t","cutoff","c","A_negative","A_compact_size")


def norm_squared(values, other=None):
    if other is None:
        return sum((arb(float(v.real))**2+arb(float(v.imag))**2 for v in values),arb(0))
    return sum(((arb(float(v.real))-arb(float(w.real)))**2 +
                (arb(float(v.imag))-arb(float(w.imag)))**2
                for v,w in zip(values,other,strict=True)),arb(0))


def fix_profiles(full_result):
    rows=[]
    for row in full_result["rows"]:
        frozen=[]
        for group in row["groups"]:
            selected=next((q for q in group["attempts"] if q.get("verified_lower",0.)>0), None)
            proposals=[{k:selected[k] for k in PROFILE_FIELDS}] if selected else []
            frozen.append(dict(group=group["group"],proposals=proposals))
        rows.append({**{k:v for k,v in row.items() if k not in {
            "groups","verified_lower_relative","verified_lower"}}, "groups":frozen})
    return {"rows":rows}


def main(prepared, primary, physical, output):
    start=timer()
    ctx.prec=96
    local=prepared.parent/"reciprocal"
    local.mkdir(parents=True,exist_ok=True)
    private_projected=local/"projected_private.npz"
    obs={}
    with np.load(prepared) as f:
        arrays=dict(f)
    for i,freq in enumerate(FREQUENCIES):
        for t in TARGETS:
            key=f"f{i}_{t}"
            y=arrays[key]
            yp=(y+y.T)/2
            yn=norm_squared(y[USED]).sqrt()
            perturbation=norm_squared(y[USED],yp[USED]).sqrt()
            full_perturbation=norm_squared(y[MASK],yp[MASK]).sqrt()
            projection_norm=norm_squared(yp[USED]).sqrt()
            obs[(float(freq),t)]=dict(original_norm=yn,projected_norm=projection_norm,
                                     perturbation=perturbation,full_perturbation=full_perturbation)
            arrays[key]=yp
    np.savez_compressed(private_projected,**arrays)
    projected=verify(private_projected,fix_profiles(primary))
    records={row["id"]:row for row in physical["cases"]}
    original={(r["case_id"],r["cap"],r["arm"]):r for r in primary["rows"]}
    results=[]
    for row in projected["rows"]:
        meta=records[row["case_id"]]
        audit=obs[(row["frequency_GHz"],row["target"])]
        unprojected=audit["original_norm"]
        pert=audit["perturbation"]
        bound=arb(row["verified_lower"])
        r={
            **{k:v for k,v in row.items() if k not in {"field_norm","numeric_bound_relative"}},
            "primary_verified_lower_relative":original[(row["case_id"],row["cap"],row["arm"])]["verified_lower_relative"],
            "projected_norm_in_primary_units_upper":hi(audit["projected_norm"]),
            "projected_certified_lower_over_primary_norm":lo(bound/unprojected) if bound>0 else 0.,
            "projected_minus_original_field_norm_upper":hi(pert),
            "projected_minus_original_relative_upper":hi(pert/unprojected),
            "projected_minus_original_field_norm_all_552_upper":hi(audit["full_perturbation"]),
            "primary_two_calibrator_difference_relative":meta["two_reference_columns_difference_relative"],
            "primary_drift_plus_held_reference_relative":meta["drift_plus_holdout_over_scatter"],
        }
        r["margin_relative_after_both_primary_proxies"]=r["projected_certified_lower_over_primary_norm"]-max(
            r["primary_two_calibrator_difference_relative"], r["primary_drift_plus_held_reference_relative"])
        results.append(r)
    counters=[]
    for cap in (.25,.5,.75,1.):
        for arm in ("real","positive_optical","free_optical"):
            group=[r for r in results if r["cap"]==cap and r["arm"]==arm]
            counters.append(dict(cap=cap,arm=arm,count=len(group),
                                 primary_positive=sum(r["primary_verified_lower_relative"]>0 for r in group),
                                 projected_positive=sum(r["verified_lower"]>0 for r in group),
                                 projected_exceeds_both_primary_proxies=sum(
                                     r["margin_relative_after_both_primary_proxies"]>0 for r in group)))
    output.write_text(json.dumps(dict(
        projection="ordinary transpose average (Y+Y.T)/2 before splitting by the same frozen port masks",
        selection="only original interval-passing candidate per case/cap/arm/group is reused; no search",
        proxy_denominator="original calibrated selected data norm, original calibration proxies unchanged",
        profiles_source="verified_scan.json", rows=results, counts=counters,
        costs=projected["costs"], verification_wall_seconds=projected["wall_seconds"],
        **elapsed(start)),indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared",type=Path,default=PREPARED)
    parser.add_argument("--output",type=Path,default=HERE/"reciprocal_verified.json")
    args=parser.parse_args()
    primary=json.loads((HERE/"verified_scan.json").read_text())
    physical=json.loads((HERE/"physical_audit.json").read_text())
    with threadpool_limits(1):
        main(args.prepared,primary,physical,args.output)
