"""Independent bare VIE plus explicit 25-port loaded-circuit point replay."""
from frequency_classes import HERE, elapsed, fro, material, now, np, rd, up, write_json, verify_midpoint
from probe import Circuit, R, at_load
from scipy.linalg import lu_factor, lu_solve
import json


def evaluate(side,band,ix,full_x):
    t0=now()
    w=Circuit(side,band,ix)
    local=full_x[ix]
    a=np.eye(len(ix))-local[:,None]*w.g
    rhs=local[:,None]*w.b
    f=lu_factor(a,check_finite=False)
    j=lu_solve(f,rhs,check_finite=False)
    state_defect=fro(a@j-rhs)/max(fro(rhs),1e-100)
    loaded=w.z0-1j*w.c@j
    s=np.eye(25)-2*R*np.linalg.solve(loaded+R*np.eye(25),np.eye(25))
    s0=np.eye(25)-2*R*np.linalg.solve(w.z0+R*np.eye(25),np.eye(25))
    observed=(at_load(s,1.)-at_load(s0,1.))[8:16,:8]/np.sqrt(8)
    return observed,j,s,dict(verified_center_replay_numeric_defect_Fro=state_defect,
                             bare_material_LU=1,bare_material_RHS=25,circuit_solves=2,**elapsed(t0))


def main():
    total=now()
    rows=[]
    arrays={}
    for side in (32,96):
        t,name,r,obs,ix,roi,a,b=material(side)
        data=json.loads((HERE/f"frequency_classes_{side}.json").read_text())
        with np.load(HERE/f"frequency_classes_{side}.npz",allow_pickle=False) as raw:
            for band in (400,460):
                suffix=f"{side}_{band}"
                ref=data["frequency"][str(band)]
                ps=[]
                packets=[]
                for tag,x in (("plus",a),("minus",b)):
                    p,j,s,stats=evaluate(side,band,ix,x)
                    guard=next(q for label,q in zip(("plus","minus"),ref["open_vs_frozen"]) if label==tag)
                    delta=fro(p-raw[f"p_{band}_{tag}"])
                    if delta>=guard["enclosure_compatibility_slack"]:
                        raise ArithmeticError("Bare wave / loaded port replay disagrees with certified frozen coefficients")
                    arrays[f"response_{suffix}_{tag}"]=p
                    arrays[f"bare_currents_{suffix}_{tag}"]=j
                    arrays[f"scattering_{suffix}_{tag}"]=s
                    stats.update(contrast_state=tag,raw_difference_vs_frozen_norm=delta,
                                 agreed_with_frozen=True,guard=guard["enclosure_compatibility_slack"])
                    packets.append(stats)
                    ps.append(p)
                arr=ps[0]-ps[1]
                d=np.linalg.norm(arr)/(2*r)
                result=dict(side=side,frequency_MHz=band,source_file=name,readout_radius=r,
                            radius_source_id=obs,center_separation_independent_over_two_r=d,
                            incident_W=1.,readout_complex_channels=64,results=packets)
                if band==400:
                    mid=raw["midpoint_400"]
                    result.update(common_midpoint_audit_over_r=[float(np.linalg.norm(p-mid)/r) for p in ps],
                                  strict_original_common_midpoint_bound=ref["midpoint_residual_upper_over_r"])
                    if max(result["common_midpoint_audit_over_r"])>=1:
                        raise RuntimeError("Independent midpoint replay failed")
                    arrays[f"midpoint_{suffix}"]=mid
                rows.append(result)
    with (HERE/"frequency_classes_independent.npz").open("wb") as f:
        np.savez_compressed(f,**arrays)
    write_json(HERE/"frequency_classes_independent.json",dict(rows=rows,**elapsed(total)))
    print([(a["side"],a["frequency_MHz"],a["center_separation_independent_over_two_r"]) for a in rows],flush=True)


if __name__=="__main__":
    main()
