"""Held-out native 96² response and independent 192² pulse-grid sensitivity."""
from pathlib import Path
import json

from ports import (A,B,INCIDENCE,PRIOR,PortEngine,cut_id,dump,elapsed,load,now,project,save)
from rank import cut_power
import numpy as np
from model import DiskGreen, network_terms, physical_solve
from threadpoolctl import threadpool_limits

OUT=Path(__file__).resolve().parent/"native"
PRETRAIN=Path(__file__).resolve().parent/"fixed_pair"


def material_pair(case):
    p=PRIOR/"results"/"extended"/f"{case}.json"
    card=json.loads(p.read_text())
    one=load(p.with_suffix(".npz"),("material","roi","patch"))
    zero=load((PRIOR/card["zero_file"]).with_suffix(".npz"),("material",))["material"].ravel()
    return zero,one["material"].ravel(),one["roi"].ravel().astype(bool),one["patch"],card


def fine_bare_response(x,side):
    started=now()
    disk=DiskGreen(side,400)
    n=network_terms(disk,8,48)
    excitation=n["port_excitation"]
    readout=n["port_readout"]
    # The fine operator and every coefficient are independently recreated.
    action=lambda a: disk.kernel.apply(a[:,None])[:,0] if a.ndim==1 else disk.kernel.apply(a)
    j,cost=physical_solve(x,action,excitation,disk.kernel.diagonal)
    z_load=n["antenna"]-1j*readout@j
    eye=np.eye(24)
    r=50.
    s0=eye-2*r*np.linalg.solve(n["matrix"],eye)
    s=eye-2*r*np.linalg.solve(z_load+r*eye,eye)
    return s-s0,s,dict(**cost,**elapsed(started))


def main():
    started=now()
    tr=load(PRETRAIN/"cut_scores.npz")
    # The eight 32-side data rows are not read. Choices use only the seven
    # historical-pair competition scores and predate native response replay.
    learn={key:int(np.argmax(np.median(v,axis=0))) for key,v in tr.items()}
    learn["exact_train_maximin"]=int(np.argmax(np.min(tr["mixed_exact"],axis=0)))
    learn.update(A=cut_id(A),B=cut_id(B))
    dump(OUT/"freeze.json",dict(cases=["n96_c01","n96_c02"],learned_from="seven n32 historical witness pairs only",
            frozen_actions=learn,scope="No new acquisition. Numerical replay of old independently constructed native 96² witnesses and a one-case material-grid sensitivity check at 192²."))
    engine=PortEngine(96)
    rows=[]
    for case in ("n96_c01","n96_c02"):
        start=now()
        x0,x1,roi,patch,card=material_pair(case)
        old=load(PRIOR/"ports"/"public"/"cases"/f"{case}.npz",("Y","radius","observation_id"))
        r=float(old["radius"])
        z=[]
        independent=[]
        for x in (x0,x1):
            st=engine.state(x)
            test,total,cost=engine.independent_response(x)
            e=float(np.linalg.norm(st["delta"]-test)/np.linalg.norm(test))
            assert e<1e-8
            z.append(st["delta"])
            independent.append(dict(relative_loaded_vs_bare=e,source_count=24,state_defect=st["state_defect"],**cost))
        delta=z[1]-z[0]
        shots=np.sqrt(cut_power(delta))/(2*r)
        rec=dict(case=case,full_observation_id=str(old["observation_id"].item()),n=96**2,frequency_MHz=400,
             radius=r,ROI_area_m2=float(roi.sum()*(2/96)**2),patch_area_m2=float(len(patch)*(2/96)**2),
             ROI_mean_difference=float((x1-x0)[roi].mean()),
             source_available_W_sum=1.,pulses=8,complex_reads=64,
             published_strict_A_patch_upper_over_r=card["upper_residual_over_r"],
             A_original_RMS_contract_ratio=[float(np.linalg.norm(project(m,A)-old["Y"])/r) for m in z],
             loaded_and_bare=independent, best_allcuts=int(np.argmax(shots)),best_allcuts_ratio=float(shots.max()),
             frozen_action_ratios={k:float(shots[i]) for k,i in learn.items()})
        arrays=dict(x0=x0,x1=x1,roi=roi,patch=patch,native_zero=z[0],native_high=z[1],native_ratios_all_cuts=shots)
        if case=="n96_c01":
            fine=[]
            for x in (x0,x1):
                # Piecewise-constant extension of each 96² domain cell to four subcells.
                hi=np.repeat(np.repeat(x.reshape(96,96),2,0),2,1).ravel()
                v,s,cost=fine_bare_response(hi,192)
                fine.append((v,cost))
            dh=fine[1][0]-fine[0][0]
            ratios=np.sqrt(cut_power(dh))/(2*r)
            rec["fine_grid"]=dict(side=192,cost=[c for z,c in fine],
                existing_A_relative_to_original_numerical_data_only=[float(np.linalg.norm(project(z,A)-old["Y"])/r) for z,c in fine],
                frozen_action_ratios={k:float(ratios[i]) for k,i in learn.items()},
                native_fine_full_matrix_relative_change=[float(np.linalg.norm(f-z)/np.linalg.norm(z)) for (f,c),z in zip(fine,z)],
                scope="Independent 192² discretization of the same piecewise constant material. The numerical cross-grid residual does not certify the 96² data model nor a continuum witness.")
            arrays.update(fine_zero=fine[0][0],fine_high=fine[1][0],fine_ratios_all_cuts=ratios)
        rows.append(dict(**rec,**elapsed(start)))
        save(OUT/f"{case}.npz",**arrays)
        dump(OUT/"report.json",dict(native=rows,counters=dict(engine.counters),**elapsed(started)))
        print(case,rec["frozen_action_ratios"],rec.get("fine_grid",{}).get("existing_A_relative_to_original_numerical_data_only"),flush=True)


if __name__=="__main__":
    with threadpool_limits(limits=1,user_api="blas"):
        main()
