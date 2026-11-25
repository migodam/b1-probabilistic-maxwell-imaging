"""Offline-only comparison of sealed map masks and synthetic truth."""
import csv
from collections import defaultdict
import json
import os

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np

from common import FREEZE, HERE, save_json


RESULTS=HERE/"results"
B=FREEZE["model"]["background_upper"]
SHAPES=("solid","ring","double","austria")
ACQ=json.loads((HERE/"public"/"acquisition.json").read_text())["observations"]


def read(path):
    return np.load(RESULTS/path,allow_pickle=False)


def run():
    records=[]
    cond={}
    primary_weak_witness={}
    for key,meta in ACQ.items():
        n,shape=meta["side"],meta["shape"]
        yvac=meta["scene"]=="vacuum"
        path=RESULTS/f"witnesses_{n}_{shape}.json"
        witnesses=json.loads(path.read_text())["witnesses"][key] if path.exists() else {}
        strict_file=RESULTS/f"{key}__b{B:g}__primary_sector.npz"
        strict=(read(strict_file.name) if strict_file.exists() else None)
        with np.load(HERE/"offline"/f"{key}.npz",allow_pickle=False) as z:
            truth,xy,cells,p = z["truth"],z["centre"],z["cells"],z["profile"]
            clean,born=z["clean"],z["born"]
        with np.load(HERE/"public"/f"{key}.npz",allow_pickle=False) as obs:
            Y,r=obs["Y"],float(obs["r"])
        node=int(xy[1]*n+xy[0])
        born_error=np.linalg.norm(born-clean)/np.linalg.norm(clean)
        truth_interior=truth[cells]/p
        assert np.all(truth_interior>=3.75) and np.all(truth_interior<=4.25)
        active=np.zeros(len(truth),dtype=bool)
        active[cells]=True
        assert np.linalg.norm(Y-clean)<r
        assert np.all(truth[~active]<=B)
        for b in ((0.,B) if yvac else (B,)):
            controls=read(f"{key}__b{b:g}__primary.npz")
            neu=controls["margin_neumann"]
            if strict is None:
                sector=np.full(n*n,-np.inf)
                rank=np.full(n*n,-1)
            else:
                sector,rank=strict["margin_sector"],strict["negative_index"]
            union=(neu>0)|(controls["margin_coercive"]>0)|(sector>0)
            assert not union[node], (key,b,"false exclusion")
            state=np.where(union,0,2).astype(np.uint8)
            geometry_valid=bool(np.all(truth[~active]<=b))
            assert geometry_valid
            rank_bp=np.argsort(-controls["osm"].ravel())
            witness=witnesses.get(str(b),{}).get("witness")
            wm={}
            if witness:
                wi=read(witness["file"])
                wmat,wx,wc,wp=wi["material"],wi["centre"],wi["active_cells"],wi["geometric_profile"]
                k=int(wx[1]*n+wx[0])
                assert not union[k]
                state[k]=1
                outside=np.delete(wmat,wc)
                assert np.all(outside>=0) and np.all(outside<=b)
                assert np.all(wmat[wc]>=3.75*wp) and np.all(wmat[wc]<=4.25*wp)
                if not yvac:
                    primary_weak_witness[key]=(k,witness["file"])
                wm=dict(witness_position=[int(wx[0]),int(wx[1])],
                        witness_at_true_centre=bool(k==node),
                        witness_residual_ratio=witness["residual_ratio_upper"],
                        witness_rel_material=float(np.linalg.norm(wmat-truth)/np.linalg.norm(truth)),
                        witness_rel_interior=float(np.linalg.norm((wmat-truth)[active])/np.linalg.norm(truth[active])),
                        witness_exterior_rmse=float(np.mean(((wmat-truth)[~active])**2)**.5))
            peak_idx=int(rank_bp[0])
            pos=(peak_idx%n,peak_idx//n)
            distance=float(np.linalg.norm(np.asarray(pos)-xy))
            row=dict(obs_id=key,grid=n,shape=shape,actual_scene=meta["scene"],background_upper=b,
                     complete_sector_file=bool(strict is not None),n_centres=n*n,
                     excluded=int(np.sum(state==0)),witness=int(np.sum(state==1)),
                     unresolved=int(np.sum(state==2)),true_centre=[int(xy[0]),int(xy[1])],
                     true_centre_retained=int(state[node]!=0),true_family_interval_verified=geometry_valid,
                     sector_excluded=int(np.sum(sector>0)),
                     sector_best_indefinite=int(np.sum((sector>0)&(rank>0))),
                     directional_neumann_excluded=int(np.sum(neu>0)),
                     directional_plus_coercivity_excluded=int(np.sum((neu>0)|(controls["margin_coercive"]>0))),
                     osm_top1_xy=list(pos),osm_location_error_cells=distance,
                     osm_true_rank_1_based=int(np.flatnonzero(rank_bp==node)[0]+1),
                     born_error=float(born_error),**wm)
            records.append(row)
            np.savez_compressed(RESULTS/f"{key}__b{b:g}__score_map.npz",
                                status=state.reshape(n,n),true_centre=xy,
                                sector=sector.reshape(n,n),neumann=neu.reshape(n,n))
        if not yvac:
            cond[key]=(node,xy,n,shape)
    lib={}
    for key,(node,xy,n,shape) in cond.items():
        if n!=32:
            continue
        arm={}
        for new in ("library_all","library"):
            bundle=[RESULTS/f"{s}__{key}__b{B:g}__{new}_sector.npz" for s in SHAPES]
            if not all(p.exists() for p in bundle):
                continue
            stacks=np.stack([read(path.name)["margin_sector"] for path in bundle])
            dead=(stacks>0)
            union=dead.all(axis=0)
            assert not union[node]
            possible=np.flatnonzero(~dead.all(axis=1))
            ternary=np.where(union,0,2)
            witness_file=None
            if len(possible)==1 and key in primary_weak_witness:
                idx,witness_file=primary_weak_witness[key]
                assert not dead[possible[0],idx]
                ternary[idx]=1
            arm[new]=dict(excluded=int(union.sum()),witness=int(np.sum(ternary==1)),
                          unresolved=int(np.sum(ternary==2)),witness_file=witness_file,
                          surviving_shapes=[SHAPES[i] for i in possible],
                          wrong_shapes_fully_excluded=[s for i,s in enumerate(SHAPES)
                                                       if s!=shape and dead[i].all()],
                          per_shape_excluded=dict(zip(SHAPES,map(int,dead.sum(axis=1)),strict=True)))
            np.savez_compressed(RESULTS/f"{key}__b{B:g}__{new}_union.npz",
                                status=ternary.reshape(n,n),shapes=SHAPES,
                                shape_margins=stacks.reshape(len(SHAPES),n,n),true_centre=xy)
        if arm:
            lib[key]=arm
    broad={}
    for key,(node,xy,n,shape) in cond.items():
        path=RESULTS/f"{key}__b{B:g}__broad_sector.npz"
        if path.exists():
            a=read(path.name)
            marked=a["margin_sector"]>0
            assert not marked[node]
            pp=RESULTS/f"{key}__b{B:g}__broad_psd_sector.npz"
            part=RESULTS/f"{key}__b{B:g}__broad_psd_interior_sector.npz"
            npsd_partial=0 if part.exists() else None
            if part.exists():
                npsd_partial=int(np.sum(read(part.name)["margin_sector"]>0))
            witness=int(key in primary_weak_witness)
            if witness:
                idx,witness_file=primary_weak_witness[key]
                assert not marked[idx]
            broad[key]=dict(denominator=n*n,excluded=int(marked.sum()),witness=witness,
                            unresolved=n*n-int(marked.sum())-witness,
                            indefinite_best=int(((a["negative_index"]>0)&marked).sum()),
                            psd_excluded=int(np.sum(read(pp.name)["margin_sector"]>0)) if pp.exists() else None,
                            psd_interior_excluded=npsd_partial,
                            psd_interior_denominator=5184 if part.exists() else None,
                            state_negative_indices_at_exclusions=sorted(set(
                                map(int,a["negative_index"][marked]))),
                            verified_margin_available=int(np.isfinite(a["margin_sector"]).sum()),
                            true_centre_retained=True)
    groups=defaultdict(list)
    for row in records:
        groups[f"{row['grid']}__actual_{row['actual_scene']}__b{row['background_upper']:g}"].append(row)
    aggregates={}
    for name,data in groups.items():
        values=lambda k:[r[k] for r in data]
        aggreg=dict(independent_objects=len(data),all_charged_centres=sum(values("n_centres")),
                    excluded=sum(values("excluded")),witnesses=sum(values("witness")),
                    unresolved=sum(values("unresolved")),true_centre_retained=sum(values("true_centre_retained")),
                    median_excluded=float(np.median(values("excluded"))),
                    excluded_range=[min(values("excluded")),max(values("excluded"))],
                    neumann_excluded=sum(values("directional_neumann_excluded")),
                    sector_excluded=sum(values("sector_excluded")),
                    best_proof_indefinite=sum(values("sector_best_indefinite")),
                    osm_top1_correct=sum(v==1 for v in values("osm_true_rank_1_based")),
                    osm_top5_covers=sum(v<=5 for v in values("osm_true_rank_1_based")),
                    osm_median_error_cells=float(np.median(values("osm_location_error_cells"))),
                    born_error_range=[min(values("born_error")),max(values("born_error"))],
                    sector_files_complete=all(values("complete_sector_file")))
        witnesses=[r for r in data if r["witness"]]
        if witnesses:
            aggreg["witness_relative_residual_max"]=max(r["witness_residual_ratio"] for r in witnesses)
            aggreg["witness_true_centre"]=sum(r["witness_at_true_centre"] for r in witnesses)
            aggreg["witness_median_relative_l2"]=float(np.median([r["witness_rel_material"] for r in witnesses]))
            aggreg["witness_median_internal_relative_l2"]=float(np.median([r["witness_rel_interior"] for r in witnesses]))
            aggreg["witness_median_exterior_rmse"]=float(np.median([r["witness_exterior_rmse"] for r in witnesses]))
        aggregates[name]=aggreg
    save_json(RESULTS/"summary.json",dict(aggregates=aggregates,rows=records,library=lib,broad=broad))
    keys=list(dict.fromkeys(k for row in records for k in row))
    with (RESULTS/"centres_summary.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=keys)
        writer.writeheader()
        writer.writerows(records)
    for filename,rows in (
        ("broad_centres.csv", [dict(obs_id=k,**v) for k,v in broad.items()]),
        ("library_per_object.csv", [dict(obs_id=k,**entry["library_all"],
                                       **entry["library_all"]["per_shape_excluded"])
                                for k,entry in lib.items() if "library_all" in entry]),
    ):
        h=list(dict.fromkeys(k for r in rows for k in r))
        with (RESULTS/filename).open("w",newline="") as stream:
            w=csv.DictWriter(stream,fieldnames=h)
            w.writeheader()
            w.writerows(rows)
    print(json.dumps(aggregates,indent=2),flush=True)


if __name__=="__main__":
    run()
