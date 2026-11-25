"""Paired noise sensitivities on four frozen 32² objects and the two native objects."""
import argparse
import json
import os

for variable in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    os.environ[variable]="1"

import numpy as np
from threadpoolctl import threadpool_limits

from certified_map import (ARITH, ReceiverGeometry, active_coercivity, apply_bounds, fnorm,
                           local_inverse, local_models)
from common import FREEZE, HERE, Geometry, clock, elapsed, geom_path, pattern, save_json
from sector_map import SectorGeometry, local_certificate, margin_distances, response_reduction
from witness import forecast, certify_witness


STRESS=HERE/"public"/"stress"
OUT=HERE/"results"/"stress"
FRACTIONS=(.01,.03,.10)
B=FREEZE["model"]["background_upper"]
PANEL=("M32-D01","M32-R01","M32-P01","M32-A01","M96-R01","M96-P01")
LO,HI=FREEZE["material_classes"]["primary_core_range"]


def prepare():
    STRESS.mkdir(exist_ok=True)
    acq=json.loads((HERE/"public"/"acquisition.json").read_text())["observations"]
    manifest={}
    for oid,m in acq.items():
        if oid.split("__")[1] not in PANEL:
            continue
        old=np.load(HERE/"public"/f"{oid}.npz",allow_pickle=False)
        sealed=np.load(HERE/"offline"/f"{oid}.npz",allow_pickle=False)
        clean=sealed["clean"]
        unit=old["Y"]-clean
        unit/=fnorm(unit)
        for fraction in FRACTIONS:
            label=oid.replace("__n03",f"__n{int(fraction*100):02d}")
            radius=(fraction+FREEZE["model"]["noise_radius_arithmetic_slack_fraction"])*fnorm(clean)
            Y=old["Y"] if fraction==.03 else clean+fraction*fnorm(clean)*unit
            r=float(old["r"]) if fraction==.03 else radius
            np.savez_compressed(STRESS/f"{label}.npz",Y=Y,r=r)
            manifest[label]=dict(old_obs=oid,noise=fraction,**m)
    save_json(STRESS/"manifest.json",manifest)


def verify_old_witness(g,key,m,b,Y,r,diff,sector):
    rows=json.loads((HERE/"results"/f"witnesses_{g.side}_{m['shape']}.json").read_text())["witnesses"]
    record=rows[m["old_obs"]][str(b)]["witness"]
    if record is None:
        return None
    from pathlib import Path
    source=np.load(HERE/"results"/record["file"],allow_pickle=False)
    allowance=record["checked_residual_upper"]+diff
    if allowance<r:
        return dict(centre=record["centre"],ratio=allowance/r,method="physical witness replay via triangle",
                    file=record["file"])
    if sector is None:
        return None
    x0=source["material"]
    cells,p=source["active_cells"],source["geometric_profile"]
    ref=float(np.median(np.delete(x0,cells)))
    physical=forecast(g,x0,cells,p,Y,r,b,extra_slack=(g.side==96))
    if physical is None:
        return None
    x,j,raw=physical
    if raw>=r:
        return None
    checked=certify_witness(sector,x,cells,j,Y,r,ref)
    if checked is None:
        return None
    path=OUT/f"{key}_witness.npz"
    np.savez_compressed(path,material=x,current=j,active_cells=cells,centre=record["centre"],
                        geometric_profile=p)
    return dict(centre=record["centre"],ratio=checked["residual_ratio_upper"],
                method="fresh local full-wave witness",file=path.name)


def primary32(g,shape,observations):
    geometry=SectorGeometry(g,B)
    ports=ReceiverGeometry(g)
    operator=json.loads((HERE/"results"/"operator_bounds.json").read_text())["32"]
    lambda_=operator["full_operator_norm_upper"]
    off,p=pattern(g.side,shape)
    xy=off+g.side//2
    ga=g.restrict(xy[:,1]*g.side+xy[:,0])
    floor,theta=active_coercivity(ga,p*HI,operator["kernel_entry_error"])
    snapshot={key:np.full(g.n,-np.inf) for key in observations}
    scales=(.5,1.,.25,2.,.125,4.,.0625,.03125,.015625)
    for ids,cells,p in local_models(g,shape):
        for w in (.125,.5,1.,2.,8.):
            keys=[key for key,(_,_,b) in observations.items() if b]
            work=np.any([snapshot[key][ids]<=0 for key in keys],axis=0)
            if not np.any(work):
                break
            local=local_certificate(geometry,cells[0],p,LO,HI,w)
            if local is None:
                continue
            work_ids,work_cells=ids[work],cells[work]
            for lower in range(0,len(work_ids),128):
                q=slice(lower,lower+128)
                reduced=response_reduction(geometry,work_cells[q],local)
                for key in keys:
                    Y,r,_=observations[key]
                    a,_,_=margin_distances(reduced,Y,r,scales)
                    subset=work_ids[q]
                    snapshot[key][subset]=np.maximum(snapshot[key][subset],a)
        local=local_inverse(g,cells[0],p,LO,HI,floor,theta,operator)
        for i in range(0,len(ids),128):
            q=slice(i,i+128)
            at=cells[q]
            BA=g.B[at]
            SA=g.S[:,at].transpose(1,0,2)
            J0=local["T"]@BA
            PA=SA@local["T"]
            F0=SA@J0
            for key,(Y,r,b) in observations.items():
                if b:
                    continue
                _,proof=apply_bounds(ports,local,Y,r,BA,SA,J0,PA,F0,0.,lambda_)
                snapshot[key][ids[q]]=np.maximum(*proof)
    return geometry,snapshot


def run(side,shape,fresh_native_witness=False):
    OUT.mkdir(parents=True,exist_ok=True)
    started=clock()
    g=Geometry(geom_path(side))
    manifest=json.loads((STRESS/"manifest.json").read_text())
    selected={key:m for key,m in manifest.items() if m["side"]==side and m["shape"]==shape}
    data={}
    for key,m in selected.items():
        op=np.load(STRESS/f"{key}.npz",allow_pickle=False)
        b=B if m["scene"]=="weak" else 0.
        data[key]=(op["Y"],float(op["r"]),b)
    geometry,screen=(primary32(g,shape,data) if side==32 else
                     (SectorGeometry(g,B) if fresh_native_witness else None,{}))
    online=elapsed(started)
    rows={}
    for key,m in selected.items():
        Y,r,b=data[key]
        old=np.load(HERE/"public"/f"{m['old_obs']}.npz",allow_pickle=False)
        delta=float(fnorm(Y-old["Y"]))
        diff=delta+ARITH*(delta+fnorm(Y)+fnorm(old["Y"])) if delta else 0.
        original=np.load(HERE/"results"/f"{m['old_obs']}__b{b:g}__score_map.npz",allow_pickle=False)
        bounds=np.maximum(original["sector"],original["neumann"]).ravel()+float(old["r"])
        if side==96:
            margin=bounds-diff-r-ARITH*(1.+r+diff)
            method="Lipschitz transfer of existing 3-percent distance lower"
        else:
            margin=screen[key]
            method="fresh strict all-centre replay"
        excluded=margin>0
        status=np.where(excluded,0,2).astype(np.uint8)
        witness=verify_old_witness(g,key,m,b,Y,r,diff,geometry)
        if witness:
            x,y=witness["centre"]
            assert not excluded[y*g.side+x]
            status[y*g.side+x]=1
        truth=np.load(HERE/"offline"/f"{m['old_obs']}.npz",allow_pickle=False)
        x,y=truth["centre"]
        assert not excluded[y*g.side+x]
        rows[key]=dict(n_centres=g.n,excluded=int(excluded.sum()),witness=witness,
                       unresolved=int((status==2).sum()),actual_scene=m["scene"],
                       inference_background=b,noise=m["noise"],true_centre_retained=True,method=method)
        np.savez_compressed(OUT/f"{key}.npz",status=status.reshape(side,side),
                            strict_margin=margin.reshape(side,side))
    save_json(OUT/f"stress_{side}_{shape}.json",
              dict(rows=rows,online_before_replay=online,all=elapsed(started)))
    print(side,shape,{k.split("__")[2]+"-"+k.split("__")[-1]:[v["excluded"],bool(v["witness"])]
                      for k,v in rows.items()},elapsed(started),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--prepare",action="store_true")
    parser.add_argument("side",type=int,nargs="?")
    parser.add_argument("shape",nargs="?")
    parser.add_argument("--fresh-native-witness",action="store_true")
    args=parser.parse_args()
    with threadpool_limits(limits=1,user_api="blas"):
        if args.prepare:
            prepare()
        else:
            run(args.side,args.shape,args.fresh_native_witness)
