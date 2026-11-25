"""Material intervals at the selected position and one-site peak contracts."""
import json
import os

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
from threadpoolctl import threadpool_limits

from certified_map import ARITH, ReceiverGeometry, add_bound, fnorm
from common import FREEZE, HERE, Geometry, clock, elapsed, geom_path, pattern, place, save_json
from sector_map import SectorGeometry, local_certificate, margin_distances, response_reduction


WEIGHTS = (.03125,.0625,.125,.25,.5,1.,2.,8.)
SCALES = tuple(2.**np.arange(-12,6))


def test(sector, cells, p, lo, hi, Y, r):
    best, rank = -np.inf, None
    for w in WEIGHTS:
        local = local_certificate(sector,cells,p,lo,hi,w)
        if local is None:
            continue
        if sector.g.side==96 and local[-1]>=15:
            continue
        reduction=response_reduction(sector,cells[None,:],local)
        value,nu,_=margin_distances(reduction,Y,r,SCALES)
        if value[0] > best:
            best=float(value[0])
            rank=int(nu[0])
        if best > 0:
            break
    return (None if not np.isfinite(best) else best),rank


def status(margin):
    return "EXCLUDED" if margin is not None and margin>0 else "UNRESOLVED"


def global_peak(ports,Y,r,Lambda):
    nrm=fnorm(Y)
    W=(Y/nrm)[None,:,:]
    direction=ports.sensor_direction(W)[0]
    first=ports.background_linear(W)[0]
    trials=[]
    for p in (0.015625,0.03125,0.046875,0.0625,0.078125):
        uv=1-p*Lambda
        need=add_bound(r*fnorm(W))
        if uv<=0:
            break
        error=p**2*direction*Lambda*ports.bnorm/uv
        margin=(W.conj()*Y).sum().real - p*first - error - need - ARITH*(1+nrm+error)
        trials.append(dict(cap=p,status=status(margin),distance_margin=float(margin)))
    successes=[v["cap"] for v in trials if v["status"]=="EXCLUDED"]
    return dict(excluded_uniform_global_caps=trials,max_excluded=max(successes,default=None))


def run(side,shape):
    started=clock()
    g=Geometry(geom_path(side))
    ports=ReceiverGeometry(g)
    operator=json.loads((HERE/"results"/"operator_bounds.json").read_text())[str(side)]
    Lambda=operator["full_operator_norm_upper"]
    b=FREEZE["model"]["background_upper"]
    geometry=SectorGeometry(g,b)
    witnesses=json.loads((HERE/"results"/f"witnesses_{side}_{shape}.json").read_text())["witnesses"]
    output={}
    for oid,row in witnesses.items():
        one=clock()
        if "__weak__" not in oid:
            continue
        registration=row[str(b)]["witness"]
        data=np.load(HERE/"public"/f"{oid}.npz",allow_pickle=False)
        Y,r=data["Y"],float(data["r"])
        if registration is not None:
            center=registration["centre"]
            mat=np.load(HERE/"results"/registration["file"],allow_pickle=False)["material"]
        else:
            other=np.load(HERE/"results"/f"{oid}__b{b:g}__primary.npz",allow_pickle=False)
            center=g.xy[np.argmin(other["reference_error"])].tolist()
            mat=None
        offset,p=pattern(side,shape)
        a,p=place(side,center,offset,p)
        entries=[]
        for lo,hi in FREEZE["material_classes"]["class_partition"]:
            value,rank=test(geometry,a,p,lo,hi,Y,r)
            state=status(value)
            if mat is not None and np.all(mat[a]>=lo*p) and np.all(mat[a]<=hi*p):
                assert state!="EXCLUDED"
                state="WITNESS_FOUND"
            entries.append(dict(interval=[lo,hi],state=state,margin=value,negative_index=rank))
        peak_lower=[]
        for cap in (2.5,3.,3.25,3.5,3.75,4.0):
            hi=np.minimum(6.,cap/p)
            value,rank=test(geometry,a,p,2.,hi,Y,r)
            peak_lower.append(dict(cap=cap,state=status(value),margin=value,negative_index=rank))
        peak_upper=[]
        for floor in (4.25,4.5,5.,5.5,6.):
            trials=0
            if mat is not None and np.max(mat[a])>=floor:
                state="WITNESS_FOUND"
            elif side==96:
                state="UNRESOLVED"
            else:
                state="EXCLUDED"
                for j in np.flatnonzero(6*p>=floor):
                    low=np.full(len(p),2.)
                    low[j]=floor/p[j]
                    value,rank=test(geometry,a,p,low,6.,Y,r)
                    trials+=1
                    if status(value)!="EXCLUDED":
                        state="UNRESOLVED"
                        break
            peak_upper.append(dict(floor=floor,state=state,site_boxes_tested=trials))
        rp=[v for v in entries if v["state"]!="EXCLUDED"]
        shape_interval=(min(v["interval"][0] for v in rp),max(v["interval"][1] for v in rp)) if rp else None
        exlo=max([2.]+[v["cap"] for v in peak_lower if v["state"]=="EXCLUDED"])
        exhi=min([6.]+[v["floor"] for v in peak_upper if v["state"]=="EXCLUDED"])
        whole_peak=global_peak(ports,Y,r,Lambda)
        high=np.full(g.n,2,dtype=np.uint8)
        if mat is not None:
            high[mat>=FREEZE["material_classes"]["independent_peak_threshold"]]=1
        np.savez_compressed(HERE/"results"/f"{oid}__global_peak_map.npz",state=high)
        output[oid]=dict(centre=center,centre_has_physical_witness=registration is not None,
                        amplitude_bins=entries,surviving_bin_hull=shape_interval,
                        local_independent_peak_low_boxes=peak_lower,
                        local_independent_peak_high_union=peak_upper,
                        local_peak_outer_interval=[exlo,exhi],
                        local_peak_witness=None if mat is None else float(np.max(mat[a])),
                        global_peak=whole_peak,global_single_cell_high_witnesses=int(np.sum(high==1)),
                        **elapsed(one))
        print("contrast",oid,shape_interval,[exlo,exhi],
              "witness x_peak",output[oid]["local_peak_witness"],
              "global lower",whole_peak["max_excluded"],elapsed(one),flush=True)
    save_json(HERE/"results"/f"contrast_{side}_{shape}.json",dict(results=output,elapsed=elapsed(started)))


if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument("side",type=int)
    p.add_argument("shape")
    a=p.parse_args()
    with threadpool_limits(limits=1,user_api="blas"):
        run(a.side,a.shape)
