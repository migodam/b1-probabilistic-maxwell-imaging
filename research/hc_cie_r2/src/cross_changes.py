"""Known-background differential imaging on native SRAP 96-to-32 data contracts."""
import argparse
from itertools import combinations
import os
from time import process_time,perf_counter
from pathlib import Path

from interaction import np, json, LAB, CFG, NORM, prepare, forward, atoms, select_pair, shortlist, threadpool_limits
from model_gap import original_sources,plan,select_ports,write_csv
from scipy.optimize import least_squares

OUT=LAB/"results/cross_changes"
C=json.loads((LAB/"configs/cross_changes.json").read_text())


def generate(workspace):
    _,data=original_sources(workspace)
    op=data.SyntheticOperators(96,C["layouts"])
    mpath=LAB/"results/model_gap"
    rows=[]
    for p,index,fa in plan(data):
        fine_m=data.random_material(op.fine_centers,index,fa)
        known=data.area_average(fine_m,96,32)
        for dd,dist in enumerate(C["changes"]):
            wall,cpu=perf_counter(),process_time()
            rng=np.random.default_rng(np.random.SeedSequence([C["seed"],index,dd]))
            pair=select_pair(rng,len(known),dd==0)
            change=np.zeros_like(known);change[pair]=C["amplitudes"]
            target_f=fine_m+np.repeat(np.repeat(change.reshape(32,32),3,0),3,1).ravel()
            newfine,rec=data.history.fine_solve(target_f,op.fine_e,op.fine_s,op.green.apply,op.green.diagonal)
            for arr in C["layouts"]:
                case0=data.synthetic_name("confirmation",p,arr)
                key=case0+"__"+dist
                with np.load(mpath/f"public/{case0}.npz") as d:
                    b,s=d["B"],d["S"]
                    fine_ref=d["cross_clean"]
                    coarse_ref=d["matched_clean"]
                with np.load(mpath/"public/coarse_green.npz") as z:G=z["G"]
                model=dict(E0=b,S=s,G=G)
                matched=forward(model,known+change)-coarse_ref
                cross=select_ports(newfine,op,arr)-fine_ref
                n=rng.normal(size=cross.shape)+1j*rng.normal(size=cross.shape);n/=NORM(n)
                np.savez_compressed(OUT/f"public/{key}.npz",b=b,s=s,known=known,
                    matched_clean=matched,matched_noisy=matched+.1*NORM(matched)*n,
                    cross_clean=cross,cross_noisy=cross+.1*NORM(cross)*n)
                np.savez_compressed(OUT/f"scoring/{key}.npz",pair=pair,steps=change[pair])
                rows.append(dict(case=key,parent=p,distance=dist,array=arr,
                    mismatch_over_change=NORM(cross-matched)/NORM(cross),mismatch_over_10pct=NORM(cross-matched)/(.1*NORM(cross)),
                    nominal_reference_gap=NORM(coarse_ref-fine_ref)/NORM(fine_ref)))
            rows[-1].update({"shared_generation_cpu_s":process_time()-cpu,"shared_generation_wall_s":perf_counter()-wall,
                             **rec,"coarse_lu":2,"coarse_rhs":20})
    write_csv(OUT/"reference_model_gap.csv",rows)


def compare(bg,model,y,method):
    ids=shortlist(bg,y)
    if method=="music_closed":pairs=[tuple(sorted(ids[:2]))]
    else:
        scans=[]
        for p in combinations(ids,2):
            p=tuple(sorted(p))
            q=np.linalg.qr(atoms(bg,p),mode="reduced")[0]
            scans.append((-NORM(q.conj().T@y.ravel()),p))
        pairs=[p for _,p in sorted(scans)[:CFG["final_pairs"]]]
    best=np.inf;cnt=0
    for p in pairs:
        ix=np.array(p);o=bg.o[:,ix];u=bg.u[ix];k=bg.k[np.ix_(ix,ix)]
        v=np.linalg.lstsq(o,y,rcond=None)[0];e=u+k@v
        den=np.sum(abs(e)**2,axis=1)
        d=np.divide(np.sum((e.conj()*v).real,axis=1),den,out=np.zeros(2),where=den>0)
        lb,ub=-bg.m[ix],20-bg.m[ix]
        d=np.clip(d,lb,ub)
        def residual(change):
            r=o@np.linalg.solve(np.eye(2)-change[:,None]*k,change[:,None]*u)-y
            return np.r_[r.real.ravel(),r.imag.ravel()]/NORM(y)
        if method=="interaction_fit":
            fit=least_squares(residual,d,bounds=(lb,ub),ftol=1e-11,gtol=1e-11,xtol=1e-11,max_nfev=80)
            d=fit.x;cnt+=fit.nfev
        err=NORM(residual(d))
        if err<best:best=err;selected=ix;steps=d
    checked=bg.m.copy();checked[selected]+=steps
    return dict(i=int(selected[0]),j=int(selected[1]),di=steps[0],dj=steps[1],nf=cnt,
                reduced=best,port=NORM(forward(model,checked)-bg.y-y)/NORM(y),
                full_lu=1,full_rhs=y.shape[1])


def solve():
    rows=[];setups=[]
    with np.load(LAB/"results/model_gap/public/coarse_green.npz") as a:g=a["G"]
    for src in sorted((OUT/"public").glob("*.npz")):
        parent,layout,distance=src.stem.split("__")
        with np.load(src) as raw:
            data=dict(G=g,E0=raw["b"],S=raw["s"])
            wall,cpu=perf_counter(),process_time();bg=prepare(data,raw["known"])
            setups.append(dict(case=src.stem,wall_s=perf_counter()-wall,cpu_s=process_time()-cpu,lu=1,
                               rhs=len(bg.m)+bg.y.shape[0]+bg.y.shape[1]))
            for condition in C["conditions"]:
                y=raw[condition]
                for mode in C["arms"]:
                    wall,cpu=perf_counter(),process_time()
                    info=compare(bg,data,y,mode)
                    rows.append(dict(case=src.stem,parent=parent,layout=layout,distance=distance,
                                     condition=condition,method=mode,**info,
                                     cpu_s=process_time()-cpu,wall_s=perf_counter()-wall))
    write_csv(OUT/"setup.csv",setups);write_csv(OUT/"predictions.csv",rows)
    (OUT/"seal.json").write_text(json.dumps(dict(rows=len(rows),truth_access=False))+"\n")


def score():
    import csv
    with (OUT/"predictions.csv").open() as f:rows=list(csv.DictReader(f))
    for row in rows:
        with np.load(OUT/("scoring/"+row["case"]+".npz")) as s:p,d=s["pair"],s["steps"]
        actual=np.zeros(1024);actual[p]=d
        predicted=np.zeros(1024);ix=np.array([int(row["i"]),int(row["j"])])
        predicted[ix]=[float(row["di"]),float(row["dj"])]
        row.update(exact_support=int(np.array_equal(p,ix)),correct_sites=len(set(p)&set(ix)),
                   change_error=NORM(predicted-actual)/NORM(actual))
    write_csv(OUT/"offline.csv",rows)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("workspace",type=Path);args=parser.parse_args()
    for p in ["public","scoring"]:(OUT/p).mkdir(exist_ok=True,parents=True)
    with threadpool_limits(1):
        generate(args.workspace.resolve());solve();score()
