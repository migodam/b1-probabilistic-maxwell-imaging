"""Two-site changes in a known reciprocal Maxwell background; study entry point."""
import argparse
import csv
import json
import os
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from time import perf_counter, process_time

for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","VECLIB_MAXIMUM_THREADS","MKL_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.optimize import least_squares
from threadpoolctl import threadpool_limits

LAB = Path(__file__).resolve().parents[1]
OUT = LAB / "results/interaction"
CFG = json.loads((LAB / "configs/interaction.json").read_text())
NORM = np.linalg.norm


def historical():
    override = os.environ.get("HC_CIE_R2_INPUT_ROOT") or os.environ.get("MAXWELL_WORKSPACE")
    roots = [Path(override).expanduser().resolve()] if override else [LAB, *LAB.parents]
    for root in roots:
        for relative in ("research/research_max_struct_r13",
                         "research/delegated/max_struct_r1/r13_worktree/research/research_max_struct_r13"):
            candidate = root / relative
            if (candidate / "data_public/fresh48/r13_fresh_0000.npz").is_file():
                return candidate
    raise FileNotFoundError("Extract the replay inputs into a checkout or set HC_CIE_R2_INPUT_ROOT")


def public_background(number):
    key = f"r13_fresh_{number:04d}"
    study = historical()
    with np.load(study / f"data_public/fresh48/{key}.npz") as p:
        data = {k: p[k] for k in ("G","S","E0","lower","upper","delta","Y")}
    with np.load(study / f"RESULTS/confirmation_core48/{key}__cie/primary.npz") as p:
        m = p["m"].reshape(-1)
    return data, m


@dataclass
class Background:
    m: np.ndarray
    u: np.ndarray
    o: np.ndarray
    k: np.ndarray
    y: np.ndarray


def forward(data, m):
    g,b,s = (data[k] for k in ("G","E0","S"))
    a = np.eye(m.size) - m[:,None]*g
    j = lu_solve(lu_factor(a), m[:,None]*b)
    return s @ j


def prepare(data, m):
    g,b,s = (data[k] for k in ("G","E0","S"))
    lu = lu_factor(np.eye(m.size) - m[:,None]*g)
    j = lu_solve(lu, m[:,None]*b)
    return Background(m, b + g@j, lu_solve(lu, s.T, trans=1).T,
                      lu_solve(lu, g.T, trans=1).T, s@j)


def write_csv(name, rows):
    with (OUT / name).open("w",newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def select_pair(rng,n,adjacent):
    side = int(np.sqrt(n))
    i = int(rng.integers(n))
    yy,xx = np.indices((side,side))
    distance = np.abs(yy.ravel()-i//side) + np.abs(xx.ravel()-i%side)
    allowed = np.flatnonzero(distance == 1) if adjacent else np.flatnonzero(distance >= side//2)
    return np.sort([i,rng.choice(allowed)])


def generate():
    (OUT / "public").mkdir(parents=True,exist_ok=True)
    (OUT / "scoring").mkdir(exist_ok=True)
    rows = []
    for number in CFG["parents"]:
        wall,cpu = perf_counter(),process_time()
        data, m = public_background(number)
        y_ref = forward(data,m)
        for q,distance in enumerate(CFG["distances"]):
            rng = np.random.default_rng(np.random.SeedSequence([CFG["seed"],number,q]))
            pair = select_pair(rng,len(m),q==0)
            mid = (float(data["upper"])+float(data["lower"]))/2
            step = np.array(CFG["change_amplitude"]) * np.where(m[pair]<mid,1.,-1.)
            target = m.copy()
            target[pair] += step
            y1 = forward(data,target)
            delta = y1-y_ref
            unit = rng.standard_normal(delta.shape) + 1j*rng.standard_normal(delta.shape)
            unit /= NORM(unit)
            keys = {}
            for noise,norm in zip(CFG["noise"],[0.,.1*NORM(delta),.1*NORM(y1)]):
                keys["data_"+noise] = delta+norm*unit
                keys["radius_"+noise] = norm
            name = f"r13_{number:04d}_{distance}"
            np.savez_compressed(OUT/"public"/(name+".npz"),reference=y_ref,**keys)
            np.savez_compressed(OUT/"scoring"/(name+".npz"),pair=pair,step=step,clean_difference=delta)
        rows.append(dict(parent=number,full_lu=3,rhs=3*y_ref.shape[1],wall_s=perf_counter()-wall,cpu_s=process_time()-cpu))
    write_csv("generation_cost.csv",rows)


def atoms(bg,pair):
    i,j = pair
    a = np.outer(bg.o[:,i],bg.u[i])
    b = np.outer(bg.o[:,j],bg.u[j])
    c = np.outer(bg.o[:,i],bg.u[j])+np.outer(bg.o[:,j],bg.u[i])
    return np.stack([a.ravel(),b.ravel(),c.ravel()],axis=1)


def shortlist(bg,y):
    left,_,right = np.linalg.svd(y,full_matrices=False)
    obs = np.sum(abs(left[:,:2].conj().T@bg.o)**2,axis=0) / np.sum(abs(bg.o)**2,axis=0)
    src = np.sum(abs(bg.u@right[:2].conj().T)**2,axis=1) / np.sum(abs(bg.u)**2,axis=1)
    return np.argsort(-(obs+src),kind="stable")[:CFG["shortlist"]]


def fit(bg, ix, y, lower, upper):
    pair = np.array(ix)
    k = bg.k[np.ix_(pair,pair)]
    a = atoms(bg,pair)
    yy = y.ravel()
    def residual(d):
        t = np.linalg.solve(np.eye(2)-d[:,None]*k,np.diag(d))
        z = a @ np.array([t[0,0],t[1,1],(t[0,1]+t[1,0])/2])-yy
        return np.r_[z.real,z.imag]
    bounds = (lower-bg.m[pair],upper-bg.m[pair])
    # Data normalization makes the stopping convention independent of port units.
    scale = NORM(yy)
    result = least_squares(lambda d: residual(d)/scale,[0.,0.],bounds=bounds,
                           ftol=1e-11,xtol=1e-11,gtol=1e-11,max_nfev=80)
    return NORM(residual(result.x)),result.x,result.nfev,result.njev


def predict(bg, data, y, mode):
    startw,startc = perf_counter(), process_time()
    ranking = shortlist(bg,y)
    tested = 0
    if mode == "music":
        pairs = [tuple(sorted(ranking[:2]))]
    else:
        screened = []
        for pair in combinations(ranking,2):
            pair = tuple(sorted(pair))
            basis = atoms(bg,pair)[:,:2 if mode=="diagonal" else 3]
            q = np.linalg.qr(basis,mode="reduced")[0]
            projected = q.conj().T @ y.ravel()
            screened.append((-np.vdot(projected,projected).real,pair))
        pairs = [p for _,p in sorted(screened)[:CFG["final_pairs"]]]
        tested = len(screened)
    best_value = float("inf")
    nfev = njev = 0
    for p in pairs:
        value,d,n,j=fit(bg,p,y,float(data["lower"]),float(data["upper"]))
        nfev += n
        njev += j
        if value < best_value:
            best_value,best_pair,best_d=value,p,d
    mat = bg.m.copy()
    mat[list(best_pair)] += best_d
    audited = forward(data,mat)-bg.y
    return dict(i=int(best_pair[0]),j=int(best_pair[1]),di=float(best_d[0]),dj=float(best_d[1]),
                residual=float(NORM(audited-y)),relative_data=float(NORM(audited-y)/NORM(y)),
                reduced_audit_difference=float(abs(best_value-NORM(audited-y))),screened=tested,
                legal_fits=len(pairs),nfev=int(nfev),njev=int(njev),full_audit_lu=1,
                full_audit_rhs=y.shape[1],cpu_s=process_time()-startc,wall_s=perf_counter()-startw)


def solve():
    results = []
    setups = []
    for number in CFG["parents"]:
        wall,cpu = perf_counter(),process_time()
        data,m = public_background(number)
        bg = prepare(data,m)
        setups.append(dict(parent=number,lu=1,rhs=m.size+bg.y.shape[0]+bg.y.shape[1],
                           cpu_s=process_time()-cpu,wall_s=perf_counter()-wall))
        for distance in CFG["distances"]:
            name = f"r13_{number:04d}_{distance}"
            with np.load(OUT/"public"/(name+".npz")) as obs:
                for noise in CFG["noise"]:
                    y = obs["data_"+noise]
                    for mode in CFG["methods"]:
                        p = predict(bg,data,y,mode)
                        results.append(dict(case=name,parent=number,distance=distance,noise=noise,method=mode,**p))
    write_csv("setup_cost.csv",setups)
    write_csv("predictions.csv",results)
    (OUT/"prediction_seal.json").write_text(json.dumps(dict(predictions=len(results),parents=len(setups),
        scoring_directory_used=False,selection="full differential data residual on legal fitted pair"),indent=2)+"\n")


def offline():
    if not (OUT/"prediction_seal.json").is_file():
        raise RuntimeError("Run solve before scoring")
    with (OUT/"predictions.csv").open() as f:
        predictions = list(csv.DictReader(f))
    rows=[]
    interactions=[]
    for number in CFG["parents"]:
        data,m = public_background(number)
        bg = prepare(data,m)
        for distance in CFG["distances"]:
            case = f"r13_{number:04d}_{distance}"
            with np.load(OUT/"scoring"/(case+".npz")) as truth:
                pair,step = truth["pair"],truth["step"]
            ground_change = np.zeros_like(m)
            ground_change[pair] = step
            a = atoms(bg,pair)
            psi = a[:,2]-a[:,:2]@np.linalg.lstsq(a[:,:2],a[:,2],rcond=None)[0]
            eps = float(NORM(psi))
            k = bg.k[np.ix_(pair,pair)]
            t = np.linalg.solve(np.eye(2)-step[:,None]*k,np.diag(step))
            c = t[0,1]
            field_area = NORM(np.outer(bg.u[pair[0]],bg.u[pair[1]])-np.outer(bg.u[pair[1]],bg.u[pair[0]]))/np.sqrt(2.)
            with np.load(OUT/"public"/(case+".npz")) as obs:
                for noise in CFG["noise"]:
                    y = obs["data_"+noise]
                    radius = float(obs["radius_"+noise])
                    c_hat = np.vdot(psi,y.ravel())/eps**2
                    interactions.append(dict(case=case,parent=number,distance=distance,noise=noise,
                        i=pair[0],j=pair[1],di=step[0],dj=step[1],psi_perp=eps,
                        wedge=float(field_area),coupling_abs=abs(c),curvature_per_cell=abs(c)*field_area,
                        visible_joint_norm=eps*abs(c),noise_radius=radius,
                        relative_cross=eps*abs(c)/NORM(obs["data_clean"]),
                        coefficient_error=abs(c_hat-c),conditional_error_bound=radius/eps))
            for row in (r for r in predictions if r["case"]==case):
                selected = np.array([int(row["i"]),int(row["j"])])
                dif = np.zeros_like(m)
                dif[selected] = [float(row["di"]),float(row["dj"])]
                rows.append({**row,"correct_sites":len(set(selected)&set(pair)),
                    "exact_support":int(np.array_equal(pair,selected)),
                    "change_relative_l2":NORM(dif-ground_change)/NORM(ground_change),
                    "material_relative_l2":NORM(dif-ground_change)/NORM(m+ground_change),
                    "foreground_rmse":NORM((dif-ground_change)[pair])/np.sqrt(2)})
    write_csv("offline.csv",rows)
    write_csv("port_visibility.csv",interactions)
    folds=[]
    # Grouped reproducibility of the frozen algorithm, with no labels used to refit the solver.
    for n in CFG["noise"]:
        for method in CFG["methods"]:
            for f in range(4):
                subset=[r for r in rows if r["noise"]==n and r["method"]==method and
                        CFG["parents"].index(int(r["parent"]))%4==f]
                folds.append(dict(fold=f,noise=n,method=method,parents=len(subset)//2,
                    support_correct=sum(int(r["exact_support"]) for r in subset),conditions=len(subset),
                    mean_change_error=float(np.mean([r["change_relative_l2"] for r in subset]))))
    write_csv("folds.csv",folds)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("phase",choices=["generate","solve","score","all"])
    phase = p.parse_args().phase
    with threadpool_limits(1):
        if phase in ("generate","all"): generate()
        if phase in ("solve","all"): solve()
        if phase in ("score","all"): offline()
