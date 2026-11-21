"""Geometry-only receiver selection and matched material-cap checks."""

import argparse
from collections import Counter
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

from geometry import HERE, SOURCE, U_GRID, norm, read_json

import numpy as np
from scipy.linalg import eigh, svd
from scipy.sparse.linalg import LinearOperator, eigsh
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(SOURCE/"sector"))
import sector


SEED = 20261122
DEG_TOL = 1e-9
DEPENDENCE_TOL = 1e-10
CANDIDATES = ("A_uniform", "B_adjacent", "C_random", "D_hidden_mode")
RING_AND_AUSTRIA = ("r10_twod_0001","austria_eps4_400MHz")
EPS = np.finfo(float).eps


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n",encoding="utf-8")


def get_archives():
    with np.load(SOURCE/"public_vie_inputs.npz", allow_pickle=False) as infile:
        keys = [k.removesuffix("__G_key") for k in infile.files if k.endswith("__G_key")]
        greens = {key: infile[key] for key in sorted({str(infile[x+"__G_key"].item()) for x in keys})}
        arrays = {x: {"G": greens[str(infile[x+"__G_key"].item())], "E0": infile[x+"__E0"],
                      "S": infile[x+"__S"]} for x in keys}
    geometry = read_json(SOURCE/"reflection/results/SELECTION.json")
    layout = {case["id"]: case["layout"] for case in geometry["view_manifest"]}
    receivers = {name: np.array(geometry["geometries"][layout[name]]["receivers"]) for name in keys}
    return arrays, receivers


def selected_space(H, s, costs):
    _, singular, z = svd(s, full_matrices=False, check_finite=False)
    rank = int(np.sum(singular > max(s.shape)*EPS*singular[0]))
    rows = z[:rank].conj().T

    def projection(v):
        return v-rows@(rows.conj().T@v)

    def mv(v):
        return projection(H@projection(v))

    rng = np.random.default_rng(18333+len(H))
    seed = projection(rng.standard_normal(len(H))+1j*rng.standard_normal(len(H)))
    count = 8
    op = LinearOperator(H.shape, matvec=mv, dtype=complex)
    while True:
        values, vecs = eigsh(op, k=count, ncv=max(24, 2*count+4), which="LA", v0=seed, tol=1e-12)
        order = np.argsort(values)
        values, vecs = values[order], vecs[:, order]
        top = values[-1]
        active = np.abs(values-top) <= DEG_TOL*max(1.,abs(top))
        costs["projected_eigensolver_calls"] += 1
        if not np.all(active):
            break
        count *= 2
        if count > 64:
            raise ArithmeticError("Top eigenspace tie exceeds the configured basis window")
    modes = np.linalg.qr(projection(vecs[:,active]))[0]
    v = projection(vecs[:,-1])
    v /= norm(v)
    top = float(np.vdot(v,H@v).real)
    costs["receiver_svd"] += 1
    meta = {"hidden_u_star": 1/top if top > 0 else None, "lambda_max": top,
            "near_null_receiver_norm": norm(s@v), "projection_eigen_residual": norm(mv(v)-top*v),
            "top_tied_eigenspace_size": modes.shape[1], "numerical_receiver_rank": rank,
            "sigma_min": float(singular[-1]), "sigma_ratio": float(singular[-1]/singular[0])}
    return meta, rows, modes, v


def arc_coordinates(receiver):
    spread = np.linalg.svd(receiver-receiver.mean(axis=0),compute_uv=False)
    if spread[-1] <= 1e-11*spread[0]:
        direction = receiver[-1] - receiver[0]
        coordinate = receiver @ (direction / norm(direction))
        return coordinate, "longitudinal line coordinate in physical units"
    angles = np.unwrap(np.arctan2(receiver[:,1],receiver[:,0]))
    return np.rad2deg(angles),"azimuth in degrees, unwrapped in original receiver order"


def uniform_choice(remaining, coords):
    values = coords[remaining]
    boundaries = np.linspace(values.min(),values.max(),5)
    remaining = set(remaining)
    picks, detail = [], []
    for i, (a,b) in enumerate(zip(boundaries[:-1],boundaries[1:],strict=True)):
        target = (a+b)/2
        pool = [j for j in remaining if a <= coords[j] <= b]
        fallback = not pool
        j = min(pool or remaining,key=lambda idx: (abs(coords[idx]-target),idx))
        detail.append({"quarter": i+1,"edges":[float(a),float(b)],"target":float(target),
                       "selected":int(j),"nearest_unused_outside_quarter":fallback})
        picks.append(j)
        remaining.remove(j)
    return picks,detail


def design(S, H, coords, time_count):
    base = list(range(len(S)//4))
    remaining = list(range(len(base),len(S)))
    abscissa, units = arc_coordinates(coords)
    uniform, traceA = uniform_choice(remaining,abscissa)
    random = np.random.default_rng(SEED).choice(remaining,4,replace=False).tolist()
    selected = {"A_uniform": uniform, "B_adjacent": remaining[:4], "C_random": random}
    traceD, intermediate = [], []
    base_meta,_,_,_ = selected_space(H,S[base],time_count)
    used = base.copy()
    for _ in range(4):
        state, rows, modes, _ = selected_space(H,S[used],time_count)
        pending = [j for j in range(len(S)) if j not in used]
        trial = S[pending]
        null = trial-(trial@rows)@rows.conj().T
        nullnorm = np.linalg.norm(null,axis=1)
        eligible = nullnorm > DEPENDENCE_TOL * np.linalg.norm(trial,axis=1)
        numerator = np.sqrt(np.mean(abs(trial @ modes)**2,axis=1))
        score = np.divide(numerator,nullnorm,out=np.zeros(len(trial)),where=eligible)
        top = float(np.max(score[eligible]))
        ties = np.flatnonzero(eligible & (abs(score-top) <= 1e-12*max(1.,top)))
        winner = int(min(pending[j] for j in ties))
        traceD.append({"state":state,"winner":winner,"scores": {
                        str(j):float(z) if legal else None
                        for j,z,legal in zip(pending,score,eligible,strict=True)},
                       "normalizer": {str(j):float(v) for j,v in zip(pending,nullnorm,strict=True)},
                       "candidates_skipped_near_dependent":[j for j,ok in zip(pending,eligible,strict=True) if not ok],
                       "tied_candidates":[pending[j] for j in ties]})
        used.append(winner)
        intermediate.append(winner)
    selected["D_hidden_mode"] = intermediate
    chosen = {}
    vectors = {}
    for arm in CANDIDATES:
        indices = base+selected[arm]
        state,_,_,v = selected_space(H,S[indices],time_count)
        chosen[arm] = {"added_receiver_ids":selected[arm],"receiver_order":indices,**state}
        vectors[arm] = v
    return {
        "baseline_receiver_ids":base,"baseline":base_meta,"physical_abscissa":abscissa.tolist(),
        "physical_abscissa_definition":units,"actual_receiver_coordinates":coords.tolist(),
        "uniform_quartile_attempts":traceA,"greedy_steps":traceD,"arms":chosen,
    },vectors


def check_negative_vector(S,G,v,u,mu):
    p,n = S.shape
    vnorm = norm(v)
    center = norm(S@v)**2+mu*(vnorm**2-u*np.vdot(v,G@v).real)
    gamma = 32*(n+p+1)*EPS
    upper = center+gamma*(2*vnorm)**2*((2*norm(S))**2+mu*(1+u*2*norm(G)))
    return float(upper)


def scan(arrays,vectors,plan):
    start,cpu=perf_counter(),process_time()
    saved_radii = read_json(SOURCE / "observation_radii.json")
    spectral_cache = {}
    case_map={}
    for key,array in arrays.items():
        with np.load(SOURCE/"results"/f"{key}.npz",allow_pickle=False) as f:
            y=f["Y_archive"]
        budgets = saved_radii[key]
        table={"radii":budgets,"arms":{}}
        for method in CANDIDATES:
            arm=plan[key]["arms"][method]
            ids=arm["receiver_order"]
            prob=sector.Sector(key,y[ids],{**array,"S":array["S"][ids]},budgets,spectral_cache)
            points=[]
            for upper in U_GRID:
                gram_store={}
                trials=[]
                for mu in sector.WEIGHTS:
                    ktest=check_negative_vector(prob.S,prob.G,vectors[key][method],upper,mu)
                    if ktest<0.:
                        trials.append({"u":upper,"mu":mu,"status":"negative_K_Rayleigh_upper","upper_test":ktest})
                        prob.cost["negative_rayleigh_verified_without_factor"]+=1
                    else:
                        result=prob.numerical(upper,mu,gram_store)
                        trials.append(result)
                valid=[v for v in trials if v["status"]=="numerical_only"]
                best=max(valid,key=lambda v:v["positive_gram_bound_squared"],default=None)
                points.append({"u":upper,"trials":trials,"best_positive_gram":best})
            max_caps={r:max((p["u"] for p in points if p["best_positive_gram"]
                        and p["best_positive_gram"]["positive_gram_excludes"][r]),default=None)
                      for r in budgets}
            table["arms"][method]={"grid":points,"largest_numerically_excluded_sampled_cap":max_caps,
                                   "cost":dict(prob.cost)}
            print(key,method,arm["hidden_u_star"],max_caps,flush=True)
        case_map[key]=table
    return case_map,{"wall_seconds":perf_counter()-start,"cpu_seconds":process_time()-cpu}


def audit_ring_admissibility(arrays, data):
    """Explain the ring's first unrealized sampled cap without changing mu."""
    start, cpu = perf_counter(), process_time()
    case = RING_AND_AUSTRIA[0]
    g = arrays[case]["G"]
    h = (g+g.conj().T)/2
    index = data["geometry"][case]["arms"]["D_hidden_mode"]["receiver_order"]
    S = arrays[case]["S"][index]
    samples = []
    for mu in sector.WEIGHTS:
        K = S.conj().T@S+mu*(np.eye(len(g))-3*h)
        lowest,basis=eigh(K,subset_by_index=(0,0),check_finite=False)
        upper = check_negative_vector(S,g,basis[:,0],3.,mu)
        samples.append({"mu":mu,"numeric_eig_min":float(lowest[0]),
                        "negative_Rayleigh_upper":upper,"strictly_not_positive_definite": upper<0})
    return {"case":case,"arm":"D_hidden_mode","u":3.,"samples":samples,
            "wall_seconds":perf_counter()-start,"cpu_seconds":process_time()-cpu}


def audit(arrays,data):
    start,cpu=perf_counter(),process_time()
    certified={}
    spectral_cache={}
    for case, row in data["data_replay"].items():
        by_method=row["arms"]
        a=by_method["A_uniform"]["largest_numerically_excluded_sampled_cap"]["frozen_observation_ball"]
        d=by_method["D_hidden_mode"]["largest_numerically_excluded_sampled_cap"]["frozen_observation_ball"]
        if case in RING_AND_AUSTRIA:
            choices=("A_uniform","D_hidden_mode")
        elif (d or 0)>(a or 0):
            choices=("D_hidden_mode",)
        else:
            choices=()
        if not choices:
            continue
        with np.load(SOURCE/"results"/f"{case}.npz",allow_pickle=False) as f:
            Y=f["Y_archive"]
        for method in choices:
            arm=by_method[method]
            positive=arm["largest_numerically_excluded_sampled_cap"]
            cap=max((x for x in positive.values() if x is not None),default=None)
            if cap is None:
                certified[f"{case}__{method}"]={"status":"no_positive_sample_to_audit"}
                continue
            point=next(p["best_positive_gram"] for p in arm["grid"] if p["u"]==cap)
            index=data["geometry"][case]["arms"][method]["receiver_order"]
            array = arrays[case]
            model=sector.Sector(case,Y[index],{**array,"S":array["S"][index]},row["radii"],spectral_cache)
            print("Arb audit",case,method,cap,point["mu"],flush=True)
            try:
                result=sector.audit(model,cap,point["mu"],project_positive_gram=True)
            except (ArithmeticError,sector.LinAlgError) as exc:
                result={"status":"audit_failed","reason":str(exc)}
            certified[f"{case}__{method}"]=result
    return {"entries":certified,"wall_seconds":perf_counter()-start,"cpu_seconds":process_time()-cpu}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--phase",choices=("grid","audit","all"),default="all")
    args=p.parse_args()
    with threadpool_limits(limits=1,user_api="blas"):
        arrays, receivers=get_archives()
        output=HERE/"design.json"
        if args.phase in ("grid","all"):
            start,cpu=perf_counter(),process_time()
            freeze,vectors,cost={}, {}, Counter()
            for key, arr in arrays.items():
                g=arr["G"]
                freeze[key],vectors[key]=design(arr["S"],(g+g.conj().T)/2,receivers[key],cost)
            data={
                "input":"Five frozen synthetic archives, same physical sources and original common absolute joint ball; only receiving rows are selected.",
                "source_measurement_use_at_geometry_freeze":0,"material_GT_reads":0,"maxwell_forward_solves":0,
                "geometry_rule": "prefix floor(original receivers/4), add four using physical quartile centers, adjacency, fixed seeded randomness, or largest normalized visibility of current hidden top eigenspace",
                "rng":SEED,"degenerate_subspace_tolerance":DEG_TOL,"near_dependence_cut":DEPENDENCE_TOL,
                "geometry":freeze,"caps":U_GRID,"fixed_mu":sector.WEIGHTS,
                "audit_precommit": "Audit the maximal positive-spectral passing grid point for uniform and guided in R10 and Austria; also guidance on any other layout where its maximum strictly exceeds uniform under the archived frozen budget.",
                "geometry_cost":{**dict(cost),"wall_seconds":perf_counter()-start,"cpu_seconds":process_time()-cpu},
            }
            dump(HERE/"DESIGN_FREEZE.json",data)
            data["data_replay"],data["data_replay_cost"]=scan(arrays,vectors,freeze)
            dump(output,data)
        else:
            data=read_json(output)
        if args.phase in ("audit","all"):
            data["strict_arb_audit"]=audit(arrays,data)
            data["strict_admissibility_check"]=audit_ring_admissibility(arrays,data)
            dump(output,data)


if __name__=="__main__":
    main()
