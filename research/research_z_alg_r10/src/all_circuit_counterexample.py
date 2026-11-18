"""Check every singleton and pair among the six-cell cofactor circuits."""

from itertools import combinations
import json
from pathlib import Path

import numpy as np
from flint import acb_mat,arb
from scipy.optimize import brentq

from bounded_certificate import E,GAMMA,conjdot,dot,feasible_witness,matrix,row,square
from joint_errors import Case
from pivot_certificate import cell_tables, determinant

ROOT = Path(__file__).resolve().parents[1]/"results"
TRIPLES = list(combinations(range(6),3))
PAIRS = list(combinations(range(len(TRIPLES)),2))
BOX = (.5,30.)


def load(data,name):
    H,R = (data[f"{name}/{key}"] for key in ("H","R"))
    B,host,y = (data[f"{name}/active/{key}"] for key in ("B","y_host","y"))
    return H,R,B,host,y


def circuit_design(case,q):
    U = np.eye(len(q))-case.H*q
    C = case.N@U
    cols = []
    intercepts = []
    for triple in TRIPLES:
        ids = list(triple)
        weights = np.array([(-1)**j*np.linalg.det(np.delete(C[:,ids],j,axis=1)) for j in range(3)])
        cols.append(U[:,ids]@weights)
        intercepts.append((q[ids]*case.R[ids])@weights)
    A = case.B_inv @ np.column_stack(cols)
    c = case.host@A + np.asarray(intercepts)
    residual = case.y@A-c
    design = np.vstack((.01*c,2e-4*A))
    return design,residual


def floating_circuits(case,q):
    design,residual = circuit_design(case,q)
    scales = np.linalg.norm(design,axis=0)
    if np.any(scales==0):
        return dict(usable=False,reason="zero raw noise-design column")
    design /= scales
    residual /= scales
    singles = abs(residual)**2
    pair_indices = np.asarray(PAIRS)
    vectors = np.moveaxis(design[:,pair_indices],1,0)
    _, R = np.linalg.qr(vectors)
    pair_gaps = np.linalg.svd(R,compute_uv=False)
    ratios = pair_gaps[:,1]/pair_gaps[:,0]
    if np.min(ratios)<1e-12:
        return dict(usable=False,reason="nearly dependent pair",
                    minimum_pair_singular_ratio=float(np.min(ratios)),
                    near_dependent_pairs=pair_indices[ratios<1e-12].tolist())
    rhs = residual[pair_indices].conj()[...,None]
    scores = np.sum(abs(np.linalg.solve(R.conj().transpose(0,2,1),rhs))**2,axis=(1,2))
    i,j = int(np.argmax(singles)),int(np.argmax(scores))
    return dict(usable=True,maximum_scalar=float(singles[i]),maximum_pair=float(scores[j]),
                maximizing_scalar=i,maximizing_pair=PAIRS[j],
                minimum_pair_singular_ratio=float(np.min(ratios)),
                maximum_from_fixed_chart=float(max(singles[t] for t,tri in enumerate(TRIPLES)
                                                   if set(case.pivots).issubset(tri))))


def projection_diagnostics(case,q,score):
    design,_ = circuit_design(case,q)
    design /= np.linalg.norm(design,axis=0)
    U,s,_ = np.linalg.svd(design,full_matrices=False)
    X = U[:,:4].conj().T @ design
    X /= np.linalg.norm(X,axis=0)
    centers = np.moveaxis(X[:,np.asarray(PAIRS)],1,0)
    pair_basis = np.linalg.qr(centers,mode="reduced")[0]
    mean_single = X@X.conj().T / len(TRIPLES)
    mean_pair = np.einsum("air,ajr->ij",pair_basis,pair_basis.conj()) / len(PAIRS)
    eig_single = np.linalg.eigvalsh(mean_single)
    eig_pair = np.linalg.eigvalsh(mean_pair)
    return dict(physical_noise_dimension=4,retained_relative_singular_values=(s[:4]/s[0]).tolist(),
                discarded_relative_singular_value=float(s[4]/s[0]),
                mean_single_projection_eigenvalues=eig_single.tolist(),
                mean_pair_projection_eigenvalues=eig_pair.tolist(),
                eta_single_lower_from_min_mean_projection=float(eig_single[0]),
                eta_pair_lower_from_min_mean_projection=float(eig_pair[0]),
                eta_single_upper_from_observation=score["maximum_scalar"]/score["forward_score"],
                eta_pair_upper_from_observation=score["maximum_pair"]/score["forward_score"])


def forward_box_witness(inputs,q):
    # The direct profile's cleared polynomial uses exact bordered determinants.
    H,R,B,host,observed = inputs
    H,R,B,host,observed = matrix(H),list(row(R).entries()),matrix(B),list(row(host).entries()),list(row(observed).entries())
    k,m = B.nrows(),B.ncols()
    A = [[int(i==j)-H[i,j]*q[j] for j in range(k)] for i in range(k)]
    D = acb_mat(A).det()
    tail = [-R[j]*q[j] for j in range(k)]
    Q = [acb_mat([A[i]+[B[i,a]] for i in range(k)]+[tail+[host[a]]]).det() for a in range(m)]
    residual = [D*observed[a]-Q[a] for a in range(m)]
    c = (E/GAMMA)**2
    wedge = sum((square(observed[a]*Q[b]-observed[b]*Q[a])
                 for a,b in combinations(range(m),2)),arb(0))
    psi = (wedge+c*sum((square(a) for a in residual),arb(0))
           -E**2*(sum((square(a) for a in Q),arb(0))+c*square(D)))
    regular = square(D)>0
    normalizer = E**2*(sum((square(a) for a in Q),arb(0))+c*square(D))
    normalized = 1+psi/normalizer if regular and normalizer>0 else None
    return dict(cleared_profile_enclosure=str(psi),resolvent_determinant_squared=str(square(D)),
                normalized_joint_score=str(normalized) if normalized is not None else None,
                normalized_joint_score_lower=str(normalized.lower()) if normalized is not None else None,
                regular=bool(regular),outside=bool(regular and psi>0))


def rigorous_circuits(inputs,q,halfwidth=0.):
    H,R,B,host,observed = inputs
    base = cell_tables(H,R,B,host,observed,BOX)
    cells = []
    qballs = [arb(float(qi),halfwidth) for qi in q]
    for i,qi in enumerate(qballs):
        lo,_,hi = base[i]
        t = (qi-arb(1)/2)/(arb(59)/2)
        cells.append(([a+t*(b-a) for a,b in zip(lo[0],hi[0])],
                      [a+t*(b-a) for a,b in zip(lo[1],hi[1])],
                      lo[2]+t*(hi[2]-lo[2]),lo[3]+t*(hi[3]-lo[3])))
    channels,residual = [],[]
    for ids in TRIPLES:
        cols = [cells[i] for i in ids]
        null_columns = [entry[0] for entry in cols]
        weights = [(-1)**j*determinant(null_columns[:j]+null_columns[j+1:]) for j in range(3)]
        alpha = [dot(weights,terms) for terms in zip(*(entry[1] for entry in cols))]
        c = dot(weights,[entry[2] for entry in cols])
        channels.append([GAMMA*c,*[E*a for a in alpha]])
        residual.append(dot(weights,[entry[3] for entry in cols]))
    gram_diags = [sum((square(t) for t in x),arb(0)) for x in channels]
    residual_sq = [square(t) for t in residual]
    margins = [g-r for g,r in zip(gram_diags,residual_sq)]
    singles = []
    for j,margin in enumerate(margins):
        g = gram_diags[j]
        singles.append(dict(triple=TRIPLES[j],gram_diagonal=str(g),psd_diagonal=str(margin),
                            score=str(residual_sq[j]/g) if g>0 else None,
                            strict_pass=bool(margin>0)))
    pairs = []
    pair_upper_bounds = []
    for i,j in PAIRS:
        a,b = gram_diags[i],gram_diags[j]
        gij = conjdot(channels[i],channels[j])
        gram_det = a*b-square(gij)
        gij_shifted = gij-residual[i].conjugate()*residual[j]
        psd_det = margins[i]*margins[j]-square(gij_shifted)
        numerator = b*residual_sq[i]+a*residual_sq[j]-2*(residual[i]*gij*residual[j].conjugate()).real
        score = numerator/gram_det if gram_det>0 else None
        if score is not None:
            pair_upper_bounds.append(score.upper())
        pairs.append(dict(circuit_indices=[i,j],gram_determinant=str(gram_det),
                          psd_determinant=str(psd_det),score=str(score) if score is not None else None,
                          strict_pass=bool(margins[i]>0 and margins[j]>0 and psd_det>0)))
    full = (forward_box_witness(inputs,qballs) if halfwidth else
            feasible_witness(q,H,R,B,host,observed))
    complete = all(s["strict_pass"] for s in singles) and all(p["strict_pass"] for p in pairs)
    bounds = dict(max_single_score_upper=str(max((r/g).upper() for r,g in zip(residual_sq,gram_diags) if g>0)),
                  max_pair_score_upper=str(max(pair_upper_bounds)) if pair_upper_bounds else None)
    return dict(precision_bits=128,strict_singletons=len([s for s in singles if s["strict_pass"]]),
                strict_pairs=len([p for p in pairs if p["strict_pass"]]),
                singular_or_uncertain_pairs=len([p for p in pairs if p["score"] is None]),
                all_singletons_and_pairs_pass=complete, full_forward_witness=full,
                strict_counterexample=bool(complete and full["outside"]),bounds=bounds,
                singletons=singles,pairs=pairs)


def main():
    old = json.loads((ROOT/"joint_errors.json").read_text())["cases"]
    records,exact = [],None
    for cohort in ("natural","confirm"):
        with np.load(ROOT/f"{cohort}_inputs.npz") as data:
            for name in data["case_ids"]:
                case = Case(data,name)
                inputs = load(data,name)
                for point in old[name]["diagnostic_points"]:
                    q = np.array(point["q"])
                    floating = floating_circuits(case,q)
                    rec = dict(cohort=cohort,case=str(name),q=q.tolist(),axis=point["axis"],sign=point["sign"],
                               step=point["step"],fixed_pivots=case.pivots,
                               best_fixed_pair=point["maximum_pair"],forward_score=point["forward_score"],
                               **floating)
                    records.append(rec)
                    if exact is None and floating["usable"] and floating["maximum_pair"]<1:
                        certificate = rigorous_circuits(inputs,q)
                        if certificate["strict_counterexample"]:
                            exact = dict(rec,selected_target=1.2,certificate=certificate)
    out = dict(material_order_zero_based=True,noise_scales=[2e-4,.01],triples=TRIPLES,pair_indices=PAIRS,
               selection="First strict example in stored natural/confirm case order and signed direction order. Points were diagnosed from stored truth and active readings; not a blind benchmark. A 1.01 target is attempted only if no stored 1.2 point gives a strict example.",
               n_points=len(records),n_numeric_pair_pass=sum(r["usable"] and r["maximum_pair"]<1 for r in records),
               n_numeric_single_pass=sum(r["usable"] and r["maximum_scalar"]<1 for r in records),
               n_numeric_usable=sum(r["usable"] for r in records),
               min_numeric_pair_singular_ratio=min(
                   r["minimum_pair_singular_ratio"] for r in records if r["usable"]),
               n_fixed_chart_pair_pass=sum(r["best_fixed_pair"]<1 for r in records))
    if exact is None:
        # A fallback uses the first stored path only, fixed independently of the all-circuit scan.
        rec = records[0]
        name = rec["case"]
        with np.load(ROOT/f"{rec['cohort']}_inputs.npz") as data:
            case = Case(data,name)
            direction = np.eye(6)[rec["axis"]]*rec["sign"]
            m = data[f"{name}/truth"]
            step = brentq(lambda t:case.profile(m+t*direction)-1.01,0,rec["step"],xtol=1e-12)
            q = m+step*direction
            cert = rigorous_circuits(load(data,name),q)
            exact = dict(rec,q=q.tolist(),step=step,forward_score=case.profile(q),
                         selected_target=1.01,certificate=cert)
            exact.update(floating_circuits(case,q))
            exact["best_fixed_pair"] = case.evaluate(q,"fallback")["maximum_pair"]
    with np.load(ROOT/f"{exact['cohort']}_inputs.npz") as data:
        case = Case(data,exact["case"])
        center = np.array(exact["q"])
        inputs = load(data,exact["case"])
        exact["acquisition_archive"] = f"results/{exact['cohort']}_inputs.npz"
        exact["input_complex_real_imag"] = {
            key:np.stack((value.real,value.imag),axis=-1).tolist()
            for key,value in zip(("H","R","active_B","active_host","active_observed"),inputs)
        }
        exact["projection_diagnostics"] = projection_diagnostics(case,center,exact)
        attempts = []
        exact["box_certificate"] = None
        for depth in (10,14,18,22):
            radius = 2.**-depth
            certificate = rigorous_circuits(inputs,center,halfwidth=radius)
            keys = ("strict_counterexample","strict_singletons","strict_pairs",
                    "singular_or_uncertain_pairs","full_forward_witness","bounds")
            attempts.append(dict(halfwidth=radius,**{key:certificate[key] for key in keys}))
            if certificate["strict_counterexample"]:
                exact["box_certificate"] = dict(halfwidth=radius,axis_bounds=np.c_[center-radius,center+radius].tolist(),
                                                certificate=certificate)
                break
        exact["dyadic_box_attempts"] = attempts
    out.update(example=exact,records=records)
    (ROOT/"all_circuit_counterexample.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print({k:out[k] for k in ("n_points","n_numeric_single_pass","n_numeric_pair_pass","n_numeric_usable")})
    print("example",exact["case"],exact["axis"],exact["sign"],exact["selected_target"],
          exact.get("maximum_scalar"),exact.get("maximum_pair"),
          exact["certificate"]["bounds"],exact["certificate"]["full_forward_witness"])


if __name__ == "__main__":
    main()
