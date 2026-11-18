"""Six-dimensional direct SIVIA controls for R10's three-plus-one readings."""

import argparse
from itertools import combinations
import json
import math
import time

import numpy as np
from flint import acb, acb_mat, arb, arb_mat

from all_pairs import distance
from bounded_certificate import E, GAMMA, ROOT, ball, dot, feasible_witness, matrix, row, square
from pivot_certificate import cell_tables, cover, halves, tensor


BOX = (.5, 30.)
CASES = (("natural", "fresh6-00"), ("natural", "fresh6-06"),
         ("confirm", "confirm6-05"), ("confirm", "confirm6-00"))


def choose_case(cohort, name):
    """Freeze an alternate-centred strip from the previously archived R10 projections."""
    source = ROOT/"results"
    data = np.load(source/f"{cohort}_inputs.npz")
    summary = json.loads((source/f"{cohort}_all_pairs.json").read_text())
    mask_data = np.load(source/f"{cohort}_all_pairs_masks.npz")
    key = f"{name}/active"
    record = summary["records"][key]
    z = np.asarray(record["witness_vectors"])
    eligible = [j for j in range(1, len(z)) if record["old_witnesses"][j]["inside"]
                and record["new_witnesses"][j]["outside"]]
    provably_feasible = not eligible
    j = max(eligible, key=record["all_charts_linf_radii"].__getitem__) if eligible else 1
    depth = int(summary["depth"])
    resolution = 2**depth
    step = (BOX[1]-BOX[0])/resolution
    pairs = summary["pairs"]
    gaps = [distance(np.argwhere(np.unpackbits(bits).reshape(resolution,resolution)),
                     pair, z[j], step) for bits,pair in zip(mask_data[key],pairs)]
    i_pair = int(np.argmax(gaps))
    p = pairs[i_pair] if eligible else data[f"{name}/pivots"].astype(int).tolist()
    radius = math.floor(gaps[i_pair]*512)/1024 if eligible else .25
    lo = np.full(6, BOX[0])
    hi = np.full(6, BOX[1])
    lo[p] = np.maximum(lo[p], z[j,p]-radius)
    hi[p] = np.minimum(hi[p], z[j,p]+radius)
    small_lo = np.maximum(lo,z[j]-1/32)
    small_hi = np.minimum(hi,z[j]+1/32)
    H, R = data[f"{name}/H"], data[f"{name}/R"]
    B, hy, y = (data[f"{key}/{k}"] for k in ("B","y_host","y"))
    info = dict(cohort=cohort,sample=name,channels=B.shape[1],pivots=p,
                old_witness_id=j, witness=z[j].tolist(), saved_guaranteed_radius=gaps[i_pair],
                region_half_width=radius, region_each_axis=np.c_[lo,hi].tolist(),
                positive_task=not provably_feasible, old_witness=record["old_witnesses"][j],
                new_witness=record["new_witnesses"][j])
    return (H,R,B,hy,y), (lo,hi), (small_lo,small_hi), info


def restricted_cell_table(inputs, box):
    H,R,B,hy,y = inputs
    lo,hi = box
    table = cell_tables(H,R,B,hy,y,BOX)

    def interpolate(left,right,t):
        return ((1-t)*left+t*right)

    for j,(lower,upper) in enumerate(zip(lo,hi)):
        if (lower,upper) == BOX:
            continue
        first, _, last = table[j]
        one, two = [arb(float(x)) for x in (lo[j],hi[j])]
        table[j] = []
        for q in (one,(one+two)/2,two):
            t = (q-arb(BOX[0]))/arb(BOX[1]-BOX[0])
            table[j].append(
                ([interpolate(a,b,t) for a,b in zip(first[0],last[0])],
                 [interpolate(a,b,t) for a,b in zip(first[1],last[1])],
                 interpolate(first[2],last[2],t), interpolate(first[3],last[3],t)))
    return table


def projective_proof(inputs, box, pivots, depth=8):
    """Apply the R10 two-coordinate Bernstein cover to exactly the same strip."""
    tic, wall = time.process_time(), time.perf_counter()
    table = restricted_cell_table(inputs,box)
    pixels, tested = cover(table,pivots,depth)
    return dict(excluded_all=len(pixels)==0,unresolved_boxes=len(pixels),visited_boxes=tested,
                unresolved_relative_measure=len(pixels)/2**(2*depth), depth=depth,
                cpu_seconds=time.process_time()-tic,wall_seconds=time.perf_counter()-wall)


def projective_chart_fallback(inputs, box, original, depth=8):
    """Try the other R10 charts in lexicographic order on the same fixed region."""
    tic,wall = time.process_time(),time.perf_counter()
    table = restricted_cell_table(inputs,box)
    count = len(table)
    coefficients = {triple:tensor(table,list(triple)) for triple in combinations(range(count),3)}
    support = [np.ones(2**depth,dtype=bool) for _ in range(count)]
    pair_masks = {}
    steps = []
    for pair in combinations(range(count),2):
        if set(pair) == set(original):
            continue
        pixels,visits = cover(table,list(pair),depth,coefficients)
        steps.append(dict(pivots=pair,visited_boxes=visits,unresolved_boxes=len(pixels)))
        if not len(pixels):
            return dict(excluded_all=True,reason="one alternative chart excludes the whole strip",
                        decisive_pair=pair,charts=steps,cpu_seconds=time.process_time()-tic,
                        wall_seconds=time.perf_counter()-wall)
        kept = np.zeros((2**depth,2**depth),dtype=bool)
        kept[tuple(pixels.T)] = True
        pair_masks[pair] = kept
    while True:
        previous = sum(np.count_nonzero(m) for m in support)
        for (i,j),kept in pair_masks.items():
            available = kept & support[i][:,None] & support[j][None,:]
            support[i] &= np.any(available,axis=1)
            support[j] &= np.any(available,axis=0)
        after = [int(np.count_nonzero(v)) for v in support]
        if 0 in after or sum(after)==previous:
            break
    return dict(excluded_all=0 in after,reason="intersect pairwise support on shared coordinate bins",
                surviving_bins=after,charts=steps,cpu_seconds=time.process_time()-tic,
                wall_seconds=time.perf_counter()-wall)


def determinant_controls(inputs, box):
    """Evaluate the six-dimensional exact gain-profile polynomial in Bernstein form."""
    H,R,B,hy,y = inputs
    k,m = B.shape
    H, R, Y, obs = matrix(H), list(row(R).entries()), matrix(B), list(row(y).entries())
    host = list(row(hy).entries())
    gamma_ratio = (E/GAMMA)**2
    pairs = list(combinations(range(m),2))
    lo,hi = box
    q = [(arb(float(a)),(arb(float(a))+arb(float(b)))/2,arb(float(b))) for a,b in zip(lo,hi)]
    c = np.empty((3,)*k,dtype=object)
    for ids in np.ndindex(c.shape):
        z = [q[j][idx] for j,idx in enumerate(ids)]
        A = [[int(i==j)-H[i,j]*z[j] for j in range(k)] for i in range(k)]
        D = acb_mat(A).det()
        # det([[A,b],[-RZ,hy]]) = D*y_model. Each variable occurs only in one column.
        tail = [-R[j]*z[j] for j in range(k)]
        Ynum = [acb_mat([A[i]+[Y[i,a]] for i in range(k)]+[tail+[host[a]]]).det()
                for a in range(m)]
        residual = [D*obs[a]-Ynum[a] for a in range(m)]
        # Exact joint gain/additive profile when D != 0. If the singular
        # forward equations are consistent, all bordered determinants vanish.
        wedge = sum((square(obs[a]*Ynum[b]-obs[b]*Ynum[a]) for a,b in pairs),arb(0))
        c[ids] = (wedge+gamma_ratio*sum((square(a) for a in residual),arb(0))
                  -E**2*sum((square(a) for a in Ynum),arb(0))
                  -E**2*gamma_ratio*square(D))
    for axis in range(k):
        g = np.moveaxis(c,axis,0)
        g[1] = 2*g[1]-(g[0]+g[2])/2
    peak = max(float(abs(a).upper()) for a in c.flat)
    scale = arb(2)**(45-math.frexp(peak)[1]) if peak else arb(1)
    out = [int((a*scale).lower().floor().unique_fmpz()) for a in c.flat]
    if np.max(np.abs(out)) >= 2**48:
        raise ValueError("Fixed-point control scaling overflow")
    return np.asarray(out,dtype=np.int64).reshape(c.shape)


def sieve_polynomial(inputs, box, cpu_limit, depth_limit=64):
    """Direct six-dimensional Bernstein bisection of the profiled forward residual."""
    tic, wall = time.process_time(), time.perf_counter()
    coefficient = determinant_controls(inputs,box)
    build = time.process_time()-tic
    width0 = box[1]-box[0]
    stack = [(coefficient,width0,0)]
    visited = excluded = leaf_count = 0
    leaf_measure = excluded_measure = 0.
    max_depth = 0
    while stack and time.process_time()-tic < cpu_limit:
        C, widths, depth = stack.pop()
        visited += 1
        max_depth = max(max_depth,depth)
        if np.all(C>0):
            excluded += 1
            excluded_measure += 2.**-depth
            continue
        if depth >= depth_limit:
            leaf_count += 1
            leaf_measure += 2.**-depth
            continue
        # The largest adjacent control difference measures variation in each axis.
        # Favor actual width as a fallback in nearly constant directions.
        variation = [np.max(abs(np.diff(C,axis=i))) for i in range(len(widths))]
        axis = int(np.argmax(variation)) if max(variation)>0 else int(np.argmax(widths))
        children = halves(C,axis)
        new_widths = widths.copy()
        new_widths[axis] /= 2
        # Resolve the positive-nearest sibling first when children differ.
        for child in sorted(children,key=lambda branch: np.min(branch)):
            stack.append((child,new_widths,depth+1))
    count = len(stack)+leaf_count
    return dict(excluded_all=count==0,visited_boxes=visited,excluded_boxes=excluded,
                unresolved_boxes=count,unresolved_at_depth_limit=leaf_count,
                unresolved_relative_measure=leaf_measure+sum(2.**-b[2] for b in stack),
                excluded_relative_measure=excluded_measure, deepest_visited=max_depth,
                build_cpu_seconds=build,cpu_seconds=time.process_time()-tic,
                wall_seconds=time.perf_counter()-wall,cpu_limit_seconds=cpu_limit)


class ForwardBoxes:
    """Centered rigorous Neumann bound for the original current and receiver equations."""
    def __init__(self, inputs):
        H,R,B,hy,y = inputs
        self.H = matrix(H)
        self.R = list(row(R).entries())
        self.B = matrix(B)
        self.host = row(hy)
        self.observed = list(row(y).entries())
        self.observed_float = y
        self.k = B.shape[0]
        self.calls = self.contractions = 0

    def test(self, lo, hi):
        self.calls += 1
        lo,hi = [arb(float(v)) for v in lo],[arb(float(v)) for v in hi]
        z = [(a+b)/2 for a,b in zip(lo,hi)]
        h = [(b-a)/2 for a,b in zip(lo,hi)]
        k = self.k
        A = acb_mat([[int(i==j)-self.H[i,j]*z[j] for j in range(k)] for i in range(k)])
        try:
            K = A.inv()
        except ZeroDivisionError:
            return False
        C = K*self.H
        abs_C_h = [[abs(C[i,j])*h[j] for j in range(k)] for i in range(k)]
        G = np.array([[float(v.upper()) for v in r] for r in abs_C_h])
        # Approximate positive weights are only proposals. Arb verifies contraction.
        try:
            weights = np.linalg.solve(np.eye(k)-G,np.ones(k))
        except np.linalg.LinAlgError:
            return False
        if not np.all(weights>0):
            return False
        weights = [arb(float(w)) for w in weights]
        ratios = [sum((g*w for g,w in zip(r,weights)),arb(0))/wi for r,wi in zip(abs_C_h,weights)]
        beta = max(t.upper() for t in ratios)
        if not beta < 1:
            return False
        self.contractions += 1

        Et = K*self.B
        zR = acb_mat([[self.R[i]*z[i] for i in range(k)]])
        yc = zR*Et+self.host
        yc_mid = np.array([complex(float(v.real.mid()),float(v.imag.mid())) for v in yc.entries()])
        r0 = self.observed_float-yc_mid
        # Precision direction from the center's exact-gain-profile metric.
        dual = r0-yc_mid*(np.vdot(yc_mid,r0)/(np.vdot(yc_mid,yc_mid).real+.02**2))
        peak = np.linalg.norm(dual)
        if not np.isfinite(peak) or not peak>0:
            return False
        q = [ball(a.conjugate()/peak) for a in dual]
        lnorm2 = sum((square(a) for a in q),arb(0))
        field = list((Et*acb_mat([[a] for a in q])).entries())
        adjoint = [r+p for r,p in zip(self.R,(zR*C).entries())]
        flow = [sum((g*abs(v) for g,v in zip(r,field)),arb(0))/wi
                for r,wi in zip(abs_C_h,weights)]
        relative_remainder = max(t.upper() for t in flow)/(1-beta)
        deriv = [p*v for p,v in zip(adjoint,field)]
        tail = relative_remainder*sum((abs(p)*hi*wi for p,hi,wi in zip(adjoint,h,weights)),arb(0))
        y_center = dot(q,yc.entries())
        real_radius = sum((abs(a.real)*b for a,b in zip(deriv,h)),arb(0))+tail
        imag_radius = sum((abs(a.imag)*b for a,b in zip(deriv,h)),arb(0))+tail
        enclosure = y_center + acb(arb(0,real_radius.upper()),arb(0,imag_radius.upper()))
        measured = dot(q,self.observed)
        return bool(square(measured-enclosure) > E**2*lnorm2 + GAMMA**2*square(enclosure))


def sieve_forward(inputs, box, cpu_limit, depth_limit=64):
    tic,wall = time.process_time(),time.perf_counter()
    tester = ForwardBoxes(inputs)
    lo,hi = box
    stack = [(lo,hi,0)]
    omitted = excluded = 0
    omitted_measure = excluded_measure = 0.
    max_depth=0
    while stack and time.process_time()-tic < cpu_limit:
        a,b,depth = stack.pop()
        max_depth=max(depth,max_depth)
        if tester.test(a,b):
            excluded += 1
            excluded_measure += 2.**-depth
            continue
        if depth >= depth_limit:
            omitted += 1
            omitted_measure += 2.**-depth
            continue
        axis = int(np.argmax(b-a))
        center = (a[axis]+b[axis])/2
        alo,bhi = a.copy(),b.copy()
        alo[axis]=bhi[axis]=center
        stack.extend(((alo,b,depth+1),(a,bhi,depth+1)))
    count = len(stack)+omitted
    return dict(excluded_all=count==0,visited_boxes=tester.calls,resolvent_contracted_boxes=tester.contractions,
                excluded_boxes=excluded,unresolved_boxes=count,unresolved_at_depth_limit=omitted,
                unresolved_relative_measure=omitted_measure+sum(2.**-b[2] for b in stack),
                excluded_relative_measure=excluded_measure,deepest_visited=max_depth,
                cpu_seconds=time.process_time()-tic,wall_seconds=time.perf_counter()-wall,
                cpu_limit_seconds=cpu_limit)


def unselected_tasks():
    """Define every region from archived old candidate arrays alone."""
    studies = (("confirm","confirm6-00"), ("confirm","confirm6-05"),
               ("natural","fresh6-00"), ("natural","fresh6-02"))
    source = ROOT/"results"
    tasks = []
    for cohort,name in studies:
        with np.load(source/f"{cohort}_inputs.npz") as data:
            witnesses = data[f"{name}/clean_witnesses"]
        if not len(witnesses):
            continue
        spans = np.ptp(witnesses,axis=0)
        pair = np.argsort(-spans,kind="stable")[:2].tolist() if len(witnesses)>1 else [0,1]
        for j,z in enumerate(witnesses):
            lo = np.full(6,BOX[0])
            hi = np.full(6,BOX[1])
            lo[pair] = np.maximum(BOX[0],z[pair]-.5)
            hi[pair] = np.minimum(BOX[1],z[pair]+.5)
            tasks.append(dict(cohort=cohort,sample=name,old_witness_id=j,old_witness=z.tolist(),
                              old_candidate_count=len(witnesses),old_coordinate_spans=spans.tolist(),
                              pivots=pair,region_each_axis=np.c_[lo,hi].tolist()))
    return tasks


def run_unselected(args):
    report = dict(precision_bits=128,noise_scales=[.0002,.01],box=BOX,strip_halfwidth=.5,
                  scope="All archived clean_witnesses in confirm6-00, confirm6-05, fresh6-00 and fresh6-02.",
                  region_rule="Freeze all strips before opening active observations. Take the two largest old candidate spans, stable by index, or indices 0 and 1 for a single candidate. Center one strip at each old candidate and keep the other four coordinates unrestricted.",
                  witness_shortcut="The same rigorous Arb point witness check is shared by both methods. A feasible witness inside a strip establishes that full exclusion is impossible, so neither exclusion algorithm runs.",
                  cpu_limit_seconds_direct=args.seconds, r10_projection_depth=8,
                  r10_fallback="If the old-span chart is inconclusive, try every other pair in lexicographic order before intersecting their surviving support. Fallback time is recorded separately and in the total.",
                  regions=[])
    path = ROOT/"results"/"direct_box_unselected.json"
    tasks = unselected_tasks()
    for meta in tasks:
        name = meta["sample"]
        if args.case and name not in args.case:
            continue
        tic,wall = time.process_time(),time.perf_counter()
        with np.load(ROOT/"results"/f"{meta['cohort']}_inputs.npz") as data:
            H,R = (data[f"{name}/{key}"] for key in ("H","R"))
            B,hy,y = (data[f"{name}/active/{key}"] for key in ("B","y_host","y"))
        box = np.array(meta["region_each_axis"]).T
        z = np.array(meta["old_witness"])
        check_tic,check_wall = time.process_time(),time.perf_counter()
        meta["old_point_witness"] = feasible_witness(z,H,R,B[:,:-1],hy[:-1],y[:-1])
        meta["active_point_witness"] = feasible_witness(z,H,R,B,hy,y)
        meta["shared_witness_cpu_seconds"] = time.process_time()-check_tic
        meta["shared_witness_wall_seconds"] = time.perf_counter()-check_wall
        if meta["active_point_witness"]["inside"]:
            for label in ("r10_pivot_bernstein","direct_full_bernstein"):
                meta[label] = dict(status="feasible_point_shared_shortcut",
                                   executed=False,excluded_all=False,visited_boxes=0,cpu_seconds=0.,wall_seconds=0.)
        else:
            inputs = H,R,B,hy,y
            meta["r10_pivot_bernstein"] = projective_proof(inputs,box,meta["pivots"],depth=8)
            meta["direct_full_bernstein"] = sieve_polynomial(inputs,box,args.seconds)
            for label in ("r10_pivot_bernstein","direct_full_bernstein"):
                r = meta[label]
                r["executed"] = True
                r["status"] = ("excluded_all" if r["excluded_all"] else
                               "partial_exclusion" if r["unresolved_relative_measure"] < 1. else "unresolved")
            if meta["r10_pivot_bernstein"]["excluded_all"] and meta["direct_full_bernstein"]["excluded_all"]:
                a,b = (meta[l]["cpu_seconds"] for l in ("r10_pivot_bernstein","direct_full_bernstein"))
                w = meta["shared_witness_cpu_seconds"]
                meta["matched_proof_cpu_ratio_direct_over_r10"] = b/a
                meta["matched_proof_cpu_ratio_including_shared_witness"] = (b+w)/(a+w)
            if not meta["r10_pivot_bernstein"]["excluded_all"]:
                fallback = projective_chart_fallback(inputs,box,meta["pivots"],depth=8)
                fallback["total_cpu_seconds_including_first_chart"] = (
                    fallback["cpu_seconds"]+meta["r10_pivot_bernstein"]["cpu_seconds"])
                fallback["total_wall_seconds_including_first_chart"] = (
                    fallback["wall_seconds"]+meta["r10_pivot_bernstein"]["wall_seconds"])
                meta["r10_chart_fallback"] = fallback
        meta["combined_cpu_seconds"] = time.process_time()-tic
        meta["combined_wall_seconds"] = time.perf_counter()-wall
        print(name,meta["old_witness_id"],
              [(label,meta[label]["status"],round(meta[label]["cpu_seconds"],6))
               for label in ("r10_pivot_bernstein","direct_full_bernstein")],flush=True)
        report["regions"].append(meta)
        path.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")


def main(args):
    report = dict(precision_bits=128,noise_model="Shared complex gain 1/100 and additive complex vector 1/5000 on one Euclidean ellipsoid.",
                  region_rule="Iff old alternative is feasible and excluded by active observations, take the highest-radius saved pair and round half the radius down to 1/1024. Remaining four coordinates span [0.5,30]. Otherwise run a known-feasible negative control.",
                  central_reference_half_width=1/32,regions=[])
    path=ROOT/"results"/"direct_box_control.json"
    for cohort,name in CASES:
        inputs,box,small,meta=choose_case(cohort,name)
        if args.case and name not in args.case:
            continue
        meta["r10_pivot_bernstein"] = projective_proof(inputs,box,meta["pivots"])
        print(name,"R10",meta["r10_pivot_bernstein"],flush=True)
        for label,fn in (("direct_forward_neumann",sieve_forward),("direct_full_bernstein",sieve_polynomial)):
            meta[label] = fn(inputs,box,args.seconds)
            meta[label]["reference_small_6d_box"] = fn(inputs,small,1.,16)
            print(name,label,meta[label],flush=True)
        report["regions"].append(meta)
        path.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--seconds",type=float)
    cli.add_argument("--case",action="append")
    cli.add_argument("--unselected",action="store_true",
                     help="Emit direct_box_unselected.json for strips chosen from old witnesses only")
    args = cli.parse_args()
    if args.seconds is None:
        args.seconds = 8. if args.unselected else 4.
    if args.unselected:
        run_unselected(args)
    else:
        main(args)
