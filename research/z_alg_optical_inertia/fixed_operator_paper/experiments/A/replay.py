"""Independent one-sided reassembly from original arrays and frozen frames."""
import argparse
from collections import Counter
import json
from time import perf_counter

from common import PUBLIC, REGION, RESULTS, SEAL, digest, load, write_json
from baseline_neumann import replay_neumann
from bounds import (Acquisition, acb, acb_mat, arb, ball, bt, bdiag, ctx, frob, gersh,
                    high, low, midpoint, propose_frame, rows, select, shift, strict_pd, validate_frame, weight_vector)
import numpy as np
from scipy.linalg import eigvalsh, hadamard, qr

ctx.prec = 192
ctx.threads = 1


def original(acq, ell, upp, weight, optical):
    j, e = acq.Jb, acq.Eb
    l, u, w, v = ([arb(float(x)) for x in p] for p in (ell, upp, weight, optical))
    cross = [acb(-z*(a+b)/2, -q/2) for a,b,z,q in zip(l,u,w,v,strict=True)]
    joint = bt(j)*bdiag(cross,e)
    return bt(acq.Ob)*acq.Ob + bt(j)*bdiag(w,j) + joint + bt(joint) + bt(e)*bdiag(
        [a*b*z for a,b,z in zip(l,u,w,strict=True)], e)


def replay_terms(matrix, record, frames, radius):
    lower = arb(0)
    states = record["states"]
    gap_minimum = None
    for term in record["terms"]:
        threshold = term["threshold"]
        frame = frames[term["frame"]]
        assert frame.shape == (matrix.nrows(), states+term["j"])
        gap = gersh(bt(ball(frame))*shift(matrix, states, threshold)*ball(frame))
        if not gap > 0:
            return False
        gap_minimum = gap if gap_minimum is None else min(gap, gap_minimum)
        lower += arb(threshold)*term["term_multiplicity"]
    return bool(lower > arb(radius)**2)


def independent_psd(acq, ell, upp, rec, frames):
    identity, H = np.eye(acq.s), hadamard(acq.s, dtype=np.int8)
    total = arb(0)
    for row in rec["proof"]["dual_terms"]:
        w, o = np.array(row["weights"]), np.array(row["optical"])
        frame = frames[row["frame"]]
        if np.any(w < 0) or frame.shape != (acq.n,acq.n):
            return False
        qform = original(acq, ell, upp, w, o)
        native = list(range(acq.n))
        source = list(range(acq.n,acq.n+acq.s))
        K = select(qform, native, native)
        gap, _ = validate_frame(K, frame)
        if not gap > 0:
            return False
        L = select(qform, native, source)
        Schur = select(qform, source, source)-bt(L)*K.solve(L)
        if row["source_mode"] is None:
            value = sum((Schur[j,j].real for j in range(acq.s)), arb(0))
        else:
            mat, scale = (identity, 1) if row["source_basis"] == "I" else (H, acq.s)
            q = row["source_mode"]
            direction = ball(mat[:,q:q+1])
            value = (bt(direction)*Schur*direction)[0,0].real/scale
        total += max(arb(0),value.lower())
    return bool(total > arb(acq.radius)**2)


class TraceElimination:
    def __init__(self, acq):
        T = np.vstack((np.eye(acq.n)[REGION], acq.G[REGION]))
        d = len(T)
        _, _, permutation = qr(T, pivoting=True)
        pivot, free = permutation[:d], permutation[d:]
        n, s = acq.n, acq.s
        exactH = ball(T[:,pivot]).solve(ball(np.eye(d)))
        inner = -(exactH*ball(T[:,free]))
        U = acb_mat(n+s,d+s)
        W = acb_mat(n+s,n-d)
        for slot, k in enumerate(pivot):
            for j in range(d):
                U[int(k),j] = exactH[slot,j]
            for j in range(n-d):
                W[int(k),j] = inner[slot,j]
        for j,k in enumerate(free):
            W[int(k),j] = 1
        for j in range(s):
            U[n+j,d+j] = 1
        self.U, self.W = U,W
        self.positive_cache = {}
        self.n = n
        self.d = d

    def reduce(self, A, key):
        H = bt(self.W)*A*self.W
        if key not in self.positive_cache:
            self.positive_cache[key] = strict_pd(H)[0]
        if not self.positive_cache[key]:
            raise ArithmeticError("replay: exact eliminated null part not positive")
        Q = bt(self.W)*A*self.U
        return bt(self.U)*A*self.U-bt(Q)*H.solve(Q)


def leaf_cover(leaves, start, end):
    """Verify recursively that the output covers the actual requested closed box."""
    by_path = {leaf["path"]: leaf for leaf in leaves}
    if len(by_path) != len(leaves):
        return False
    visited = set()
    maxlen = max(map(len, by_path), default=0)

    def children(prefix):
        return [leaf for path, leaf in by_path.items() if path.startswith(prefix)]

    def walk(path, l, u):
        if path in by_path:
            leaf = by_path[path]
            if not np.array_equal(leaf["box"],(l,u)) or "verified" not in leaf:
                return False
            visited.add(path)
            return True
        if len(path) >= maxlen:
            return False
        left, right = children(path+"L"), children(path+"R")
        if not left or not right:
            return False
        left_max = np.max([leaf["box"][1] for leaf in left], axis=0)
        right_min = np.min([leaf["box"][0] for leaf in right], axis=0)
        changes = np.flatnonzero(left_max < u)
        if len(changes) != 1:
            return False
        index = changes[0]
        cut = left_max[index]
        if not l[index] < cut < u[index] or right_min[index] != cut:
            return False
        expected_l, expected_u = l.copy(), u.copy()
        expected_l[index] = expected_u[index] = cut
        if not np.array_equal(left_max,expected_u) or not np.array_equal(right_min,expected_l):
            return False
        return walk(path+"L",l,expected_u) and walk(path+"R",expected_l,u)
    return walk("root", start,end) and len(visited)==len(leaves)


def replay_chart(acq, matrix, exact_leaf, selected, frames):
    d = 2*len(REGION)
    free = [i for i in range(d) if i not in selected]
    states = len(free)
    assert states == exact_leaf["verified"]["states"]
    lower, upper = [[arb(float(v)) for v in edge] for edge in exact_leaf["box"]]
    center, radius = [(a+b)/2 for a,b in zip(lower,upper)], [(b-a)/2 for a,b in zip(lower,upper)]
    full = acb_mat(d+acq.s, states+acq.s)
    for j,i in enumerate(free):
        full[i,j] = 1
    for j in range(acq.s):
        full[d+j,states+j] = 1
    field = {}
    for i, c in zip(selected,center,strict=True):
        direction = acb_mat(1,states+acq.s)
        direction[0,free.index(len(REGION)+i)] = 1
        for j in range(acq.s):
            direction[0,states+j] = acq.Bb[int(REGION[i]),j]
        field[i] = direction
        for j in range(states+acq.s):
            full[i,j] = c*direction[0,j]
    square_lower = arb(0)
    curvature = select(matrix, selected, selected)
    quadratic_positive = strict_pd(curvature)[0]
    for term in exact_leaf["verified"]["terms"]:
        F = ball(frames[term["frame"]])
        if F.nrows() != states+acq.s or F.ncols()!=states+term["j"]:
            return False
        D = [field[i]*F for i in selected]
        Ft = full*F
        shifted_total = shift(matrix,d,term["threshold"])
        center_pullback = bt(Ft)*shifted_total*Ft
        derivative = []
        for i, di in zip(selected,D,strict=True):
            b = bt(Ft)*select(shifted_total, range(d+acq.s), [i])*di
            derivative.append(b+bt(b))
        quad = arb(0) if quadratic_positive else frob(curvature)*sum(
            (abs(delta)**2*frob(di)**2 for delta,di in zip(radius,D,strict=True)),arb(0))
        for i in range(F.ncols()):
            margin = center_pullback[i,i].real-quad
            for j in range(F.ncols()):
                if i!=j:
                    margin -= abs(center_pullback[i,j])
                margin -= sum((delta*abs(dq[i,j]) for delta,dq in zip(radius,derivative,strict=True)),arb(0))
            if not margin > 0:
                return False
        square_lower += arb(term["threshold"])*term["term_multiplicity"]
    return bool(square_lower > arb(acq.radius)**2)


def prove_endpoint_divergence(acq, ell, upp):
    intervals = [(arb(float(a)),arb(float(u))) for a,u in zip(ell,upp,strict=True)]
    center = [(a+u)/2 for a,u in intervals]
    widths = [(u-a)/2 for a,u in intervals]
    A = ball(np.eye(acq.n))-bdiag(center,acq.Gb)
    T = A.transpose().solve(acq.Gb.transpose()).transpose()
    work = bdiag(widths,T.transpose()).transpose()
    for power in (2,4,8,16):
        work = work*work
        total = abs(sum((work[i,i] for i in range(acq.n)),acb(0)))/acq.n
        if total > 1:
            return {"power":power,"rho_lower":low((total.log()/power).exp()),
                    "legal_material_endpoint":upp.tolist()}
    return None


def inertia_count(matrix):
    spectrum = eigvalsh(midpoint(matrix))
    negative = int(np.count_nonzero(spectrum < 0.))
    positive = int(np.count_nonzero(spectrum > 0.))
    exact = []
    for sign,count in ((1,positive),(-1,negative)):
        if not count:
            exact.append(0)
            continue
        sided = sign*matrix
        frame,_ = propose_frame(midpoint(sided),count)
        if frame is None or not validate_frame(sided,frame)[0] > 0:
            return None
        exact.append(count)
    if sum(exact)==matrix.nrows():
        return exact[0], exact[1], 0
    return None


def check_pair(a, raw, pair, *, audit_inertia=False):
    path = raw / f"pair_{pair+1:02d}.json"
    result = json.loads(path.read_text())
    frames = load(raw / f"pair_{pair+1:02d}_arrays.npz")
    checks, failures, divergence = Counter(), [], []
    seconds = Counter()
    signed_counts = []
    current = None
    start = perf_counter()
    for row in result["rows"]:
        target,bg,k,slot = row["box_positions"]
        if current != (target,bg):
            current = target,bg
            acq = Acquisition(a["G"],a["S"],a["B"],a["Y"][pair,target,bg],float(a["radii"][pair,target,bg]))
            elimination = TraceElimination(acq)
        ell,upp = a["boxes"][pair,bg,k,slot]
        if row["witness"]["audit"]["proved"]:
            t0 = perf_counter()
            checks["witness"] += 1
            if not acq.strict_witness(frames[row["witness"]["material"]],ell,upp)["proved"]:
                failures.append([*row["box_positions"],"witness"])
            seconds["witness"] += perf_counter()-t0
        if row["methods"]["A"].get("positive_endpoint_perturbation_spectral_radius_float",0.) > 1.05 and len(divergence)<8:
            receipt = prove_endpoint_divergence(acq,ell,upp)
            if receipt:
                divergence.append(dict(box_positions=row["box_positions"],**receipt))
        for method in "ABCDE":
            if not row["methods"][method]["excluded"]:
                continue
            checks[method] += 1
            t0 = perf_counter()
            rec = row["methods"][method]
            try:
                if method=="A":
                    check = replay_neumann(acq,ell,upp,rec,frames[rec["receiver_direction"]])
                elif method=="B":
                    check = independent_psd(acq,ell,upp,rec,frames)
                elif method in "CD":
                    proof = rec["proof"]
                    w,o = weight_vector(acq.n,*proof["weights"])
                    matrix = original(acq,ell,upp,w,o)
                    if method=="D":
                        matrix = elimination.reduce(matrix,proof["weights"][0])
                    check = replay_terms(matrix,proof,frames,acq.radius)
                    if audit_inertia:
                        cutoff = proof["states"]
                        counts = inertia_count(select(matrix,range(cutoff),range(cutoff)))
                        if counts is not None and method=="D":
                            counts = counts[0]+proof["positive_eliminated_states"],counts[1],counts[2]
                        signed_counts.append(dict(position=row["box_positions"], method=method,
                                                  positive_negative_zero=counts))
                elif rec["stage"]=="inherited_D":
                    check = row["methods"]["D"]["excluded"]
                else:
                    successful = [attempt for attempt in rec["attempts"] if attempt["excluded"]]
                    assert len(successful)==1
                    proof = successful[0]
                    selected = tuple(proof["trace_rows_within_joint_bank"])
                    els,ups = ell[REGION[list(selected)]],upp[REGION[list(selected)]]
                    check = leaf_cover(proof["leaves"],els,ups)
                    w,o = weight_vector(acq.n,*proof["searched_weights"],skipped=selected)
                    reduced = elimination.reduce(original(acq,ell,upp,w,o),proof["searched_weights"][0])
                    check &= all(replay_chart(acq,reduced,leaf,selected,frames) for leaf in proof["leaves"])
                    checks["E_new"] += 1
            except (ArithmeticError, ValueError, ZeroDivisionError, AssertionError) as error:
                failures.append(dict(row=row["box_positions"],method=method,exception=str(error)))
            else:
                if not check:
                    failures.append(dict(row=row["box_positions"],method=method))
            seconds[method] += perf_counter()-t0
    return dict(pair=pair+1, pass_all=not failures, checked=dict(checks), failures=failures,
                confirmed_positive_endpoint_divergence_examples=divergence,
                seconds=perf_counter()-start, seconds_by_method=dict(seconds), precision_bits=ctx.prec,
                signed_state_inertia=signed_counts,
                replay="reads sealed original arrays, weights and saved frames; reassembles full row quadratic; independently eliminates T-null complement")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--begin",type=int,default=1)
    parser.add_argument("--end",type=int,default=24)
    parser.add_argument("--audit-inertia",action="store_true")
    args = parser.parse_args()
    seal = json.loads(SEAL.read_text())
    assert digest(PUBLIC)==seal["files"][PUBLIC.name]
    a = load(PUBLIC)
    raw = RESULTS / "raw"
    destination = RESULTS / ("replay_inertia" if args.audit_inertia else "replay")
    for pair in range(args.begin-1,args.end):
        output = destination / f"pair_{pair+1:02d}.json"
        if output.exists():
            continue
        row = check_pair(a,raw,pair,audit_inertia=args.audit_inertia)
        write_json(output,row)
        print(f"pair {pair+1:02d}: {row['checked']}, pass {row['pass_all']}, "
              f"{row['seconds']:.2f} s",flush=True)


if __name__=="__main__":
    main()
