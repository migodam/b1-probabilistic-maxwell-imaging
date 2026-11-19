"""Certify the whole two-cell box from shared broadband observations."""
import argparse
from collections import Counter, deque
import json
from time import perf_counter, process_time

import numpy as np
from flint import acb_mat, arb
from scipy.stats import chi2

from geometry_controls import coherent_mirror_block, CENTER, Mirror, V
from geometry_experiment import Block, FullTM, ROOT, SIGMA, material_rms
from geometry_bounds import centers, cball, frobenius, lower, upper

OUT = ROOT/"results"
BAND = (-.08,.012,.08)
ARMS = {
    "symmetric_broadband": (False,BAND),
    "asymmetric_broadband": (True,BAND),
    "asymmetric_same_frequency": (True,(0.,0.,0.)),
}
BOX = (.5,30.)
GRID_DEPTH = 32
VISIT_CAP = 90_000


def model_blocks(arms):
    equal = Mirror([0,0])
    old = coherent_mirror_block(equal,0.,1.)
    asymmetric, offsets = arms
    blocks = [old]
    for d in offsets:
        block = coherent_mirror_block(equal,d,1/3)
        if asymmetric:
            _,b,*_ = equal.wave.kernels(d)
            # One left line source per pair: power is 1 on each exposure
            # before the same sqrt(1/3) duration/energy weighting.
            block = Block(block.green,b[:,::2]*np.sqrt(1/3),block.receiver,block.host,block.active)
        blocks.append(block)
    return blocks


def acquire():
    material = CENTER+.6*V
    rng = np.random.default_rng(14_250819)
    noise = SIGMA/np.sqrt(2)*(rng.normal(size=24)+1j*rng.normal(size=24))
    epsilon = SIGMA*np.sqrt(chi2.ppf(.99,2*len(noise))/2)
    archive = dict(epsilon=epsilon)
    for name, arm in ARMS.items():
        model = FullTM(model_blocks(arm))
        y,*_ = model.evaluate(material)
        archive[name] = y+noise
    np.savez_compressed(OUT/"geometry-mirror-box-reads.npz",**archive)


def setup(blocks, readings):
    midpoint = np.full(2,np.mean(BOX))
    center = list(centers(blocks,midpoint))
    pieces = []
    start = 0
    for c in center:
        size = c.y.nrows()
        demeaned = [c.y[i,0]-cball(y) for i,y in enumerate(readings[start:start+size])]
        coefficients = [[c.a[r,i]*c.e[j,s] for i in range(2) for j in range(2)]
                        for s in range(c.e.ncols()) for r in range(c.a.nrows())]
        pieces.append((c.k,coefficients,demeaned))
        start += size
    return midpoint[0], pieces


def domain(lo,hi,midpoint):
    return [arb((a+b)/2-midpoint,(b-a)/2) for a,b in zip(lo,hi,strict=True)]


def residual_bounds(x,pieces):
    z1,z2 = x
    for k,coefficients,center_residual in pieces:
        a,d = 1-k[0,0]*z1,1-k[1,1]*z2
        b,c = k[0,1]*z1,k[1,0]*z2
        det = a*d-b*c
        if not abs(det)>0:
            yield None
            continue
        numerator = (z1*d,b*z2,c*z1,z2*a)
        yield [(sum(u*v for u,v in zip(row,numerator,strict=True))+base*det)/det
               for row,base in zip(coefficients,center_residual,strict=True)]


def interval_lower_joint(residual, y, eps, gamma):
    """Optional shared-gain lower bound at a material box, with gain disk enlarged."""
    n = len(y)
    penalty = arb(eps)/arb(gamma)
    y_norm_squared = sum(abs(z)**2 for z in y)+penalty**2
    bp = sum(y[i].conjugate()*residual[i] for i in range(n))/y_norm_squared
    norm_perp_squared = sum(abs(residual[i]-bp*y[i])**2 for i in range(n)) + abs(penalty*bp)**2
    prediction_norm_squared = sum(abs(y[i]+residual[i])**2 for i in range(n))+penalty**2
    return norm_perp_squared*y_norm_squared/prediction_norm_squared


def certify(blocks, readings, epsilon, gamma=0.):
    start,cpu = perf_counter(),process_time()
    midpoint, pieces = setup(blocks,readings)
    ya = [cball(x) for x in readings]
    budget = arb(epsilon)**2
    raw_budget = budget
    if gamma:
        # Feasibility bounds the uncalibrated residual independently of the box.
        obs_norm = sum(abs(y)**2 for y in ya).sqrt()
        raw_budget = (arb(epsilon)+arb(gamma)*(obs_norm+arb(epsilon))/(1-arb(gamma)))**2
    todo = deque([((BOX[0],BOX[0]),(BOX[1],BOX[1]),0)])
    keep = []
    counts = Counter()
    while todo and counts["visited"]<VISIT_CAP:
        lo,hi,depth = todo.popleft()
        counts["visited"] += 1
        pieces_used = list(residual_bounds(domain(lo,hi,midpoint),pieces))
        counts["block_bounds"] += len(pieces_used)
        bounds = [abs(z)**2 for r in pieces_used if r is not None for z in r]
        reject = sum(bounds)>raw_budget
        if gamma and not reject and all(r is not None for r in pieces_used):
            reject = interval_lower_joint(sum(pieces_used,[]),ya,epsilon,gamma)>budget
        if reject:
            counts["excluded"] += 1
        elif depth == GRID_DEPTH:
            keep.append((lo,hi,depth))
        else:
            axis = int(hi[1]-lo[1] > hi[0]-lo[0])
            center = (lo[axis]+hi[axis])/2
            hi_left, lo_right = list(hi), list(lo)
            hi_left[axis], lo_right[axis] = center, center
            todo.append((tuple(lo_right),hi,depth+1))
            todo.append((lo,tuple(hi_left),depth+1))
    counts["terminal_kept"] = len(keep)
    counts["unprocessed_kept"] = len(todo)
    keep.extend(todo)
    if keep:
        vlo = np.min([q[0] for q in keep],axis=0)
        vhi = np.max([q[1] for q in keep],axis=0)
        # Outward upper bound on Euclidean RMS diameter via the hull.
        gap = [arb(h)-arb(l) for h,l in zip(vhi,vlo,strict=True)]
        diameter = upper((sum(w*w for w in gap)/2).sqrt())
        halfspace = ("m1>m2" if all(q[0][0]>q[1][1] for q in keep) else
                     "m2>m1" if all(q[0][1]>q[1][0] for q in keep) else "unresolved")
    else:
        vlo, vhi, diameter, halfspace = [],[],0.,"empty"
    return dict(cells=keep,hull=(list(vlo),list(vhi)),rms_diameter_upper=diameter,
                cert_direction=halfspace, counts=dict(counts), cpu_s=process_time()-cpu, wall_s=perf_counter()-start)


def analyze(shared_gain):
    with np.load(OUT/"geometry-mirror-box-reads.npz",allow_pickle=False) as f:
        epsilon = float(f["epsilon"])
        gamma = 0.
        if shared_gain:
            epsilon *= np.sqrt(2)
            gamma = np.sqrt(2)*4e-4
        result = dict(epsilon=epsilon,gamma=gamma,full_box=BOX,depth=GRID_DEPTH,
                      center_inverse_factorizations_per_arm=4,new_total_energy_per_feed_pair=1.,
                      arms={})
        for name,arm in ARMS.items():
            row = certify(model_blocks(arm),f[name],epsilon,gamma)
            row.update(offsets=arm[1],power_each_new_pair=[1/3]*3)
            result["arms"][name] = row
            print(name,row["counts"],row["hull"],row["rms_diameter_upper"],row["cert_direction"],flush=True)
    tag = "-gain" if shared_gain else ""
    save_outer_cover(result, tag)


def save_outer_cover(result, tag):
    cells = {}
    for name,row in result["arms"].items():
        boxes = row.pop("cells")
        cells[name] = np.array([[lo,hi] for lo,hi,_ in boxes]).reshape(-1,2,2)
        cells[name+"_depth"] = np.array([depth for _,_,depth in boxes],dtype=np.uint8)
    result["cells_file"] = f"geometry-mirror-box{tag}-cells.npz"
    np.savez_compressed(OUT/result["cells_file"], **cells)
    (OUT/f"geometry-mirror-box{tag}.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")


def fixed_score(blocks, point, observed, epsilon, gamma, gain):
    c = list(centers(blocks,point))
    y = acb_mat([[v.y[i,0]] for v in c for i in range(v.y.nrows())])
    observed = acb_mat([[cball(value)] for value in observed])
    residual = y*cball(gain)-observed
    additive = frobenius(residual)/arb(epsilon)
    penalty = abs(cball(gain)-1)/arb(gamma) if gamma else arb(0)
    return upper((additive**2+penalty**2).sqrt())


def score(shared_gain):
    tag = "-gain" if shared_gain else ""
    sealed = json.loads((OUT/f"geometry-mirror-box{tag}.json").read_text())
    material = CENTER+.6*V
    partner = material[::-1]
    with (np.load(OUT/"geometry-mirror-box-reads.npz",allow_pickle=False) as f,
          np.load(OUT/sealed["cells_file"],allow_pickle=False) as outer):
        gap = arb(float(material[0]))-arb(float(material[1]))
        scored = dict(pair_rms=material_rms(material,partner),pair_rms_lower=lower(gap),
                      truth=material.tolist(),other=partner.tolist(),arms={})
        for name,arm in ARMS.items():
            boxes = outer[name]
            pts = [bool(np.any(np.all((boxes[:,0,:]<=x)&(x<=boxes[:,1,:]),axis=1)))
                   for x in (material,partner)]
            blocks = model_blocks(arm)
            wave = FullTM(blocks)
            eps, gamma = sealed["epsilon"],sealed["gamma"]
            observed = f[name]
            profile = []
            validated_bound = []
            for x in (material,partner):
                y = wave.evaluate(x)[0]
                g = 1. if not gamma else (np.vdot(y,observed)+(eps/gamma)**2)/(np.vdot(y,y).real+(eps/gamma)**2)
                q = np.linalg.norm(g*y-observed)/eps
                profile.append(float(np.hypot(q,abs(g-1)/gamma)) if gamma else float(q))
                validated_bound.append(fixed_score(blocks,x,observed,eps,gamma,g))
            scored["arms"][name] = dict(truth_in_cover=pts[0],exchange_in_cover=pts[1],
                                        profile_scores=profile,fixed_gain_witness_score_upper=validated_bound)
        (OUT/f"geometry-mirror-box{tag}-scores.json").write_text(json.dumps(scored,indent=2,allow_nan=False)+"\n")
        print(scored,flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage",choices=("acquire","certify","score"))
    p.add_argument("--shared-gain",action="store_true")
    args = p.parse_args()
    if args.stage=="acquire":
        acquire()
    elif args.stage=="certify":
        analyze(args.shared_gain)
    else:
        score(args.shared_gain)
