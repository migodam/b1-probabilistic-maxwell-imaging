"""Material-coordinate sets from physical linear combinations of the measured waves."""
import json
from pathlib import Path
import os
import time

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
from scipy.linalg import null_space

ROOT = Path(__file__).resolve().parents[1]
BOX = (.5, 12.)
CENTER = sum(BOX) / 2
HALF = (BOX[1] - BOX[0]) / 2
EPS = 2e-4
GAIN = .01
DEPTH = 13


def norm2(array):
    return np.sum(abs(array)**2, axis=-1)


def square_polynomial(c):
    """Squared norm of vector coefficients in 1,x,p,xp, over real x,p."""
    out = np.zeros((len(c), 3, 3))
    for a, b in np.ndindex(2, 2):
        for e, f in np.ndindex(2, 2):
            out[:, a+e, b+f] += np.einsum("jk,jk->j", c[:, a, b].conj(), c[:, e, f]).real
    return out


def pack_runs(mask, depth=DEPTH):
    cuts = np.flatnonzero(np.diff(np.r_[False, mask, False].astype(int)))
    return (BOX[0]+(BOX[1]-BOX[0])*cuts.reshape(-1, 2)/2**depth).tolist()


def coefficients(H, R, B, host, reads):
    inv = np.linalg.pinv(B)
    U = np.stack([np.eye(len(R))-CENTER*H, -HALF*H])
    X = np.moveaxis(inv @ U, 1, 2)  # degree, cell, measured waveform
    predict = X @ host + np.stack([CENTER*R, HALF*R])
    residual = X @ reads - predict
    null = null_space(B.T)
    nU = (null[:, 0] @ U) if null.shape[1] == 1 else None
    return X, predict, residual, nU


def pair_bilinears(A, nu, i):
    # X_i(q)*b_j(p) - X_j(p)*b_i(q). Leading axis in the output indexes j.
    other = np.arange(nu.shape[1]) != i
    return (np.einsum("ac,bj->jabc", A[:, i], nu[:, other])
            - np.einsum("bjc,a->jabc", A[:, other], nu[:, i]))


def pair_mask(X, predicted, residual, nu, i, eps, gamma, depth):
    rc = pair_bilinears(residual[..., None], nu, i)
    noise = np.concatenate((eps*pair_bilinears(X, nu, i),
                            gamma*pair_bilinears(predicted[..., None], nu, i)), axis=-1)
    denominator = np.maximum(np.max(abs(rc), axis=(1,2,3)), np.max(abs(noise), axis=(1,2,3)))
    denominator = np.maximum(denominator, np.finfo(float).tiny)
    rc = rc/denominator[:, None, None, None]
    noise = noise/denominator[:, None, None, None]
    sqr, sqv = square_polynomial(rc), square_polynomial(noise)
    poly = sqr-sqv
    # If some x in a bin is feasible, at its center |r| <= ||v||+h*L;
    # in p the two derivative norms attain their upper bounds at endpoints.
    lbound = np.maximum(abs(rc[:,1,0,0]+rc[:,1,1,0]), abs(rc[:,1,0,0]-rc[:,1,1,0]))
    lbound += np.sqrt(np.maximum(norm2(noise[:,1,0]+noise[:,1,1]),norm2(noise[:,1,0]-noise[:,1,1])))
    pad = 256*np.finfo(float).eps*np.sum(abs(sqr)+abs(sqv), axis=(1,2))
    keep = np.array([0])
    for level in range(1, depth+1):
        keep = (2*keep[:, None] + np.array([0,1])).ravel()
        if not len(keep):
            break
        h = 2.**-level
        x = -1+(2*keep+1)*h
        A0 = poly[:,0,None,:]+x[None,:,None]*(poly[:,1,None,:]+x[None,:,None]*poly[:,2,None,:])
        C, D, E = np.moveaxis(A0, 2, 0)
        lowest = C+E-abs(D)
        with np.errstate(divide="ignore", invalid="ignore"):
            stationary = np.clip(-D/(2*E), -1., 1.)
        interior = (E>0) & np.isfinite(stationary)
        safe_p = np.where(interior, stationary, 0.)
        lowest = np.where(interior, np.minimum(lowest, C+safe_p*(D+safe_p*E)), lowest)
        vn = sqv[:,0,None,:]+x[None,:,None]*(sqv[:,1,None,:]+x[None,:,None]*sqv[:,2,None,:])
        v_max = np.sqrt(np.maximum(0.,vn[:,:,0]+vn[:,:,2]+abs(vn[:,:,1])+pad[:,None]))
        delta = h*lbound[:,None]
        keep = keep[np.all(lowest <= delta*(2*v_max+delta)+pad[:,None], axis=0)]
    mask = np.zeros(2**depth, dtype=bool)
    mask[keep] = True
    return mask


def full_mask(X, predicted, residual, i, eps, gamma, depth):
    x = -1+(2*np.arange(2**depth)+1)*2.**-depth
    r0, r1 = residual[:, i]
    n = np.column_stack((eps*X[:, i], gamma*predicted[:, i]))
    mid_r = abs(r0+r1*x)
    mid_v = np.sqrt(np.maximum(0., norm2(n[0])+2*x*np.vdot(n[0],n[1]).real+x*x*norm2(n[1])))
    return mid_r <= mid_v+2.**-depth*(abs(r1)+np.linalg.norm(n[1]))+1e-13*(abs(r0)+abs(r1)+np.linalg.norm(n))


def old_sets(H,R,B,hy,yd,eps,gamma):
    X,p,r,n = coefficients(H,R,B,hy,yd)
    return np.array([pair_mask(X,p,r,n,i,eps,gamma,DEPTH) for i in range(len(R))])


def extended_sets(H,R,B,hy,yd,eps,gamma,old):
    X,p,r,_ = coefficients(H,R,B,hy,yd)
    return old & np.array([full_mask(X,p,r,i,eps,gamma,DEPTH) for i in range(len(R))])


def stats(mask, truth):
    bins = np.clip(((truth-BOX[0])/(BOX[1]-BOX[0])*2**DEPTH).astype(int),0,2**DEPTH-1)
    intervals = [pack_runs(m) for m in mask]
    widths = [max(b for a,b in v)-min(a for a,b in v) if v else None for v in intervals]
    second = [sum(b-a for a,b in v) for v in intervals]
    return dict(complete_truth_in=bool(np.all(mask[np.arange(len(bins)),bins])),
                coords_in=int(sum(mask[np.arange(len(bins)),bins])),
                median_hull_width=None if None in widths else float(np.median(widths)),
                max_hull_width=None if None in widths else max(widths),
                mean_total_length=float(np.mean(second)), intervals=intervals)


def source_choices(D, rng):
    full = D["observable_modes"]
    B = D["B_old"]
    P = D["q"]@full
    n = null_space(B.T)[:,0]
    a = n@P
    old_part = np.linalg.pinv(B)@P
    x = np.linalg.solve(np.eye(P.shape[1])+old_part.conj().T@old_part, a.conj())
    x /= np.linalg.norm(x)
    matched = a.conj()/np.linalg.norm(a)
    random = rng.normal(size=len(a))+1j*rng.normal(size=len(a))
    random /= np.linalg.norm(random)
    ctl = D["p_old"][:,0]
    def describe(z):
        return dict(gap_gain=float(abs(a@z)),
                    old_leakage=float(np.linalg.norm(old_part@z)),
                    inverse_noise_gain=float(np.hypot(1.,np.linalg.norm(old_part@z))/abs(a@z)))
    return dict(gap=full@matched, noise=full@x, random=full@random),dict(
        gap=describe(matched), noise=describe(x), random=describe(random),
        old_repeat_gain=float(abs(n@D["q"]@ctl)),noise_gap_alignment=float(abs(np.vdot(x,matched))))


def noisy_vector(rng, size):
    v = rng.normal(size=size)+1j*rng.normal(size=size)
    return v/np.linalg.norm(v)


def main():
    rows, geometry = [], {}
    start = time.process_time()
    for k in (12,24):
        D = np.load(ROOT/"results"/f"lift_k{k}.npz")
        H,R,B,D0 = D["H"],D["R"],D["B_old"],D["y_host"]
        old_feeds = D["p_old"]
        host_old = D0@old_feeds
        sources, details = source_choices(D, np.random.default_rng(8871+k))
        geometry[k] = details
        truth = np.vstack((D["pair_materials"],D["heldout_materials"]))
        for sample,(z,delta) in enumerate(zip(truth,D["independent_all_source_data"])):
            rng = np.random.default_rng(394760+100*k+sample)
            total = D0+delta
            noise = noisy_vector(rng, k)*EPS*.55
            gerr = .8*GAIN*np.exp(2j*np.pi*rng.uniform())
            data = (1+gerr)*(total@old_feeds)+noise[:-1]
            # One fixed bounded receiver gain is shared by old and new readings.
            before = old_sets(H,R,B,host_old,data,EPS,GAIN)
            assumed_exact_gain = old_sets(H,R,B,host_old,data,EPS,0.)
            row = dict(k=k, sample=sample, truth=z.tolist(), gain_error=[gerr.real,gerr.imag],
                       joint_error_budget_used=float(np.hypot(.55,.8)),
                       old=stats(before,z), old_ignore_gain=stats(assumed_exact_gain,z))
            for name,a in sources.items():
                source_matrix = np.column_stack((B,D["q"]@a))
                h = np.r_[host_old,D0@a]
                y = np.r_[data,(1+gerr)*(total@a)+noise[-1]]
                m = extended_sets(H,R,source_matrix,h,y,EPS,GAIN,before)
                row[name] = stats(m,z)
            rows.append(row)
        print(k,details,[(key,sum(v[key]['complete_truth_in'] for v in rows if v['k']==k),
                           np.median([v[key]['mean_total_length'] for v in rows if v['k']==k]))
                           for key in ('old','old_ignore_gain','gap','noise','random')],flush=True)
    out=dict(box=BOX, joint_additive_scale=EPS, joint_gain_scale=GAIN, material_bin_depth=DEPTH,
             split="First two states are the designed pair. Other twelve were not used to design the array.",
             cpu_seconds=time.process_time()-start,source_metrics=geometry, rows=rows)
    (ROOT/"results"/"noise_sets.json").write_text(json.dumps(out,ensure_ascii=False,separators=(',',':'))+'\n')


if __name__ == '__main__':
    main()

