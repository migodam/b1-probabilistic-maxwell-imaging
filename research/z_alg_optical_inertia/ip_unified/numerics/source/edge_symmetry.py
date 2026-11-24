"""Prove the geometric edge-weight inertia using exact reflection-sector enclosures."""
import argparse
from pathlib import Path
import sys

from flint import acb, arb, ctx
import numpy as np
from scipy.linalg import cholesky, eigh
from threadpoolctl import threadpool_limits

from edge_core import HERE, MULTIPLIERS, clock, configure, herm, load_json, profile, save_json, times
from native import high, low

sys.path.insert(0,str(HERE.parents[2]/"synthetic"))
from nonvacuity import ExactFourierGeometry, binary_frob_ceil, gamma


def charge_norm(array):
    return binary_frob_ceil(array)


def max_coeff(array):
    return charge_norm(array)/arb(array.size).sqrt()


def exact_symmetric_kernel(fine, n):
    p = fine.pad
    modes = range(1-n,n)
    exact = fine.exact
    exact_sym, asymmetry = {}, arb(0)
    values = np.zeros((p,p),complex)
    for dy in modes:
        for dx in modes:
            idx = (dy%p)*p+(dx%p)
            sdy,sx = abs(dy),abs(dx)
            key = sdy,sx
            if key not in exact_sym:
                sym = sum((exact[(sy*sdy%p)*p+(tx*sx%p)]
                           for sy,tx in ((1,1),(1,-1),(-1,1),(-1,-1))),0)/4
                mid = complex(sym)
                delta = abs(sym - acb(mid.real,mid.imag)).upper()
                exact_sym[key] = sym,mid,delta
            sym,mid,_ = exact_sym[key]
            values[dy%p,dx%p] = mid
            asymmetry = max(asymmetry,abs(exact[idx]-sym).upper())
    return values, asymmetry, max(v[2] for v in exact_sym.values())


def reflection_data(g):
    start = clock()
    side, p = g.side,g.pad
    half = side//2
    local = np.arange(half)
    yy,xx = np.meshgrid(local,local,indexing="ij")
    orbit = tuple((side-1-yy if sy else yy,side-1-xx if sx else xx)
                  for sy,sx in ((0,0),(0,1),(1,0),(1,1)))
    orbit = tuple((y.ravel(),x.ravel()) for y,x in orbit)
    basis_rows = tuple(y*side+x for y,x in orbit)
    exact = ExactFourierGeometry(HERE.parents[2]/"synthetic/public/fine_geometry.npz")
    kernel, asymmetry, rounding = exact_symmetric_kernel(exact,side)
    del exact
    ix = np.arange(side)
    rec = np.arange(g.sensor.shape[0])
    receiver_defects = []
    for sy,sx in ((0,1),(1,0),(1,1)):
        y = side-1-ix if sy else ix
        x = side-1-ix if sx else ix
        permutation = (y[:,None]*side+x[None,:]).ravel()
        if sy and sx:
            angles = rec+len(rec)//2
        elif sy:
            angles = -rec
        else:
            angles = len(rec)//2-rec
        sq = sum((arb(float(a.real)-float(b.real))**2 + arb(float(a.imag)-float(b.imag))**2
                  for a,b in zip(g.sensor[:,permutation].flat,g.sensor[angles%len(rec),:].flat,strict=True)),arb(0))
        # Explicitly cover binary64 subtraction even for nearly identical ports.
        size = arb(g.sensor.size).sqrt()
        receiver_defects.append((sq.sqrt()+4*gamma(4)*size*max_coeff(g.sensor)).upper())
    largest_defect = max(receiver_defects)
    sensors = [g.sensor[:,ids] for ids in basis_rows]
    kernel_blocks = []
    ay,ax = orbit[0]
    for by,bx in orbit:
        kernel_blocks.append(kernel[(ay[:,None]-by[None,:])%p,(ax[:,None]-bx[None,:])%p])
    return sensors,kernel_blocks,basis_rows[0],(asymmetry,rounding,largest_defect),times(*start)


def sector_matrices(g, w, shared, signs):
    """Build one of four parity blocks and bound its representation defect."""
    sensors,kernels,active,defects,_ = shared
    _,kernel_error,_ = defects
    k = len(active)
    r = sensors[0].shape[0]
    G = kernels[0].copy()
    S = sensors[0].copy()
    for sign,g_block,s_block in zip(signs[1:],kernels[1:],sensors[1:],strict=True):
        G += sign*g_block
        S += sign*s_block
    S *= .5
    omega = np.asarray(w.real)[g.blocks[active]]/w.scale
    z = (omega+1j*w.optical[0]/w.scale)/2
    delta_weights = max(abs(arb(float(v))/arb(w.scale)-arb(float(v/w.scale))).upper() for v in set(w.real))
    delta_z = delta_weights/2 + abs(arb(w.optical[0])/arb(w.scale)-arb(w.optical[0]/w.scale))/2
    nG,nS = charge_norm(G),charge_norm(S)
    green_error = k*(4*kernel_error+8*gamma(18)*max(max_coeff(x) for x in kernels))
    port_error = arb(k*r).sqrt()*(4*gamma(18)*max(max_coeff(x) for x in sensors))
    z_norm = max_coeff(z)
    ss_error = (gamma(8*r+32)*4*nS*nS+2*nS*port_error+port_error**2)/arb(w.scale)
    dg_error = 2*((z_norm+delta_z)*green_error+delta_z*nG)
    round_guard = gamma(40)*(arb(k).sqrt()*max_coeff(omega)+nS*nS/arb(w.scale)+4*z_norm*nG)
    guard = (arb(k).sqrt()*delta_weights+ss_error+dg_error+round_guard).upper()
    A = S.conj().T@S/w.scale
    pg = z[:,None]*G
    A -= pg
    A -= pg.conj().T
    A.flat[::k+1] += omega
    A = herm(A)
    return A,guard


def positive_outside_few(A, guard, asymmetry, gap):
    """Verify PSD after a numerical low-rank update and a declared gap."""
    start = clock()
    n = len(A)
    eigen,q = eigh(A,subset_by_value=(-np.inf,0.),check_finite=False,driver="evr")
    W = q*np.sqrt(-2*eigen)
    Wn = charge_norm(W) if W.shape[1] else arb(0)
    L_init = A.copy()
    if W.shape[1]:
        L_init += W@W.conj().T
    L_init.flat[::n+1] -= gap
    L = cholesky(L_init,lower=True,check_finite=False)
    approx = L_init-L@L.conj().T
    defect = charge_norm(approx)
    rounding = gamma(8*n+32)*(4*charge_norm(L)**2+4*Wn**2+charge_norm(A)+arb(n).sqrt())
    # Both this reduced block and the off-block symmetry mismatch are covered.
    positive = arb(gap)-defect-rounding-guard-asymmetry
    if not positive>0:
        raise ArithmeticError("Shifted reduced Cholesky residual did not establish the positive gap.")
    return dict(necessary_max_rank=int(W.shape[1]),negative_eigenvalues=eigen.tolist(),
                positive_gap_lower=low(positive),operand_model_error_upper=high(guard),
                dense_binary64_cholesky_defect_upper=high(defect),
                dense_binary64_guard_upper=high(rounding),**times(*start))


def main(factors, verify_sources, path):
    start = clock()
    geo,native,data,_ = configure("u100_best")
    from native import Geometry
    g = Geometry(geo.spec,geo.sensor,geo.incident,block_side=4)
    del geo
    work = reflection_data(g)
    asym, _, sensor = work[3]
    mu, optical = native.w.scale,native.w.optical[0]
    out = load_json(path) if path.exists() else dict(cases=[],reflection_preparation=work[-1])
    out["cases"] = [row for row in out["cases"] if row["edge_multiplier"] not in factors]
    numerical = load_json(HERE/"edge_core.json")
    for factor in factors:
        started = clock()
        w = profile(mu,optical,factor)
        z_max = max(abs(arb(x)/arb(mu)/2).upper() for x in w.real)+arb(1)
        asymmetry = (2*z_max*g.n*asym+2*charge_norm(g.sensor)*sensor/arb(mu)).upper()
        row = dict(edge_multiplier=factor,proof_normalization=mu,declared_gap=1./256,
                   full_state_reflection_error_upper=high(asymmetry),
                   kernel_symmetrization_max_upper=high(asym),
                   sensor_reflection_error_F_upper=high(sensor),
                   sector_proofs=[])
        for sx,sy in ((1,1),(-1,1),(1,-1),(-1,-1)):
            ps = clock()
            A,guard = sector_matrices(g,w,work,(1,sx,sy,sx*sy))
            pc = positive_outside_few(A,guard,asymmetry,row["declared_gap"])
            row["sector_proofs"].append(dict(reflection_signs=[sx,sy],assembly_and_proof=times(*ps),**pc))
            del A
            print("sector",factor,sx,sy,pc["necessary_max_rank"],pc["positive_gap_lower"],flush=True)
        upper = sum(part["necessary_max_rank"] for part in row["sector_proofs"])
        count = next(case for case in numerical["cases"] if case["edge_multiplier"]==factor)
        lower = count["outward_negative_count_at_least"]
        row.update(inertia_upper=upper,inertia_lower=lower)
        if upper==lower:
            row["exact_inertia"] = upper
            row["positive_gap_lower"] = min(part["positive_gap_lower"] for part in row["sector_proofs"])
        row.update(times(*started))
        print("inertia",factor,lower,upper,row["wall_seconds"],flush=True)
        out["cases"].append(row)
        save_json(path,{**out,**times(*start)})
        if verify_sources and upper==lower:
            from edge_core import analyze_row, reusable_stationarity
            from replay import NativeState
            source_start = clock()
            before = g.cost.copy()
            state = dict(exact_nu=upper,negative_gap_lower=count["outward_negative_gap_lower"],
                         positive_gap_lower=row["positive_gap_lower"])
            engine = NativeState(g,w)
            forms,defects,indices,setup = reusable_stationarity(engine,{32:g.incident})
            row["outward_rows"] = [
                dict(id=meta["id"],radius=meta["padded_radius"],
                     **analyze_row(y,state,mu,forms,defects,indices[32],meta["padded_radius"])) for meta,y in data
            ]
            row["source_proof"] = dict(stationary_setup=setup,green_columns=dict(g.cost-before),**times(*source_start))
            print("source distance",factor,max(r["distance_upper"]/r["radius"] for r in row["outward_rows"]),flush=True)
            save_json(path,{**out,**times(*start)})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--factor", action="append", type=float, choices=MULTIPLIERS)
    parser.add_argument("--source",action="store_true")
    parser.add_argument("--output",type=Path,default=HERE/"edge_symmetry_verified.json")
    args = parser.parse_args()
    ctx.prec,ctx.threads=112,1
    with threadpool_limits(1):
        main(args.factor or list(MULTIPLIERS),args.source,args.output)
