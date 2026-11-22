"""Compare scalar sector and optical multipliers on the frozen EIGOR port masks.

The common Fourier envelope lowers the native operator on the 40x40 grid.
Each numerical finalist is rebuilt with outward-rounded Arb arithmetic.
"""
import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import eigh
from scipy.linalg.blas import dgemm
from threadpoolctl import threadpool_limits

from physical import FREQUENCIES, GROUPS, HERE, PREPARED, SIDE, TARGETS, USED, elapsed, timer

STUDY = HERE.parents[1]
sys.path.insert(0, str(STUDY / "fine"))
sys.path.insert(0, str(STUDY / "core"))
from fourier_bound import gram, grid, h
from verify_fourier import ExactGeometry, hi, lo, sub
from certificate import ball, frobenius_squared, gershgorin, star

CAPS = (.25, .5, .75, 1.)
RATIOS_OF_CAP = (-1., -.5, -.25, -.125, 0., .125, .25, .5, 1., 2., 4., 8.)
MU_DYADIC_POWERS = tuple(range(-22, 5, 2))
CUTOFFS = (7, 10)
ARM_TESTS = {"real": lambda r: r == 0, "positive_optical": lambda r: r > 0,
             "free_optical": lambda r: True}
REC = 24


def coefficient_spec(symbol, q, cap, t):
    z = 1-cap*symbol.real+t*symbol.imag
    outside = np.ones(symbol.shape, bool)
    outside[q[:, 0], q[:, 1]] = False
    base = np.floor(64*.975*np.min(z[outside]))/64
    if base <= 0:
        return None
    return float(base), z[q[:, 0],q[:, 1]]-base


def inverse_gram_numeric(full_gram, coeff, base):
    n = len(coeff)
    active = np.flatnonzero(coeff != 0)
    all_ = np.r_[active,n+np.arange(2*REC)]
    g = full_gram[np.ix_(all_,all_)]
    d = coeff[active]
    r = len(d)
    small = g[:r,:r]/base
    small.flat[::r+1] += 1/d
    ev,v = eigh(small, check_finite=False)
    if min(abs(ev)) <= 1.e-8:
        return None
    ng = int(np.count_nonzero(d < 0) - np.count_nonzero(ev < 0))
    if not 0 <= ng <= r:
        return None
    cross = g[:r,r:]
    coord = v.conj().T@cross
    inv = h(g[r:,r:]/base - (coord.conj().T/ev) @ coord / base**2)
    return inv, ng, r


def numeric_group(inv, base_nu, y_by_target, rx, tx, cap, ratio, scales, cost):
    ob = np.ix_(rx,rx)
    source = np.ix_(REC+tx,REC+tx)
    mixed = np.ix_(rx,REC+tx)
    z = (cap+1j*ratio)/2
    d,v = eigh(inv[ob],check_finite=False)
    cost["receiver_eigenproblems"] += 1
    ds = [v.conj().T@(y-z*inv[mixed]) for y in y_by_target]
    bb = abs(z)**2*inv[source]
    candidates = []
    for mu in scales:
        dr=d+mu
        if min(abs(dr)) < 1e-11:
            continue
        nu=base_nu-int(np.count_nonzero(dr<0))
        if not 0 <= nu < len(tx):
            continue
        for obj, field in enumerate(ds):
            free = h(mu*((field.conj().T/dr)@field - bb))
            eig = eigh(free,eigvals_only=True,check_finite=False)
            cost["source_eigenproblems"] += 1
            chosen = np.maximum(eig[:len(eig)-nu],0.)
            floor = np.sqrt(chosen.sum())
            if floor>1e-9:
                candidates.append((obj,dict(mu=float(mu),nu=nu,
                                            numeric_lower=float(floor))))
    return candidates


def scan(path, audits):
    started=timer()
    details, by_frequency = [], []
    records={(row["frequency_GHz"],row["target"]):row for row in audits["cases"]}
    with np.load(path) as prepared:
        for fi,freq in enumerate(FREQUENCIES):
            t0=timer()
            spectrum,s,b=(prepared[f"f{fi}_{name}"] for name in ("spectrum","S","B"))
            y=[prepared[f"f{fi}_{name}"] for name in TARGETS]
            cutoffs={}
            for cutoff in CUTOFFS:
                q=grid(cutoff,len(spectrum))
                cutoffs[cutoff]=(q,gram(s,b,q,len(spectrum)))
            scales=[np.linalg.norm(s[rx],ord=2)**2*2.**np.array(MU_DYADIC_POWERS)
                    for rx,tx in GROUPS]
            frequency_cost=Counter()
            winners={}
            for cap in CAPS:
                for multiplier in RATIOS_OF_CAP:
                    t=float(multiplier*cap)
                    for cutoff,(q,full_g) in cutoffs.items():
                        config=coefficient_spec(spectrum,q,cap,t)
                        if config is None:
                            frequency_cost["omitted_symbol_nonpositive"]+=1
                            continue
                        c,d=config
                        proposal=inverse_gram_numeric(full_g,d,c)
                        frequency_cost["small_A_eigenproblems"]+=1
                        if proposal is None:
                            frequency_cost["numeric_A_unresolved"]+=1
                            continue
                        cinv,base_nu,dimension=proposal
                        for group, (rx,tx) in enumerate(GROUPS):
                            if base_nu >= len(rx)+len(tx):
                                frequency_cost["A_inertia_exhausted_profiles"]+=1
                                continue
                            yy=[a[np.ix_(rx,tx)] for a in y]
                            for obj,row in numeric_group(cinv,base_nu,yy,rx,tx,cap,t,scales[group],frequency_cost):
                                for arm,predicate in ARM_TESTS.items():
                                    if not predicate(t):
                                        continue
                                    key=(obj,cap,group,arm)
                                    pool=winners.setdefault(key,[])
                                    pool.append(dict(**row,t=t,cutoff=cutoff,c=c,
                                                     A_negative=base_nu,A_compact_size=dimension))
                                    pool.sort(key=lambda a:a["numeric_lower"],reverse=True)
                                    del pool[3:]
            for obj,target in enumerate(TARGETS):
                case=records[(freq,target)]
                for cap in CAPS:
                    for arm in ARM_TESTS:
                        record=dict(case_id=case["id"],frequency_GHz=float(freq),target=target,cap=cap,arm=arm,
                                    field_norm=case["field_norm"], groups=[])
                        for group in range(len(GROUPS)):
                            pool=winners.get((obj,cap,group,arm),[])
                            record["groups"].append(dict(group=group,proposals=pool))
                        best_counts=[entry["proposals"][0]["numeric_lower"] if entry["proposals"] else 0.
                                     for entry in record["groups"]]
                        record["numeric_bound_relative"]=float(np.linalg.norm(best_counts)/case["field_norm"])
                        details.append(record)
            row=dict(frequency_GHz=float(freq),**frequency_cost,**elapsed(t0))
            by_frequency.append(row)
            print(freq, dict(frequency_cost), "sec", round(row["wall_seconds"],2), flush=True)
    return dict(caps=CAPS,ratio_of_cap=RATIOS_OF_CAP,mu_dyadic_powers=MU_DYADIC_POWERS,cutoffs=CUTOFFS,
                grid_side=SIDE,rows=details,by_frequency=by_frequency,**elapsed(started))


def rigorously_count_inertia(a, approximate):
    values,basis=eigh(h(approximate),check_finite=False)
    vp=ball(basis)
    gram_min=gershgorin(star(vp)*vp,upper=False)
    if not gram_min>0:
        raise ArithmeticError("Eigenbasis congruence is unresolved.")
    transform=star(vp)*a*vp
    negative=0
    for i,value in enumerate(values):
        off=sum((abs(transform[i,j]).upper() for j in range(len(values)) if i!=j),arb(0))
        center=transform[i,i].real
        if value<0 and center+off<0:
            negative+=1
        elif value>=0 and center-off>0:
            continue
        else:
            raise ArithmeticError("A verified Gershgorin interval crosses zero.")
    return negative


def center(matrix):
    return np.array([complex(a) for a in matrix.entries()]).reshape(matrix.nrows(),matrix.ncols())


def roundoff(k):
    factor=arb(k)*arb(2)**-53
    return factor/(1-factor)


def binary_complex_product(a,b):
    return ((dgemm(1.,a.real,b.real)-dgemm(1.,a.imag,b.imag))+
            1j*(dgemm(1.,a.real,b.imag)+dgemm(1.,a.imag,b.real)))


def norm_upper(a):
    # Squaring a midpoint-radius complex interval as a scalar product can
    # enclose a negative number. Positive individual upper endpoints avoid it.
    return sum((abs(v).upper()**2 for v in a.entries()),arb(0)).sqrt().upper()


def inverse_certificate_fast(a, rhs):
    approximate=center(a)
    d,v=eigh(h(approximate),check_finite=False)
    vd=v*d
    vstar=v.conj().T
    vb=ball(v)
    size=len(d)
    vnorm2=frobenius_squared(vb).upper()
    dmax=arb(float(max(abs(d))))
    guard=roundoff(4*size+32)
    # A conservative underflow allowance for more scalar operations than
    # present in a real dot product and the surrounding complex sums.
    underflow=arb((16*size+128)*size*size)*arb(2)**-1022
    gram=binary_complex_product(vstar,v)
    grem=ball(gram)
    for i in range(size):
        grem[i,i]-=1
    orthogonal_defect=(norm_upper(grem)+4*guard*vnorm2+underflow).upper()
    # The comparison matrix is the exact V diag(d) V*, rather than the
    # rounded center product. Its inertia follows from congruence.
    product=binary_complex_product(vd,vstar)
    scaled_guard=roundoff(3)*dmax*vnorm2+underflow
    reconstructed=(4*guard*(dmax*vnorm2+scaled_guard) + scaled_guard+underflow)
    epsilon=(norm_upper(a-ball(product))+reconstructed).upper()
    gap=((1-orthogonal_defect)*arb(float(min(abs(d))))-epsilon).lower()
    if not (orthogonal_defect<1 and gap>0):
        raise ArithmeticError("Fast binary64 inertia residual did not resolve the spectral gap.")
    xx=ball(v@((vstar@center(rhs))/d[:,None]))
    err=rhs-a*xx
    residual=norm_upper(err)
    shield=(residual**2/gap).upper()
    approx=star(rhs)*xx+star(xx)*err
    return approx,shield,int(np.count_nonzero(d<0)),dict(
        compact_gap_lower=lo(gap), compact_decomposition_error_upper=hi(epsilon),
        compact_basis_defect_upper=hi(orthogonal_defect),
        compact_rhs_defect_upper=hi(residual))


def inverse_gram_verified(exact, symbol, q, cap, t, c):
    started=timer()
    u,t,base=map(arb,(cap,t,c))
    kept=np.zeros(symbol.shape,bool)
    kept[q[:,0],q[:,1]]=True
    slack=arb(100)
    for g in symbol[~kept]:
        outside=1-u*arb(float(g.real))+t*arb(float(g.imag))-base
        if not outside>0:
            raise ArithmeticError("An omitted mode has no strict lower slack.")
        slack=min(slack,outside.lower())
    keep,d=[],[]
    for index,(yy,xx) in enumerate(q):
        g=symbol[yy,xx]
        coeff=1-u*arb(float(g.real))+t*arb(float(g.imag))-base
        if coeff==0:
            continue
        if not (coeff>0 or coeff<0):
            raise ArithmeticError("A retained Fourier coefficient straddles zero.")
        keep.append(index)
        d.append(coeff)
    ix=range(2*REC)
    small=sub(exact.unit,keep,keep)/base
    for i,dd in enumerate(d):
        small[i,i]+=1/dd
    cross=sub(exact.up,keep,ix)
    try:
        correction,shield,small_negative,method=inverse_certificate_fast(small,cross)
        method["inverse_verifier"]="binary64_reconstruction_with_outward_defect"
    except ArithmeticError:
        small_negative=rigorously_count_inertia(small,center(small))
        solved=small.solve(cross)
        correction=star(cross)*solved
        shield=arb(0)
        method={"inverse_verifier":"direct_Arb_congruence_and_solve"}
    negative=int(sum(x<0 for x in d)-small_negative)
    if not 0<=negative<=len(keep):
        raise ArithmeticError("Invalid A inertia after exact congruence.")
    inverse_gram=exact.g0/base-correction/base**2
    delta=(shield/base**2).upper()
    expansion=acb(arb(0,delta),arb(0,delta))
    for i in ix:
        for j in ix:
            inverse_gram[i,j]+=expansion
    return inverse_gram,negative, dict(state_dimension=SIDE**2,proof_dimension=len(keep),
                                       A_negative=negative,omitted_symbol_slack_lower=lo(slack),
                                       **method,**elapsed(started))


def rigorous_group(inv, base_nu, y, rx, tx, cap, ratio, mu):
    q=range(len(tx))
    rmat=sub(inv,rx,rx)
    bb=sub(inv,REC+tx,REC+tx)
    c=sub(inv,rx,REC+tx)
    z=acb(arb(cap),arb(ratio))/2
    scal=arb(mu)
    for i in q:
        rmat[i,i]+=scal
    nsmall=rigorously_count_inertia(rmat,center(rmat))
    nu=base_nu-nsmall
    if not 0<=nu<len(tx):
        raise ArithmeticError("Verified K inertia exhausts the source group.")
    defect=ball(y)-c*z
    squared=z.real**2+z.imag**2
    hf=(star(defect)*rmat.solve(defect)-bb*squared)*scal
    _,modes=eigh(h(center(hf)),check_finite=False)
    values=[]
    for k in range(nu+1,len(tx)+1):
        v=ball(modes[:,-k:])
        floor=gershgorin(star(v)*hf*v,upper=False)
        if not floor>0:
            break
        gram_ceil=gershgorin(star(v)*v,upper=True)
        values.append((floor/gram_ceil).lower())
    return sum(values,arb(0)),dict(nu=nu,retained_positive_source_directions=len(values))


def verify(path, scan_report, previous=None, frequency_filter=None):
    ctx.prec,ctx.threads=96,3
    started=timer()
    results,costs=[],[]
    total_prefix_wall,total_prefix_cpu=0.,0.
    if previous:
        results,costs=previous["rows"],previous["costs"]
        total_prefix_wall=previous["wall_seconds"]
        total_prefix_cpu=previous["cpu_seconds"]
    completed = {row["frequency_GHz"] for row in costs}
    with np.load(path) as prepared:
        for fi,frequency in enumerate(FREQUENCIES):
            if frequency in completed or (frequency_filter is not None and frequency not in frequency_filter):
                continue
            fstart=timer()
            selected=[row for row in scan_report["rows"] if row["frequency_GHz"]==frequency]
            spectrum,s,b=(prepared[f"f{fi}_{name}"] for name in ("spectrum","S","B"))
            cost=Counter()
            geometry, base_cache, group_cache={}, {}, {}
            profile_proofs=[]
            for row in selected:
                y=prepared[f"f{fi}_{row['target']}"]
                exact_norm=sum((arb(float(z.real))**2+arb(float(z.imag))**2 for z in y[USED]),arb(0)).sqrt()
                total=arb(0)
                group_records=[]
                for item,(rx,tx) in zip(row["groups"],GROUPS,strict=True):
                    out=dict(group=item["group"],attempts=[])
                    for prop in item["proposals"]:
                        cutoff=prop["cutoff"]
                        cap,t,c,mu=row["cap"],prop["t"],prop["c"],prop["mu"]
                        state_key=cutoff,cap,t,c
                        try:
                            if state_key not in base_cache:
                                if cutoff not in geometry:
                                    q=grid(cutoff,len(spectrum))
                                    geometry[cutoff]=(ExactGeometry(s,b,q,len(spectrum),cutoff),q)
                                    cost["exact_geometry_builds"]+=1
                                exact,q=geometry[cutoff]
                                cost["exact_A_attempts"]+=1
                                try:
                                    base_cache[state_key]=inverse_gram_verified(exact,spectrum,q,cap,t,c)
                                    iv,nu,summary=base_cache[state_key]
                                    profile_proofs.append(dict(cutoff=cutoff,cap=cap,t=t,c=c,**summary))
                                except (ArithmeticError,ZeroDivisionError,ValueError) as err:
                                    base_cache[state_key]=str(err)
                            verified=base_cache[state_key]
                            if isinstance(verified,str):
                                raise ArithmeticError(verified)
                            iv,nu,summary=verified
                            case_key=state_key,mu,item["group"],row["target"]
                            if case_key not in group_cache:
                                cost["exact_stationary_group_solves"]+=1
                                group_cache[case_key] = rigorous_group(iv,nu,y[np.ix_(rx,tx)],rx,tx,cap,t,mu)
                            score,details=group_cache[case_key]
                            out["attempts"].append(dict(**prop,verified_lower=lo(score.sqrt()) if score>0 else 0.,
                                                        **{k:v for k,v in details.items() if k!="nu"},
                                                        verified_inertia=details["nu"]))
                            if score>0:
                                total+=score
                                break
                        except (ArithmeticError,ZeroDivisionError,ValueError) as err:
                            out["attempts"].append(dict(**prop,unresolved=str(err)))
                    group_records.append(out)
                floor=lo(total.sqrt()) if total>0 else 0.
                ratio=lo(total.sqrt()/exact_norm) if total>0 else 0.
                results.append(dict(**{key:value for key,value in row.items() if key!="groups"},
                                    groups=group_records,verified_lower=floor,verified_lower_relative=ratio))
            crow=dict(frequency_GHz=float(frequency),proofs=profile_proofs,**cost,**elapsed(fstart))
            costs.append(crow)
            print("verified",frequency,"A attempts",cost["exact_A_attempts"],
                  "positive",sum(r["verified_lower"]>0 for r in results if r["frequency_GHz"]==frequency),
                  "seconds",round(crow["wall_seconds"],2),flush=True)
            # Keep a resumable evidence prefix while proceeding by frequency.
            cost=elapsed(started)
            checkpoint = "proof_progress.json" if frequency_filter is None else "proof_selected_progress.json"
            (path.parent/checkpoint).write_text(json.dumps(
                dict(grid_side=SIDE,rows=results,costs=costs,
                     wall_seconds=cost["wall_seconds"]+total_prefix_wall,
                     cpu_seconds=cost["cpu_seconds"]+total_prefix_cpu),indent=2)+"\n")
    cost=elapsed(started)
    return dict(grid_side=SIDE,rows=results,costs=costs,interval_threads=ctx.threads,
                interval_precision=ctx.prec,passed_prefix_frequencies=sorted(completed),
                wall_seconds=cost["wall_seconds"]+total_prefix_wall,
                cpu_seconds=cost["cpu_seconds"]+total_prefix_cpu)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared",type=Path,default=PREPARED)
    parser.add_argument("--physical",type=Path,default=HERE/"physical_audit.json")
    parser.add_argument("--verify",action="store_true")
    parser.add_argument("--resume",type=Path)
    parser.add_argument("--frequency",type=float,action="append",choices=FREQUENCIES.tolist())
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    with threadpool_limits(1):
        if args.verify:
            proposed=json.loads((HERE/"numeric_scan.json").read_text())
            prev=json.loads(args.resume.read_text()) if args.resume else None
            data=verify(args.prepared,proposed,prev,args.frequency)
            output=args.output or (HERE/"verified_scan.json" if args.frequency is None
                                   else args.prepared.parent/"selected_verification.json")
        else:
            data=scan(args.prepared,json.loads(args.physical.read_text()))
            output=args.output or HERE/"numeric_scan.json"
        output.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
