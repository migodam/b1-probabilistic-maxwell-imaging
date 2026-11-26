"""Run frozen, GT-blind profile/ROI scans and unit-normalized A1 energy requirements."""
import argparse
import json
import re
import time

import numpy as np
from scipy.linalg import eigh, norm, solve
from scipy.ndimage import distance_transform_edt
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

from sector import (HERE, PUBLIC, UPPER, Green, Profile, adj, blocks, dark_compression,
                    dark_trial, herm, surrogate_uniform, tail_from_blocks)

OPTICAL = (0.0, 0.25, 1.0, 4.0, 16.0, 64.0, 256.0)
WEIGHT = (0.002, 0.012, 0.06, 0.3)
CAPS = (0.0, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 3.0)
TIME = time.perf_counter


def profile(upper, s, w=1.0, optical=0.0, lower=0.0, receive=1.0):
    v = np.ones_like(upper)
    return Profile(v*lower, upper, v*w, v*optical, np.ones(len(s))*receive)


def spectra(a):
    return eigh(a, eigvals_only=True, check_finite=False)


def dark_profile(trial, scattered, s, p):
    h = dark_compression(trial, scattered, s, p)
    e = spectra(h)
    return int((e < -1e-8).sum()), float(e[min(len(e), len(s))-1]), e


def optimize_pixels(j, q, s, upper, start, roi):
    """Direct per-pixel negative-index proposal, with a smooth spatial step prior on weights."""
    n = len(j)
    safe = max(1, s.shape[0])
    test = (8, 16) if len(j.T) >= 16 else (4, 8)
    optical_bound = 256.0
    w0 = np.where(roi, 2.0, 0.7)
    w0 = w0 / w0.mean()
    x0 = np.r_[w0, np.full(n, start / optical_bound)]
    calls = 0

    def objective(x, index):
        nonlocal calls
        calls += 1
        weights = x[:n] / np.mean(x[:n])
        optical = optical_bound * x[n:]
        p = Profile(np.zeros(n), upper, weights, optical, np.ones(len(s)))
        k = dark_compression(j, q, s, p)
        eig, vectors = eigh(k, subset_by_index=[index-1, index-1],
                            driver="evr", check_finite=False)
        current = j @ vectors[:, 0]
        field = q @ vectors[:, 0]
        local = current.conj() * field
        derivative = abs(current)**2-upper*local.real
        dw = (derivative - (derivative @ weights)/n) / np.mean(x[:n])
        regularizer = 0.005/safe
        grad = -np.r_[dw, optical_bound * local.imag] + regularizer*(x-x0)
        value = -eig[0] + 0.5*regularizer*norm(x-x0)**2
        return value, grad

    best = None
    for index in test:
        run = minimize(objective, x0, args=(index,), jac=True, method="L-BFGS-B",
                       bounds=[(0.035, 4.0)]*n + [(-1.0, 1.0)]*n,
                       options={"maxiter":28, "maxfun":45, "ftol":1e-8, "gtol":1e-6})
        w = run.x[:n] / np.mean(run.x[:n])
        lam = optical_bound * run.x[n:]
        p = Profile(np.zeros(n), upper, w, lam, np.ones(len(s)))
        nu, edge, _ = dark_profile(j, q, s, p)
        row = (nu, -edge, p, str(run.message), int(run.nit))
        if best is None or row[:2] < best[:2]:
            best = row
    return best[2], dict(calls=calls, stop=best[3], iterations_for_selected=best[4],
                         receiver_weight="identity; nonnegative reweight leaves exact receiver-null trials unchanged")


def extended_regions(mask):
    if not np.any(mask):
        return {}
    side = mask.shape[0]
    d = distance_transform_edt(~mask)
    order = np.argsort(d.ravel(), kind="stable")
    out = {}
    for fraction in (0.50, 0.75, 0.90, 0.97):
        region = np.zeros(mask.size, bool)
        region[order[:int(np.ceil(mask.size*fraction))]] = True
        out[f"expanded_{int(fraction*100)}"] = region
    return out


def energy_requirements(k, rhs, h0, radius, reference):
    e, v = eigh(k, driver="evd", check_finite=False)
    source = adj(v) @ rhs
    base = max(1e-15, np.max(abs(e)))
    shifts = base*np.geomspace(1e-4, 16.0, 64)
    margin = []
    for shift in shifts:
        val = e + shift
        nu = int((val < 0).sum())
        if nu >= h0.shape[0] or np.min(np.abs(val)) < 1e-13:
            continue
        schur = spectra(herm(h0 - adj(source) @ (source/val[:, None])))[::-1]
        tail2 = float(np.maximum(schur[nu:], 0).sum())
        margin.append(dict(shift=float(shift), nu=nu, tail2=tail2,
                           energy_required_upper=(tail2-radius**2)/shift,
                           energy_required_over_incident=(tail2-radius**2)/(shift*reference)))
    best = max(margin, key=lambda r: r["energy_required_upper"], default=None)
    return dict(base_negative=int((e < 0).sum()), shifted=margin, best=best,
                reference_incident_fro2=reference)


def numeric_support(green, s, b, y, radius, roi, store, key, obstruction):
    active = np.flatnonzero(~roi)
    if not len(active):
        return dict(active=0, best=dict(distance=float(norm(y)), nu=0,
                                        margin=float(norm(y)-radius), exact_zero_forward=True))
    if len(active) > 1200:
        return dict(active=int(len(active)), best=None, reason="support system exceeds dense audit limit")
    p0 = len(active)
    min_nu=obstruction["profile_with_smallest_dark_count"]["negative_lower_numeric"]
    if min_nu >= y.shape[1] and key != "roi_01":
        return dict(active=int(p0), best=None, reason="all trial profiles have at least s negative directions on exterior support")
    g = green.matrix(active, active)
    ref = float(norm(b[active])**2)
    scale = float(norm(s[:, active], 2)**2)
    profile_rows = []
    energy = None
    verified_candidate = None
    se,be=s[:,active],b[active]
    h0=herm(adj(y)@y)
    for optical in OPTICAL:
        name="real_uniform" if optical==0 else f"optical_uniform_{optical:g}"
        neg=next((a["negative_lower_numeric"] for a in obstruction["trials"] if a["method"]==name),0)
        if neg>=y.shape[1] and optical!=0:
            continue
        mix=(UPPER+1j*optical)/2
        state0=np.eye(p0)-mix*g-mix.conjugate()*adj(g)
        de,v=eigh(herm(state0),driver="evd",check_finite=False)
        base_nu=int((de<0).sum())
        pivot_gap=float(min(abs(de)))
        receiver=se@v
        incident=adj(v)@(mix*be)
        lhs=adj(receiver)@y
        for i, weight in enumerate(WEIGHT):
            if neg>=y.shape[1] and i!=2:
                continue
            w=scale*weight
            p=profile(np.full(p0,UPPER),s,w,optical*w)
            if optical==0 and i==2:
                k,l,c=blocks(g,se,be,y,p)
                energy=energy_requirements(k,l,c,radius,ref)
            if pivot_gap<1e-11:
                k,l,c=blocks(g,se,be,y,p)
                distance,nu,gap,vals=tail_from_blocks(k,l,c)
                extra=dict(gap=gap,pivot="original dense active state")
            else:
                inv=1/(w*de)
                mids=herm(np.eye(len(s))+(receiver*inv)@adj(receiver))
                aux=spectra(mids)
                nu=int(base_nu+np.sum(aux>0)-len(s))
                distance,vals=0.,[]
                if nu<y.shape[1] and min(abs(aux))>1e-12:
                    proj=lhs+w*incident
                    sol=inv[:,None]*proj
                    corr=receiver@sol
                    h=herm(h0-adj(proj)@sol+adj(corr)@solve(mids,corr,assume_a="her"))
                    vals=spectra(h)[::-1].tolist()
                    distance=float(np.maximum(vals[max(0,nu):],0).sum()**0.5)
                extra=dict(core_pivot_gap=pivot_gap,receiver_pivot_gap=float(min(abs(aux))),
                           pivot="base profile plus rank-receivers update")
            profile_rows.append(dict(weight=float(w),ratio=optical,nu=nu,
                                     distance=distance,margin=distance-radius,eigenvalues=vals,**extra))
            if distance>radius and (verified_candidate is None or distance>verified_candidate[0]):
                verified_candidate=distance,p
    best = max(profile_rows, key=lambda x: (x["distance"], -x["nu"]))
    if verified_candidate is not None:
        store[key+"_active"] = active
        store[key+"_weight"] = verified_candidate[1].weight
        store[key+"_optical"] = verified_candidate[1].optical
    return dict(active=int(len(active)), best=best, profiles=profile_rows, energy_condition=energy)


def support_optimization(trial, scattered, s, roi, store, key, use_optimizer):
    upper = np.where(roi, 0.0, UPPER)
    if trial.shape[1] == 0:
        result=dict(method="real_uniform",negative_lower_numeric=0,edge_at_s=None)
        return dict(trials=[result],profile_with_smallest_dark_count=result)
    trials = []
    configurations = {}
    mass=herm(adj(trial)@trial)
    propagation=adj(trial)@scattered
    realpart=herm(propagation)
    opticalpart=(propagation-adj(propagation))/(2j)
    rec=s@trial
    residual=herm(adj(rec)@rec)
    for ratio in OPTICAL:
        p = profile(upper, s, optical=ratio)
        e=spectra(residual+mass-UPPER*realpart+ratio*opticalpart)
        count=int((e < -1e-8).sum())
        edge=float(e[min(len(e),len(s))-1])
        name = "real_uniform" if ratio == 0 else f"optical_uniform_{ratio:g}"
        configurations[name] = p
        trials.append(dict(method=name, negative_lower_numeric=count, edge_at_s=edge))
    for level in (0.05, 0.25, 4.0):
        for ratio in (0.0, 1.0, 16.0, 64.0):
            w = np.where(roi, 1.0, level)
            p = profile(upper, s, w=w, optical=ratio*w)
            e=spectra(residual+level*(mass-UPPER*realpart+ratio*opticalpart))
            count=int((e < -1e-8).sum())
            edge=float(e[min(len(e),len(s))-1])
            name = f"roi_weight_{level:g}_optical_{ratio:g}"
            configurations[name] = p
            trials.append(dict(method=name, negative_lower_numeric=count, edge_at_s=edge))
    if use_optimizer and np.any(roi) and trial.shape[1]>=16:
        best = min(trials, key=lambda x: (x["negative_lower_numeric"], -x["edge_at_s"]))
        start = float(np.median(configurations[best["method"]].optical))
        pix, cost = optimize_pixels(trial, scattered, s, upper, start, roi)
        count, edge, _ = dark_profile(trial, scattered, s, pix)
        name = "optimized_pixel"
        configurations[name] = pix
        trials.append(dict(method=name, negative_lower_numeric=count, edge_at_s=edge, optimization=cost))
    sel = min(trials, key=lambda x: (x["negative_lower_numeric"], -x["edge_at_s"]))
    for name in {"real_uniform", sel["method"], "optimized_pixel"} & set(configurations):
        p = configurations[name]
        k = dark_compression(trial, scattered, s, p)
        e, v = eigh(k, subset_by_index=(0,min(len(s),len(k))-1), check_finite=False)
        store[key+"_"+name+"_currents"] = trial @ v
        store[key+"_"+name+"_upper"] = upper
        store[key+"_"+name+"_weight"] = p.weight
        store[key+"_"+name+"_optical"] = p.optical
    return dict(trials=trials, profile_with_smallest_dark_count=sel)


def signed_updates(green, s, b, y, roi, scale):
    active = np.flatnonzero(roi)
    if not len(active) or len(active) > 256 or green.side != 32:
        return None
    clock = TIME()
    all_cells = np.arange(green.n)
    g = green.matrix(all_cells, all_cells)
    p = profile(np.full(green.n, UPPER), s, scale, scale)
    base, lbase, h0 = blocks(g, s, b, y, p)
    e, v = eigh(base, driver="evd", check_finite=False)
    nu0 = int((e < 0).sum())
    traces = np.vstack((np.eye(green.n)[active], g[active]))
    vt = adj(v) @ adj(traces)
    pretrace = herm(adj(vt) @ (vt/e[:, None]))
    cached = TIME()-clock
    rows = []
    for cap in CAPS:
        a = (UPPER-cap)*scale/2
        c = TIME()
        rhs = lbase.copy()
        rhs[active] -= a*b[active]
        proj = adj(v) @ rhs
        small = pretrace.copy()
        m = len(active)
        small[np.arange(m), m+np.arange(m)] += 1/a
        small[m+np.arange(m), np.arange(m)] += 1/a
        eig = spectra(herm(small))
        nu = int(nu0+np.sum(eig>0)-m)
        if nu < y.shape[1]:
            routed = adj(vt) @ (proj / e[:, None])
            schur = herm(h0 - adj(proj) @ (proj/e[:, None]) + adj(routed) @ solve(small, routed, assume_a="her"))
            es = spectra(schur)[::-1]
            dist = float(np.maximum(es[max(0,nu):],0).sum()**0.5)
        else:
            dist = 0.
        after_update = TIME()-c
        c = TIME()
        endpoints = np.full(green.n, UPPER)
        endpoints[active] = cap
        direct = tail_from_blocks(*blocks(g,s,b,y,profile(endpoints,s,scale,scale)))
        rows.append(dict(cap=cap, nu_lowrank=nu, nu_direct=direct[1], proof_size=2*m,
                         distance_lowrank=dist, distance_direct=direct[0], update_seconds=after_update,
                         direct_seconds=TIME()-c, proof="numerical equality"))
    return dict(base_nu=nu0, preparation_seconds=cached, rounds=rows)


def run_case(side, cid, with_images, geometry_only=False, optimize=False, matrix=False, image_version="published-stop"):
    data = np.load(PUBLIC/"cases"/f"{cid}.npz")
    op_id = str(data["operator_id"])
    op = PUBLIC/"operators"/f"{op_id}.npz"
    green = Green.read(op, side)
    raw = np.load(op)
    s, b, y = raw["S"], raw["E0"], data["Y"]
    radius = float(data["radius"])
    scale = float(norm(s,2)**2)
    clock, cpu = TIME(),time.process_time()
    out = dict(case=cid, operator=op_id, side=side, sources=b.shape[1], receivers=s.shape[0],
               radius=radius, noise_over_data=radius/float(norm(y)), shape_prior="whole field 0 <= chi <= 4.5",
               sections={})
    section = out["sections"]
    store = {}
    all_cells = np.arange(green.n)
    trial, scattered = dark_trial(green, s, all_cells, rank=104)
    section["global_at_prior_cap_dark"] = support_optimization(trial, scattered, s,
                    np.zeros(green.n,bool),store,"global_prior",False)
    if not geometry_only:
        physical = {}
        mass=herm(adj(trial)@trial)
        propagation=adj(trial)@scattered
        realpart=herm(propagation)
        opticalpart=(propagation-adj(propagation))/(2j)
        for cap in CAPS[1:]:
            choices = []
            for ratio in (0.0,0.25,1.,4.,16.):
                bound=int((spectra(mass-cap*realpart+ratio*opticalpart)<-1e-8).sum())
                if bound >= y.shape[1]:
                    choices.append(dict(weight=None,optical=ratio,distance=0.,nu_surrogate=None,
                                        compressed_negative_at_least_numeric=bound,h=[],margin=-radius))
                    continue
                for w in WEIGHT:
                    p = profile(np.full(green.n,cap), s, w=w*scale, optical=ratio*w*scale)
                    d, nu, eig, meta = surrogate_uniform(green,s,b,y,p)
                    choices.append(dict(weight=w*scale, optical=ratio, distance=d, nu_surrogate=nu, h=eig,
                                        margin=d-radius, **meta))
            physical[str(cap)] = sorted(choices, key=lambda r:r["distance"],reverse=True)
        section["global_cap"] = physical
    image = HERE.parent/"imaging"/"results"/cid
    image = image/("stationary/baseline.npz" if image_version=="stationary" else "baseline.npz")
    out["image_version"]=image_version
    if with_images and image.exists():
        source = np.load(image)
        masks = {f"roi_{j+1:02d}": m.ravel() for j,m in enumerate(source["rois"])}
        masks.update({f"background_{j+1:02d}": m.ravel() for j,m in enumerate(source["background_tiles"])})
        masks["roi_union"] = np.any(source["rois"],axis=0).ravel()
        for j,m in enumerate(source["background_tiles"]):
            masks[f"all_except_background_{j+1:02d}"] = ~m.ravel()
        yy,xx=np.indices((side,side))
        masks["left_half"]=(xx<side//2).ravel()
        masks["right_half"]=(xx>=side//2).ravel()
        masks["upper_half"]=(yy<side//2).ravel()
        masks["lower_half"]=(yy>=side//2).ravel()
        masks.update(extended_regions(source["roi"]))
        section["roi"] = {}
        for key, roi in masks.items():
            # For lower-cap hypotheses, exact exterior support keeps the dark
            # compression identical for every cap on the indicated ROI.
            local_start = TIME()
            active = np.flatnonzero(~roi)
            tr, sc = dark_trial(green,s,active,rank=104)
            record = dict(roi_cells=int(roi.sum()), obstruction=support_optimization(
                          tr, sc, s, roi, store,key,optimize and key=="roi_01" and side==32))
            if not geometry_only and (key.startswith(("expanded","all_except","roi_union")) or key.endswith("_half") or matrix):
                record["zero_roi_test"] = numeric_support(green,s,b,y,radius,roi,store,key,record["obstruction"])
            if matrix and key=="roi_01" and cid in ("n32_c01","n32_c05"):
                record["signed_updates"] = signed_updates(green,s,b,y,roi,scale*0.06)
            record["wall_seconds"] = TIME()-local_start
            section["roi"][key] = record
            store[key+"_mask"] = roi.reshape(side,side)
    out.update(wall_seconds=TIME()-clock, cpu_seconds=time.process_time()-cpu)
    return out,store


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append")
    parser.add_argument("--suffix",default="main")
    parser.add_argument("--with-images",action="store_true")
    parser.add_argument("--geometry-only",action="store_true")
    parser.add_argument("--optimize",action="store_true")
    parser.add_argument("--matrix",action="store_true")
    parser.add_argument("--image-version",choices=("published-stop","stationary"),default="published-stop")
    args=parser.parse_args()
    cases=args.case or [x.stem for x in sorted((PUBLIC/"cases").glob("n*_c*.npz"))
                         if re.match(r"n(?:32|96)_c\d+$",x.stem)]
    if not cases:
        raise FileNotFoundError("Public material observations are not available yet")
    for cid in cases:
        side=int(cid.split("_")[0][1:])
        out,store=run_case(side,cid,args.with_images,args.geometry_only,args.optimize,args.matrix,args.image_version)
        prefix=HERE/"results"/f"{cid}_{args.suffix}"
        prefix.with_suffix(".json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
        np.savez_compressed(prefix.with_suffix(".npz"), **store)
        print(cid, out["wall_seconds"], "ROI", len(out["sections"].get("roi",[])), flush=True)


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
