"""Joint real-voxel readout of per-frequency Fourier currents on archived R15 data."""

from argparse import ArgumentParser
from collections import Counter
import importlib
import json
import os
from pathlib import Path
import sys
from time import perf_counter, process_time
from types import ModuleType

for name in ("OPENBLAS_NUM_THREADS","VECLIB_MAXIMUM_THREADS","OMP_NUM_THREADS"):
    os.environ[name] = "1"

import numpy as np
from scipy.fft import ifft2
from scipy.linalg import cho_factor, cho_solve, lu_factor, lu_solve, qr, svd
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT/"results"
PACKS = ("holdout","holdout_debye")
METHODS = ("fixed","fft_joint","full_gn")
DESIGN_COUNTS = ((1,8,8),(2,4,8),(4,2,8))
N = 24
TRADEOFF = .32
FFT_HALFWIDTH = 3


def read_json(path):
    return json.loads(path.read_text())


def write_json(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n")


def r15_modules(path):
    candidates = [path] if path else [ROOT.parent/"research_z_alg_r15",
        ROOT.parents[2]/"r15_publication/research/research_z_alg_r15"]
    for location in candidates:
        location = location.resolve()
        if (location/"src/ndt.py").exists() and (location/"results/holdout_observations.npz").exists():
            break
    else:
        raise FileNotFoundError("R15 source and observations are required. Pass --r15 with the archived run directory.")
    package = ModuleType("r15_reference")
    package.__path__ = [str(location/"src")]
    sys.modules[package.__name__] = package
    return location,tuple(importlib.import_module(f"{package.__name__}.{mod}") for mod in ("ndt","designs","inversion"))


def fft_bank():
    a,b = np.meshgrid(np.arange(-FFT_HALFWIDTH,FFT_HALFWIDTH+1),
                      np.arange(-FFT_HALFWIDTH,FFT_HALFWIDTH+1),indexing="ij")
    impulses = np.zeros((N,N,a.size),dtype=complex)
    impulses[a.ravel()%N,b.ravel()%N,np.arange(a.size)]=1
    return ifft2(impulses,axes=(0,1),norm="ortho").reshape(N*N,-1)


def power(z):
    return float(np.vdot(z,z).real)


def prep(models, y, sigma, dictionary, bright_soft=False):
    t0,cpu0 = perf_counter(),process_time()
    blocks = []
    for model,data in zip(models,y,strict=True):
        left,s,right = svd(model.Sh,full_matrices=False,check_finite=False)
        transformed = left.conj().T@data
        bright_ok = (s>=.13*s[0]) & (np.linalg.norm(transformed,axis=1)/np.sqrt(data.shape[1])>=2.75*sigma)
        hit = np.flatnonzero(bright_ok)
        rank = int(hit[-1]+1) if hit.size else 0
        v = right[:rank].conj().T
        anchor = v@(transformed[:rank,:]/s[:rank,None])
        blocks.append(dict(anchor=anchor,anchor_field=model.Eh+model.Gh@anchor,bright=rank,v=v,
                           beta=float((TRADEOFF*s[0])**2),target=data-model.Sh@anchor))
    anchor_cost = dict(wall=perf_counter()-t0,cpu=process_time()-cpu0)
    t0,cpu0 = perf_counter(),process_time()
    for model,b in zip(models,blocks,strict=True):
        v = b.pop("v")
        q,se,_ = svd(dictionary-v@(v.conj().T@dictionary),full_matrices=False,check_finite=False)
        q = q[:,se>1e-9]
        if bright_soft:
            q = np.column_stack((v,q[:,:dictionary.shape[1]-v.shape[1]]))
            if q.shape[1]!=dictionary.shape[1]:
                raise ValueError("The soft-bright control must preserve current dimension")
        sq = model.Sh@q
        b.update(q=q,gq=model.Gh@q,sqs=sq.conj().T@sq,projected=sq.conj().T@b["target"])
    return blocks,anchor_cost,dict(wall=perf_counter()-t0,cpu=process_time()-cpu0)


def correction(block,multiplier,material):
    z = multiplier*material[:,None]
    a = block["q"]-z*block["gq"]
    target = z*block["anchor_field"]-block["anchor"]
    beta = block["beta"]
    normal = block["sqs"]+beta*a.conj().T@a
    rhs = block["projected"]+beta*a.conj().T@target
    solution = cho_solve(cho_factor(normal,check_finite=False),rhs,check_finite=False)
    kkt = np.linalg.norm(normal@solution-rhs)/max(np.linalg.norm(rhs),1e-30)
    return block["anchor"]+block["q"]@solution,float(kkt)


def joint_radial_fit(currents,fields,beta,factors,previous,penalty,weight,bounds):
    d = np.zeros_like(previous)
    b = np.zeros_like(previous)
    for j,e,w,f in zip(currents,fields,beta,factors,strict=True):
        z = f*e
        d += w*np.sum(abs(z)**2,axis=1)
        b += w*np.sum((z.conj()*j).real,axis=1)
    prior = np.divide(b,d,out=previous.copy(),where=d>0)
    initial = np.clip(prior,*bounds)
    d *= weight
    b *= weight

    def objective(m):
        val,grad = penalty.evaluate(m)
        return float(m@(.5*d*m-b)+val), d*m-b+grad

    ans = minimize(objective,initial,jac=True,method="L-BFGS-B",bounds=[bounds]*len(previous),
                   options=dict(maxiter=45,gtol=1e-9,ftol=1e-12,maxcor=10))
    material = min((ans.x,previous,initial),key=lambda x:objective(x)[0])
    return material,int(ans.nfev)


def current_solver(model,y,blocks,factors,reference,penalty,bounds,method,iterations):
    cpu0,wall0 = process_time(),perf_counter()
    m = np.zeros(N*N)
    b = [c["beta"] for c in blocks]
    weight = 1/(reference**2*y.size)
    counts = Counter()
    hist = []
    stop = "iteration_cap"
    for it in range(iterations):
        if method=="fixed":
            currents = [c["anchor"] for c in blocks]
            kkt = None
        else:
            inner = [correction(z,f,m) for z,f in zip(blocks,factors,strict=True)]
            currents = [row[0] for row in inner]
            kkt = max(row[1] for row in inner)
            counts.update(reduced_factorizations=len(currents),
                          reduced_complex_rhs=sum(j.shape[1] for j in currents))
        fields = [h.Eh+h.Gh@j for h,j in zip(model,currents,strict=True)]
        m,evals = joint_radial_fit(currents,fields,b,factors,m,penalty,weight,bounds)
        counts["material_objective_evaluations"] += evals
        sv = sum(w*power(j-f*m[:,None]*e) for j,e,w,f in zip(currents,fields,b,factors,strict=True))
        dv = sum(power(h.Sh@j-data) for h,j,data in zip(model,currents,y,strict=True))
        value = .5*weight*(sv+dv)+penalty.evaluate(m)[0]
        hist.append(dict(iteration=it+1,value=value,data_cost=.5*weight*dv,state_cost=.5*weight*sv,
                         current_normal_relative_residual=kkt))
        if method=="fixed":
            stop = "material_fit"
            break
        if len(hist)>1 and abs(hist[-2]["value"]-value)<1e-10:
            stop = "objective_tolerance"
            break
    return m,dict(wall=perf_counter()-wall0,cpu=process_time()-cpu0,history=hist,
                  online_counts=dict(counts),stop=stop)


class FullForwardCounter:
    def __init__(self,inversion):
        self.forward = inversion.fullwave
        self.counts = Counter()

    def __call__(self, models, m, *, jacobian=False, material_factors=None):
        self.counts["requested_full_forward_batches"] += 1
        if np.any(m):
            self.counts["full_factorizations"] += len(models)
            self.counts["full_complex_rhs"] += sum(
                h.Eh.shape[1]+(h.Sh.shape[0] if jacobian else 0) for h in models)
            self.counts["analytic_receiver_transpose_rhs"] += sum(h.Sh.shape[0] for h in models) if jacobian else 0
        return self.forward(models,m,jacobian=jacobian,material_factors=material_factors)


def full_gn(inversion,models,y,factors,reference,penalty,iterations):
    callback = FullForwardCounter(inversion)
    inversion.fullwave = callback
    t0,c0 = perf_counter(),process_time()
    try:
        m,history,_ = inversion.projected_gn(models,y,reference,penalty,
                                             iterations=iterations,material_factors=factors)
    finally:
        inversion.fullwave = callback.forward
    return m,dict(wall=perf_counter()-t0,cpu=process_time()-c0,online_counts=dict(callback.counts),history=history)


def tangent_audit(models,bank,factors,material,vectors):
    """Apply K without assembling its spectrum, and check its diagonal explicitly."""
    wall,cpu = perf_counter(),process_time()
    n = material.size
    whole = np.zeros(n)
    diagonal = np.zeros(n)
    explicit_diagonal = np.zeros(n)
    capture = np.zeros(vectors.shape[1])
    action = np.zeros_like(vectors)
    block_measures = []
    coeff_rows = []
    fwd = rhs = 0
    for h,b,f in zip(models,bank,factors,strict=True):
        if np.any(material):
            mat = -f*material[:,None]*h.Gh
            mat.flat[::n+1]+=1
            J = lu_solve(lu_factor(mat,check_finite=False),f*material[:,None]*h.Eh,check_finite=False)
            e = h.Eh+h.Gh@J
            fwd += 1
            rhs += h.Eh.shape[1]
        else:
            e = h.Eh
        arm = b["q"]-(f*material)[:,None]*b["gq"]
        s = qr(arm,mode="economic",check_finite=False)[0]
        e = np.sqrt(b["beta"])*f*e
        g = np.sum(abs(e)**2,axis=1)
        leverage = np.sum(abs(s)**2,axis=1)
        whole += g
        diagonal += g*leverage
        r = np.einsum("ij,is->ijs",s.conj(),e)
        explicit_diagonal += np.sum(abs(r)**2,axis=(1,2))
        for j,d in enumerate(vectors.T):
            t = e*d[:,None]
            proj = s.conj().T@t
            capture[j] += power(proj)
            action[:,j] += np.sum((e.conj()*(s@proj)).real,axis=1)
        coeff_rows.append(g)
        block_measures.append(float(leverage.sum()))
    lit = whole>0
    tr = np.sum(diagonal[lit]/whole[lit])
    tr_directional = np.sum(explicit_diagonal[lit]/whole[lit])
    all_power = whole@vectors**2
    by_operator = np.sum(vectors*action,axis=0)/all_power
    direct = capture/all_power
    assert np.max(abs(direct-by_operator))<1e-12
    assert abs(tr-tr_directional)<1e-11
    return dict(trace=float(tr),explicit_voxel_trace=float(tr_directional),
                diagonal_error=float(np.max(abs(diagonal[lit]-explicit_diagonal[lit])/whole[lit])),
                quadratic_error=float(np.max(abs(direct-by_operator))),
                random_direction_retained_fractions=direct.tolist(),
                per_frequency_projector_trace=block_measures,
                frequency_illumination_fraction_min=np.min(np.stack(coeff_rows)[:,lit]/whole[lit],axis=1).tolist(),
                frequency_illumination_fraction_max=np.max(np.stack(coeff_rows)[:,lit]/whole[lit],axis=1).tolist(),
                cpu=process_time()-cpu,wall=perf_counter()-wall,full_factorizations=fwd,full_complex_rhs=rhs)


def data_and_models(directory,ndt,designs):
    plans = tuple(designs.make(*k) for k in DESIGN_COUNTS)
    if any(a.shots!=8 or a.samples!=64 or a.shot_energy!=.125 for a in plans):
        raise ValueError("Unexpected R15 acquisition schedule")
    cold = {}
    for pack in PACKS:
        with np.load(directory/"results"/f"{pack}_observations.npz",allow_pickle=False) as data:
            observations = {k:data[k] for k in ("case_ids","frequencies_ghz","clean","noise","reference_amplitude",
                                                "sigma_2pct","sources_m","receivers_m","dielectric_factor","background_epsilon")}
        assert np.array_equal(observations["sources_m"],ndt.SOURCES)
        assert np.array_equal(observations["receivers_m"],ndt.RECEIVERS)
        hostfn = ndt.host if pack=="holdout" else ndt.host_debye
        calibrate = np.ones(len(observations["frequencies_ghz"])) if pack=="holdout" else ndt.dielectric_factor(
            observations["frequencies_ghz"])
        assert np.array_equal(calibrate,observations["dielectric_factor"])
        assert len(observations["case_ids"])==9
        for f in sorted({v for x in plans for v in x.freqs}):
            wall,cpu = perf_counter(),process_time()
            hostfn(N,f)
            cold[f"{pack}/{f:g}"] = dict(cpu=process_time()-cpu,wall=perf_counter()-wall,
                                         full_factorizations=1,complex_rhs=N*N+len(ndt.SOURCES)+len(ndt.RECEIVERS))
        for plan in plans:
            models = []
            factors = []
            for f in plan.freqs:
                h = hostfn(N,f)
                models.append(ndt.Host(h.Eh[:,plan.tx],h.Sh[plan.rx,:],h.Gh,h.xy))
                idx = int(np.flatnonzero(observations["frequencies_ghz"]==f)[0])
                factors.append(observations["dielectric_factor"][idx])
            yield pack,plan,observations,models,factors,cold


def fit(directory,modules,args):
    ndt,designs,inversion = modules
    result = {}
    arrays = {}
    dictionary = fft_bank()
    penalty = inversion.MaterialPenalty(tv=.0055,amplitude=.0035)
    tangent_vectors = np.random.default_rng(20261120).standard_normal((N*N,4))
    for pack,plan,ob,models,factors,cold in data_and_models(directory,ndt,designs):
        cases = {}
        freq_info = dict(freqs=plan.freqs,tx=plan.tx,rx=plan.rx,shots=plan.shots,
                         samples=plan.samples,energy=1.,shot_energy=plan.shot_energy)
        reference = float(ob["reference_amplitude"])
        for idx,cid in enumerate(ob["case_ids"]):
            if args.case and cid not in args.case:
                continue
            cid = str(cid)
            y,sigma = designs.read_arm(ob,idx,plan)
            assert y.size == plan.samples
            bank,anchor_setup,fft_setup = prep(models,y,sigma,dictionary)
            record = dict(sigma=sigma,bright_ranks=[x["bright"] for x in bank],
                          fft_ranks=[int(x["q"].shape[1]) for x in bank],
                          anchor_setup=anchor_setup,fft_specific_setup=fft_setup,outputs={})
            for method in METHODS:
                if method=="full_gn":
                    m,info = full_gn(inversion,models,y,factors,reference,penalty,args.iterations)
                else:
                    m,info = current_solver(models,y,bank,factors,reference,penalty,ndt.BOUNDS,method,args.iterations)
                    extra = fft_setup if method=="fft_joint" else {"cpu":0.,"wall":0.}
                    info["cpu"] += anchor_setup["cpu"]+extra["cpu"]
                    info["wall"] += anchor_setup["wall"]+extra["wall"]
                observer = FullForwardCounter(inversion)
                estimate = observer(models,m,material_factors=factors)
                info["full_residual_noise_rms"] = float(np.linalg.norm(estimate-y)/(sigma*np.sqrt(y.size)))
                info["full_residual_relative"] = float(np.linalg.norm(estimate-y)/np.linalg.norm(y))
                info["validation_counts"] = dict(observer.counts)
                key = f"{pack}/{plan.name}/{cid}/{method}"
                arrays[key] = m
                record["outputs"][method] = dict(array=key,**info)
            record["weighted_trace"] = {stage:tangent_audit(models,bank,factors,m,tangent_vectors)
                for stage,m in (("background",np.zeros(N*N)),("full_gn_fit",arrays[f"{pack}/{plan.name}/{cid}/full_gn"]))}
            cases[cid] = record
            print(pack,plan.name,cid,[(m,round(rec["full_residual_noise_rms"],3),round(rec["cpu"],3))
                 for m,rec in record["outputs"].items()],flush=True)
        result[f"{pack}/{plan.name}"] = dict(**freq_info,cases=cases)
    np.savez_compressed(RESULTS/"broadband_current_materials.npz",**arrays)
    write_json(RESULTS/"broadband_current.json",dict(settings=dict(iterations=args.iterations,tradeoff=TRADEOFF,
        fft_halfwidth=FFT_HALFWIDTH,support_threshold=.20,gain_fixed=1.,inversion_side=N,data_pack="original 48-by-48 archived",
        exact_per_frequency_host=True,case_subset=args.case), designs=result,cold_host_setup=cold))


def material_score(m,truth):
    signal = abs(truth)>=.2
    found = abs(m)>=.2
    overlap = (signal&found)&(truth*m>0)
    return dict(RMSE=float(np.linalg.norm(m-truth)/np.sqrt(m.size)),
                signed_support_IoU=float(overlap.sum()/(signal|found).sum()) if signal.any() else None,
                false_support_pixels=int((found&~signal).sum()), detected=bool(found.any()),
                maximum_abs_contrast=float(np.max(abs(m))))


def score(directory,modules):
    ndt,designs,inversion=modules
    study=read_json(RESULTS/"broadband_current.json")
    scored={}
    for pack,plan,ob,models,factors,_ in data_and_models(directory,ndt,designs):
        with np.load(directory/"results"/f"{pack}.npz",allow_pickle=False) as archive:
            true={str(cid):m.ravel() for cid,m in zip(archive["case_ids"],archive["truth"],strict=True)}
        rows=study["designs"][f"{pack}/{plan.name}"]["cases"]
        with np.load(RESULTS/"broadband_current_materials.npz",allow_pickle=False) as m:
            for idx,cid in enumerate(ob["case_ids"]):
                if cid not in rows:
                    continue
                truth = true[cid]
                y,sigma=designs.read_arm(ob,idx,plan)
                observer=FullForwardCounter(inversion)
                predicted=observer(models,truth,material_factors=factors)
                scored[f"{pack}/{plan.name}/{cid}"] = dict(is_empty=bool(not np.any(truth)),
                  truth_residual_noise_rms=float(np.linalg.norm(predicted-y)/(sigma*np.sqrt(y.size))),
                  truth_evaluation_cost=dict(observer.counts),
                  methods={key:material_score(m[r["array"]],truth) for key,r in rows[cid]["outputs"].items()})
    write_json(RESULTS/"broadband_current_scores.json",scored)


def ablation(directory,modules,args):
    """Change one current-step assumption on the existing two- and four-frequency data."""
    ndt,designs,inversion=modules
    dictionary,penalty=fft_bank(),inversion.MaterialPenalty(tv=.0055,amplitude=.0035)
    records,arrays={},{}
    for pack,plan,ob,models,factors,_ in data_and_models(directory,ndt,designs):
        if len(plan.freqs)==1:
            continue
        for idx,cid in enumerate(ob["case_ids"]):
            if args.case and cid not in args.case:
                continue
            y,sigma=designs.read_arm(ob,idx,plan)
            bank,t1,t2=prep(models,y,sigma,dictionary,bright_soft=args.stage=="soft")
            old_beta=[b["beta"] for b in bank]
            if args.stage=="flat":
                average=float(np.mean(old_beta))
                for b in bank:
                    b["beta"]=average
            m,meta=current_solver(models,y,bank,factors,float(ob["reference_amplitude"]),penalty,
                                   ndt.BOUNDS,"fft_joint",args.iterations)
            meta["cpu"] += t1["cpu"]+t2["cpu"]
            meta["wall"] += t1["wall"]+t2["wall"]
            observer=FullForwardCounter(inversion)
            output=observer(models,m,material_factors=factors)
            key=f"{pack}/{plan.name}/{cid}"
            arrays[key]=m
            records[key]=dict(**meta,former_beta=old_beta,new_beta=[b["beta"] for b in bank],
                              bright_ranks=[b["bright"] for b in bank],
                              total_current_ranks=[int(b["q"].shape[1]) for b in bank],
                              physical_noise_ratio=float(np.linalg.norm(output-y)/(sigma*np.sqrt(y.size))))
    path=RESULTS/f"broadband_current_{args.stage}_materials.npz"
    np.savez_compressed(path,**arrays)
    # Scoring is offline after all control materials are saved.
    with np.load(path,allow_pickle=False) as m:
        for pack in PACKS:
            with np.load(directory/"results"/f"{pack}.npz",allow_pickle=False) as file:
                truth={str(cid):z.ravel() for cid,z in zip(file["case_ids"],file["truth"],strict=True)}
            for key,r in records.items():
                p,_,case=key.split("/")
                if p==pack:
                    r["metrics"]=material_score(m[key],truth[case])
    write_json(RESULTS/f"broadband_current_{args.stage}.json",records)


if __name__=="__main__":
    p=ArgumentParser(description=__doc__)
    p.add_argument("stage",choices=("fit","score","flat","soft"))
    p.add_argument("--r15",type=Path)
    p.add_argument("--iterations",type=int,default=24)
    p.add_argument("--case",action="append")
    args=p.parse_args()
    directory,modules=r15_modules(args.r15)
    if args.stage=="fit":
        fit(directory,modules,args)
    elif args.stage=="score":
        score(directory,modules)
    else:
        ablation(directory,modules,args)
