"""Alternating current correction and full-pixel material recovery on archived reads."""
from argparse import ArgumentParser
from collections import Counter
import json
import os
from pathlib import Path
from time import perf_counter, process_time

for name in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    os.environ[name]="1"

import numpy as np
from scipy.fft import ifft2
from scipy.linalg import cho_factor, cho_solve, lu_factor, lu_solve, svd
from scipy.optimize import minimize
from scipy.stats import chi2

from cohorts import load_cohort, load_truth
from score_materials import material_metrics
from voxel_gn import Prior, Wave, voxel_gn

OUT = Path(__file__).resolve().parents[1]/"results"
METHODS = ("fixed_current", "fft_dark", "fft_bright", "full_current", "voxel_only")


def fourier_modes(shape, halfwidth):
    a,b = np.meshgrid(np.arange(-halfwidth,halfwidth+1),
                      np.arange(-halfwidth,halfwidth+1),indexing="ij")
    impulses = np.zeros((*shape,a.size),dtype=complex)
    impulses[a.ravel()%shape[0],b.ravel()%shape[1],np.arange(a.size)] = 1
    return ifft2(impulses,axes=(0,1),norm="ortho").reshape(np.prod(shape),-1)


def norm2(a):
    return float(np.vdot(a,a).real)


def rms(a):
    return float(np.linalg.norm(a)/np.sqrt(a.size))


class CurrentStep:
    def __init__(self, model, case, method, halfwidth, tradeoff):
        setup = perf_counter()
        self.G,self.S,self.E0 = [model[k] for k in ("G","S","E0")]
        self.Y = case["Y"]
        self.method = method
        self.counts = Counter()
        left,s,right = svd(self.S,full_matrices=False,check_finite=False)
        z = left.conj().T@self.Y
        detectable = np.flatnonzero((s>=.05*s[0]) & (np.linalg.norm(z,axis=1)/np.sqrt(z.shape[1])>=2.75*case["sigma"]))
        self.bright = int(detectable[-1]+1) if detectable.size else 0
        v = right[:self.bright].conj().T
        self.anchor = v@(z[:self.bright]/s[:self.bright,None])
        self.field_anchor = self.E0+self.G@self.anchor
        self.data_rhs = self.Y-self.S@self.anchor
        self.beta = (tradeoff*s[0])**2
        self.rank = 0
        if method in ("fft_dark","fft_bright","fft_dark_pixel"):
            low = fourier_modes(model["shape"],halfwidth)
            q,se,_ = svd(low-v@(v.conj().T@low),full_matrices=False,check_finite=False)
            q = q[:,se>1e-9]
            self.q = np.column_stack((v,q)) if method=="fft_bright" else q
            self.rank = self.q.shape[1]
            self.Gq = self.G@self.q
            obs = self.S@self.q
            self.normal_data = obs.conj().T@obs
            self.observed_projection = obs.conj().T@self.data_rhs
        elif method=="full_current":
            self.eye = np.eye(len(self.G),dtype=complex)
            self.rank = len(self.G)
        self.time_setup = perf_counter()-setup

    def update(self, material):
        if self.method=="fixed_current":
            return self.anchor
        if self.method=="full_current":
            A = self.eye - material[:,None]*self.G
            factor = lu_factor(A,check_finite=False)
            physical = lu_solve(factor,material[:,None]*self.E0,check_finite=False)
            transfer = lu_solve(factor,self.S.T,trans=1,check_finite=False).T
            normal = transfer@transfer.conj().T
            normal.flat[::len(normal)+1] += self.beta
            dual = cho_solve(cho_factor(normal,check_finite=False),self.Y-self.S@physical,check_finite=False)
            corrected = lu_solve(factor,transfer.conj().T@dual,check_finite=False)
            self.counts.update(full_factorizations=1,full_complex_rhs=2*self.E0.shape[1]+len(self.S),
                               small_factorizations=1,small_complex_rhs=self.E0.shape[1])
            return physical+corrected
        lin = self.q-material[:,None]*self.Gq
        goal = material[:,None]*self.field_anchor-self.anchor
        normal = self.normal_data+self.beta*(lin.conj().T@lin)
        rhs = self.observed_projection+self.beta*(lin.conj().T@goal)
        weights = cho_solve(cho_factor(normal,check_finite=False),rhs,check_finite=False)
        self.counts.update(small_factorizations=1,small_complex_rhs=self.E0.shape[1])
        return self.anchor+self.q@weights

    def stationarity(self, material, current):
        field = self.E0+self.G@current
        defect = current-material[:,None]*field
        row = self.S.conj().T@(self.S@current-self.Y)
        grad = row+self.beta*(defect-self.G.conj().T@(material[:,None]*defect))
        denom = np.linalg.norm(row)+self.beta*np.linalg.norm(defect-self.G.conj().T@(material[:,None]*defect))
        if self.method=="fixed_current":
            return None
        return float(np.linalg.norm(grad if self.method=="full_current" else self.q.conj().T@grad)
                     /max(denom,1e-30))


def material_update(field, current, previous, prior, scale, beta, bounds, pixel):
    energy = np.sum(abs(field)**2,axis=1)
    response = np.real(np.sum(field.conj()*current,axis=1))
    ratio = np.divide(response,energy,out=previous.copy(),where=energy>0)
    start = np.clip(ratio,*bounds)
    if pixel:
        return start,dict(evaluations=0,status="pixel_closed_form",gradient_mapping=None)
    a = beta*energy/scale
    b = beta*response/scale
    n = len(previous)

    def fun(x):
        reg,edge_gradient,_ = prior.evaluate(x)
        grad = a*x-b + np.bincount(prior.left,edge_gradient,minlength=n) - np.bincount(
            prior.right,edge_gradient,minlength=n)
        # Drop the fixed term beta ||J||² / (2||Y||²).
        return float(np.dot(x,.5*a*x-b)+reg),grad

    answer = minimize(fun,start,method="L-BFGS-B",jac=True,bounds=[bounds]*n,
                      options={"maxiter":65,"ftol":1e-13,"gtol":3e-9,"maxcor":10})
    chosen = answer.x if fun(answer.x)[0]<=min(fun(start)[0],fun(previous)[0]) else min(
        (previous,start),key=lambda x:fun(x)[0])
    grad = fun(chosen)[1]
    projected = chosen-np.clip(chosen-grad,*bounds)
    return chosen,dict(evaluations=answer.nfev,status=int(answer.status),
                       gradient_mapping=float(np.linalg.norm(projected,np.inf)))


def alternating(model,case,method,iterations,halfwidth,tradeoff):
    wall,cpu = perf_counter(),process_time()
    prior = Prior(model["shape"])
    pixel = method=="fft_dark_pixel"
    step = CurrentStep(model,case,method,halfwidth,tradeoff)
    m = np.clip(np.ones(len(step.G)),*model["bounds"])
    scale = norm2(case["Y"])
    current_time = material_time = 0.
    trace, saved = [], {}
    before = float("inf")
    monotonic = True
    stop = "iteration_cap"
    design_check = None
    for it in range(iterations):
        t = perf_counter()
        current = step.update(m)
        current_time += perf_counter()-t
        if it==0:
            design_check = step.stationarity(m,current)
        t = perf_counter()
        field = step.E0+step.G@current
        m, info = material_update(field,current,m,prior,scale,step.beta,model["bounds"],pixel)
        material_time += perf_counter()-t
        data = norm2(step.S@current-step.Y)
        state = norm2(current-m[:,None]*field)
        reg = 0. if pixel else prior.evaluate(m)[0]
        obj = .5*(data+step.beta*state)/scale+reg
        monotonic &= obj <= before+1e-12
        trace.append(dict(iteration=it+1,objective=obj,rel_auxiliary_data=float(np.sqrt(data/scale)),
                          rel_constitutive=float(np.sqrt(state/max(norm2(current),1e-300))),
                          wall=perf_counter()-wall,**info))
        if it+1 in (1,4,8,16):
            saved[str(it+1)] = m.copy()
        if before-obj < 1e-11 and obj<=before+1e-12:
            stop = "gain_tolerance"
            break
        before = obj
    online = dict(wall=perf_counter()-wall,cpu=process_time()-cpu)
    # Validation only. Its full solve was not used to accept a current/material update.
    wave = Wave(model)
    final_field = wave.state(m)[0]
    diagnostic = float(np.linalg.norm(final_field-case["Y"])/(case["sigma"]*np.sqrt(case["Y"].size)))
    costs = {k:step.counts[k] for k in ("full_factorizations","full_complex_rhs","small_factorizations","small_complex_rhs")}
    result = dict(method=method,cohort_id=case["id"],current_dimension=step.rank,
                  free_real_material_coordinates=len(m),dark_or_bright_cut=step.bright,
                  state_tradeoff=tradeoff,lambda_state=float(np.sqrt(step.beta)),stop=stop,
                  monotone_objective=bool(monotonic),first_current_projected_kkt_relative=design_check,
                  residual_over_sigma_rms=diagnostic,near_noise_threshold=float(
                      np.sqrt(chi2.ppf(.99,2*case["Y"].size)/(2*case["Y"].size))),
                  online_cost=dict(**online,setup_wall=step.time_setup,current_wall=current_time,
                                   material_wall=material_time,**costs),
                  final_validation=dict(full_factorizations=wave.lu_count,full_complex_rhs=wave.rhs),
                  active_count=int(np.sum((m<=model["bounds"][0]+1e-12)|(m>=model["bounds"][1]-1e-12))),
                  trace=trace)
    return m,result,saved


def fit(a):
    metadata = {k:v for k,v in vars(a).items() if k!="stage"}
    rows, outputs = [], {}
    for family in a.families:
        model,cases=load_cohort(family)
        for case in cases:
            if a.cases and case["id"] not in a.cases:
                continue
            for method in a.methods:
                if method=="voxel_only":
                    m,gn = voxel_gn(model,case,12+72*a.iterations,a.iterations)
                    result = dict(method=method,free_real_material_coordinates=len(m),current_dimension=0,
                                  near_noise_threshold=float(np.sqrt(chi2.ppf(.99,2*case["Y"].size)/(2*case["Y"].size))),
                                  **gn)
                    saved={}
                else:
                    m,result,saved=alternating(model,case,method,a.iterations,a.halfwidth,a.tradeoff)
                key = case["id"]+"/"+method
                outputs.update({key+"/"+it:val for it,val in saved.items()})
                outputs[key]=m
                rows.append(dict(cohort=family,case=case["id"],array=key,shape=model["shape"],**result))
                print(case["id"],method,round(result["residual_over_sigma_rms"],4),
                      round(result["online_cost"]["cpu"],3),len(result["trace"]),flush=True)
    np.savez_compressed(OUT/f"alternating{a.tag}.npz",**outputs)
    (OUT/f"alternating{a.tag}.json").write_text(
        json.dumps(dict(settings=metadata,cases=rows),indent=2,allow_nan=False)+"\n")


def score(a):
    stored = json.loads((OUT/f"alternating{a.tag}.json").read_text())
    scores = {}
    with np.load(OUT/f"alternating{a.tag}.npz",allow_pickle=False) as x:
        truth_cache={}
        for r in stored["cases"]:
            ci=(r["cohort"],r["case"])
            if ci not in truth_cache:
                truth_cache[ci]=load_truth(*ci).ravel()
            known,estimated=truth_cache[ci],x[r["array"]]
            support = max(float(np.max(known))/2,.02)
            mask,guess = known>support,estimated>support
            metrics = material_metrics(estimated,known,shape=tuple(r["shape"]))
            truth_norm = np.linalg.norm(known)
            checkpoints = {k.rsplit("/",1)[1]:float(np.linalg.norm(x[k]-known)/truth_norm)
                           for k in x.files if k.startswith(r["array"]+"/")}
            scores[r["array"]] = dict(
                **{k:metrics[k] for k in ("relativeL2","RMSE","amplitude_NRMSE","fixed_support_IoU")},
                half_peak_support_IoU=float(np.sum(mask&guess)/max(1,np.sum(mask|guess))),
                checkpoint_relativeL2=checkpoints)
            print(r["array"],round(metrics["relativeL2"],4),round(metrics["amplitude_NRMSE"],4),
                  round(metrics["fixed_support_IoU"],4),flush=True)
    (OUT/f"alternating{a.tag}_scores.json").write_text(json.dumps(scores,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    p = ArgumentParser(description=__doc__)
    p.add_argument("stage",choices=("fit","score"))
    p.add_argument("--families",nargs="+",default=["r10_twod"])
    p.add_argument("--cases",nargs="+")
    p.add_argument("--methods",nargs="+",choices=(*METHODS,"fft_dark_pixel"),default=list(METHODS))
    p.add_argument("--iterations",type=int,default=28)
    p.add_argument("--halfwidth",type=int,default=3)
    p.add_argument("--tradeoff",type=float,default=.32)
    p.add_argument("--tag",default="")
    a=p.parse_args()
    if a.stage=="fit":
        fit(a)
    else:
        score(a)
