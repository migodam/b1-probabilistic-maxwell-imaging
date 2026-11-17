"""Read a homogeneous patch of a nonlinear TM VIE using a rank-deficient feed bank."""
import argparse
from pathlib import Path
import json
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import least_squares
from scipy.special import j0
from scipy.stats import qmc

from geometry_controls import physics_module, power, ROOT, whitening


K = 12
D = 10
Q0 = 12.
# Nested connected patches, from two central sites outward over the top two rows.
ORDER = np.array([5,6,2,1,0,4,8])
SIZES = range(1,8)
SIGMAS = [0., 1e-4, 1e-3]
LABELS = ["nominal_optimal", "random_valid", "projected_uniform"]
HET = [0., .1, .5, 1.]
SAMPLES = 128


def geometry(physics):
    grid = np.array([[x,y] for y in [-.22,0.,.22] for x in [-.33,-.11,.11,.33]])
    rng = np.random.default_rng(7111802)
    active = grid + rng.uniform(-.016,.016,grid.shape)
    theta = np.arange(5)*2*np.pi/5+.34
    host = .59*np.column_stack([np.cos(theta), np.sin(theta)])
    m0 = np.r_[np.zeros(K),[2.,2.5,3.,1.8,2.2]]
    recv = [[.21,.91]]
    az = .27+np.arange(D)*2*np.pi/D
    feed = 1.03*np.column_stack([np.cos(az),np.sin(az)])
    centers = np.vstack([active,host])
    m = physics.Model(centers=centers,sources=feed,receivers=recv,radius=.039,k=2*np.pi)
    kernel = physics.ReducedKernel(m,m0,active=np.arange(K))
    gram = j0(2*np.pi*np.linalg.norm(feed[:,None]-feed[None,:],axis=-1))
    return m,m0,kernel,gram


def normed(x):
    return x / np.linalg.norm(x) * np.exp(-1j*np.angle(x[np.argmax(abs(x))]))


def forward(model, host, z):
    m = host.copy()
    m[:K] = z
    return model.forward(m, independent=True)["y"][0]


def choose_mode(P, vrow, synthesis, gram):
    # Maximize nominal |v p| per sqrt(sum of the two scaled shot variances).
    eigmetric = synthesis.conj().T @ gram @ synthesis
    W = whitening(eigmetric)
    return normed(P @ W @ (vrow @ P @ W).conj())


def choose_two(P, vrow, basis_a, basis_h, gram):
    F, G = basis_a@P,basis_h@P
    metric = F.conj().T@gram@F + Q0**2*(G.conj().T@gram@G)
    W = whitening(metric)
    return normed(P @ W @ (vrow @ P @ W).conj())


def total_ratio(y, a, tau, Rp):
    reads=y@a
    denom = Rp + tau[1]*reads[...,1]
    return tau[0]*reads[...,0]/denom, denom


def store_two(raw, key, coeff, arrays):
    for name, arr in arrays.items():
        raw[f"{key}_{name}"] = arr
    raw[f"{key}_mode"] = coeff


def fit_same_data(kernel, idx, true_read, a):
    outside = np.setdiff1d(np.arange(K), idx)
    def material(x):
        z=np.full(K,x[0])
        z[outside]=x[1:]
        return z
    def funjac(x, derivative=False):
        z = material(x)
        L = np.eye(K)-z[:,None]*kernel.K
        j = np.linalg.solve(L,z[:,None]*(kernel.q@a))
        if derivative:
            v = np.linalg.solve(L.T,kernel.L.T).T[0]
            Jcells = (v[:,None]*(kernel.q@a+kernel.K@j)).T
            J = np.column_stack([Jcells[:,idx].sum(axis=1),Jcells[:,outside]])
            return np.vstack([J.real,J.imag])
        err = (kernel.L@j)[0]-true_read
        return np.r_[err.real,err.imag]
    res=least_squares(funjac,np.full(len(outside)+1,Q0),jac=lambda x:funjac(x,True),
                      method="trf",bounds=(.5,30.),ftol=1e-13,xtol=1e-13,gtol=1e-14,max_nfev=300)
    return [res.x[0],np.linalg.norm(res.fun),res.nfev,res.status]


def first_pair(model, host, kernel, idx, a, tau, p, mode_label):
    """Compare a +1 patch change after fitting every unknown exterior cell."""
    outside = np.setdiff1d(np.arange(K), idx)
    z0 = np.full(K, Q0)
    y0 = forward(model, host, z0) @ a
    def fit(q1, fixed_outside):
        z1 = z0.copy()
        z1[idx] = q1
        free = outside
        if fixed_outside:
            z1[outside[0]] += 4.
            free = outside[1:]
        def residual(x):
            z1[free] = x
            y1 = kernel.evaluate(z1)[0] @ a
            diff = y1-y0
            return np.r_[diff.real,diff.imag]
        fit = least_squares(residual,z0[free],bounds=(.5,30.),gtol=1e-15,ftol=1e-15,
                            xtol=1e-15,max_nfev=1500)
        z1[free]=fit.x
        y1 = forward(model,host,z1) @ a
        gap = np.linalg.norm(y1-y0)
        return dict(z0=z0.tolist(),z1=z1.tolist(), data0=np.stack((y0.real,y0.imag),axis=-1).tolist(),
                    data1=np.stack((y1.real,y1.imag),axis=-1).tolist(), full_data_gap=float(gap),
                    max_nuisance_change=float(max(abs(z1[outside]-z0[outside]))), nfev=fit.nfev)
    pair=fit(Q0+1.,False)
    nuisance=fit(Q0,True)
    if mode_label=="nominal_optimal":
        centered=y0-kernel.y0[0]@a
        visibility=kernel.L[0]@p+tau[1]*centered[1]
        bound=abs(visibility)/np.linalg.norm(tau*[1.,Q0+1.])
        pair["analytic_lower_bound_for_q_plus_one"]=float(bound)
    return dict(patch_shifted=pair, patch_fixed_nuisance_altered=nuisance)


def experiment(physics):
    start = time.perf_counter()
    cpu = time.process_time()
    model,host,kernel,gram = geometry(physics)
    B,H,R = kernel.q,kernel.K,kernel.L[0]
    C = null_space(B.conj().T).conj().T
    # B has full column rank; its unique LS feed realizes the orthogonal projection.
    drive = np.linalg.pinv(B)
    zref=np.full(K,Q0)
    fixed=forward(model,host,zref)
    vrow=np.linalg.solve((np.eye(K)-zref[:,None]*H).T,R)
    rng=np.random.default_rng(181102)
    samples=qmc.Sobol(K+1,scramble=True,seed=2026111811).random_base2(7)
    ztarget=2.+26.*samples[:,0]
    zothers=.5+29.5*samples[:,1:]
    ref_host=kernel.y0[0]
    error=(rng.standard_normal((SAMPLES,2))+1j*rng.standard_normal((SAMPLES,2)))/np.sqrt(2)
    raw=dict(centers=model.centers,receivers=model.receivers,sources=model.sources,gram=gram,
             H=H,B=B,R=R,host=host,target=ztarget,outside_samples=zothers,
             raw_unit_noise=error,labels=LABELS,sigmas=SIGMAS,heterogeneity_peak=HET,
             patch_order=ORDER,nominal=Q0)
    summary=dict(k=K,d=D,samples=SAMPLES,box=[.5,30.],target_range=[2.,28.],
                 nominal=Q0,cell_order_zero_based=ORDER.tolist(),cond_B=float(np.linalg.cond(B)), rows={})
    for r in SIZES:
        idx=ORDER[:r]
        E=np.eye(K)[:,idx]
        fixed_cond=C @ ((np.eye(K)-Q0*H)@E)
        simultaneous=np.vstack([C@E,C@H@E])
        one=null_space(fixed_cond)
        two=null_space(simultaneous)
        row=dict(dim_one=one.shape[1],dim_two=two.shape[1],
                 one_singular_values=np.linalg.svd(fixed_cond,compute_uv=False).tolist(),
                 two_singular_values=np.linalg.svd(simultaneous,compute_uv=False).tolist())
        summary["rows"][str(r)]=row
        # Identical nuisance and target for all controls at this patch size.
        ytrue=np.empty((SAMPLES,D),complex)
        yfixed=np.empty_like(ytrue)
        varied_z=zothers.copy()
        varied_z[:,idx]=ztarget[:,None]
        pattern=np.cos(np.arange(r)*1.618+.31)
        pattern-=pattern.mean()
        pattern/=max(1e-20,np.max(abs(pattern)))
        yhet=np.empty((len(HET),SAMPLES,D),complex)
        for s,z in enumerate(varied_z):
            ytrue[s]=forward(model,host,z)-ref_host
            qmat=z.copy()
            qmat[idx]=Q0
            yfixed[s]=forward(model,host,qmat)
            for h_i,h in enumerate(HET):
                if h==0 or r==1:
                    yhet[h_i,s]=ytrue[s]
                else:
                    zhet=z.copy()
                    zhet[idx]+=h*pattern
                    yhet[h_i,s]=forward(model,host,zhet)-ref_host
        key=f"r{r}"
        raw.update({key+"_z":varied_z,key+"_ytrue":ytrue,key+"_yfixed":yfixed,
                    key+"_heterogeneity_pattern":pattern,key+"_yhetero":yhet})
        if one.shape[1]:
            candidates=E@one
            targetmap=(np.eye(K)-Q0*H)
            p=choose_mode(candidates,vrow,drive@targetmap@candidates,gram)
            one_label="reachable"
        else:
            # Negative control: choose best LS approximation to the nominal desired B a form.
            p=E@np.linalg.svd(fixed_cond,full_matrices=True)[2][-1].conj()
            one_label="least_unreachable_component"
        b=(np.eye(K)-Q0*H)@p
        a=drive@b
        miss=np.linalg.norm(B@a-b)/np.linalg.norm(b)
        a/=np.sqrt(power(a,gram))
        drift=(yfixed-fixed)@a
        raw.update({key+"_one_mode":p,key+"_one_drive":a,key+"_one_null_drift":drift})
        row.update(one_mode=one_label,one_defect=float(miss),one_drift_max=float(max(abs(drift))),
                   one_drift_median=float(np.median(abs(drift))))
        row["two"]={}
        modes={}
        if two.shape[1]:
            P = E@two
            modes["nominal_optimal"]=choose_two(P,vrow,drive,drive@H,gram)
            modes["random_valid"]=normed(P@(rng.normal(size=P.shape[1])+1j*rng.normal(size=P.shape[1])))
        modes["projected_uniform"]=E@(np.ones(r)/np.sqrt(r))
        for label,p in modes.items():
            targets=np.column_stack([p,H@p])
            a=drive@targets
            defects=np.linalg.norm(B@a-targets,axis=0)/np.linalg.norm(targets,axis=0)
            tau=np.array([np.sqrt(power(x,gram)) for x in a.T])
            a/=tau
            Rp=R@p
            values=np.empty((len(SIGMAS),SAMPLES),complex)
            reads=ytrue@a
            denom=Rp+tau[1]*reads[:,1]
            for j,sigma in enumerate(SIGMAS):
                noisy=reads+sigma*error
                values[j]=tau[0]*noisy[:,0]/(Rp+tau[1]*noisy[:,1])
            het_values,_=total_ratio(yhet,a,tau,Rp)
            clean_err = abs(values[0]-ztarget)
            # The specified real material estimate is Re(complex quotient), without clipping.
            med_noise=np.median(abs(values.real-ztarget),axis=1)
            p95_noise=np.quantile(abs(values.real-ztarget),.95,axis=1)
            detuning=np.median(abs(het_values.real-ztarget),axis=1)
            nom_noise_norm=np.linalg.norm(tau*np.array([1.,Q0]))/abs(vrow@p)
            record=dict(defects=defects.tolist(),tau=tau.tolist(),nominal_noise_gain=float(nom_noise_norm),
                        min_abs_denom=float(min(abs(denom))),median_abs_denom=float(np.median(abs(denom))),
                        feed_norm=np.linalg.norm(a,axis=0).tolist(),
                        median_clean_complex_error=float(np.median(clean_err)),max_clean_complex_error=float(max(clean_err)),
                        median_real_error_by_noise=med_noise.tolist(),p95_real_error_by_noise=p95_noise.tolist(),
                        median_real_error_by_heterogeneity=detuning.tolist())
            store_two(raw,key+"_"+label,p,dict(drives=a,scale=tau,readout=reads,
                      denominator=denom,estimates=values,heterogeneity_estimates=het_values))
            if label=="nominal_optimal":
                htargets=varied_z.copy()
                htargets[:,idx]+=HET[-1]*pattern
                v=np.array([np.linalg.solve((np.eye(K)-z[:,None]*H).T,R) for z in htargets])
                weights=v[:,idx]*p[idx]
                av=np.sum(htargets[:,idx]*weights,axis=1)/weights.sum(axis=1)
                record["heterogeneous_weighted_quotient_max_difference"]=float(max(abs(av-het_values[-1])))
                outside=(het_values[-1].real<np.min(htargets[:,idx],axis=1))|(het_values[-1].real>np.max(htargets[:,idx],axis=1))
                record["outside_patch_range_count_at_peak_1"]=int(np.sum(outside))
                raw[key+"_"+label+"_heterogeneous_weights_peak_1"]=weights
            if r==6 and label in ("nominal_optimal","projected_uniform"):
                fitstart=time.perf_counter()
                # Identical two complex observations and box, single nominal start for 24 cases.
                trf=np.array([fit_same_data(kernel,idx,y,a) for y in reads[:24]])
                record.update(trf_same_data_24=dict(
                    median_patch_error=float(np.median(abs(trf[:,0]-ztarget[:24]))),
                    max_patch_error=float(max(abs(trf[:,0]-ztarget[:24]))),
                    max_residual=float(max(trf[:,1])),median_nfev=float(np.median(trf[:,2])),
                    count_hit_limit=int(np.sum(trf[:,3]==0)),wall_seconds=time.perf_counter()-fitstart))
                raw[key+"_"+label+"_trf_24"]=trf
                record["nominal_material_pair"]=first_pair(model,host,kernel,idx,a,tau,p,label)
            row["two"][label]=record
    np.savez_compressed(ROOT/"results/geometry_patch_raw.npz",**raw)
    summary["full_forward_states"]=model.counts["INDEPENDENT_FORWARD"]
    make_plot(raw,summary)
    summary["wall_seconds_with_plots"]=time.perf_counter()-start
    summary["cpu_seconds_with_plots"]=time.process_time()-cpu
    (ROOT/"results/geometry_patch.json").write_text(json.dumps(summary,indent=2)+"\n")


def make_plot(raw,summary):
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})
    fig, axes=plt.subplots(1,3,figsize=(14.4,4.1),layout="constrained")
    r=list(SIZES)
    rows=summary["rows"]
    axes[0].plot(r,[rows[str(i)]["dim_one"] for i in r],"o-",color="#1a797b",label="One prescribed value")
    axes[0].plot(r,[rows[str(i)]["dim_two"] for i in r],"s-",color="#ae5528",label="Two-shot readout")
    axes[0].set(xlabel="Number of cells sharing material",ylabel="Reachable mode dimension",title="A  Ten feeds for twelve active cells",xticks=r)
    axes[0].legend(frameon=False,fontsize=9)
    for key,clr,lbl in zip(LABELS,["#007181","#8a83a4","#bb651d"],["Nominal optimum","Random reachable","Projected uniform"]):
        sel=[i for i in r if key in rows[str(i)]["two"]]
        axes[1].plot(sel,[rows[str(i)]["two"][key]["median_clean_complex_error"] for i in sel],"o",color=clr,label=lbl,ls="--" if key=="random_valid" else "-")
        v=abs(raw["r6_"+key+"_estimates"][1].real-raw["target"])
        axes[2].plot(np.sort(v),np.arange(1,len(v)+1)/len(v),color=clr,label=lbl,ls="--" if key=="random_valid" else "-")
    axes[1].set(yscale="log",xlabel="Number of cells sharing material",ylabel="Median complex material error",title="B  Noiseless quotient",xticks=r)
    axes[1].legend(frameon=False,fontsize=8)
    axes[2].set(xscale="log",xlabel="Absolute error in real material value",ylabel="Empirical cumulative fraction",title="C  Six-cell patch, noise RMS = 0.0001")
    for ax in axes: ax.grid(True,alpha=.18)
    for ext in ["pdf","png"]:fig.savefig(ROOT/f"figures/geometry_patch.{ext}",dpi=200)
    plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--physics",type=Path)
    experiment(physics_module(parser.parse_args().physics))
