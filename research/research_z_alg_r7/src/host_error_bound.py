"""Propagate the measured scalar host uncertainty to a material-uniform pencil bound."""
import json
import time
import numpy as np
from numpy.polynomial.polynomial import polyval
from patch_pencil import focused_polynomial, noisy_components
from pencil_experiment import CENTER, SCALE, ROOT, REGIONS, ROUNDING, host_model


EPS = 2e-4
PRIOR = [.97, 1.03]
BINS = 512
PATCH = np.array(REGIONS["square"])


def opnorm(a):
    return np.linalg.norm(a,ord=2)


def calibration_range(model,host,cal_feed,reading,estimate):
    C=np.flatnonzero(host)
    L=np.diag(np.sqrt(host[C]))
    K=L@model.G[np.ix_(C,C)]@L
    u=model.S[:,C]@L
    v=L@model.E0[C]@cal_feed
    edges=np.linspace(*PRIOR,BINS+1)
    accepted=[]
    max_lip=0.
    min_denom=1.
    for lo,hi in zip(edges[:-1],edges[1:]):
        center=(lo+hi)/2
        width=(hi-lo)/2
        Z=np.linalg.solve(np.eye(len(C))-center*K,np.eye(len(C)))
        mid=(center*u@Z@v).item()
        denom=1-width*opnorm(Z@K)
        min_denom=min(min_denom,denom)
        if denom<=0:
            accepted.append([lo,hi])
            continue
        lip=opnorm(u@Z)*np.linalg.norm(Z@v)/denom**2
        max_lip=max(max_lip,lip)
        if abs(mid-reading)<=EPS+width*lip+ROUNDING:
            accepted.append([lo,hi])
    if not accepted:
        raise ValueError("Scale prior is inconsistent with the calibration observation.")
    merged=[]
    for lo,hi in accepted:
        if merged and lo==merged[-1][1]:
            merged[-1][1]=hi
        else:
            merged.append([lo,hi])
    delta=max(abs(x-estimate) for span in merged for x in span)
    return dict(components=merged,radius=delta,bins_retained=len(accepted),
                max_slope_envelope=max_lip,min_bin_resolvent_margin=min_denom)


def fitted_fields(model,nominal,estimate,feeds,coords,template,delta):
    host=nominal.m0
    A=nominal.active
    C=np.flatnonzero(host)
    L=np.diag(np.sqrt(host[C]))
    D=estimate*host
    Z=np.linalg.solve(np.eye(model.n)-D[:,None]*model.G,np.eye(model.n))
    full_G=model.G@Z
    full_V=model.S@Z
    Et=np.linalg.solve(np.eye(model.n)-model.G*D,model.E0)
    H=full_G[np.ix_(A,A)]
    gamma=(H-H.conj().T)/(2j)
    eig,U=np.linalg.eigh(gamma)
    if min(eig)<=0:
        raise ValueError("Positive-definite flux Gram required for this bound")
    W=(U/np.sqrt(eig))@U.conj().T
    sqrt_gamma=(U*np.sqrt(eig))@U.conj().T
    direct=opnorm(L@full_G[np.ix_(C,C)]@L)
    left=opnorm(L@full_G[np.ix_(C,A)]@W)
    right=opnorm(W@full_G[np.ix_(A,C)]@L)
    receive=opnorm(full_V[:,C]@L)
    active_receive=opnorm(full_V[:,A]@W)
    kmax=direct+left*right
    vmax=receive+active_receive*right
    rho=kmax*delta
    factor=delta*vmax/(1.-rho) if rho<1. else None
    # Two coefficients of p(t) produce degree-two q(t) p(t).
    prompt=full_G[np.ix_(C,A[PATCH])]@template
    f=Et[C]@feeds@coords
    f[:,:template.shape[1]]+=CENTER*prompt
    f[:,1:template.shape[1]+1]+=SCALE*prompt
    Lf=L@f
    bounds=dict(delta=delta,kmax=kmax,vmax=vmax,rho=rho,neumann_margin=1-rho,c=factor,
                min_gram_eigenvalue=min(eig),max_gram_eigenvalue=max(eig),
                direct_GCC=direct,weighted_GCA=left,weighted_GAC=right,
                weighted_VC=receive,weighted_VA=active_receive)
    return H,full_V[:,A][0],sqrt_gamma,Lf,bounds


def interval_stats(components,truth):
    values=[[CENTER+SCALE*l,CENTER+SCALE*h] for l,h in components]
    return dict(components=values,covered=any(l<=truth<=h for l,h in values),
                whole_box=(values==[[1.,15.]]),length=sum(h-l for l,h in values),
                width=values[-1][1]-values[0][0] if values else None)


def summarize(rows,key):
    states=[row[key] for row in rows]
    widths=[s["width"] for s in states if s["width"] is not None]
    return dict(n=len(states),covered=sum(s["covered"] for s in states),
                empty=sum(s["width"] is None for s in states),whole_box=sum(s["whole_box"] for s in states),
                width_min=min(widths,default=None),width_median=float(np.median(widths)) if widths else None,
                width_max=max(widths,default=None),median_union_length=float(np.median([s["length"] for s in states])))


def run():
    start=time.perf_counter()
    cpu=time.process_time()
    folder=ROOT/"results"
    cal=json.loads((folder/"geometry_calibrated.json").read_text())
    obs=np.load(folder/"geometry_calibrated_raw.npz")
    sample=np.load(folder/"pencil_observations.npz")
    rows=[]
    layouts={}
    for layout in ["primary","confirmation"]:
        model,kernel,_,_=host_model(layout)
        estimate=cal["layouts"][layout]["calibration"]["fitted_scale"]
        lut=calibration_range(model,kernel.m0,obs[layout+"_cal_feed"],obs[layout+"_cal_reading"].item(),estimate)
        prefix=layout+"_one_calibration"
        F=obs[prefix+"_feeds"]
        coords=obs[prefix+"_coords"]
        P=obs[prefix+"_templates"]
        design=dict(coord=coords,templates=P)
        H,R,root_g,Lf,constants=fitted_fields(model,kernel,estimate,F,coords,P,lut["radius"])
        c=constants["c"]
        if c is not None:
            rhs_coeff=np.sqrt(2)*np.vstack([(EPS+ROUNDING)*coords,c*Lf])
        layout_rows=[]
        for i,x in enumerate(sample[layout+"_square_truth_parameters"]):
            # Truth labels enter below for scoring, never for scale/calibration or interval generation.
            name=f"{layout}_square_{i:02d}_pencil_one_calibration"
            d=obs[name+"_response"]
            g=focused_polynomial(d,R[PATCH],design,CENTER,SCALE)
            nominal=noisy_components(g,coords,EPS+ROUNDING)
            outer=noisy_components(g,rhs_coeff,1.) if c is not None else [[-1.,1.]]
            q=x[0]
            t=(q-CENTER)/SCALE
            a=polyval(t,coords.T)
            sys_res=abs(obs[name+"_model_error"]@a)
            projected_total=abs(polyval(t,g))
            field_norm=np.linalg.norm(polyval(t,Lf.T))
            systematic=c*field_norm if c is not None else None
            noise_term=(EPS+ROUNDING)*np.linalg.norm(a)
            system_bound_pass=sys_res<=systematic+1e-12 if c is not None else True
            combined_bound_pass=projected_total<=noise_term+systematic+1e-12 if c is not None else True
            z=np.full(len(H),q)
            z[np.setdiff1d(np.arange(len(H)),PATCH)]=x[1:]
            T=np.linalg.solve(np.eye(len(H))-z[:,None]*H,np.diag(z))
            normT=opnorm(root_g@T@root_g)
            row=dict(layout=layout,case=i,truth=q,nominal=interval_stats(nominal,q),
                     outer=interval_stats(outer,q),projected_systematic=sys_res,
                     systematic_envelope=systematic,sys_over_bound=sys_res/systematic if systematic else None,
                     projected_total=projected_total,noise_envelope=noise_term,system_bound_pass=bool(system_bound_pass),
                     combined_bound_pass=bool(combined_bound_pass),weighted_T_norm=normT)
            rows.append(row)
            layout_rows.append(row)
        # Simulated hidden scale enters here only for scoring the LUT's data-only set.
        true_scale=cal["simulated_host_scale"]
        layouts[layout]=dict(estimate=estimate,lut=lut,constants=constants,
            scale_covered=any(l<=true_scale<=u for l,u in lut["components"]),
            nominal=summarize(layout_rows,"nominal"),outer=summarize(layout_rows,"outer"),
            max_sys_over_bound=max((r["sys_over_bound"] for r in layout_rows if r["sys_over_bound"] is not None),default=None),
            max_weighted_T_norm=max(r["weighted_T_norm"] for r in layout_rows))
    result=dict(prior=PRIOR,calibration_bins=BINS,calibration_noise_radius=EPS,target_noise_radius=EPS,
        numerical_margin=ROUNDING,outer_rule="sqrt(2) * norm([eps alpha; c L f_C]); physical scalar real host family, float kernels.",
        layouts=layouts,nominal=summarize(rows,"nominal"),outer=summarize(rows,"outer"),
        system_projected_bound_pass=sum(row["system_bound_pass"] for row in rows),
        total_projected_bound_pass=sum(row["combined_bound_pass"] for row in rows),
        lut_24x24_center_solves=2*BINS,new_full_loaded_solves=0,scoring_12x12_T_solves=len(rows),
        rows=rows,wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu)
    (folder/"host_error_bound.json").write_text(json.dumps(result,indent=2)+"\n")


if __name__=="__main__":run()
