"""Native upper witnesses for the *relaxed* supported T-operator ball, not physical media."""
import json
import numpy as np
from scipy.linalg import eigh,svd
from flint import arb,ctx
from threadpoolctl import threadpool_limits
from branches import expand,as_branch
from native_box import HERE,clock,elapsed,ball,star,gershgorin,hi,frobenius_squared,herm
from site_driver import DEFAULT_PUBLIC

def norm(x):
    return frobenius_squared(x).sqrt()

def spectral_upper(x):
    xx=star(x)*x
    mat=np.array([[complex(v) for v in row] for row in xx.tolist()])
    q=ball(eigh(herm(mat),check_finite=False)[1])
    return (gershgorin(star(q)*xx*q,upper=True)/gershgorin(star(q)*q,upper=False)).sqrt().upper()

def clip(x):
    u,d,v=svd(x,full_matrices=False,check_finite=False)
    return (u*np.minimum(d,.99995))@v

def project_data(x,k,obs,r):
    e=obs-k*x
    if np.linalg.norm(e)<=r: return x
    kk=k*k
    def noise(t):
        return np.linalg.norm(e/(1+t*kk))
    hi_t=1.
    while noise(hi_t)>r and hi_t<1.e60:hi_t*=16.
    low=0.
    for _ in range(36):
        mid=(low+hi_t)/2
        if noise(mid)>r:low=mid
        else:hi_t=mid
    return x+hi_t*k*e/(1+hi_t*kk)

def candidate(k,obs,r):
    x=np.zeros_like(obs);z=np.zeros_like(obs)
    for step in range(120):
        x=clip(z)
        if np.linalg.norm(k*x-obs)<.9995*r:break
        y=project_data(2*x-z,k,obs,.998*r)
        z+=y-x
    return x,float(np.linalg.norm(k*x-obs)),step+1

class Controls:
    def __init__(self,s,b,C):
        u,self.d,vh=svd(s,full_matrices=False,check_finite=False)
        ub,self.e,wh=svd(b,full_matrices=False,check_finite=False)
        self.u,self.w=u,wh.conj().T
        self.v=vh.conj().T; self.p=ub
        self.s,self.b=s,b;self.C=C
        self.k=C*np.outer(self.d,self.e)
        self.form=None
    def verify(self,x,y):
        if self.form is None:
            v,p=ball(self.v),ball(self.p)
            self.form=ball(self.s)*v,star(p)*ball(self.b),spectral_upper(v)*spectral_upper(p)
        left,right,vp=self.form
        x=ball(x)
        bound=(vp*spectral_upper(x)).upper()
        e=norm(ball(y)-arb(self.C)*left*x*right).upper()
        return hi(bound),hi(e)

def main():
    ctx.prec=112;ctx.threads=1
    start=clock()
    online=json.loads((HERE/"site_online.json").read_text())
    manifest=json.loads((DEFAULT_PUBLIC/"manifest.json").read_text())
    geometry=np.load(DEFAULT_PUBLIC/"geometry96.npz")
    trial,cache=[],{}
    witnesses={}
    for f in manifest["cases"]:
        obs=f["case"];ap=f["array"]; old=online["predictions"][obs]
        with np.load(DEFAULT_PUBLIC/f["measurement_file"]) as raw:
            y,r=raw["Y"],float(raw["delta"])
        for item in old["tile_scores"]:
            p=item["primary"]
            if p["rigor"]=="strict_native_translated_tile_port_bound":continue
            coords=tuple(item["tile"]);key=ap,coords
            if key not in cache:
                caps=np.zeros((4,4));caps[coords]=4.
                ix=expand(as_branch("one",caps),96)>0
                C=next(s["primary"]["current_norm_factor_upper"] for s in old["tile_scores"] if
                       "current_norm_factor_upper" in s["primary"])
                cache[key]=Controls(geometry[ap+"__S"][:,ix],geometry[ap+"__B"][ix],C)
            ctl=cache[key]
            m=ctl.u.conj().T@y@ctl.w
            x,err,it=candidate(ctl.k,m,r)
            row=dict(observation=obs,tile=coords,original_decision=p["decision"],radius=r,iterations=it,
                     proposed_relaxed_residual=err,current_gain_C=ctl.C,strict_energy_only_fit=False)
            if err<r:
                upper,dist=ctl.verify(x,y)
                row.update(T_norm_over_C_upper=upper,relaxed_data_distance_upper=dist,
                           strict_energy_only_fit=(upper<=1. and dist<=r))
                if row["strict_energy_only_fit"]:
                    key=f"{obs}__tile{coords[0]}_{coords[1]}"
                    witnesses[key]=x
                    row["witness_key"]=key
            trial.append(row)
    np.savez_compressed(HERE/"contraction_witnesses.npz",**witnesses)
    out=dict(scope="Exploratory, after compact control result. Feasible arbitrary bounded transfer, no material assertion.",
        rows=trial,strict_count=sum(a["strict_energy_only_fit"] for a in trial),
        strict_fit_for_originally_excluded=sum(a["strict_energy_only_fit"] and a["original_decision"]=="PRUNE" for a in trial),
        **elapsed(start))
    (HERE/"contraction_control.json").write_text(json.dumps(out,indent=2)+"\n")
    print(len(trial),out["strict_count"],out["strict_fit_for_originally_excluded"],elapsed(start),flush=True)
if __name__=="__main__":
    with threadpool_limits(1): main()
