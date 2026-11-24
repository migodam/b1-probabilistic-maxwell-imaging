"""A rejected hard support also implies a positive lower bound on outside contrast."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits
from flint import arb,ctx
from branches import expand,primary_menu
from development import inputs,SLIT,SHELL
from native_box import HERE,OPTICS,ball,herm,star,lo,hi,gershgorin,frobenius_squared,clock,elapsed

def norm(z):
    return frobenius_squared(z).sqrt()

def outside(geometry,branch,y,r,cert):
    spec,s,b=geometry
    mask=expand(branch,96)>0
    _,v=eigh(herm(s[:,mask]@s[:,mask].conj().T),check_finite=False)
    _,w=eigh(herm(b[mask].conj().T@b[mask]),check_finite=False)
    v=ball(v[:,:cert["left_modes"]]); w=ball(w[:,:cert["right_modes"]])
    vn=gershgorin(star(v)*v,upper=True).sqrt().upper()
    wn=gershgorin(star(w)*w,upper=True).sqrt().upper()
    c=arb(cert["current_norm_factor_upper"])
    ylow=norm(star(v)*ball(y)*w).lower()
    a_in=norm(ball(b[mask])*w).upper()
    a_out=norm(ball(b[~mask])*w).upper()
    s_in=norm(star(v)*ball(s[:,mask])).upper()
    s_out=norm(star(v)*ball(s[:,~mask])).upper()
    # Compression of a unitary-diagonalized circulant convolution.
    gmax=max((arb(float(x.real))**2+arb(float(x.imag))**2).upper() for x in spec.ravel()).sqrt().upper()
    excess=(ylow-vn*wn*arb(r)-s_in*c*a_in).lower()
    if excess<=0:
        return dict(contrast_lower=0.,current_lower=0.,status="unresolved")
    ell=(excess/(s_out+s_in*c*gmax)).lower()
    aa=(a_out+gmax*c*a_in).upper()
    bb=(gmax+gmax*gmax*c).upper()
    lower=(ell/(aa+bb*ell)).lower()
    return dict(contrast_lower=lo(lower),current_lower=lo(ell),projected_signal_excess_lower=lo(excess),
                G_norm_upper=hi(gmax),outside_total_field_intercept_upper=hi(aa),
                outside_total_field_slope_upper=hi(bb),left_norm_upper=hi(vn),right_norm_upper=hi(wn),
                Q_bounds_status="one native coercivity certificate plus outward norm arithmetic",
                status="positive" if lower>0 else "unresolved")

def main():
    ctx.prec=128;ctx.threads=1
    start=clock()
    base=json.loads((HERE/"fast_development.json").read_text())["cases"]
    with np.load(OPTICS/"synthetic/public/fine_geometry.npz") as a:
        geom=a["spectrum"],a["PW180_10x20__S"],a["PW180_10x20__B"]
    result=[]
    for name in (SLIT,SHELL):
        Y,r,_,_=inputs(name)
        with np.load(HERE.parent/"diagnosis/trajectories"/name/"step_0064.npz") as a:
            menu=primary_menu(a["m"],a["J"])["branches"][1:3]
        for bran,cert in zip(menu,base[name],strict=True):
            assert cert["branch"]==bran.name
            t=clock()
            row=outside(geom,bran,Y,r,cert)
            result.append(dict(observation=name,branch=bran.name,noise_radius=r,**row,**elapsed(t)))
    (HERE/"outside_contrast.json").write_text(json.dumps(dict(rows=result,**elapsed(start)),indent=2)+"\n")
    for r in result:
        print(r["observation"],r["branch"],r["current_lower"],r["contrast_lower"],flush=True)
if __name__=="__main__":
    with threadpool_limits(1):main()
