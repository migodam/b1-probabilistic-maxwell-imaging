"""Energy-only floor and one-sided whole-field small-contrast control."""
import json

import numpy as np
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits

from sector import HERE,PUBLIC,Green
from verify import arb,ball,fhi,lower,norm
from flint import ctx


def data_energy(s,y,radius):
    left,sigma,_=np.linalg.svd(s,full_matrices=False)
    yd=left.conj().T@y
    unseen=np.linalg.norm(y-left@yd)**2
    r2=radius**2
    if unseen>=r2:
        return None
    if np.linalg.norm(y)<=radius:
        return 0.
    l2=sigma**2
    def tail(z):
        return np.linalg.norm((z/(l2+z))[:,None]*yd)**2+unseen-r2
    top=1.
    while tail(top)<0:
        top*=2
    lam=brentq(tail,0,top,xtol=1e-18)
    return float(np.linalg.norm((sigma/(l2+lam))[:,None]*yd)**2)


def main():
    ctx.prec=112
    out=[]
    for filename in sorted((HERE/"results").glob("n*_main.json")):
        obs=json.loads(filename.read_text())
        case=obs["case"]
        side=obs["side"]
        arrays=np.load(filename.with_suffix(".npz"))
        data=np.load(PUBLIC/"cases"/f"{case}.npz")
        loc=PUBLIC/"operators"/f"{obs['operator']}.npz"
        raw=np.load(loc)
        g=Green.read(loc,side)
        s,b,y,r=raw["S"],raw["E0"],data["Y"],float(data["radius"])
        sf,bf,yf=norm(ball(s)),norm(ball(b)),norm(ball(y))
        spectrum=max(abs(z).upper() for z in ball(g.spectrum).entries())
        gnorm=spectrum+(norm(ball(g.left))*norm(ball(g.right))).upper()
        # Use a distinct lower bound for ||Y|| in the reverse triangle.
        yb=ball(y)
        ylower=sum((z.real**2+z.imag**2 for z in yb.entries()),arb(0)).sqrt().lower()
        a=(ylower-arb(r)-arb("1e-8")*yf).lower()
        elementary=(a/(sf*bf+a*gnorm)).lower() if a>0 else arb(0)
        row=dict(case=case,roi={},spectral_green_norm_upper=fhi(gnorm),
                 strict_global_peak_floor_elementary_neumann=lower(elementary/arb(2)))
        for key, p in obs["sections"].get("roi",{}).items():
            energy=p.get("zero_roi_test",{}).get("energy_condition")
            if not energy:
                continue
            e=energy["best"]
            active=np.flatnonzero(~arrays[key+"_mask"].ravel())
            amin=data_energy(s[:,active],y,r)
            row["roi"][key]=dict(energy_screen_critical=e["energy_required_upper"],
                         energy_screen_critical_over_incident=e["energy_required_over_incident"],
                         energy_only_minimum_given_data_numeric=amin,
                         energy_only_minimum_over_incident=None if amin is None
                             else amin/energy["reference_incident_fro2"],
                         shifted_negative_index=e["nu"],initial_negative_index=energy["base_negative"],
                         spectral_shift=e["shift"],scope="positive only under an external certified total-current limit")
        out.append(row)
    (HERE/"results"/"energy_controls.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    with threadpool_limits(1):
        main()
