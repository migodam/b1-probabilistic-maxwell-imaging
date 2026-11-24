"""Validate the real directional derivatives and the unchanged Dual-CIE branch."""

from pathlib import Path
import sys
import json

import numpy as np
from threadpoolctl import threadpool_limits

from cross_engine import FrozenLaw, PixelLaw, blocks, iterate, mature, norm_squared

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/"diagnosis"))
from diagnose import CASES, Slotted, get_case, read_npz, write


def model_for(name):
    case,pub = get_case(name)
    model = mature.legacy.PaperModel(pub,case,mature.legacy.Clock(float("inf")))
    return model,case,pub


def derivatives():
    model,case,pub = model_for(Slotted)
    sample = read_npz(HERE.parent/"diagnosis/trajectories"/Slotted/"step_0064.npz")
    J,E = sample["J"],sample["E"]
    rng = np.random.default_rng(20261124)
    dJ = rng.normal(size=J.shape)+1j*rng.normal(size=J.shape)
    dJ *= np.linalg.norm(J)/np.linalg.norm(dJ)
    dE = model.green.apply(dJ)
    cap = np.linspace(.25,20.,model.n)
    kinds = {
        "diag":FrozenLaw(J,E,model.E,cap,blocks(model.shape,1),diagonal=True),
        "joint":FrozenLaw(J,E,model.E,cap,blocks(model.shape,1)),
        "regional4_passive_sector":FrozenLaw(J,E,model.E,cap,blocks(model.shape,2),passive=True,sector_weight=.31),
        "regional16_sector":FrozenLaw(J,E,model.E,cap,blocks(model.shape,4),sector_weight=.31),
        "pixel":PixelLaw(model.E,E,False),
        "pixel_passive":PixelLaw(model.E,E,True),
    }
    out={}
    for name,law in kinds.items():
        direct,adj_rhs,opt,sec = law.gradient_parts(J,E)
        g = direct+model.green.adjoint(adj_rhs)
        slope = float(np.vdot(g,dJ).real)
        weight=getattr(law,"sector_weight",0.)
        terms=law.line(J,E,dJ,dE)
        eps=1.e-5
        def value(opt_sec):
            return opt_sec[0]+weight*opt_sec[1]
        def direct_value(t):
            return value(law.gradient_parts(J+t*dJ,E+t*dE)[2:])
        finite=(direct_value(eps)-direct_value(-eps))/(2*eps)
        line_error=max(abs(value(law.line_value(terms,t))-direct_value(t)) for t in (-.04,0.,.025))
        rel=abs(slope-finite)/max(abs(slope),abs(finite),1e-12)
        if rel>3e-6 or line_error>1.e-10*max(1.,opt+weight*sec):
            raise ArithmeticError(f"Law derivative mismatch: {name} {rel}, {line_error}.")
        out[name]=dict(analytic=slope,finite_difference=finite,relative_error=rel,
                       path_absolute_error=line_error)
    return out


def baseline(name,n=32):
    a,case,pub = model_for(name)
    b,_,_ = model_for(name)
    m0=np.zeros(a.n)
    old=mature.run_cie_arm(a,pub["Y"],case,"dual",m0.copy(),float(pub["delta"]))
    new=iterate(b,pub["Y"],case,"dual",m0.copy(),20.)
    errors={"m":0.,"J":0.,"E":0.}
    try:
        for _ in range(n):
            m,j,stage,event=next(old)
            m2,j2,s2,e2=next(new)
            if (j,stage,event)!=(j2,s2,e2):
                raise ArithmeticError(f"Baseline stages differ up to step {j}.")
            terms = {"m":(m,m2), "J":(old.gi_frame.f_locals["J"],new.gi_frame.f_locals["J"]),
                     "E":(old.gi_frame.f_locals["Q"],new.gi_frame.f_locals["Q"])}
            for k,(left,right) in terms.items():
                errors[k]=max(errors[k],float(np.max(abs(left-right))))
    finally:
        old.close();new.close()
    if max(errors.values())>1.e-12:
        raise ArithmeticError(f"Cloned Dual-CIE baseline differs: {name} {errors}.")
    return dict(steps=n,maximum_difference=errors,counts_match=dict(a.clock.counts)==dict(b.clock.counts))


def varying_box():
    model,case,pub = model_for(Slotted)
    cap = np.repeat([0.,.25,.75,3.],model.n//4)
    m0 = np.ones(model.n)
    gen = iterate(model,pub["Y"],case,"regional4",m0,cap,strength=.25,cadence=8)
    excess = 0.
    try:
        for count in range(32):
            m,_,_,_=next(gen)
            excess=max(excess,float(np.max(m-cap)),float(-m.min()))
    finally:
        gen.close()
    if excess>1.e-12:
        raise ArithmeticError(f"Spatial box exceeded by {excess}.")
    return dict(steps=32,maximum_violation=excess,
                corrections=model.clock.counts["cross_corrections_accepted"])


def main():
    with threadpool_limits(1):
        out=dict(derivatives=derivatives(),baseline={name:baseline(name) for name in CASES},
                 spatial_box=varying_box())
        write(HERE/"checks.json",out)
        print(json.dumps(out,indent=2))


if __name__=="__main__":
    main()
