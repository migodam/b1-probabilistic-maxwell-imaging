"""Finite normal-memory escape in the exact shared-material projection geometry.
This is an analytic two-channel constitutive toy, not a Maxwell recovery test.
"""
import json, math
from pathlib import Path
import numpy as np

LO, HI = 0.1, 3.0
thlo, thhi=np.arctan(LO),np.arctan(HI)

def pc(z):
    B=z@z.T
    vals,vecs=np.linalg.eigh(B)
    choices=[thlo,thhi]
    # Rank-one approximation directions, modulo sign.
    for u in vecs.T:
        if u[1]<0:u=-u
        th=np.arctan2(u[0],u[1])
        if thlo <= th <= thhi: choices.append(th)
    objs=[]
    for th in choices:
        u=np.array([np.sin(th),np.cos(th)])
        p=np.outer(u,u@z)
        objs.append(float(np.linalg.norm(z-p)**2))
    i=int(np.argmin(objs));th=choices[i]
    u=np.array([np.sin(th),np.cos(th)])
    return np.outer(u,u@z),float(np.tan(th))

def run(kind):
    theta=np.arctan(1.0) if kind=='interior' else thlo
    u=np.array([np.sin(theta),np.cos(theta)])
    up=np.array([np.cos(theta),-np.sin(theta)])
    v=np.array([3.2,0.0])
    w=np.array([0.0,1.0]) if kind=='interior' else np.array([-.1875,np.sqrt(1-.1875**2)])
    C=np.outer(u,v);D=np.outer(up,w);A=C+D
    a=float(v@v);b=float(w@w);g=float(v@w)
    if kind=='interior':tc=np.sqrt(a/b)
    else:
        ct=1/np.tan(thhi-thlo)
        tc=(-g*ct+np.sqrt((g*ct)**2+a*b))/b
    def pa(z):return z+(float(np.vdot(D,A-z).real)/float(np.vdot(D,D).real))*D
    z=A.copy();rows=[]
    for k in range(500):
        c,m=pc(z);p=pa(2*c-z);zp=z+p-c
        rows.append(dict(k=k,material=m,shadow_gap=float(np.linalg.norm(p-c)),
                         ray_distance=float(np.linalg.norm(z-(C+(k+1)*D))),
                         C_shadow_change=float(np.linalg.norm(c-C))))
        z=zp
    first=next(r['k'] for r in rows if r['C_shadow_change']>1e-8)
    ap=A.copy()
    for k in range(500):ap=pa(pc(ap)[0])
    # Skip exactly the unchanged-shadow prefix: same start as ordinary DR.
    kskip=int(np.floor(tc)) # t from 1 through floor(tc) remain before event.
    zfast=C+(1+kskip)*D
    for k in range(500-kskip):
        c,m=pc(zfast);zfast=zfast+pa(2*c-zfast)-c
    return dict(kind=kind,theta=theta,a=a,b=b,c=g,threshold=tc,
                predicted_first_switch_k=int(np.floor(tc)),observed_first_switch_k=first,
                initial_AP_gap=float(np.linalg.norm(pc(A)[0]-A)),
                final_AP_gap=float(np.linalg.norm(pc(ap)[0]-ap)),
                final_DR_shadow_gap=rows[-1]['shadow_gap'],
                final_DR_material=rows[-1]['material'],
                skipped_projection_pairs=kskip,fast_forward_state_difference=float(np.linalg.norm(z-zfast)),
                trace_first_ten=rows[:10])

out=dict(scope='Analytic constitutive rank-one toy, not Maxwell material recovery; all fields real with two shared illuminations.',cases=[run('interior'),run('lower_bound')])
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
