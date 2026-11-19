"""A continuum-of-frequency ambiguity witness for the analytic pulse TM model.

The frequency is an interval here. All special functions at reference points,
jets, full-wave solves and powers are evaluated in Arb complex balls.
"""
import argparse
import json
from math import inf, nextafter
from time import process_time, perf_counter

from flint import arb, acb, acb_mat, ctx

from data14 import load_cases, cn_radius, ROOT

ctx.prec=112


def pt(x):
    return arb(float(x))


def up(x):
    return nextafter(float(x.upper()),inf)


def lower(x):
    return nextafter(float(x.lower()),-inf)


def norm(v):
    return sum((abs(z).upper()**2 for z in v.entries()),arb(0)).upper().sqrt().upper()


def rowsum(M):
    return max(sum((abs(M[i,j]).upper() for j in range(M.ncols())),arb(0)).upper()
               for i in range(M.nrows()))


def dist(a,b):
    return sum(((pt(x)-pt(y))**2 for x,y in zip(a,b)),arb(0)).sqrt()


def outward_disk(z,r):
    e=arb(0,r.upper())
    return z+acb(e,e)


def hankels(x0,c,w):
    """Order 0/1 recurrences bound the entire real interval via Gronwall."""
    mid=x0*(1+c)
    a=acb(mid.bessel_j(0),mid.bessel_y(0))
    b=acb(mid.bessel_j(1),mid.bessel_y(1))
    if not w.is_zero():
        dx=abs(x0*w).upper()
        normA=1+1/(mid-dx).lower()
        high=max(abs(a).upper(),abs(b).upper())*(normA*dx).exp()
        a=outward_disk(a,dx*high)
        b=outward_disk(b,dx*high*normA)
    return a,b


def multiply(a,b):
    return (a[0]*b[0],a[1]*b[0]+a[0]*b[1],a[2]*b[0]+2*a[1]*b[1]+a[0]*b[2])


class Physical:
    def __init__(self,case):
        g=case.geometry;p=g['sources'];x=g['centers'];rx=g['receivers'][0]
        k=pt(g['k'])
        self.self_x=k*pt(g['radius'])
        self.green=[[k*dist(a,b) if i!=j else None for j,b in enumerate(x)] for i,a in enumerate(x)]
        self.source=[[k*dist(a,b) for b in p] for a in x]
        self.receiver=[[k*dist(rx,b) for b in x]]
        self.ports=[[k*dist(a,b) if i!=j else None for j,b in enumerate(p)] for i,a in enumerate(p)]
        self.solves=self.factor_attempts=self.rhs=self.kernel_calls=self.bessel_pairs=0

    def kernels(self,c,w):
        self.kernel_calls+=1
        u=1+c+arb(0,w.upper())
        t=self.self_x;xr=t*u
        a,b=hankels(t,c,w)
        scalar=acb(0,arb.pi()/2)
        coeff=(scalar*xr*b.real,scalar*t*xr*a.real,
               scalar*t*t*(a.real-xr*b.real))
        diag=(scalar*xr*b-1,scalar*t*xr*a,scalar*t*t*(a-xr*b))
        sqrt=u.sqrt()
        amp=(acb(0,sqrt/4),acb(0,1/(8*sqrt)),acb(0,-1/(16*u*sqrt)))
        cache={}
        def h0(x):
            # Cache only exact reused geometry entries; use stable object identity.
            key=id(x)
            if key not in cache:
                H,J=hankels(x,c,w)
                cache[key]=(H,-x*J,x*x*(J/(x*u)-H))
            return cache[key]
        def matrices(ds,mode):
            values=[]
            for line in ds:
                row=[]
                for x in line:
                    if x is None:
                        v=diag if mode=="external" else (arb(1),arb(0),arb(0))
                    elif mode=="gram":
                        v=[z.real for z in h0(x)]
                    else:
                        v=multiply(coeff if mode=="external" else amp,h0(x))
                    row.append(v)
                values.append(row)
            return tuple(acb_mat([[v[n] for v in row] for row in values]) for n in range(3))
        G=matrices(self.green,'external');B=matrices(self.source,'source');S=matrices(self.receiver,'external')
        P=matrices(self.ports,'gram')[0]
        self.bessel_pairs+=1+len(cache)
        return G,B,S,P

    def response(self,material,kernels,second):
        G,B,S,_=kernels
        inv_m=acb_mat([[1/pt(m) if i==j else 0 for j in range(len(material))]
                       for i,m in enumerate(material)])
        self.factor_attempts+=1
        inverse=(inv_m-G[0]).inv()
        if not all(v.is_finite() for v in inverse.entries()):
            raise ValueError("forward inverse not enclosed")
        self.solves+=1
        j=inverse*B[0]
        jp=inverse*(B[1]+G[1]*j)
        self.rhs+=2*B[0].ncols()
        if second:
            jpp=inverse*(B[2]+G[2]*j+2*G[1]*jp)
            self.rhs+=B[0].ncols()
            return S[2]*j+2*S[1]*jp+S[0]*jpp
        return S[0]*j,S[1]*j+S[0]*jp


def pnorm(row,Pi):
    return abs((row*Pi*row.conjugate().transpose())[0,0]).upper().sqrt().upper()


def interval_limit(phy,c,w,m,g):
    center=phy.kernels(c,arb(0))
    interval=phy.kernels(c,w)
    a,b=[phy.response(x,center,False) for x in m]
    dd=[phy.response(x,interval,True) for x in m]
    y=g[0]*a[0]-g[1]*b[0]
    yp=g[0]*a[1]-g[1]*b[1]
    ypp=g[0]*dd[0]-g[1]*dd[1]
    ends=[y-w*yp,y+w*yp]
    remainder=w*w*norm(ypp)/2
    fixed=max(norm(v) for v in ends)+remainder
    P,Pw=center[3],interval[3]
    if any(not acb_mat([[P[i,j] for j in range(n)] for i in range(n)]).det().real>0
           for n in range(1,P.nrows()+1)):
        raise ValueError("port Gram not verified positive at midpoint")
    inv=P.inv()
    norm_inv=rowsum(inv)
    eta=norm_inv*rowsum(Pw-P)
    if not eta<1: raise ValueError("port Gram variation")
    coherent=arb(3).sqrt()*(max(pnorm(v,inv) for v in ends)+norm_inv.sqrt()*remainder)/(1-eta).sqrt()
    center_coherent=arb(3).sqrt()*(y*inv*y.conjugate().transpose())[0,0].real.sqrt()
    return fixed,coherent,center_coherent


def prove(case,wit,b=.02,coherent=True,max_depth=13,ratio=1.):
    wall,cpu=perf_counter(),process_time()
    phy=Physical(case);m=wit['materials']
    g=[1+acb(pt(x),pt(y)) for x,y in wit['gains_zeta_ri']]
    z=phy.kernels(arb(0),arb(0))
    s=[]
    for x,gg in zip(m,g):
        pre=phy.response(x,z,False)[0]
        obs=acb_mat([[acb(pt(v.real),pt(v.imag)) for v in case.old]])
        s.append(up(abs((gg-1)/pt(case.gain_radius))**2+norm(gg*pre-obs)**2/pt(case.additive_radius)**2))
    distance=sum(((pt(x)-pt(y))**2 for x,y in zip(*m)),arb(0))/len(m[0])
    half=arb(str(b));limit=2*pt(cn_radius(case.sigma,1 if coherent else len(case.old)))
    target=arb(str(ratio))*limit
    pending=[(arb(0),half,0)];kept=[];unknown=[];reasons={};max_mid=0.
    while pending:
        c,w,d=pending.pop();failure="width"
        try:
            fix,coh,mid=interval_limit(phy,c,w,m,g)
            bound=coh if coherent else fix
            max_mid=max(max_mid,lower(mid/limit) if coherent else 0.)
            if bound<target:
                kept.append([float(c.mid()),up(w),up(fix),up(coh)]);continue
            if coherent and mid>target:
                unknown.append([float(c.mid()),up(w),"point exceeds requested gap bound"]);break
        except (ValueError,ZeroDivisionError) as e:
            failure=str(e)
        reasons[failure]=reasons.get(failure,0)+1
        if d>=max_depth:
            unknown.append([float(c.mid()),up(w),failure]);continue
        pending.extend(((c-w/2,w/2,d+1),(c+w/2,w/2,d+1)))
    ok=not (pending or unknown) and max(s)<1
    return dict(case_id=case.case_id,coherent=coherent,half_band=b,materials=m,gains_zeta_ri=wit['gains_zeta_ri'],
        rms_lower=lower(distance.sqrt()),analytic_old_scores_upper=s, threshold_2r=lower(limit),
        certified_continuum=ok,db_coverage="real interval in relative frequency",
        terminal_intervals=kept,unresolved=unknown,pending_after_contradiction=len(pending),
        upper_ratio=(max(v[3 if coherent else 2] for v in kept)/lower(limit) if ok else None),target_ratio=ratio,
        max_point_coherent_ratio_lower=max_mid,retries=reasons,interval_solves=phy.solves,
        successful_complex_rhs=phy.rhs,forward_inverse_attempts=phy.factor_attempts,operator_jet_calls=phy.kernel_calls,
        hankel_order_pairs=phy.bessel_pairs,wall_s=perf_counter()-wall,cpu_s=process_time()-cpu)


def main():
    p=argparse.ArgumentParser();p.add_argument('--case',action='append',required=True)
    p.add_argument('--known',action='store_true');p.add_argument('--fixed',action='store_true')
    p.add_argument('--same-gain',action='store_true')
    p.add_argument('--provided',help='Saved witness table, or materials and gains keyed by case id')
    p.add_argument('--ratio',type=float,default=1.)
    p.add_argument('--band',type=float,default=.02);p.add_argument('--tag',default='')
    a=p.parse_args();r=ROOT/'results'
    candidates=json.loads((r/'band_probe.json').read_text())
    menu=json.loads((r/'pivot_menu_witnesses.json').read_text())['cases']
    cases={x.case_id:x for x in load_cases()}
    folder=r/'continuous_band';folder.mkdir(exist_ok=True)
    for name in a.case:
        if a.provided:
            table=json.loads((r/a.provided).read_text())
            entry=table.get('cases',table)[name]
            wit=entry.get('best',entry)
        elif a.known:
            source=json.loads((r/'pivot_menu_gain_controls.json').read_text())['cases'][name]
            wit=source['modes']['known_unit']['best']
        elif a.same_gain:
            source=json.loads((r/'pivot_menu_gain_slack.json').read_text())['cases'][name]
            wit=source['modes']['common']['best']
        elif a.fixed:
            wit=menu[name]['best']
        else:
            options=[v for v in candidates[name] if v['coherent']<.96 and v['rms']>.2]
            wit=max(options,key=lambda v:v['rms'])
        if wit is None:continue
        row=prove(cases[name],wit,a.band,not a.fixed,ratio=a.ratio)
        (folder/(name+a.tag+'.json')).write_text(json.dumps(row,allow_nan=False)+'\n')
        print(name,row['certified_continuum'],row['rms_lower'],row['upper_ratio'],len(row['terminal_intervals']),
              row['retries'],row['wall_s'],flush=True)

if __name__=='__main__':
    main()
