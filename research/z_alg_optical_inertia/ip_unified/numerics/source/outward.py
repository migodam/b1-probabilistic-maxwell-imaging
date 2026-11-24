"""An independent outward lower endpoint from the full-box local defect.
No negative Ritz vector or positive eigenbasis is used in this calculation.
"""
import argparse,json
from pathlib import Path
from time import perf_counter,process_time
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits
from flint import arb,acb,acb_mat,ctx
from profile import HERE,configure,FourierEnvelope,herm,save_json
from native import ball,star,low,high,frobenius_squared

def norm(a):
    return frobenius_squared(a).sqrt().upper()

def eye(n):
    return acb_mat(n,n,[int(i==j) for i in range(n) for j in range(n)])

class FrameFloor(FourierEnvelope):
    def select(self,w,q):
        for cut in [4,6,8,10,12]:
            rec=self.numerical_matrix(w,cut,1.,0.,self.floor)
            if rec is None:continue
            mat,k,modes,iso,level=rec
            need=max(0,k-q)
            if not need or eigh(herm(mat),eigvals_only=True,subset_by_index=(need-1,need-1),check_finite=False)[0]<-1.e-8:
                return level,1.,0.,self.floor,modes,iso,need
        raise ArithmeticError('independent inertia proposal unavailable')

def frame_floor(geo, w, omega, eta,q):
    g,gb=geo.spec,geo.spec_ball
    u=arb(w.cap); om=arb(omega)
    # The disk contains every row multiplier x_i between 0 and u. Jensen's
    # inequality for the compressed Fourier operator is used before this bound.
    def shifted(z):
        b=2*abs(1-u*z/2)**2+u*u*abs(z)**2/2
        return acb(z.real+om*b/u,z.imag)
    geo.spec=g+omega/w.cap*(2*abs(1-w.cap*g/2)**2+(w.cap**2/2)*abs(g)**2)
    geo.spec_ball=[[shifted(z) for z in row] for row in gb]
    try:
        env=FrameFloor(geo);env.floor=eta
        return env.verify(w,q)
    finally:geo.spec,geo.spec_ball=g,gb

def lower_columns(nat,jtrial,hfloat,y,eta,omega,q,r):
    geo=nat.g; mu=arb(nat.w.scale);u=arb(nat.w.cap); om=arb(omega)*mu
    p,s=y.shape
    h,v=eigh(herm(hfloat));h,v=h[::-1],v[:,::-1]
    vv=ball(v); jb=ball(jtrial); yy=ball(y)
    rhs=star(geo.sensor_ball)*yy/mu + acb(arb(nat.w.cap)/2,arb(nat.w.optical[0])/(2*mu))*geo.incident_ball
    kj=nat.exact_action(jb)
    tt=kj-rhs
    hh=star(yy)*yy+mu*(star(jb)*(kj-rhs)-star(rhs)*jb)
    # H_eta for the binary trial with its full nonstationarity, scaled K/scale.
    hv=star(vv)*hh*vv
    gram=star(vv)*vv-eye(s)
    for i in range(s):
        hv[i,i]-=arb(float(h[i]))
    norm_slack=norm(hv)+(1+max(abs(float(x)) for x in h))*norm(gram)
    stat=(mu/arb(eta))*norm(tt*vv)**2
    err=(norm_slack+stat).upper()
    jj=jb*vv
    ee=(geo.incident_ball+geo.exact_green(jb))*vv
    mis=[arb(0) for _ in range(s)]
    end=[arb(0) for _ in range(s)]
    for i in range(jj.nrows()):
        for k in range(s):
            j,e=jj[i,k],ee[i,k]
            cj,ce=complex(j),complex(e)
            x=0 if abs(ce)==0 else min(float(u),max(0.,(ce.conjugate()*cj).real/abs(ce)**2))
            # For any (rounded) complex test b, |j-xe|^2 >=
            # min_{v in {0,u}} 2 Re b*(j-ve) - |b|^2.
            b=acb(cj-x*ce);sq=abs(b)**2
            left=2*(b.conjugate()*j).real-sq
            right=left-2*u*(b.conjugate()*e).real
            mis[k]+=max(arb(0),min(left.lower(),right.lower()))
            end[k]+=max((abs(j)**2).upper(),(abs(j-u*e)**2).upper())
    # Evaluate proposed eigenvalue lower endpoints directly. Increasing j is
    # not required for correctness; each uses its own j-dimensional trial.
    out=[]
    for j in range(q,s):
        hs=[arb(float(h[k]))-err for k in range(j+1)]
        lo=float(h[j]-float(err))-1.e-12
        hi=min(float(h[j-1]) if j else 1.e100,float(h[j]+float(om*mis[j])))
        for _ in range(32):
            t=(hi+lo)/2
            if t<=lo:break
            tt=arb(t)
            if any(x<=tt for x in hs[:-1]):
                hi=t;continue
            rho=sum((om*end[k]/(hs[k]-tt) for k in range(j)),arb(0))
            if hs[j]+om*mis[j]/(1+rho)>tt:lo=t
            else:hi=t
        # Check rather than rely on bisection initialization.
        t=arb(lo)
        if any(x<=t for x in hs[:-1]):
            out.append(0.);continue
        rho=sum((om*end[k]/(hs[k]-t) for k in range(j)),arb(0))
        if not (hs[j]+om*mis[j]/(1+rho)>t and abs(t)<=max(abs(z) for z in h)+1):
            out.append(0.);continue
        out.append(max(0.,lo))
    dd=sum((arb(x) for x in out),arb(0)).sqrt()
    return dict(distance_lower=max(0.,low(dd)),ratio_lower=max(0.,low(dd/arb(r))),
                residual_eigenvalue_lowers=out,source_basis_and_center_form_error=high(norm_slack),
                stationarity_deduction=high(stat),full_template_defect_lower=[max(0.,low(z)) for z in mis],
                endpoint_energy_upper=[high(z) for z in end])

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--case',default='u100_best');ap.add_argument('--weight',type=float,default=.001)
    ap.add_argument('--count',type=int,default=13);ap.add_argument('--objects',default='123456');a=ap.parse_args()
    ctx.prec=82;ctx.threads=2
    ts,cs=perf_counter(),process_time()
    geo,nat,data,_=configure(a.case)
    c=HERE.parent/'native';sav=np.load(c/(a.case+'.npz'))
    jb,hs=sav['center'],sav['hs'];om,eta=a.weight,2.**-20
    proof=frame_floor(geo,nat.w,om,eta,a.count)
    print('state done',proof['upper_count'],flush=True)
    rows=[]; path=HERE/(a.case+'_strict_source.json')
    for n,(meta,y) in enumerate(data):
        if str(n+1) not in a.objects:continue
        t,cpu=perf_counter(),process_time()
        r=meta['padded_radius']
        row=lower_columns(nat,jb[n],hs[n],y,eta,om,a.count,r)
        rows.append(dict(id=meta['id'],r=r,**row,wall_seconds=perf_counter()-t,cpu_seconds=process_time()-cpu))
        save_json(path,dict(case=a.case,scaled_state_spend=om,eta=eta,state_inertia_upper=a.count,
            fourier_proof=proof,precision_bits=ctx.prec,rows=rows,wall_seconds=perf_counter()-ts,cpu_seconds=process_time()-cs))
        print(meta['id'],row['ratio_lower'],flush=True)
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
