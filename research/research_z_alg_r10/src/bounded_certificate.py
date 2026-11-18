"""Validate global projected material intervals for selected noisy R10 data with Arb."""
import json
import time
import numpy as np
from flint import acb, acb_mat, arb, ctx, __version__ as flint_version
from scipy.linalg import null_space
from noise_sets import ROOT, BOX, EPS, GAIN, source_choices, noisy_vector

ctx.prec = 128
DEPTH = 12
CENTER = arb(25)/4
HALF = arb(23)/4
# Exact binary64 inputs for the frozen operators. Noise/gain radii are exact decimals.
E = arb(1)/5000
GAMMA = arb(1)/100


def ball(z):
    z = complex(z)
    return acb(float(z.real), float(z.imag))


def row(a):
    return acb_mat([[ball(z) for z in a]])


def matrix(A):
    return acb_mat([[ball(z) for z in line] for line in A])


def conjdot(x, y):
    return sum((a.conjugate()*b for a,b in zip(x,y)), acb(0))


def dot(x, y):
    return sum((a*b for a,b in zip(x,y)), acb(0))


def square(z):
    return z.real**2+z.imag**2


def norm(x):
    value = sum((square(z) for z in x),arb(0))
    return abs(value).sqrt()


def sqpoly(coeff):
    """Input indexed (q degree,p degree,signal channels), real polynomial in q,p."""
    out = [[arb(0) for _ in range(3)] for _ in range(3)]
    for a,b in np.ndindex(2,2):
        for c,d in np.ndindex(2,2):
            out[a+c][b+d] += conjdot(coeff[a][b],coeff[c][d]).real
    return out


def make_parts(B,H,R,host,observed):
    k = len(R)
    cols = B.shape[1]
    M = matrix(B)
    if cols == k-1:
        n_float = null_space(B.T)[:,0]
        excluded = int(np.argmax(abs(n_float)))
        ix = [i for i in range(k) if i!=excluded]
        inv = matrix(B[ix]).inv()
        mrow = row(B[excluded]) * inv
        n = [acb(1) if j==excluded else -mrow[0,ix.index(j)] for j in range(k)]
    else:
        ix = list(range(k))
        inv = M.inv()
        n = None
    Ah, yd = list(row(host).entries()), list(row(observed).entries())
    rr = list(row(R).entries())
    H = matrix(H)
    X, c, r, nb = [],[],[],[]
    for i in range(k):
        T, C, Z, N = [],[],[],[]
        for degree, value in enumerate((CENTER,HALF)):
            u = [acb(int(j==i and degree==0)) - value*H[j,i] for j in range(k)]
            x=list((inv * acb_mat([[u[j]] for j in ix])).entries())
            prediction=dot(Ah,x)+value*rr[i]
            T.append(x); C.append(prediction); Z.append(dot(yd,x)-prediction)
            if n is not None:
                N.append(dot(n,u))
        X.append(T); c.append(C); r.append(Z); nb.append(N)
    return X,c,r,nb


def two_cell_data(X,c,r,n,i,j):
    rc, vc = [[None,None],[None,None]], [[None,None],[None,None]]
    for a,b in np.ndindex(2,2):
        def diff(u,v):
            return u*n[j][b]-v*n[i][a]
        rc[a][b]=[diff(r[i][a],r[j][b])]
        vc[a][b]=[E*diff(u,v) for u,v in zip(X[i][a],X[j][b])]
        vc[a][b].append(GAMMA*diff(c[i][a],c[j][b]))
    P,Q=sqpoly(rc),sqpoly(vc)
    F=[[P[a][b]-Q[a][b] for b in range(3)] for a in range(3)]
    def endpoint(k,sign):
        return [u+sign*v for u,v in zip(k[1][0],k[1][1])]
    L=norm(endpoint(rc,1)).max(norm(endpoint(rc,-1)))
    L+=norm(endpoint(vc,1)).max(norm(endpoint(vc,-1)))
    return F,Q,L


def value_in_p(f,x):
    return [f[0][j]+x*(f[1][j]+x*f[2][j]) for j in range(3)]


def eliminated(pair,x,h):
    F,Q,L=pair
    c,d,e=value_in_p(F,x)
    at_end=c+e-abs(d)
    if e>0:
        p0=-d/(2*e)
        lower=at_end if p0 < -1 or p0 > 1 else (c-d*d/(4*e)).min(at_end)
    elif e<=0:
        lower=at_end
    else:
        lower=c-abs(d)-abs(e)
    v0,v1,v2=value_in_p(Q,x)
    vmax=abs(v0+v2+abs(v1)).sqrt()
    return lower>h*L*(2*vmax+h*L)


def single_parts(X,c,r,i):
    noise=[[E*x for x in X[i][s]]+[GAMMA*c[i][s]] for s in (0,1)]
    return r[i],noise,abs(r[i][1])+norm(noise[1])


def eliminated_full(single,x,h):
    r,v,L=single
    return abs(r[0]+x*r[1]) > norm([a+x*b for a,b in zip(v[0],v[1])])+h*L


def bins_to_intervals(bins):
    segs=[]
    for t in bins:
        if segs and segs[-1][1]==t:
            segs[-1][1]=t+1
        else:
            segs.append([t,t+1])
    # Width 11.5 and dyadic endpoints are exact as binary64 here.
    return [[BOX[0]+(BOX[1]-BOX[0])*s/2**DEPTH for s in pair] for pair in segs]


def trace(old,full,i):
    X,c,r,n=old
    checks=[two_cell_data(X,c,r,n,i,j) for j in range(len(X)) if j!=i]
    sing=single_parts(*full[:3],i)
    kept=[0]
    after=[]
    counts=[0,0]
    for level in range(1,DEPTH+1):
        h=arb(1)/(2**level)
        next_kept=[]
        for prev in kept:
            for t in (2*prev,2*prev+1):
                x=-1+(2*t+1)*h
                counts[0]+=1
                if any(eliminated(ch,x,h) for ch in checks):
                    continue
                next_kept.append(t)
                if level==DEPTH and not eliminated_full(sing,x,h):
                    after.append(t)
        kept=next_kept
    counts[1]=len(kept)
    return bins_to_intervals(kept),bins_to_intervals(after),counts


def feasible_witness(z,H,R,B,hy,yd):
    # Physical row using T=(diag(1/z)-H)^(-1); verify the entire old/new data jointly.
    k=len(R)
    mat=acb_mat([[((1/ball(z[i])) if i==j else acb(0))-ball(H[i,j]) for j in range(k)] for i in range(k)])
    y=list((row(R)*mat.inv()*matrix(B)+row(hy)).entries())
    residual=[ball(c)-a for c,a in zip(yd,y)]
    g=conjdot(y,residual)/(conjdot(y,y).real+(E/GAMMA)**2)
    score=square(g)/GAMMA**2+sum((square(e-g*a)/E**2 for e,a in zip(residual,y)),arb(0))
    return dict(objective_enclosure=str(score),inside=bool(score<1),outside=bool(score>1))


def make_inputs(k):
    d=np.load(ROOT/"results"/f"lift_k{k}.npz")
    feeds=d["p_old"]; B=d["B_old"]; hy=d["y_host"]
    p,_ = source_choices(d,np.random.default_rng(8871+k))
    # Missing-direction matched beam computed from the old geometry alone.
    a=p["gap"]
    full=np.column_stack((B,d["q"]@a))
    h0=hy@feeds; hf=np.r_[h0,hy@a]
    reads=[]
    for s,delta in enumerate(d["independent_all_source_data"][:3]):
        rng=np.random.default_rng(394760+100*k+s)
        noise=noisy_vector(rng,k)*EPS*.55
        err=.8*GAIN*np.exp(2j*np.pi*rng.uniform())
        ideal=delta+hy
        reads.append((1+err)*np.r_[ideal@feeds,ideal@a]+noise)
    # Independently archived simulated readouts, not regenerated in verification mode.
    np.savez(ROOT/"results"/f"bounded_instance_{k}.npz",H=d["H"],R=d["R"],B=full,
             host=hf,observed=np.asarray(reads),z=np.vstack((d["pair_materials"],d["heldout_materials"][:1])))


def run():
    cpu=time.process_time()
    result=[]
    for k,i in ((12,5),(24,22)):
        f=np.load(ROOT/"results"/f"bounded_instance_{k}.npz")
        H,R,B,hy,reads,Z=[f[key] for key in ("H","R","B","host","observed","z")]
        for s,yd in enumerate(reads):
            old=make_parts(B[:,:-1],H,R,hy[:-1],yd[:-1])
            full=make_parts(B,H,R,hy,yd)
            pre,post,counts=trace(old,full,i)
            witnesses=[]
            for j in ([0,1] if s<2 else [2]):
                witnesses.append(dict(id=j,q=float(Z[j,i]),
                    old=feasible_witness(Z[j],H,R,B[:,:-1],hy[:-1],yd[:-1]),
                    joint=feasible_witness(Z[j],H,R,B,hy,yd)))
            record=dict(k=k,sample=s,coordinate_zero_based=i,old=pre,after_active=post,
                        q_bins_examined=counts[0],witnesses=witnesses)
            result.append(record)
            print(k,s,pre,post,[tuple((a[c]['inside'] for c in ('old','joint'))) for a in witnesses],flush=True)
    out=dict(flint=flint_version,precision_bits=ctx.prec,depth=DEPTH,model="Saved binary64 coefficients regarded as exact. The additive Euclidean vector noise and one complex receiver gain share a unit ellipsoid, radii exactly 1/5000 and 1/100.",
             cpu_seconds=time.process_time()-cpu,records=result)
    (ROOT/"results"/"bounded_certificate.json").write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    import sys
    if '--make-inputs' in sys.argv:
        for k in (12,24):
            make_inputs(k)
    run()

