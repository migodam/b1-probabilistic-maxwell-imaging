"""Amortized full-domain branch tests with a finite weak but uncertain background."""
from pathlib import Path
from time import perf_counter, process_time
import argparse
import json
import sys

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
sys.path.insert(0, str(STUDY / 'physical_revision/array'))
from screen import inputs
sys.path.insert(0, str(STUDY))
from core.certificate import ball, star, safe_float, gershgorin, frobenius_squared

SITES = ((28,30),(58,54))
CENTERS_WIDTHS = ((4,.5),(4,1),(4,2),(3,3))
BACKGROUNDS = (1/16, 1/8)
WEIGHTS = (2**-12,2**-10)
def clock(): return perf_counter(), process_time()
def elapsed(start):
    return dict(zip(('wall_s','cpu_s'),(a-b for a,b in zip(clock(),start))))
def low(x): return safe_float(x,upper=False)
def high(x): return safe_float(x,upper=True)
def mid(a): return np.array([complex(v) for v in a.entries()]).reshape(a.nrows(),a.ncols())
def sub(a,r,c): return acb_mat(len(r),len(c),[a[i,j] for i in r for j in c])
def eye(n): return acb_mat(n,n,[int(i==j) for i in range(n) for j in range(n)])
def block(a,b,c,d):
    return acb_mat(a.nrows()+c.nrows(),a.ncols()+b.ncols(),[
        v for mat1,mat2 in ((a,b),(c,d)) for i in range(mat1.nrows())
        for mat in (mat1,mat2) for v in [mat[i,j] for j in range(mat.ncols())]])
def save(name, value):
    out=HERE/name
    out.with_suffix('.tmp').write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    out.with_suffix('.tmp').replace(out)

def kernel_from_symbol(spec):
    rows=[acb.dft([acb(complex(z)) for z in row],inverse=True) for row in spec]
    return [acb.dft(list(col),inverse=True) for col in zip(*rows)]
def patch(site):
    return np.array([(i*96+j) for i in range(site[0],site[0]+8)
                     for j in range(site[1],site[1]+8)],dtype=int)

def arrays(kern,pix,model,common,base):
    s,b=model
    n,pad=b.shape[0],len(kern)
    h=len(pix)
    allidx=np.arange(n)
    dx=((pix%96)[:,None]-(allidx%96)[None,:])%pad
    dy=((pix//96)[:,None]-(allidx//96)[None,:])%pad
    # The table is indexed [x][y] after two inverse DFTs.
    gr=acb_mat(h,n,[kern[int(x)][int(y)] for x,y in zip(dx.flat,dy.flat)])
    g_on_i=sub(gr,range(h),pix.tolist())
    i_base=sub(base,pix.tolist(),range(base.ncols()))
    g_base=gr*base
    ts=block(i_base,acb_mat(h,0),g_base,acb_mat(h,0))
    tt=block(eye(h),star(g_on_i),g_on_i,gr*star(gr))
    return block(common,star(ts),ts,tt),g_on_i,sub(ball(b[pix]),range(h),range(b.shape[1]))

def get_profile(gram,incident,background,center,width,mu,max_sym,ports,sources,patch_weight=1):
    p,s,h=ports,sources,incident.nrows()
    zero=arb(0)
    g=1-arb(background)*arb(max_sym)
    if not g>0: raise ValueError('Fourier floor nonpositive')
    m=gram.nrows();ri=range(p);bi=range(p+s);ai=range(p+s,m)
    port=sub(gram,ri,range(m))
    port_gram=sub(gram,ri,ri)+arb(mu)*g*eye(p)
    full=(gram-star(port)*port_gram.solve(port))/g
    d=[arb(1)/arb(mu)]*p + [arb(background)/2]*s
    H=acb_mat(p+s,p+s,[(arb(1)/arb(mu) if i==j and i<p else zero)-d[i]*full[i,j]*d[j] for i in bi for j in bi])
    weight=arb(patch_weight)
    x=weight*arb(center)-arb(background)/2
    ah=weight*(arb(center)**2-arb(width)**2)
    aj=weight-1
    determinant=aj*ah-x*x
    if not determinant<0: raise ArithmeticError("indefinite patch update required")
    uinv=block(ah/determinant*eye(h),x/determinant*eye(h),x/determinant*eye(h),aj/determinant*eye(h))
    z=sub(full,ai,ai)+uinv
    cols=acb_mat(2*h,p+s,[full[p+s+i,j]*d[j]+(incident[i-h,j-p] if i>=h and j>=p else zero) for i in range(2*h) for j in bi])
    zn=mid(z)
    vals,vec=eigh((zn+zn.conj().T)*.5,check_finite=False)
    bv=ball(vec)
    if not gershgorin(star(bv)*bv,upper=False)>0: raise ArithmeticError('trial frame')
    transformed=star(bv)*z*bv
    neg=pos=0
    mingap=None
    for i in range(2*h):
        rad=sum((abs(transformed[i,j]) for j in range(2*h) if i!=j),arb(0))
        d0=transformed[i,i].real
        if (d0-rad)>0: pos+=1; gap=(d0-rad).lower()
        elif (d0+rad)<0: neg+=1; gap=(-d0-rad).lower()
        else: raise ArithmeticError('Unresolved Hermitian inertia')
        mingap=gap if mingap is None else min(mingap,gap)
    out=arb(mu)*(H+star(cols)*z.solve(cols))
    return out,pos-h,dict(b=background,mu=mu,center=center,halfwidth=width,local_cells=h,
                           background_floor_lower=low(g),small_dimensions=[2*h,p+s],
                           patch_weight=float(patch_weight),core_npos=pos,core_nneg=neg,congruence_gap_lower=low(mingap),
                           inertia_one_sided_state=pos-h)

def trace_lower(form,nu,readings):
    p=readings.nrows();s=readings.ncols()
    rhs=acb_mat(p+s,s,readings.entries()+eye(s).entries())
    gram=star(rhs)*form*rhs
    approx=mid(gram); lam,v=eigh((approx+approx.conj().T)*.5,check_finite=False)
    vb=ball(v); tag=star(vb)*gram*vb; metric=star(vb)*vb
    # Congruence plus a scale-relative trace test; verified eigenvalues via
    # approximate diagonalization of H itself, including orthogonality.
    orth=frobenius_squared(metric-eye(s)).sqrt()
    if not orth<1: raise ArithmeticError('eigenvector conditioning')
    diagonal=acb_mat(s,s,[float(lam[i]) if i==j else 0 for i in range(s) for j in range(s)])
    defect=frobenius_squared(gram*vb-vb*diagonal).sqrt()
    # Replace V by its unitary polar factor. The polar error is bounded in
    # operator norm by orth/(1+sqrt(1-orth)). This gives an ordered Weyl bound.
    vdist=orth/(1+(1-orth).sqrt())
    herr=frobenius_squared(gram).sqrt()
    err=defect+vdist*(herr+max(abs(arb(float(x))) for x in lam))
    squares=sum((max(arb(0),(arb(float(x))-err).lower()) for x in lam[:max(0,s-nu)]),arb(0))
    return low(squares.sqrt())

def main(args):
    ctx.prec=128
    tick=clock()
    spec, sensors, sources, data=inputs()
    S,B=sensors[32],sources[32]
    clean=[(row['id'],Y,row['padded_radius']) for row,Y in data[32][32] if row['noise_pct']==2]
    kern=kernel_from_symbol(spec)
    base=ball(np.column_stack([S.conj().T,B]))
    common=star(base)*base
    prep=elapsed(tick)
    payload=dict(physics='native binary64 inputs, exact 192² padded Fourier convolution, all 9216 cells variable',
                 base_time=prep,runs=[],controls=[],sites=[],errors=[],
                 precision=ctx.prec,source_cases=[c[0] for c in clean])
    master=[];checks={}
    for si,site in enumerate(SITES):
        start=clock()
        loc=patch(site)
        gram,gs,bs=arrays(kern,loc,(S,B),common,base)
        spent=elapsed(start)
        # Only acquisition generator sees these known values. The branch engine
        # receives their readings and radii only. The fitted branch is checked offline.
        i,j=np.indices((8,8))
        vals=(3.6+.35*np.cos((i+1)*.69)*np.sin((j+1)*.52)) if si==0 else (4.15+.32*np.sin(i*.61+j*.37))
        cv=ball(vals.reshape(-1,1).astype(complex)); gm=acb_mat(len(loc),len(loc),[cv[k,0]*gs[k,l] for k in range(len(loc)) for l in range(len(loc))])
        bm=acb_mat(len(loc),bs.ncols(),[cv[k,0]*bs[k,l] for k in range(len(loc)) for l in range(bs.ncols())])
        f=ball(S[:,loc])*(eye(len(loc))-gm).solve(bm)
        y=mid(f)
        noise=np.random.default_rng(2026112504+si)
        e=noise.normal(size=y.shape)+1j*noise.normal(size=y.shape)
        r=.020001*np.linalg.norm(y)
        y+=.02*np.linalg.norm(y)*e/np.linalg.norm(e)
        name=f'host_patch_{si+1:02d}__pw360_32x32__n02'
        rho=high(frobenius_squared(ball(y)-f).sqrt())/r
        assert rho<1
        payload['controls'].append(dict(id=name,case=si+1,reference_peak=float(vals.max()),radius=r,
                                        forward_ratio_upper=float(np.nextafter(rho,np.inf)),
                                        setup_and_validation=elapsed(start)))
        np.savez_compressed(HERE/(name+'.npz'),Y=y,r=r,active=loc,reference=vals)
        master.append((gram,bs));clean.append((name,y,r))
        payload['sites'].append(dict(site=si+1,geometry_time=spent))
    tick_all=clock()
    for si,(gram,bs) in enumerate(master):
        for b in BACKGROUNDS:
            for mu in WEIGHTS:
                for m,h in CENTERS_WIDTHS:
                    start=clock()
                    try:
                        form,nu,facts=get_profile(gram,bs,b,m,h,mu,float(np.max(spec.real)),S.shape[0],B.shape[1])
                        facts.update(site=si+1,setup=elapsed(start))
                        rows=[]
                        for name,Y,rad in clean:
                            start=clock()
                            bound=trace_lower(form,nu,ball(Y))
                            rows.append(dict(id=name,distance_lower=bound,radius=rad,
                                             ratio_lower=low(arb(bound)/arb(float(rad))),online=elapsed(start)))
                        facts['readings']=rows
                        payload['runs'].append(facts)
                        print(si+1,b,mu,m,h,nu,[round(a['ratio_lower'],2) for a in rows],flush=True)
                    except (ArithmeticError,ZeroDivisionError,ValueError) as exc:
                        payload['errors'].append(dict(site=si+1,b=b,mu=mu,center=m,halfwidth=h,error=str(exc),cost=elapsed(start)))
                    payload['profiles_and_readouts']=elapsed(tick_all)
                    save('results.json',payload)
if __name__=='__main__':
    with threadpool_limits(1): main(argparse.ArgumentParser().parse_args())
