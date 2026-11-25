"""Test the interval-width mechanism against a fitted low peak on the identical patch."""
import numpy as np
from run import (HERE, SITES, WEIGHTS, BACKGROUNDS, clock, elapsed, inputs, ball,
                 ctx, kernel_from_symbol, star, arrays, patch, get_profile, trace_lower,
                 low, high, mid, eye, acb_mat, frobenius_squared, arb, threadpool_limits, save)

def main():
    start=clock();ctx.prec=128
    spec,sensors,sources,_=inputs();s,b=sensors[32],sources[32]
    v=ball(np.column_stack((s.conj().T,b)));gram=star(v)*v
    old=kernel_from_symbol(spec)
    cases=[]
    for j,site in enumerate(SITES):
        pix=patch(site)
        gr,gs,bs=arrays(old,pix,(s,b),gram,v)
        data=np.load(HERE/('host_patch_%02d__pw360_32x32__n02.npz'%(j+1)))
        x=.625*data['reference'].flatten()/data['reference'].max()
        X=acb_mat(len(x),len(x),[float(x[i]) if i==k else 0 for i in range(len(x)) for k in range(len(x))])
        f=ball(s[:,pix])*(eye(len(x))-X*gs).solve(X*bs)
        y=mid(f);rg=np.random.default_rng(202611251171+j)
        e=rg.standard_normal(y.shape)+1j*rg.standard_normal(y.shape)
        radius=.020001*np.linalg.norm(y);y+=.02*np.linalg.norm(y)*e/np.linalg.norm(e)
        name='low_amplitude_site_%02d'%(j+1)
        np.savez_compressed(HERE/(name+'.npz'),Y=y,reference=x,r=radius,site=pix)
        physical=high(frobenius_squared(f-ball(y)).sqrt()/arb(float(radius)))
        rows=[]
        for bg in BACKGROUNDS:
            for m,h,w in [(3,3,1),(4,2,1),(4,2,1.5),(4,.5,6)]:
                for mu in WEIGHTS:
                    t=clock()
                    port,nu,info=get_profile(gr,bs,bg,m,h,mu,float(spec.real.max()),32,32,w)
                    d=trace_lower(port,nu,ball(y))
                    rows.append(dict(b=bg,center=m,h=h,patch_weight=w,mu=mu,nu=nu,ratio=low(arb(d)/arb(float(radius))),cost=elapsed(t)))
        cases.append(dict(id=name,peak=float(x.max()),known_feasible_parent_upper=physical,rows=rows))
        save('same_site.json',dict(scope='two lower amplitudes of existing native patch geometries, paired controls only; 2 percent new noisy readings; whole domain uncertain',cases=cases,total=elapsed(start)))
if __name__=='__main__':
    with threadpool_limits(1):main()
