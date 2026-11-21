"""Private analytic shapes used only by data generation and subsequent scoring."""
import itertools
from scipy.special import spherical_jn,spherical_yn
from physics import *
FAMILIES=['eccentric_shell','split_frame','separated_bodies','crossed_bars_layer']
def shape(x,idx):
    fam=idx//3;v=idx%3;phi=[.09,.43,-.34][v];theta=[-.19,.31,.24][v]
    x=x-np.array([.008*(v-1),.009*(1-v),.006*v])
    a,b,c=x.T;u=a*np.cos(phi)+b*np.sin(phi);q=-a*np.sin(phi)+b*np.cos(phi)
    y=q*np.cos(theta)+c*np.sin(theta);z=-q*np.sin(theta)+c*np.cos(theta)
    if fam==0:
        ext=(u/(.184+.006*v))**2+(y/.148)**2+(z/.157)**2<=1
        inner=((u-.029-.007*v)/.107)**2+((y+.012)/.094)**2+((z-.009)/.093)**2<1
        return (2.05+.57*v)*(ext&~inner)
    if fam==1:
        ext=(abs(u)<.177)&(abs(y)<.151)&(abs(z)<(.058+.003*v));inner=(abs(u)<.088)&(abs(y)<.078)
        slit=(u>.05)&(abs(y-.019*v)<.034)
        return (2.32+.63*v)*(ext&~inner&~slit)
    if fam==2:
        s=((u+.103)/.081)**2+((y+.032)/.093)**2+((z-.01)/.104)**2<=1
        a=((u-.092)/.069)**2+((y-.068)/.077)**2+((z+.05)/.066)**2<=1
        b=((u-.071)/.048)**2+((y+.108)/.055)**2+((z-.072)/.063)**2<=1
        return np.maximum.reduce([(1.8+.22*v)*s,(3.2+.32*v)*a,(2.6+.20*v)*b])
    rod=(abs(u)<.204)&(abs(y)<.038)&(abs(z-.049)<.030)
    rod2=(abs(y)<.179)&(abs(u+.033)<.041)&(abs(z-.037)<.037)
    layer=(abs(u)<.135)&(abs(y)<.143)&(abs(z+.03)<(.017+.002*v))
    return np.maximum.reduce([(2.25+.12*v)*rod,(2.79+.18*v)*rod2,(3.15+.4*v)*layer])
def frac(mod,idx):
    h=mod.h
    f=lambda x: (2*(np.linalg.norm(x,axis=1)<=.12)) if idx==-1 else shape(x,idx)
    return sum(f(mod.pts+h*np.asarray(d)) for d in itertools.product((-.375,-.125,.125,.375),repeat=3))/64
def mie(mu):
    k=K;r=.12;x=k*r;m=np.sqrt(3.);ns=np.arange(1,13)
    def ric(z):
        j=spherical_jn(ns,z);jp=spherical_jn(ns,z,True);y=spherical_yn(ns,z);yp=spherical_yn(ns,z,True)
        return z*j,j+z*jp,z*(j+1j*y),j+1j*y+z*(jp+1j*yp)
    px,dp,h,dh=ric(x);p,dp2,_,_=ric(m*x)
    a=(m*p*dp-px*dp2)/(m*p*dh-h*dp2);b=(p*dp-m*px*dp2)/(p*dh-m*h*dp2)
    out=np.zeros_like(mu,dtype=complex);prev=np.zeros_like(mu);pn=np.ones_like(mu)
    for n,an,bn in zip(ns,a,b):
        tau=n*mu*pn-(n+1)*prev;out+=(2*n+1)/(n*(n+1))*(an*tau+bn*pn)
        prev,pn=pn,((2*n+1)*mu*pn-(n+1)*prev)/n
    return out
def generate(grid_only=False):
    begin=time.perf_counter();cpu=time.process_time();out=ROOT/'data';out.mkdir(exist_ok=True);private=ROOT/'sealed_truth';private.mkdir(exist_ok=True)
    rows=[];validation=[]
    if grid_only:
        old=json.loads((ROOT/'development_ring/data/generation.json').read_text())
        rows=[r for r in old['rows'] if r['geom']=='hemi'];validation=[r for r in old['validation'] if r['geom'] in ('hemi','sphere_plane_wave')]
    for geom in (['plane_grid'] if grid_only else ['plane_grid','hemi']):
        coarse=Model(12,geom);fine=Model(24,geom)
        for idx in range(12):
            start=time.perf_counter();cost=dict(fine.count);m=frac(fine,idx);Y,err=fine.physical(t(m).real);Y=npy(Y)
            # Same whitened random draw on paired array layouts, per-object independent.
            rng=np.random.default_rng(710031+idx);w=rng.normal(size=Y.shape)+1j*rng.normal(size=Y.shape);noise=.1*np.linalg.norm(Y)*w/np.linalg.norm(w)
            delta=.1*np.linalg.norm(Y+noise)/np.sqrt(1.01)
            np.savez_compressed(out/f'{geom}_{idx:02d}.npz',Y=Y+noise,delta=delta)
            if geom=='plane_grid' and not grid_only:np.savez_compressed(private/f'{idx:02d}.npz',chi=m.reshape(12,2,12,2,12,2).mean((1,3,5)).flatten(),chi24=m)
            rows.append(dict(geom=geom,case=idx,delta=delta,solve_residual=err,gen_wall=time.perf_counter()-start,noise_fraction=float(np.linalg.norm(noise)/np.linalg.norm(Y)),cost=dict(fine.count-Counter(cost))))
            if idx==0:
                coarse_data,_=coarse.physical(t(frac(coarse,idx)).real);validation.append(dict(geom=geom,case='eccentric_shell',relative_12_24=float(np.linalg.norm(npy(coarse_data)-Y)/np.linalg.norm(Y))))
            write(out/'generation.json',dict(rows=rows,validation=validation,wall=time.perf_counter()-begin,cpu=time.process_time()-cpu,increments_existing_hemi=grid_only))
            print(geom,idx,rows[-1]['gen_wall'],flush=True)
        if geom=='plane_grid' and not grid_only:
            theta=np.linspace(0,np.pi,181);dirs=np.stack([np.sin(theta),theta*0,np.cos(theta)],1);eth=np.stack([np.cos(theta),theta*0,-np.sin(theta)],1);ref=1j/K*mie(np.cos(theta))[:,None]*eth
            for mod in [coarse,fine]:
                m=t(frac(mod,-1)).real;inc=np.zeros((mod.N,3),complex);inc[:,0]=np.exp(1j*K*mod.pts[:,2]);j,err=mod.solve(m,m.repeat_interleave(3)[:,None]*t(inc.reshape(-1,1)))
                far=K**2/(4*np.pi)*mod.h**3*np.exp(-1j*K*dirs@mod.pts.T)@npy(j).reshape(-1,3)
                far-=dirs*(dirs*far).sum(1)[:,None];validation.append(dict(geom='sphere_plane_wave',n=mod.n,mie_farfield_rel=float(np.linalg.norm(far-ref)/np.linalg.norm(ref)),solve_residual=err))
    write(out/'generation.json',dict(rows=rows,validation=validation,wall=time.perf_counter()-begin,cpu=time.process_time()-cpu,increments_existing_hemi=grid_only))
if __name__=='__main__':
    import sys
    generate('--grid-only' in sys.argv)
