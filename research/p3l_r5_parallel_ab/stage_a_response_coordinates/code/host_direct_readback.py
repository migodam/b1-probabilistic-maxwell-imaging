"""Analysis-only direct GL readback of already-paid transition nodes."""
from run_v3_strata import *
def main():
    cache={}
    for path in (OUT/'physics/A1_strata_v3/endpoints').glob('*.npz'):
        with np.load(path) as z:
            x=z['x'];inactive=x.copy();inactive[[0,64]]=0
            if np.any(inactive):continue
            a=x[[0,64]];cache[tuple(np.round(a,14))]=dict(z=a,f=np.concatenate((z['prediction'].real.T,z['prediction'].imag.T),axis=1).ravel()*np.sqrt(2),chi=z['chi'])
    host=min(cache.values(),key=lambda r:np.linalg.norm(r['z']-np.array([HOSTMEAN,0.])))
    out=[]
    for ni in range(3):
        name=f'S1_P0_N{ni}'
        with np.load(OUT/f'v2/datasets/{name}.npz') as z:d=z['d'];sig=float(z['sigma'])
        summaries={}
        for n in [24,48]:
            x,w=leggauss(n);nodes=[host];lw=[np.log(PMASS)-.5*np.sum((host['f']-d)**2)/sig**2]
            for a,wa in zip(CUT*x,CUT*w):
                for b,wb in zip(4*x,4*w):
                    r=cache[tuple(np.round([a,b],14))];nodes.append(r);lw.append(np.log(wa*wb)+p_log(r['z'],r['f'],d,sig))
            lz=logsumexp(lw);weights=np.exp(np.array(lw)-lz);ss=summarize(nodes,weights,d,sig,[],lz);ss['host_mass']=float(weights[0]);summaries[str(n)]=ss;report(S/f'{name}_DIRECT_GL_{n}.json',ss)
        out.append(dict(observation=name,convergence=difference(summaries['48'],summaries['24']),host_mass=summaries['48']['host_mass'],physical_calls=0,pure_stratum='omitted here for analysis only; its saved V3 contribution numerically underflowed'))
    report(S/'HOST_DIRECT_READBACK.json',out);print(json.dumps(plain(out)))
if __name__=='__main__':main()
