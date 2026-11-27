"""Two-state bounded-real upper bound for the continuum of passive loads."""
from pathlib import Path
import argparse,json
import numpy as np
from scipy.linalg import eigvalsh,cholesky
from threadpoolctl import threadpool_limits
from time import perf_counter,process_time
from passive_replay import rd,fro,opnorm,up,HERE,lower


def poslower(h,err):
    """Unitary eigendecomposition residual with explicit loss of orthogonality."""
    z,v=np.linalg.eigh((h+h.conj().T)/2)
    ve=up(fro(v.conj().T@v-np.eye(len(z)))+rd(len(z)+3)*(np.sqrt(len(z))+fro(v)**2))
    diff=up(fro(h-(v*z)@v.conj().T)+rd(len(z)+5)*(fro(h)+fro(v)**2*fro(z))+err)
    return lower(min(float(z.min())*(1-ve),float(z.min())*(1+ve))-diff)


def storage_candidate(s1,s0,gamma,margin):
    # Scale the readout once so the lower-right input block is of order unity.
    a=np.diag([s1[24,24],s0[24,24]])
    b=np.vstack([s1[24,:16],s0[24,:16]])
    c=np.column_stack([s1[:16,24],-s0[:16,24]])/gamma
    d=(s1[:16,:16]-s0[:16,:16])/gamma
    state=np.column_stack([a,b]); output=np.column_stack([c,d]);squared=output.conj().T@output
    x=np.zeros((2,2),complex)
    for i in range(400):
        cross=a.conj().T@x@b+c.conj().T@d
        pivot=np.eye(16)-b.conj().T@x@b-d.conj().T@d
        if eigvalsh(pivot)[0]<1e-9:return None
        px=a.conj().T@x@a+c.conj().T@c+cross@np.linalg.solve(pivot,cross.conj().T)
        nx=(px+px.conj().T)/2+margin*np.eye(2)
        if np.linalg.norm(nx-x)<1e-12*(1+np.linalg.norm(nx)):
            x=nx;break
        x=nx
    x=(x+x.conj().T)/2;x[np.diag_indices(2)]=x.diagonal().real
    h=np.eye(18,dtype=complex);h[:2,:2]=x
    h-=state.conj().T@x@state+squared
    # The block-diagonal A and signed C are exact selections of the stored centers;
    # d requires complex subtraction before dividing by the binary64 gamma.
    de=up(rd(6)*(fro(s1[:16,:16])+fro(s0[:16,:16])+fro(s1[:16,24])+fro(s0[:16,24]))/gamma)
    err=up(de*(2*fro(output)+de)+rd(24)*(fro(x)+np.sqrt(16)+fro(state)**2*fro(x)+fro(output)**2))
    eig=poslower(h,err); plo=poslower(x,0.)
    if min(eig,plo)<=0:return None
    return x,dict(storage_min_lower=plo,supply_min_lower=eig,arithmetic_error_upper=err,iterations=i+1,riccati_slack=margin)


def rigorous_storage(s1,s0,gamma):
    for margin in (0.01,0.1,0.001,1.0,1e-5):
        answer=storage_candidate(s1,s0,gamma,margin)
        if answer is not None:
            return answer
    return None


def modeling(s,e):
    p=up(abs(s[24,24])*(1+rd(2))+e); assert p<1
    bg,cg=fro(s[:16,24]),fro(s[24,:16]); gap=lower(1-p)
    return up(e*(1+(bg+cg+e)/gap+bg*cg/gap**2))


def run(file):
    clock=perf_counter(),process_time();q=np.load(HERE/file);s1,s0,e1,e0=(q[k] for k in ('s1','s0','e1','e0'))
    e=up(modeling(s1,float(e1))+modeling(s0,float(e0)))
    band=file.split('_')[-1].replace('.npz','')
    cover='cover_'+file[len('verified_responses_'):].rsplit('_',1)[0]+'.json'
    receipt=json.loads((HERE/cover).read_text());r=receipt['noise_radius']
    low=0.;hi=receipt['outputs'][band]['arbitrary_wave_upper_over_2r']*2*r*1.2
    best=rigorous_storage(s1,s0,hi)
    for j in range(7):
        if best:break
        hi*=1.4;best=rigorous_storage(s1,s0,hi)
    if best:
        for j in range(20):
            g=(low+hi)/2;check=rigorous_storage(s1,s0,g)
            if check:hi,best=g,check
            else:low=g
    result=dict(band_MHz=int(band),inputs=file,feasible=bool(best),model_error_upper=e)
    if best:
        x,log=best;np.save(HERE/('storage_'+file.removesuffix('.npz')+'.npy'),x)
        result.update(upper_over_2r=up((hi+e)/(2*r)),trial_center_gamma=hi,proof=log)
    result.update(wall_s=perf_counter()-clock[0],cpu_s=process_time()-clock[1])
    (HERE/('kyp_'+file.removesuffix('.npz')+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(file,result.get('upper_over_2r'),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('files',nargs='+');a=p.parse_args()
    with threadpool_limits(1):
        for f in a.files:run(Path(f).name)
