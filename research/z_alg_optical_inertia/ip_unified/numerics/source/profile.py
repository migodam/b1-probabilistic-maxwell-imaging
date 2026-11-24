"""Full-box lower bounds from the constitutive defect of a fixed current.
Search stage. All final inequalities are checked independently in outward.py.
"""
from pathlib import Path
import argparse,json,sys
from time import perf_counter,process_time
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'native'))
from replay import configure
from verify_native import FourierEnvelope
from native import ball
from screen import save_json

def herm(a):
    return (a+a.conj().T)*.5

def best_frame(geo,native,omega,eta,gr=1.):
    # The physical trial and G stay fixed. Only a lower comparison symbol changes.
    env=FourierEnvelope(geo)
    g=geo.spec
    u=native.w.cap
    b=(1+gr)*abs(1-u*g/2)**2+(1+1/gr)*(u/2)**2*abs(g)**2
    geo.spec=g+omega*b/u
    try:
        best=None
        for cut in [4,6,8,10]:
            item=env.numerical_matrix(native.w,cut,1.,0.,eta)
            if item is None:continue
            mat,neg,*_=item
            ev=eigh(herm(mat),eigvals_only=True,check_finite=False)
            count=int(neg-sum(ev< -1.e-8))
            if best is None or count<best[0]:
                best=(count,cut)
        return best if best is not None else (geo.n,None)
    finally:geo.spec=g

def spectrum(omega,hs,j,e,nu,u):
    vals,v=eigh(herm(hs));vals,v=vals[::-1],v[:,::-1]
    jj=j@v;ee=e@v
    local_arg=np.clip((ee.conj()*jj).real/np.maximum(abs(ee)**2,1.e-300),0,u)
    miss=np.sum(abs(jj-ee*local_arg)**2,axis=0)
    # Simultaneously bound the scaled norm of all stronger columns. Löwner
    # Young on the midpoint and the bounded, independent real cell variation.
    p=jj-u/2*ee
    pp=herm(p.conj().T@p);qq=herm(ee.conj().T@ee)*(u/2)**2
    ends=abs(jj)**2;ends2=abs(jj-u*ee)**2
    out=np.array(vals)
    for i in range(nu,len(vals)):
        if not i:
            out[i]=vals[i]+omega*miss[i];continue
        lo=vals[i];hi=min(vals[i-1],lo+omega*miss[i])
        for _ in range(33):
            t=(hi+lo)*.5
            if t<=lo:break
            scales=(vals[:i]-t)**-.5
            def norm(z):
                gram=z[:i,:i]*scales[:,None]*scales[None,:]
                return max(0.,eigh(gram,eigvals_only=True,subset_by_index=(i-1,i-1),check_finite=False)[0])**.5
            gram_bound=(norm(pp)+norm(qq))**2
            rho=omega*gram_bound
            if vals[i]+omega*miss[i]/(1+rho)>t:lo=t
            else:hi=t
        out[i]=lo
    return out,vals,v,miss

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--case',default='u100_best');a=ap.parse_args()
    t0,c0=perf_counter(),process_time()
    case=a.case
    geo,native,_,_=configure(case)
    data=HERE.parent/'native'
    raw=np.load(data/(case+'.npz')); js=json.loads((data/(case+'.json')).read_text())
    j,e,h=raw['center'],raw['center_field'],raw['hs']
    scale=js['profile']['scale'];u=native.w.cap
    # The small analytic choice does not need an additional Fourier solve.
    gap=js['verified_state']['positive_gap_lower'] if 'positive_gap_lower' in js['verified_state'] else js['verified_state']['fourier']['positive_gap_lower']
    old_nu=len(raw['negative_eig'])
    simple=float(gap*.95/(1+u*np.abs(geo.spec).max())**2)
    omegas=[simple,.0003,.001,.003,.01,.03,.06,.1,.18,.3]
    results=[];eta=2.**-20
    for w in omegas:
        q,cut=(old_nu,0) if w==simple else best_frame(geo,native,w,eta)
        print(case,w,q,flush=True)
        if q>=32:continue
        row=dict(omega=w,eta=eta,q=q,cut=cut,material_cap=u,objects=[])
        for jr,er,hr,entry in zip(j,e,h,js['rows']):
            r=entry['padded_radius'];lo,base,_,miss=spectrum(w*scale,hr,jr,er,q,u)
            row['objects'].append(dict(id=entry['object_id'],radius=r,base=(np.maximum(base[old_nu:],0).sum())**.5/r,
                        lower_candidate=(np.maximum(lo[q:],0).sum())**.5/r,
                        eigenvalues=lo.tolist(),old_eigenvalues=base.tolist(),exact_line_mismatch=miss.tolist()))
        results.append(row);print([round(o['lower_candidate'],4) for o in row['objects']],flush=True)
    save_json(HERE/(case+'_profile.json'),dict(case=case,procedure='floating proposals only',scales=results,wall=perf_counter()-t0,cpu=process_time()-c0))
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
