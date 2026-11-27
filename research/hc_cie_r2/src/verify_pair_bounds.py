"""Conservative finite floating-model all-pair exclusions on two nontrivial endpoints.

Assumes standard IEEE double complex BLAS error model without overflow/underflow.
Treats archived floating inputs as exact. Physical noise/model error has separate units.
"""
import json
from time import perf_counter,process_time
from interaction import np, LAB, NORM, public_background, historical, prepare, threadpool_limits
from scipy.linalg import inv

UROUND=np.finfo(float).eps/2
INFLATE=1+1e-9


def ga(k):return k*UROUND/(1-k*UROUND)


def upper(a):return INFLATE*NORM(a)


def background_errors(data,bg):
    g,s,b,m=data["G"],data["S"],data["E0"],bg.m
    n=m.size
    if not np.array_equal(g,g.T):
        raise ValueError("requires symmetric reciprocal array")
    a=np.eye(n)-m[:,None]*g
    inva=inv(a)
    ea=ga(8)*(np.sqrt(n)+max(abs(m))*upper(g))*INFLATE
    ai,an=upper(inva),upper(a)
    rho=(upper(np.eye(n)-a@inva)+(ea+ga(8*n+16)*an)*ai+ga(8*n+16)*np.sqrt(n))*INFLATE
    assert rho<1
    bound=INFLATE*ai/(1-rho)
    def solve_error(res,x,y):
        xn=upper(x)
        return INFLATE*bound*(upper(res)+ea*xn+ga(8*n+16)*(xn*an+upper(y)))
    eo=solve_error(bg.o@a-s,bg.o,s)
    eu=solve_error(a.T@bg.u-b,bg.u,b)
    y_alt=(bg.o*m)@b
    ey=INFLATE*(upper(y_alt-bg.y)+(eo+ga(8*n+16)*upper(bg.o))*upper(b)*max(abs(m)))
    return dict(inverse_norm_upper=bound,neumann_residual_upper=rho,error_O=eo,error_U=eu,error_Y0=ey)


def bound_all(data,bg):
    di=background_errors(data,bg)
    r0=data["Y"]-bg.y
    r=r0.ravel();rn=upper(r)
    er=(di["error_Y0"]+ga(2)*(upper(data["Y"])+upper(bg.y)))*INFLATE
    eo,eu=di["error_O"],di["error_U"]
    no=np.linalg.norm(bg.o,axis=0)*INFLATE;nu=np.linalg.norm(bg.u,axis=1)*INFLATE
    def atomerr(i,j):
        return INFLATE*(eo*nu[j]+eu*no[i]+eo*eu+ga(16)*no[i]*nu[j])
    ii,jj=np.triu_indices(bg.m.size,1)
    minimum=np.inf;maxmove=0.;minev=np.inf
    for k in range(0,len(ii),96):
        i,j=ii[k:k+96],jj[k:k+96]
        oi=bg.o[:,i].T[:,:,None];oj=bg.o[:,j].T[:,:,None]
        ui=bg.u[i,None,:];uj=bg.u[j,None,:]
        raw=np.stack([oi*ui,oj*uj,oi*uj+oj*ui],axis=-1).reshape(len(i),-1,3)
        scales=np.linalg.norm(raw,axis=1)
        q=raw/scales[:,None,:]
        qh=q.conj().swapaxes(1,2)
        hn=qh@q;d=qh@r
        qnorm=np.linalg.norm(q,axis=(1,2))*INFLATE
        eh=INFLATE*ga(8*r.size+16)*qnorm**2
        il=np.linalg.inv(hn)
        ilnorm=np.linalg.norm(il,axis=(1,2))*INFLATE
        den=np.linalg.norm(np.eye(3)-hn@il,axis=(1,2))*INFLATE
        den+=(eh+ga(40)*np.linalg.norm(hn,axis=(1,2))*INFLATE)*ilnorm+ga(40)*2
        lam=(1-INFLATE*den)/(INFLATE*ilnorm)
        assert min(lam)>0
        x=(il@d[:,:,None])[:,:,0];xn=np.linalg.norm(x,axis=1)*INFLATE
        qres=r-(q@x[:,:,None])[:,:,0]
        qrn=np.linalg.norm(qres,axis=1)
        erv=INFLATE*ga(40)*(rn+qnorm*xn)
        ortho=np.linalg.norm((qh@qres[:,:,None])[:,:,0],axis=1)*INFLATE
        ortho+=INFLATE*qnorm*(erv+ga(8*r.size+16)*qrn*INFLATE)
        sq=np.maximum(0.,np.maximum(0.,qrn/INFLATE-erv)**2-(ortho*INFLATE)**2/lam)
        sub0=np.sqrt(sq)/INFLATE
        ec=atomerr(i,j)+atomerr(j,i)+ga(4)*(no[i]*nu[j]+no[j]*nu[i])
        eq=INFLATE*np.linalg.norm(np.stack([atomerr(i,i),atomerr(j,j),ec],axis=1)/scales,axis=1)+ga(8)*qnorm
        mover=eq/(np.sqrt(lam)/INFLATE-eq)
        assert min(mover)>=0
        lb=np.maximum(0.,sub0-er-rn*mover*INFLATE)/INFLATE
        minimum=min(minimum,float(min(lb)));maxmove=max(maxmove,float(max(mover)))
        minev=min(minev,float(min(lam)))
    return {**di,"min_lower_distance":minimum,"radius":float(data["delta"]),
            "min_lower_distance_by_delta":minimum/float(data["delta"]),"max_subspace_gap_bound":maxmove,
            "min_normalized_gram_eigen_lower":minev,"pairs":len(ii),"R_error_upper":er}


if __name__=="__main__":
    dest=LAB/"results/support/all_pairs_verified.json"
    result={}
    with threadpool_limits(1):
        for n in [17,44]:
            wall,cpu=perf_counter(),process_time()
            data,_=public_background(n)
            with np.load(historical()/f"RESULTS/confirmation_core48/r13_fresh_{n:04d}__dual_cie/primary.npz") as f:m=f["m"].ravel()
            bg=prepare(data,m)
            out=bound_all(data,bg)
            out.update(wall_s=perf_counter()-wall,cpu_s=process_time()-cpu)
            result[f"r13_fresh_{n:04d}"]=out
            dest.write_text(json.dumps(result,indent=2)+"\n")
            print(n,out["min_lower_distance_by_delta"],flush=True)
