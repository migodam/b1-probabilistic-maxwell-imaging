"""Direct source-threshold test on the bordered constitutive form of a finite box.

The strict pass performs outward ball congruence on the selected positive
subspace. A second route uses the standard binary64 dot-product error bound.
"""
from pathlib import Path
from time import perf_counter,process_time
import sys,argparse,json
import numpy as np
from scipy.linalg import eigh
from flint import arb,acb_mat,ctx
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
POOL=HERE.parent/'pool'
sys.path.insert(0,str(POOL))
from pool_common import FullWave
from box_certificates import ball,star,diagonal,gershgorin,lower
UROUND=2.**-53
WEIGHTS=(1/32,1/128,1/8,1/512,1/2)
def clock(): return perf_counter(),process_time()
def elapsed(t): return dict(zip(('wall_seconds','cpu_seconds'),[a-b for a,b in zip(clock(),t)]))
def fn(a):
    # Finite normal IEEE inputs. Two includes hypot and positive-sum roundoff.
    return 2.*float(np.linalg.norm(a))+1e-260
def matmul_err(a,b,p): return (32*p+32)*UROUND*a*b+1e-255
def save(name,obj):
    p=HERE/name;p.with_suffix('.tmp').write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');p.with_suffix('.tmp').replace(p)

class Border:
    def __init__(self,g,s,b):
        self.n,self.p=g.shape[0],b.shape[1]
        self.s=s
        self.w=np.column_stack((np.eye(self.n),np.zeros((self.n,self.p)))).astype(complex)
        self.e=np.column_stack((g,b))
        self.norm_e=fn(self.e)
        self.sc=float(np.linalg.eigvalsh(s@s.conj().T)[-1])
        self.nw=fn(self.w)
        self.exact=None
    def measurement(self,y):
        self.v=np.column_stack((self.s,-y))
        self.d=self.v.conj().T@self.v
        self.d_error=matmul_err(fn(self.v),fn(self.v),self.v.shape[0])
    def material(self,lo,hi):
        # (J-lE)^*(J-uE) supplies the Hermitian form for each real cell.
        loww=self.w-lo[:,None]*self.e
        upperw=self.w-hi[:,None]*self.e
        def rowerr(x): return 12*UROUND*(self.nw+(1+np.max(np.abs(x)))*self.norm_e)+1e-255
        ea, eb = rowerr(lo),rowerr(hi)
        a,b=fn(loww),fn(upperw)
        tmp=loww.conj().T@upperw
        h=(tmp+tmp.conj().T)*.5
        err=2*(a*eb+b*ea+ea*eb+matmul_err(a,b,self.n)+6*UROUND*fn(tmp))
        return h,err
    def binary(self,q,vals,v,assembly):
        # Q has the proposed positive columns. Re-evaluate its residual and Gram;
        # only the displayed error bounds rely on the accuracy of BLAS products.
        dim,width=q.shape
        g=q.conj().T@q
        g.flat[::width+1]-=1
        bq=fn(q)
        eps_g=fn(g)+2*matmul_err(bq,bq,dim)+8*UROUND*(width+1)
        if eps_g>=.25: return None
        av=v@q
        qd=q*vals
        gap=float(vals[0])
        residual=fn(av-qd)+2*(matmul_err(fn(v),bq,dim)+6*UROUND*(fn(av)+fn(qd)))+2*assembly
        # ||Q|| <= sqrt(1+eps_g)<2; symmetrization is included in assembly.
        slack=gap-2*(eps_g*max(abs(vals))+2*residual+1e-248)
        return slack if np.isfinite(slack) and slack>0 else None
    def strict(self,lo,hi,y,scale,th,q,raw=False):
        if self.exact is None:
            self.exact=ball(self.w),ball(self.e)
        w,e=self.exact
        v=ball(np.column_stack((self.s,-y)))
        proj=ball(q)
        a,b=(w-diagonal([arb(float(x)) for x in endpoint])*e for endpoint in (lo,hi))
        vq=v*proj
        lq=a*proj; uq=b*proj
        cross=star(lq)*uq
        sour=acb_mat(self.p,q.shape[1],[proj[i,j] for i in range(self.n,self.n+self.p) for j in range(q.shape[1])])
        result=star(vq)*vq + arb(float(scale))*(cross+star(cross))/2-arb(float(th))*star(sour)*sour
        return lower(gershgorin(result,upper=False))
    def check(self,lo,hi,y,radius,*,verify=False,weights=WEIGHTS,jvalues=None):
        h,delta=self.material(lo,hi)
        hnorm=fn(h); dnorm=fn(self.d)
        tested=0
        for scale in [self.sc*k for k in weights]:
            a=self.d+scale*h
            for j in ((1,self.p,2) if jvalues is None else jvalues):
                if j>self.p: continue
                # margin dwarfs last-bit rounding in r^2/j; validate r^2 < j*th exactly below.
                th=float(radius*radius*(1+2**-15)/j)
                if not arb(float(radius))**2<arb(th)*j: continue
                c=a.copy()
                c.flat[self.n*(len(a)+1)::len(a)+1]-=th
                hsym=(c+c.conj().T)*.5
                round_error=2*(self.d_error+scale*delta+24*UROUND*(dnorm+abs(scale)*hnorm+abs(th)*(self.p+1))+1e-250)
                vals,vectors=eigh(hsym,check_finite=False,driver='evd')
                tested+=1
                first=self.p-j
                if vals[first]<=0: continue
                q=vectors[:,first:]
                margin=self.binary(q,vals[first:],hsym,round_error)
                if margin is None: continue
                out=dict(rejected=True,trials=tested,j=j,weight=scale,threshold=th,
                         binary_margin_lower=float(margin),assembly_error_upper=round_error)
                if verify:
                    t=clock()
                    out['arb_congruence_lower']=self.strict(lo,hi,y,scale,th,q)
                    out['independent_check']=elapsed(t)
                    if out['arb_congruence_lower']<=0: raise ArithmeticError('Independent outward congruence failed')
                return out
        return dict(rejected=False,trials=tested)
def scan(path,output):
    ctx.prec=104
    t=clock()
    with np.load(POOL/'structured_results/public/geometry.npz') as a:
        g=a['green'];s=a['sensor'][a['train_rx']];b=a['incident'][:,a['train_tx']]
    checker=Border(g,s,b)
    out=dict(setup=elapsed(t),model='stored full 24² VIE; 49 independent unknowns on shared support',weights=WEIGHTS,rows={})
    for f in sorted(path.glob('guided*.npz')):
        with np.load(f) as a:
            y,r,lo,hi=a['train'],float(a['radius']),a['lower'],a['upper']
        t=clock();checker.measurement(y)
        rows=[checker.check(l,h,y,r) for l,h in zip(lo,hi)]
        wall=elapsed(t)
        t=clock()
        checked=[checker.check(l,h,y,r,verify=True) if dec['rejected'] else dec for l,h,dec in zip(lo,hi,rows)]
        support=elapsed(t)
        out['rows'][f.stem]=dict(readout_time=wall,accepted_proof_replay=support,decisions=checked)
        print(f.stem,sum(row['rejected'] for row in rows),wall,flush=True)
    save(output,out)
if __name__=='__main__':
    with threadpool_limits(1): scan(POOL/'structured_results','development.json')
