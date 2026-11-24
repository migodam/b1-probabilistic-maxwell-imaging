"""Native one-sided Fourier certificates for fixed spatial material boxes."""

from dataclasses import dataclass
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.fft import fft2
from scipy.linalg import eigh, solve

HERE=Path(__file__).resolve().parent
OPTICS=HERE.parents[1]/'z_alg_optical_inertia'
sys.path.insert(0,str(OPTICS))
from core.certificate import ball, frobenius_squared, gershgorin, safe_float, star


def clock():
    return perf_counter(),process_time()


def elapsed(start):
    return dict(wall_seconds=perf_counter()-start[0],cpu_seconds=process_time()-start[1])


def hi(a):
    return safe_float(a,upper=True)


def lo(a):
    return safe_float(a,upper=False)


def herm(a):
    return (a+a.conj().T)*.5


def mid(a):
    return np.array([[complex(v) for v in row] for row in a.tolist()])


def block(a,rows,cols):
    rows,cols=range(*rows),range(*cols)
    return acb_mat(len(rows),len(cols),[a[i,j] for i in rows for j in cols])


def identity(side):
    a=acb_mat(side,side)
    for i in range(side):
        a[i,i]=1
    return a


@dataclass(frozen=True)
class Branch:
    """Box enforced by every later step in this branch; tile cap zero means exact air."""
    name: str
    caps: tuple[tuple[float, ...], ...]


class NativeBox:
    def __init__(self,spectrum,sensor,incident,*,root_precision=96):
        start=clock()
        self.g=np.asarray(spectrum,dtype=complex)
        self.s=np.asarray(sensor,dtype=complex)
        self.b=np.asarray(incident,dtype=complex)
        self.p=self.g.shape[0]
        self.side=int(np.sqrt(self.s.shape[1]))
        self.n=self.side**2
        self.nr,self.ns=self.s.shape[0],self.b.shape[1]
        if (self.g.shape!=(self.p,self.p) or self.b.shape!=(self.n,self.ns) or
                self.s.shape[1]!=self.n or self.p<2*self.side-1):
            raise ValueError('Require native square convolution with enough padding and common receiver rows.')
        self.z=np.column_stack((self.s.conj().T,self.b))
        self.precision=root_precision
        ctx.prec,ctx.threads=root_precision,1
        self.roots=[acb(-2*arb(k)/self.p).exp_pi_i() for k in range(self.p)]
        self.symbol_cache={}
        self.dcache={}
        self.geometry_cost=elapsed(start)

    def profile(self,c,t):
        key=c,t
        if key not in self.symbol_cache:
            d=-self.g.real+t*self.g.imag-c
            modes=np.stack(np.where(d<0),axis=1)
            if len(modes)==0:
                raise ValueError('Use a nonempty native truncation or the explicit zero-material branch.')
            strength=np.sqrt(-d[modes[:,0],modes[:,1]])
            self.symbol_cache[key]=(modes,strength)
        return self.symbol_cache[key]

    def tile_weights(self,branch,c):
        u=np.asarray(branch.caps,dtype=float)
        shape=u.shape
        if u.ndim!=2 or self.side%shape[0] or self.side%shape[1] or not np.all(np.isfinite(u)) or np.any(u<0):
            raise ValueError('Nonnegative tile caps must exactly partition the native grid.')
        widths=self.side//shape[0], self.side//shape[1]
        denominator=np.ones_like(u)
        valid=u>0
        denominator[valid]=1/u[valid]+c
        if np.any(denominator[valid]<=0):
            raise ValueError('Positive-support diagonal lost positivity.')
        w=np.zeros_like(u)
        w[valid]=1/denominator[valid]
        return u,widths,w

    def numeric_grams(self,branch,c,t):
        u,widths,w=self.tile_weights(branch,c)
        weights=np.repeat(np.repeat(w,widths[0],axis=0),widths[1],axis=1)
        modes,coeff=self.profile(c,t)
        flat=weights.ravel()[:,None]*self.z
        z=self.z.conj().T@flat
        wz=fft2(flat.reshape(self.side,self.side,-1),s=(self.p,self.p),axes=(0,1))/self.p
        e=coeff[:,None]*wz[modes[:,0],modes[:,1],:]
        delta=(modes[:,None,:]-modes[None,:,:])%self.p
        wt=fft2(weights,s=(self.p,self.p))/self.p**2
        a=coeff[:,None]*wt[delta[:,:,0],delta[:,:,1]]*coeff[None,:]
        return a,e,z

    def numeric_reduction(self,blocks,alpha,y,t):
        a,e,z=blocks
        m=z[:self.nr,:self.nr]+alpha*np.eye(self.nr)
        e0=e[:,:self.nr]
        rz=solve(herm(m),np.column_stack((z[:self.nr],e0.conj().T)),assume_a='pos',check_finite=False)
        z0=z-z[:,:self.nr]@rz[:,:self.nr+self.ns]
        e1=e-e0@rz[:,:self.nr+self.ns]
        t0=herm(np.eye(a.shape[0])-a+e0@rz[:,self.nr+self.ns:])
        ev,v=eigh(t0,check_finite=False,driver='evd')
        nu=int(np.count_nonzero(ev<0))
        if nu>=self.ns or np.min(abs(ev))<1e-9:
            return 0.,nu,(ev,v)
        coef=np.vstack((y/alpha,((1+1j*t)*.5)*np.eye(self.ns)))
        b=e1@coef
        x=v@((v.conj().T@b)/ev[:,None])
        H=herm(y.conj().T@y-alpha*(coef.conj().T@z0@coef+b.conj().T@x))
        score=float(np.sqrt(np.maximum(0,eigh(H,eigvals_only=True)[::-1][nu:]).sum()))
        return score,nu,(ev,v)

    def exact_grams(self,branch,c,t):
        """Separable Fourier sums match the exact saved binary64 finite model."""
        u,widths,_=self.tile_weights(branch,c)
        modes,_=self.profile(c,t)
        retained={tuple(v) for v in modes}
        at,ac=arb(t),arb(c)
        strength=[]
        for i,j in modes:
            g=self.g[i,j]
            h=ac+arb(float(g.real))-at*arb(float(g.imag))
            if not h>0:
                raise ArithmeticError('Unresolved retained Fourier coefficient.')
            strength.append(h.sqrt())
        omission=arb(1000)
        for i,j in np.ndindex(self.g.shape):
            if (i,j) in retained:
                continue
            g=self.g[i,j]
            slack = -arb(float(g.real))+at*arb(float(g.imag))-ac
            if not slack>=0:
                raise ArithmeticError('Unresolved omitted Fourier modal floor.')
            omission=min(omission,slack.lower())
        w=[[1/(1/arb(float(v))+ac) if v>0 else arb(0) for v in row] for row in u]
        axes=[]
        for length,nt in zip(widths,u.shape,strict=True):
            axes.append([[sum((self.roots[d*k%self.p] for k in range(j*length,(j+1)*length)),acb(0))
                          for j in range(nt)] for d in range(self.p)])
        dy=(modes[:,None,0]-modes[None,:,0])%self.p
        dx=(modes[:,None,1]-modes[None,:,1])%self.p
        atoms=[(i,j,w[i][j]) for i,j in np.ndindex(u.shape) if w[i][j]>0]
        sums={}
        denominator=self.p*self.p
        for i,j in set(zip(dy.ravel().tolist(),dx.ravel().tolist())):
            sums[i,j]=sum((v*axes[0][i][a]*axes[1][j][b] for a,b,v in atoms),acb(0))/denominator
        rank=len(modes)
        a=acb_mat(rank,rank,[strength[i]*sums[int(dy[i,j]),int(dx[i,j])]*strength[j]
                            for i in range(rank) for j in range(rank)])
        axis_y,axis_x=np.unique(modes[:,0]),np.unique(modes[:,1])
        f_y=acb_mat(len(axis_y),self.side,[self.roots[int(k)*j%self.p] for k in axis_y for j in range(self.side)])
        f_x=acb_mat(len(axis_x),self.side,[self.roots[int(k)*j%self.p] for k in axis_x for j in range(self.side)])
        pi=[(int(np.searchsorted(axis_y,k)),int(np.searchsorted(axis_x,l))) for k,l in modes]
        weighted=[]
        proj=[]
        for column in self.z.T:
            field=ball(column.reshape(self.side,self.side))
            for i,j in np.ndindex(self.side,self.side):
                field[i,j]*=w[i//widths[0]][j//widths[1]]
            weighted.extend(field.entries())
            out=f_y*field*f_x.transpose()
            proj.extend(strength[i]*out[y,x]/self.p for i,(y,x) in enumerate(pi))
        ncol=self.nr+self.ns
        wz=acb_mat(ncol,self.n,weighted).transpose()
        e=acb_mat(ncol,rank,proj).transpose()
        z=star(ball(self.z))*wz
        return a,e,z,lo(omission)

    def interval_score(self,grams,alpha,t,y,proposal):
        a,e,z,omission=grams
        n=self.nr+self.ns
        nr=self.nr
        m=block(z,(0,nr),(0,nr))+arb(alpha)*identity(nr)
        e0=block(e,(0,e.nrows()),(0,nr))
        z0=block(z,(0,nr),(0,n))
        u=m.solve(z0)
        rzz=z-star(z0)*u
        rfx=e-e0*u
        t0=identity(a.nrows())-a+e0*m.solve(star(e0))
        ev,vectors=proposal
        v=ball(vectors)
        s=star(v)*v
        if not gershgorin(s,upper=False)>0:
            raise ArithmeticError('Unresolved basis for compact congruence.')
        t1=star(v)*t0*v
        gaps=[]
        neg=0
        for i in range(t1.nrows()):
            radius=sum((abs(t1[i,j]) for j in range(t1.ncols()) if i!=j),arb(0)).upper()
            d=t1[i,i].real
            if (d-radius).lower()>0:
                gaps.append((d-radius).lower())
            elif (d+radius).upper()<0:
                neg+=1
                gaps.append(-(d+radius).upper())
            else:
                raise ArithmeticError('Unresolved sign in compact inertia.')
        if neg>=self.ns:
            return 0.,neg,dict(compact_rank=t0.nrows(),omitted_symbol_slack_lower=omission)
        inv_bound=gershgorin(s,upper=True)/min(gaps)
        matrix_y=ball(y)
        c=acb_mat(n,self.ns)
        zz=acb(1,arb(t))/2
        for i in range(nr):
            for j in range(self.ns):
                c[i,j]=matrix_y[i,j]/arb(alpha)
        for j in range(self.ns):
            c[nr+j,j]=zz
        b=rfx*c
        x=ball(vectors@((vectors.conj().T@mid(b))/ev[:,None]))
        tx=t0*x
        loss=(arb(alpha)*inv_bound*frobenius_squared(tx-b)).upper()
        h=star(matrix_y)*matrix_y-arb(alpha)*(star(c)*rzz*c+star(b)*x+star(x)*b-star(x)*tx)
        h-=loss*identity(self.ns)
        _,w=eigh(herm(mid(h)),check_finite=False)
        terms=[]
        for k in range(neg+1,self.ns+1):
            first=ball(w[:,-k:])
            value=gershgorin(star(first)*h*first,upper=False)
            if not value>0:
                break
            if not gershgorin(star(first)*first,upper=False)>0:
                raise ArithmeticError('Unresolved orthogonality of proposed source eigenvectors.')
            terms.append((value/gershgorin(star(first)*first,upper=True)).lower())
        r=sum(terms,arb(0))
        lower=max(0.,lo(r.sqrt())) if r>0 else 0.
        return lower,neg,dict(compact_rank=t0.nrows(),omitted_symbol_slack_lower=omission,
                              norm_inverse_upper=hi(inv_bound),elimination_error_upper=hi(loss),
                              source_squared_eigenvalues_lower=list(map(lo,terms)))

    def score(self,y,radius,branch,*,t_grid=(0.,.25,.5,1.),c_grid=(.9375,.5),
              alpha_grid=(2.**-17,2.**-15,2.**-13,2.**-11,2.**-9),max_proofs=2,max_rank=850):
        """Return PRUNE only on a strict interval-certified residual above the declared radius."""
        start=clock()
        y=np.asarray(y,dtype=complex)
        if y.shape!=(self.nr,self.ns) or not np.all(np.isfinite(y)) or not np.isfinite(radius) or radius<0:
            raise ValueError('Require complete declared ports and finite nonnegative Frobenius radius.')
        upper=float(np.max(branch.caps))
        if upper==0:
            b=ball(y)
            value=max(0.,lo(frobenius_squared(b).sqrt()))
            return dict(branch=branch.name,decision='PRUNE' if value>radius else 'RETAIN',radius=radius,
                        residual_lower=value,rigor='exact_zero_material',**elapsed(start))
        screen=[]
        proposals=[]
        for c_factor in c_grid:
            c=-c_factor/upper
            for t in t_grid:
                modes,_=self.profile(c,t)
                if len(modes)>max_rank:
                    screen.append(dict(c=c,t=t,rank=len(modes),status='compact_rank_budget'))
                    continue
                raw=self.numeric_grams(branch,c,t)
                for alpha in alpha_grid:
                    lower,nu,eigen=self.numeric_reduction(raw,alpha,y,t)
                    row=dict(c=c,t=t,alpha=alpha,rank=len(modes),numeric_residual=lower,nu=nu)
                    screen.append(row)
                    if lower>radius:
                        proposals.append((lower,c,t,alpha,eigen))
        verify_start=clock()
        tried=[]
        value=0.
        for _,c,t,alpha,eigen in sorted(proposals,key=lambda item:item[0],reverse=True)[:max_proofs]:
            trial=clock()
            try:
                exact=self.exact_grams(branch,c,t)
                lower,nu,meta=self.interval_score(exact,alpha,t,y,eigen)
                value=max(lower,value)
                tried.append(dict(c=c,t=t,alpha=alpha,strict_residual_lower=lower,nu=nu,
                                  interval=meta,**elapsed(trial)))
                if value>radius:
                    break
            except (ArithmeticError,ZeroDivisionError,ValueError) as exc:
                tried.append(dict(c=c,t=t,alpha=alpha,error=str(exc),**elapsed(trial)))
        return dict(branch=branch.name,decision='PRUNE' if value>radius else 'RETAIN',radius=radius,
                    residual_lower=value,rigor='strict_native_lower_envelope' if value else 'no_pruning_certificate',
                    hypothesis_max_cap=upper,native_side=self.side,screen=screen,proofs=tried,
                    interval_time=elapsed(verify_start),**elapsed(start))
