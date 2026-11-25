"""Proof primitives for exact stored-model inputs. All contractions use complex balls."""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import qr
from flint import acb, acb_mat, arb

def mat(rows):
    vals = np.asarray(rows, dtype=object)
    return acb_mat(vals.shape[0], vals.shape[1],
                   [z if isinstance(z,(acb,arb)) else acb(complex(z)) for z in vals.flat])

def adj(a):
    return a.conjugate().transpose()

def block(a, rows, cols):
    rows, cols = list(rows), list(cols)
    return acb_mat(len(rows),len(cols),[a[i,j] for i in rows for j in cols])

def arr(a):
    return np.array([[complex(a[i,j]) for j in range(a.ncols())] for i in range(a.nrows())])

def diag(seq):
    seq=list(seq)
    return acb_mat(len(seq),len(seq),[x if i==j else 0 for i,x in enumerate(seq) for j in range(len(seq))])

def assemble_border(G,S,B,V,Y,lower,upper,weight,optical):
    """Dictionary input B, measured combinations V, kernel and affine signs as in the paper."""
    e=B*V
    n=G.nrows()
    mid=diag([(arb(l)+arb(u))/2 for l,u in zip(lower,upper,strict=True)])
    w=diag(weight)
    a=diag([arb(s)*((arb(u)-arb(l))/2)**2 for s,l,u in zip(weight,lower,upper,strict=True)])
    la=diag(optical)
    ident=diag([1]*n)
    state=ident-mid*G
    k=adj(S)*S+adj(state)*w*state-adj(G)*a*G+(la*G-adj(G)*la)/acb(0,2)
    l=adj(S)*Y+(adj(state)*w*mid+adj(G)*a+acb(0,.5)*la)*e
    c=adj(Y)*Y+adj(e)*(mid*w*mid-a)*e
    return mat(np.block([[np.array(k.tolist(),object),-np.array(l.tolist(),object)],
                         [-np.array(adj(l).tolist(),object),np.array(c.tolist(),object)]]))

def local_pair_rows(G,B,V,ids):
    n=G.nrows()
    old=B*V
    t=[];e=[]
    for i in ids:
        t.extend([[int(j==i) for j in range(n)],[G[i,j] for j in range(n)]])
        e.extend([[0]*V.ncols(),[old[i,j] for j in range(V.ncols())]])
    return mat(t),mat(e)

def positive_dimension(a):
    """A posteriori, two *lower* counts. Difference is undecided, not an exact null count."""
    f=arr(a)
    ev,v=np.linalg.eigh((f+f.conj().T)/2)
    counts=[]
    for sign in (1,-1):
        inds=np.flatnonzero(sign*ev>0)
        inds=inds[np.argsort(-sign*ev[inds])]
        if not len(inds):
            counts.append(0);continue
        frame=mat(v[:,inds])
        restricted=sign*adj(frame)*a*frame
        q=len(inds)
        while q:
            bounds=[restricted[i,i].real-sum((abs(restricted[i,j]) for j in range(q) if j!=i),arb(0)) for i in range(q)]
            if all(x>0 for x in bounds):break
            q-=1
        counts.append(q)
    return counts

@dataclass
class FixedKernel:
    base:acb_mat
    trace:acb_mat
    dim:int
    eliminated:int
    pivots:list
    gram_rank_scope:str

def build_fixed_kernel(border,trace,affine,n,extra=None):
    """Exact kernel of named covectors through an interval pivot and a proved positive complement.

    Extra numerical covectors define exact stored functionals; an approximate invariant
    subspace is never asserted to be invariant. Redundant rows must be removed algebraically.
    """
    t=trace if extra is None else mat(trace.tolist()+extra.tolist())
    q=t.nrows(); s=border.nrows()-n
    col=qr(arr(t),pivoting=True,mode='economic')[2].tolist()
    pivot,free=col[:q],col[q:]
    if q>n:raise ValueError('too many independent covectors')
    kp=block(t,range(q),pivot)
    if q and free:
        wr=-kp.solve(block(t,range(q),free))
    W=acb_mat(n,len(free))
    for j,i in enumerate(free):W[i,j]=1
    for i,r in enumerate(pivot):
        for j in range(len(free)):W[r,j]=wr[i,j]
    # Retain coordinate columns indexed by the pivot, followed by the source identity.
    retained=pivot+list(range(n,n+s))
    raw=block(border,retained,retained)
    if free:
        kk=block(border,range(n),range(n))
        cc=adj(W)*kk*W
        if positive_dimension(cc)[0]!=len(free):raise ValueError('positive complement unproved')
        x=adj(W)*block(border,range(n),retained)
        raw-=adj(x)*cc.solve(x)
    f=mat([parts+aff for parts,aff in zip(block(trace,range(trace.nrows()),pivot).tolist(),affine.tolist(),strict=True)])
    return FixedKernel(raw,f,q,len(free),pivot,'interval pivot and positive physical complement')

def update_box(kernel,delta):
    return kernel.base+adj(kernel.trace)*delta*kernel.trace

def certify_inertia(a,n,r2,j,t):
    s=a.nrows()-n
    shifted=a-diag([0]*n+[t]*s)
    pos,neg=positive_dimension(shifted)
    return dict(positive_at_least=pos,negative_at_least=neg,unresolved_signs=n+s-pos-neg,
                j=j,t=float(t),decision='EXCLUDED' if pos>=n+j and arb(j)*t>=r2 else 'UNRESOLVED')
