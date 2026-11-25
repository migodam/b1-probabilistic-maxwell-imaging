"""Interval charts for affine rectangular complex equalities A(xi)x = B(xi)Z."""
from dataclasses import dataclass, field
from hashlib import blake2b
from itertools import combinations
from numbers import Integral
from time import perf_counter

from flint import acb, acb_mat, arb, ctx
import numpy as np
from scipy.linalg import qr

from fixed_trace import adj, arr, block, diag, mat


def _binary_array(a):
    a=np.array(a,dtype=np.complex128,copy=True)
    if not np.all(np.isfinite(a)):
        raise ValueError('Affine coefficients must be finite binary64 numbers.')
    a.flags.writeable=False
    return a


@dataclass(frozen=True)
class AffineEqualities:
    A0: np.ndarray
    B0: np.ndarray
    A_coeff: tuple = ()
    B_coeff: tuple = ()

    def __post_init__(self):
        for k in ('A0','B0'):
            object.__setattr__(self,k,_binary_array(getattr(self,k)))
        if self.A0.ndim!=2 or self.B0.ndim!=2 or self.A0.shape[0]!=self.B0.shape[0] or not all(self.A0.shape):
            raise ValueError('The equality and source coefficient matrices must have matching rows.')
        for k,base in (('A_coeff',self.A0),('B_coeff',self.B0)):
            terms=tuple(_binary_array(a) for a in getattr(self,k))
            if any(a.shape!=base.shape for a in terms):
                raise ValueError('The affine terms must use their respective constant matrix shape.')
            object.__setattr__(self,k,terms)
        if len(self.A_coeff)!=len(self.B_coeff):
            raise ValueError('A and B must index the same real parameters.')

    @property
    def system_hash(self):
        digest=blake2b(digest_size=20)
        for a in (self.A0,self.B0,*self.A_coeff,*self.B_coeff):
            digest.update(repr(a.shape).encode())
            digest.update(a.tobytes(order='C'))
        return digest.hexdigest()


@dataclass
class ChartLeaf:
    lower: tuple
    upper: tuple
    path: str
    status: str
    row_pivot: tuple = ()
    column_pivot: tuple = ()
    free_columns: tuple = ()
    remaining_rows: tuple = ()
    lift: acb_mat | None = None
    consistency: acb_mat | None = None
    pivot_receipt: dict = field(default_factory=dict)
    attempted_pivots: list = field(default_factory=list)
    reason: str | None = None
    system_hash: str | None = None

    def record(self):
        value=dict(vars(self))
        for name in ('lift','consistency'):
            x=value[name]
            value[name]=None if x is None else dict(shape=[x.nrows(),x.ncols()],entries=[
                [_lower(z.real),_upper(z.real),_lower(z.imag),_upper(z.imag)] for z in x.entries()])
        return value


@dataclass
class ChartCover:
    status: str
    leaves: list
    subdivisions: list
    precision_bits: int
    wall_seconds: float
    system_hash: str

    def record(self):
        return dict(status=self.status,leaves=[leaf.record() for leaf in self.leaves],
                    subdivisions=self.subdivisions,precision_bits=self.precision_bits,
                    wall_seconds=self.wall_seconds,system_hash=self.system_hash)


@dataclass(frozen=True)
class ChartPullback:
    full: acb_mat
    n: int
    consistency: acb_mat


def _lower(x):
    out=float(np.nextafter(float(x.lower()),-np.inf))
    return out if np.isfinite(out) else None


def _upper(x):
    out=float(np.nextafter(float(x.upper()),np.inf))
    return out if np.isfinite(out) else None


def _row_norm(m):
    return max((sum((abs(m[i,j]).upper() for j in range(m.ncols())),arb(0)).upper()
                for i in range(m.nrows())),default=arb(0))


def _hstack(a,b):
    return acb_mat(a.nrows(),a.ncols()+b.ncols(),
                   [z for i in range(a.nrows()) for m in (a,b) for z in (m[i,j] for j in range(m.ncols()))])


def _affine(base,coeff,lower,upper):
    out=mat(base)
    for value,a in zip((arb(l).union(arb(u)) for l,u in zip(lower,upper,strict=True)),coeff,strict=True):
        out+=acb(value)*mat(a)
    return out


def _backsolve(p,rhs,c,resolvent):
    """Enclose P^-1 B with a direct ball solve or its verified residual bound."""
    try:
        x=p.solve(rhs)
        if all(z.is_finite() for z in x.entries()):
            return x,'direct Arb solve'
    except (ArithmeticError,ValueError,ZeroDivisionError):
        pass
    with np.errstate(over='ignore',invalid='ignore'):
        numeric=arr(c)@arr(rhs)
    proposal=mat(numeric) if np.all(np.isfinite(numeric)) else acb_mat(p.nrows(),rhs.ncols())
    defect=rhs-p*proposal
    out=acb_mat(p.nrows(),rhs.ncols())
    for col in range(rhs.ncols()):
        err=max((abs(defect[i,col]).upper() for i in range(p.nrows())),default=arb(0))
        pad=arb(0,(err*resolvent).upper())
        for row in range(p.nrows()):
            out[row,col]=proposal[row,col]+acb(pad,pad)
    return out,'uniform inverse bound times column residual'


def _make_chart(a,b,lower,upper,path,rows,cols):
    n=a.ncols()
    p=block(a,rows,cols)
    try:
        approx=arr(p)
        if not np.all(np.isfinite(approx)):
            raise ValueError('midpoint is not finite binary64')
        inv=np.linalg.solve(approx,np.eye(len(cols)))
        if not np.all(np.isfinite(inv)):
            raise ValueError('inverse proposal is not finite binary64')
        c=mat(inv)
    except (np.linalg.LinAlgError,ValueError):
        return None,dict(rows=rows,columns=cols,reason='midpoint inverse proposal failed')
    eta=_row_norm(diag([1]*len(cols))-c*p)
    summary=dict(rows=rows,columns=cols,inverse_residual_upper=_upper(eta))
    if not eta<1:
        return None,dict(summary,reason='uniform inverse residual could not be bounded below one')
    bound=(_row_norm(c)/(1-eta)).upper()
    free=tuple(i for i in range(n) if i not in cols)
    remaining=tuple(i for i in range(a.nrows()) if i not in rows)
    rhs=_hstack(-block(a,rows,free),block(b,rows,range(b.ncols())))
    x,method=_backsolve(p,rhs,c,bound)
    s=b.ncols()
    lift=acb_mat(n+s,len(free)+s)
    for i,k in enumerate(free):
        lift[k,i]=1
    for i,k in enumerate(cols):
        for j in range(x.ncols()):
            lift[k,j]=x[i,j]
    for i in range(s):
        lift[n+i,len(free)+i]=1
    # The original unpivoted equations are imposed on the lifted variables.
    consistency=_hstack(block(a,remaining,range(n)),-block(b,remaining,range(s)))*lift
    proof=dict(summary,inverse_norm_upper=_upper(bound),
               condition_upper=_upper((_row_norm(p)*bound).upper()),solve_method=method)
    return ChartLeaf(lower,upper,path,'OUTWARD_CHART',rows,cols,free,remaining,lift,consistency,proof),proof


def _pivots(a,size,preferred):
    seen=set()
    aprox=arr(a)
    qr_proposal=()
    if np.all(np.isfinite(aprox)):
        try:
            columns=tuple(map(int,qr(aprox,pivoting=True,mode='economic')[2][:size]))
            rows=tuple(map(int,qr(aprox[:,columns].T,pivoting=True,mode='economic')[2][:size]))
            qr_proposal=((rows,columns),)
        except (ValueError,np.linalg.LinAlgError):
            pass
    # QR is only a proposal; direct enumeration then tries alternatives up to the budget.
    for pivots in (preferred, qr_proposal,
                   ((r,c) for r in combinations(range(a.nrows()),size)
                          for c in combinations(range(a.ncols()),size))):
        for r,c in pivots:
            key=(tuple(sorted(r)),tuple(sorted(c)))
            if key not in seen:
                seen.add(key)
                yield key


def rectangular_equality_charts(system,lower,upper,*,pivot_size=None,preferred_pivots=(),
                                max_pivots_per_leaf=16,max_depth=5,max_leaves=32,seconds=None):
    """Cover a real parameter box by verified charts and explicit unresolved leaves."""
    start=perf_counter()
    lower,upper=tuple(map(float,lower)),tuple(map(float,upper))
    if (len(lower)!=len(system.A_coeff) or len(upper)!=len(lower) or not np.all(np.isfinite((lower,upper)))
            or any(l>u for l,u in zip(lower,upper,strict=True))):
        raise ValueError('Use finite ordered closed intervals for every parameter.')
    m,n=system.A0.shape
    q=min(m,n) if pivot_size is None else pivot_size
    if (isinstance(q,bool) or int(q)!=q or not 1<=q<=min(m,n)
            or any(isinstance(x,bool) or int(x)!=x for x in (max_pivots_per_leaf,max_depth,max_leaves))
            or max_pivots_per_leaf<1 or max_depth<0 or max_leaves<1
            or seconds is not None and (not np.isfinite(seconds) or seconds<0)):
        raise ValueError('Invalid pivot or subdivision budget.')
    q=int(q)
    preferred=tuple((tuple(r),tuple(c)) for r,c in preferred_pivots)
    for r,c in preferred:
        if len(set(r))!=q or len(r)!=q or len(set(c))!=q or len(c)!=q or not (
                all(isinstance(i,Integral) and i in range(m) for i in r)
                and all(isinstance(j,Integral) and j in range(n) for j in c)):
            raise ValueError('Preferred pivots require valid distinct row and column indexes.')
    preferred=tuple((tuple(map(int,r)),tuple(map(int,c))) for r,c in preferred)
    deadline=float('inf') if seconds is None else start+seconds
    identity=system.system_hash
    leaves,splits=[],[]
    stack=[(lower,upper,'root',0)]
    while stack:
        lo,hi,path,depth=stack.pop()
        attempts=[]
        chart=None
        if perf_counter()<deadline:
            a=_affine(system.A0,system.A_coeff,lo,hi)
            b=_affine(system.B0,system.B_coeff,lo,hi)
            for r,c in _pivots(a,q,preferred):
                if len(attempts)>=max_pivots_per_leaf or perf_counter()>=deadline:
                    break
                chart,info=_make_chart(a,b,lo,hi,path,r,c)
                if chart is not None and perf_counter()>=deadline:
                    chart=None
                    info=dict(info,reason='pivot proof completed after the time budget')
                attempts.append(info)
                if chart is not None:
                    break
        if chart is not None:
            chart.attempted_pivots=attempts
            chart.system_hash=identity
            leaves.append(chart)
            continue
        # Use the widest splittable coordinate. Shared closed endpoints preserve coverage.
        widths=[(b/2-a/2,i,a/2+b/2) for i,(a,b) in enumerate(zip(lo,hi,strict=True)) if a<a/2+b/2<b]
        can_continue=perf_counter()<deadline
        if widths and can_continue and depth<max_depth and len(leaves)+len(stack)+2<=max_leaves:
            _,i,mid=max(widths,key=lambda pair:(pair[0],-pair[1]))
            left,right=list(hi),list(lo)
            left[i]=right[i]=mid
            splits.append(dict(path=path,parameter=i,cut=mid,failed_attempts=attempts))
            stack.extend(((tuple(right),hi,path+'R',depth+1),(lo,tuple(left),path+'L',depth+1)))
        else:
            reason=('time budget exhausted' if not can_continue else
                    'subdivision budget exhausted' if widths else 'no certified pivot on an unsplittable box')
            leaves.append(ChartLeaf(lo,hi,path,'UNRESOLVED',attempted_pivots=attempts,reason=reason,
                                    system_hash=identity))
    return ChartCover('OUTWARD_COVER' if all(x.status=='OUTWARD_CHART' for x in leaves) else 'UNRESOLVED',
                      leaves,splits,ctx.prec,perf_counter()-start,identity)


def evaluate_chart(system,leaf,point):
    """Independently re-evaluate the same rational chart on one parameter point."""
    point=tuple(map(float,point))
    if (leaf.status!='OUTWARD_CHART' or leaf.system_hash!=system.system_hash or len(point)!=len(leaf.lower) or not all(
            l<=t<=u for l,t,u in zip(leaf.lower,point,leaf.upper,strict=True))):
        raise ValueError('Point must belong to a verified chart leaf.')
    a=_affine(system.A0,system.A_coeff,point,point)
    b=_affine(system.B0,system.B_coeff,point,point)
    result,_=_make_chart(a,b,point,point,leaf.path,leaf.row_pivot,leaf.column_pivot)
    if result is None:
        raise ArithmeticError('Point replay needs higher precision.')
    result.system_hash=leaf.system_hash
    return result


def pullback_border(border,leaf,threshold=0.):
    """Enclose the source-threshold pullback and carry all residual equalities."""
    if leaf.status!='OUTWARD_CHART':
        raise ValueError('A pullback requires a verified chart.')
    a=border if isinstance(border,acb_mat) else mat(border)
    p=leaf.lift
    if a.nrows()!=p.nrows() or a.ncols()!=p.nrows() or not np.isfinite(threshold) or threshold<0:
        raise ValueError('Border dimensions and nonnegative threshold must match the chart.')
    states=len(leaf.free_columns)+len(leaf.column_pivot)
    shifted=a-diag([0]*states+[arb(float(threshold))]*(p.nrows()-states))
    return ChartPullback(adj(p)*shifted*p,len(leaf.free_columns),leaf.consistency)
