"""Mechanics shared by the prospective full-TM acquisition and open design."""
from dataclasses import dataclass
from hashlib import blake2b
import json
import os
from pathlib import Path
from time import perf_counter, process_time

from flint import acb, acb_mat, arb, ctx
import numpy as np
from scipy.linalg import eigh
from scipy.special import hankel1, jv

COMPONENT = Path(__file__).resolve().parent
HERE = COMPONENT/"diagnostic" if os.environ.get("ZALG_EXPERIMENT_C_DIAGNOSTIC")=="1" else COMPONENT
ROOT = COMPONENT.parents[1]
REGISTRY = HERE/'registration.json'
SEEDS = {"asymmetric":30174, "half":30175, "occluded":30176}
PORTS = {
    "asymmetric": ((17, 3, 2), (123, 2, 1), (239, 9, 5)),
    "half": ((-65, 2, 1), (0, 2, 1), (65, 2, 1)),
    "occluded": ((30, 2, 1), (60, 2, 1), (240, 2, 1)),
}
ANGLES = tuple(range(0, 360, 30))
OLD = (0, 6)
BATCHES = (1, 2, 4)
CONTRACTS = ('added', 'fixed_total')
POLICIES = ('random', 'uniform_gap', 'fisher_nuisance', 'conditional_ellipsoid')
# The same deterministic receive-noise bound per physically fired beam.
COLUMN_RADIUS = "0.00018"
WEIGHTS = (.0078125, .03125, .125, .5, 2., 8., 32., 128.)
SIDE = 4
WAVELENGTH = "0.75"
WIDTH_WAVELENGTHS = "0.70"


@dataclass(frozen=True)
class Geometry:
    layout: str
    g: np.ndarray
    s: np.ndarray
    b: np.ndarray
    ex_g: acb_mat
    ex_s: acb_mat
    ex_b: acb_mat


def clocks():
    return perf_counter(), process_time()


def elapsed(start):
    return dict(wall_seconds=perf_counter()-start[0], cpu_seconds=process_time()-start[1])


def hash_arrays(*arrays):
    out=blake2b(digest_size=20)
    for a in arrays:
        a=np.ascontiguousarray(a)
        out.update(a.dtype.str.encode())
        out.update(repr(a.shape).encode())
        out.update(a.tobytes())
    return out.hexdigest()


def hash_file(path):
    return blake2b(path.read_bytes(),digest_size=20).hexdigest()


def write_json(path,payload):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.part')
    tmp.write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    tmp.replace(path)


def read_json(path):
    return json.loads(path.read_text())


def adj(a):
    return a.conj().T


def herm(a):
    return (a+adj(a))/2


def star(a):
    return a.conjugate().transpose()


def ball(a):
    return acb_mat(a.shape[0],a.shape[1],[acb(float(z.real),float(z.imag)) for z in a.flat])


def mid(a):
    return np.fromiter((complex(z) for z in a.entries()),complex,a.nrows()*a.ncols()).reshape(a.nrows(),a.ncols())


def eye(n):
    return acb_mat(n,n,[int(i==j) for i in range(n) for j in range(n)])


def diag(v):
    return acb_mat(len(v),len(v),[v[i] if i==j else 0 for i in range(len(v)) for j in range(len(v))])


def minor(a, rows, cols):
    rows,cols=list(rows),list(cols)
    return acb_mat(len(rows),len(cols),[a[i,j] for i in rows for j in cols])


def hjoin(a,b):
    return acb_mat(a.nrows(),a.ncols()+b.ncols(),[m[i,j] for i in range(a.nrows())
                                                 for m in (a,b) for j in range(m.ncols())])


def vjoin(a,b):
    return acb_mat(a.nrows()+b.nrows(),a.ncols(),a.entries()+b.entries())


def norm2(a):
    return sum((z.real*z.real+z.imag*z.imag for z in a.entries()),arb(0))


def norm(a):
    return norm2(a).nonnegative_part().sqrt()


def low(x):
    return float(np.nextafter(float(x.lower()),-np.inf))


def high(x):
    return float(np.nextafter(float(x.upper()),np.inf))


def spectrum(a):
    """Complete verified Hermitian spectrum via polar-unitary residual bound."""
    v,u=eigh(herm(mid(a)),check_finite=False)
    proposal,d=ball(u),diag([arb(float(w)) for w in v])
    gram=norm(star(proposal)*proposal-eye(len(v)))
    if not gram<1:
        raise ArithmeticError('Hermitian eigencolumns failed the polar full-rank check.')
    # For proposal = U*H, the eigenvalues of H*D*H are those of the
    # generalized pencil (D, H^(-2)). Min-max on each sign bounds their
    # displacement from the ordered diagonal D by ||D|| times the defect.
    error=(norm(a-proposal*d*star(proposal))+gram*arb(float(max(abs(v),default=0.)))).upper()
    return [arb(float(w))-error for w in v],[arb(float(w))+error for w in v],error


def inertia(a):
    lo,hi,error=spectrum(a)
    resolved=all(x>0 or y<0 for x,y in zip(lo,hi,strict=True))
    gap=min(((x if x>0 else -y).lower() for x,y in zip(lo,hi,strict=True)),default=arb(0))
    return (dict(negative=sum(v<0 for v in hi), positive=sum(v>0 for v in lo),
                 full_count_verified=resolved, gap_lower=low(gap) if resolved else None,
                 hermitian_eigensystem_error_upper=high(error)),
            lo,hi)


def geometry(layout):
    # Coordinates in wavelengths. The volume integral uses outgoing H^(1)
    # with the exact pulse-disk self integral and disk-integrated off-diagonal.
    dx=arb(WIDTH_WAVELENGTHS)/SIDE
    k=2*arb.pi()
    q=k*dx/arb.pi().sqrt()
    coef=arb.pi()*q*q.bessel_j(1)/2

    def h(order, v):
        return acb(v.bessel_j(order),v.bessel_y(order))

    grid=[(dx*(2*(i%SIDE)+1-SIDE)/2,dx*(2*(i//SIDE)+1-SIDE)/2) for i in range(SIDE**2)]
    kc={0:acb(0,arb.pi()*q/2)*h(1,q)-1}
    for d2 in {i*i+j*j for i in range(SIDE) for j in range(SIDE)}-{0}:
        kc[d2]=acb(0,coef)*h(0,k*dx*arb(d2).sqrt())
    g=acb_mat(len(grid),len(grid),[kc[(i%SIDE-j%SIDE)**2+(i//SIDE-j//SIDE)**2]
                                    for i in range(len(grid)) for j in range(len(grid))])
    pos=[]
    for d,rr,dd in PORTS[layout]:
        t=arb(d)*arb.pi()/180
        r=arb(rr)/dd
        pos.append((r*t.cos(),r*t.sin()))
    # Equal unit integration mass across every three-receiver layout.
    s=acb_mat(len(pos),len(grid),[acb(0,coef)/arb(len(pos)).sqrt()*h(0,k*((x-u)**2+(y-v)**2).sqrt())
                                 for u,v in pos for x,y in grid])
    dirs=[(t.cos(),t.sin()) for t in (arb(d)*arb.pi()/180 for d in ANGLES)]
    b=acb_mat(len(grid),len(dirs),[acb(0,k*(x*u+y*v)).exp() for x,y in grid for u,v in dirs])
    return Geometry(layout,mid(g),mid(s),mid(b),g,s,b)


def energy(contract, size):
    """Per-beam exact energy for old and new independent shots."""
    if contract=='added':
        return arb(1),arb(1)
    if contract=='fixed_total':
        return arb(1)/2,arb(1)/size
    raise ValueError(contract)


def noise_radius(registration, scene):
    return registration.get("noise_radius_by_scene", {}).get(scene,COLUMN_RADIUS)


def radius_squared(count, column_radius=COLUMN_RADIUS):
    return count*arb(column_radius)**2


def selection_matrix(ids, amplitude):
    return acb_mat(len(ANGLES),len(ids),[amplitude if i==j else 0 for i in range(len(ANGLES)) for j in ids])


def numeric_forward(geo,x,ids):
    x=np.asarray(x)
    a=np.eye(len(x))-x[:,None]*geo.g
    b=geo.b[:,list(ids)]
    j=np.linalg.solve(a,x[:,None]*b)
    return geo.s@j


def forward_jacobian(geo,x,ids):
    a=np.eye(len(x))-x[:,None]*geo.g
    b=geo.b[:,list(ids)]
    j=np.linalg.solve(a,x[:,None]*b)
    e=b+geo.g@j
    response=np.linalg.solve(a.T,geo.s.T).T
    jac=np.einsum('ri,is->rsi',response,e).reshape(-1,len(x))
    return geo.s@j, np.vstack((jac.real,jac.imag))


def exact_forward(geo,x,ids,amplitude):
    xd=diag([arb(float(v)) for v in x])
    a=eye(len(x))-xd*geo.ex_g
    if not abs(a.det()).lower()>0:
        raise ArithmeticError('Physical state uniqueness not verified.')
    j=a.solve(xd*geo.ex_b*selection_matrix(ids,amplitude))
    return geo.ex_s*j


def require_p0():
    gate=ROOT/'results/p0.json'
    if not gate.is_file():
        raise RuntimeError(f'P0 evidence not available: {gate}')
    data=read_json(gate)
    status=data.get('status',data.get('decision',data.get('verdict')))
    explicit=data.get('passed') is True or data.get('p0_pass') is True
    accepted={'pass','passed','p0_pass','p0_passed','p0_passed_strict_checks','complete_pass'}
    if not explicit and str(status).strip().lower() not in accepted:
        raise RuntimeError(f'P0 is not recorded as passed: {status}')
    return dict(receipt=str(gate.relative_to(ROOT)),status=status)


def init_arithmetic():
    ctx.prec,ctx.threads=128,1
