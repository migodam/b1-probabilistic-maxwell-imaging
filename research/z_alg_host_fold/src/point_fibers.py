"""Count both point fibers for the analytic disk discretization, enclosing the original host solve.

The 20-root reduction follows the project's R2/R5 algorithm. Witnesses may also
be supplied by HomotopyContinuation.jl. Only the interval tests decide the count.
"""
from pathlib import Path
import argparse
import json
import time
from itertools import combinations, permutations
import numpy as np
from scipy.linalg import qr
from flint import arb, acb, acb_mat, ctx

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
ctx.prec = 160


def rball(value):
    return arb(str(value['lo'])).union(arb(str(value['hi'])))


def from_bounds(value):
    if isinstance(value, list):
        return [from_bounds(item) for item in value]
    return acb(rball(value['real']), rball(value['imag']))


def bounds(value):
    return [float(np.nextafter(float(value.lower()), -np.inf)),
            float(np.nextafter(float(value.upper()), np.inf))]


def midpoint(value):
    return complex(float(value.real.mid()), float(value.imag.mid()))


def exact(value):
    return acb(float(value.real), float(value.imag))


def identity(n):
    return acb_mat([[int(i == j) for j in range(n)] for i in range(n)])


def chart(host, obs, midpoint_B):
    B = host['B']
    _, _, pivot = qr(midpoint_B.T, pivoting=True)
    I = sorted(pivot[:3].tolist())
    J = sorted(set(range(6)) - set(I))
    inv = acb_mat([B[i] for i in I]).inv()
    delta = [exact(y)-y0 for y, y0 in zip(obs,host['y0'])]
    particular = acb_mat([delta])*inv
    changes = -acb_mat([B[i] for i in J])*inv
    W = acb_mat(4,6)
    for j,i in enumerate(I):
        W[0,i] = particular[0,j]
        for a in range(3):
            W[1+a,i] = changes[a,j]
    for a,i in enumerate(J):
        W[1+a,i] = 1
    V = W*acb_mat(host['H'])
    for i in range(6):
        V[0,i] += host['R'][0][i]
    coeff = []
    for i in range(6):
        A = [[W[a,i]*V[b,i].conjugate()-V[a,i]*W[b,i].conjugate()
              for b in range(4)] for a in range(4)]
        scale = max(abs(midpoint(x)) for row in A for x in row)
        coeff.append(acb_mat(A)/scale)
    return W,V,coeff,I


def exps(degree):
    return [(i,j,degree-i-j) for i in range(degree,-1,-1)
            for j in range(degree-i,-1,-1)]


def leading_reduction(coeff):
    # 15 maximal 4x4 minors of 6x4. Expanding determinants is done once per observation.
    basis = [e for d in range(4) for e in exps(d)]
    leading = exps(4)
    monomials = basis + leading
    zero = (0,0,0)
    result = []
    for rows in combinations(range(6),4):
        poly = {}
        for perm in permutations(range(4)):
            inv_count = sum(a>b for a,b in combinations(perm,2))
            term = {zero: acb((-1)**inv_count)}
            for i,j in zip(rows,perm):
                new = {}
                for e,v in term.items():
                    for q in range(4):
                        exp = list(e)
                        if q:
                            exp[q-1] += 1
                        exp = tuple(exp)
                        new[exp] = new.get(exp,0) + v*coeff[i][q,j]
                term = new
            for e,v in term.items():
                poly[e] = poly.get(e,0)+v
        result.append([poly[e] for e in monomials])
    T = acb_mat([r[len(basis):] for r in result])
    L = acb_mat([r[:len(basis)] for r in result])
    inverse = T.inv()
    reduction = inverse*(-L)
    return basis,leading,reduction,{'rank':15,'condition_midpoint':float(np.linalg.cond(
        np.array([[midpoint(v) for v in row] for row in T.tolist()])))}


def floating_coefficients(C):
    return np.array([[[midpoint(q) for q in row] for row in ci.tolist()] for ci in C])


def nominate(C, basis, leading, red):
    lookup = {e:i for i,e in enumerate(basis)}
    rel = {e:i for i,e in enumerate(leading)}
    L = np.array([[midpoint(q) for q in row] for row in red.tolist()])
    M = np.zeros((3,20,20),complex)
    for j in range(3):
        for col,u in enumerate(basis):
            e = list(u)
            e[j] += 1
            e = tuple(e)
            if e in lookup:
                M[j,lookup[e],col] = 1
            else:
                M[j,:,col] = L[rel[e]]
    mix = np.array([1,0.371+0.193j,-0.227+0.419j])
    _,vecs = np.linalg.eig(np.einsum('i,ijk->jk',mix,M).T)
    points = []
    for v in vecs.T:
        c = np.array([np.vdot(v,m.T@v)/np.vdot(v,v) for m in M])
        D = np.einsum('j,ijk->ik',np.r_[1.,c],C)
        d = np.linalg.lstsq(D[:,1:],-D[:,0],rcond=None)[0]
        points.append(polish(C,np.r_[c,d]))
    return points


def polish(C, x):
    x = x.copy()
    for _ in range(12):
        c,d = np.r_[1.,x[:3]],np.r_[1.,x[3:]]
        f = np.einsum('a,iab,b->i',c,C,d)
        jac = np.column_stack((np.einsum('iab,b->ia',C[:,1:,:],d),
                               np.einsum('a,iab->ib',c,C[:,:,1:])))
        dx = np.linalg.solve(jac,-f)
        x += dx
        if np.linalg.norm(dx,np.inf) < 5e-15*(1+np.linalg.norm(x,np.inf)):
            break
    return x


def values(C, point):
    c = acb_mat([[acb(1)]+point[:3]])
    d = acb_mat([[acb(1)]+point[3:]]).transpose()
    f,jac = [],[]
    for ci in C:
        left,right = c*ci, ci*d
        f.append([(left*d)[0,0]])
        jac.append([right[j,0] for j in range(1,4)]+[left[0,j] for j in range(1,4)])
    return acb_mat(f), acb_mat(jac)


def certify_root(C, x):
    symmetric = np.linalg.norm(x[3:]-x[:3].conj(),np.inf) < 1e-7
    if symmetric:
        x = np.r_[x[:3],x[:3].conj()]
    p = [exact(v) for v in x]
    xp = acb_mat([[v] for v in p])
    f,j = values(C,p)
    iv = acb_mat([[v.mid() for v in row] for row in j.tolist()]).inv()
    Y = acb_mat([[v.mid() for v in row] for row in iv.tolist()])
    for base in (1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4):
        r = base*max(1,np.linalg.norm(x,np.inf))
        offset = acb(arb(0,r),arb(0,r))
        X = [v+offset for v in p]
        _,jx = values(C,X)
        M = identity(6)-Y*jx
        # Real/imag l-infinity matrix bound of a complex-linear map.
        contraction = max(sum(v.real.abs_upper()+v.imag.abs_upper() for v in row)
                          for row in M.tolist())
        box = acb_mat([[v] for v in X])
        K = xp-Y*f+M*(box-xp)
        if contraction<1 and all(X[i].contains_interior(K[i,0]) for i in range(6)):
            return X,[K[i,0] for i in range(6)],symmetric, float(contraction.upper()),r
    raise ArithmeticError("Krawczyk failed to certify one candidate")


def classify(C, W, V, H, X, K, symmetric, legal):
    # Certify uniqueness of the d lift over the entire projected c root.
    c = [acb(1)] + X[:3]
    vec = acb_mat([c])
    rows = [(vec*ci).tolist()[0][1:] for ci in C]
    mid = np.array([[midpoint(v) for v in row] for row in rows])
    _,_,piv = qr(mid.T,pivoting=True)
    lift = acb_mat([rows[i] for i in piv[:3]])
    assert lift.det().abs_lower()>0
    if any(not X[j].conjugate().overlaps(X[3+j]) for j in range(3)):
        return {'status':'NONCONJUGATE'}
    exact_invariance = all(X[j].conjugate().contains(X[3+j]) and
                           X[3+j].contains(X[j].conjugate()) for j in range(3))
    if not (symmetric and exact_invariance):
        raise ArithmeticError('conjugacy unresolved')
    vec = acb_mat([[acb(1)]+K[:3]])
    w,v = (vec*W).tolist()[0],(vec*V).tolist()[0]
    if not all(q.abs_lower()>0 for q in v):
        raise ArithmeticError('possible dark receiver node')
    z = [(p/q).real for p,q in zip(w,v)]
    lower,upper = legal
    rec = {'material_intervals':[bounds(t) for t in z]}
    if any(t<lower or t>upper for t in z):
        return dict(rec,status='OUTSIDE')
    if not all(t>lower and t<upper for t in z):
        raise ArithmeticError('material boundary unresolved')
    A = acb_mat([[int(i==j)-z[i]*H[i][j] for j in range(6)] for i in range(6)])
    assert A.det().abs_lower()>0
    return dict(rec,status='LEGAL',material_midpoint=[float(t.mid()) for t in z],
                receiver_node_modulus_lower=min(float(q.abs_lower()) for q in v))


def closed_count(C, W, V, H, nominees, bounds_z):
    classified,roots,boxes = [],[],[]
    for x in nominees:
        X,K,sym,norm,r = certify_root(C,x)
        item = classify(C,W,V,H,X,K,sym,bounds_z)
        item.update(krawczyk_radius=r,contraction=norm,
                    candidate_re=x.real.tolist(),candidate_im=x.imag.tolist())
        classified.append(item)
        boxes.append(X[:3])
        if item['status']=='LEGAL':
            roots.append({k:v for k,v in item.items() if k.startswith('material_') or
                          k.startswith('receiver_')})
    distinct = all(any(not a.overlaps(b) for a,b in zip(A,B))
                   for A,B in combinations(boxes,2))
    if len(boxes)!=20 or not distinct:
        raise ArithmeticError("20 mutually disjoint c projections required")
    return classified, roots


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--hc',action='store_true',help='certify separate Julia homotopy nominees')
    args = parser.parse_args()
    with np.load(OUT/'host_operators.npz') as data:
        num = dict(data)
    packed = json.loads((OUT/'host_operators_interval.json').read_text())
    H = {k:from_bounds(packed[k]) for k in ('H','B','R','y0')}
    result = {'scope':'whole regular material box of analytic disk pulse discretization for each exact total observation',
              'material_box':num['material_box'].tolist(),'precision_bits':ctx.prec,'fibers':{}}
    for side in ('minus','plus'):
        obs = num['y_'+side]
        start = time.perf_counter()
        W,V,C,I = chart(H,obs,num['B'])
        cf = floating_coefficients(C)
        np.savetxt(OUT/('hc_coeff_'+side+'.txt'),np.c_[cf.real.ravel(),cf.imag.ravel()],fmt='%.17g')
        ch_end = time.perf_counter()
        basis,leading,red,diag = leading_reduction(C)
        lead_end = time.perf_counter()
        if args.hc:
            witness = np.loadtxt(OUT/('hc_roots_'+side+'.txt'))
            nominees = [polish(cf,x[:6]+1j*x[6:]) for x in witness]
        else:
            nominees = nominate(cf,basis,leading,red)
        gen_end = time.perf_counter()
        items,roots = closed_count(C,W,V,H['H'],nominees,num['material_box'])
        end = time.perf_counter()
        for root in roots:
            m = num['m0'].copy()
            m[num['active']] += root['material_midpoint']
            physical = num['S']@np.linalg.solve(np.eye(len(m))-m[:,None]*num['G'],m[:,None]*num['E0'])
            root['full_forward_midpoint_residual'] = float(np.linalg.norm(physical-obs))
        result['fibers'][side] = dict(leading=diag,pivot=I,complex_count=len(items),
            legal_count=len(roots),classifications=items,physical_roots=roots,
            wall_seconds=dict(chart=ch_end-start,leading=lead_end-ch_end,
                              nominate=gen_end-lead_end,certify=end-gen_end,total=end-start))
    target = 'point_homotopy_certificate.json' if args.hc else 'point_fibers.json'
    (OUT/target).write_text(json.dumps(result,indent=2)+'\n')
    print({s:(v['complex_count'],v['legal_count'],v['wall_seconds'])
           for s,v in result['fibers'].items()})


if __name__=='__main__':
    main()
