"""Outward Arb bounds for the exact binary64 discrete TM matrices on a full box."""

from dataclasses import dataclass

from flint import acb, acb_mat, arb, arb_mat, ctx
import numpy as np

ctx.prec = 96


def cball(value):
    return acb(float(value.real), float(value.imag))


def cmatrix(x):
    return acb_mat([[cball(value) for value in row] for row in x])


def eye(size):
    return acb_mat([[int(i == j) for j in range(size)] for i in range(size)])


def abs_matrix(x):
    return arb_mat([[abs(x[i,j]) for j in range(x.ncols())] for i in range(x.nrows())])


def cstack(matrices):
    return acb_mat([[x[i,j] for j in range(x.ncols())]
                    for x in matrices for i in range(x.nrows())])


def astack(matrices):
    return arb_mat([[x[i,j] for j in range(x.ncols())]
                    for x in matrices for i in range(x.nrows())])


def real_jacobian(x):
    return acb_mat([[part(x[i,j]) for j in range(x.ncols())]
                    for i in range(x.nrows()) for part in (lambda a: a.real, lambda a: a.imag)])


def frobenius(x):
    return sum(abs(x[i,j])**2 for i in range(x.nrows()) for j in range(x.ncols())).sqrt()


def np_mid(x):
    return np.array([[complex(x[i,j]) for j in range(x.ncols())] for i in range(x.nrows())])


def upper(x):
    return np.nextafter(float(x.upper()), np.inf)


def lower(x):
    return np.nextafter(float(x.lower()), -np.inf)


def positive_definite(x):
    """Outward complex Cholesky of a known Hermitian family."""
    n = x.ncols()
    factors = [[acb(0) for _ in range(n)] for _ in range(n)]
    for j in range(n):
        pivot = x[j,j].real - sum(abs(factors[j][k])**2 for k in range(j))
        if not pivot > 0:
            return False
        diag = pivot.sqrt()
        for i in range(j+1,n):
            cross = sum(factors[i][k] * factors[j][k].conjugate() for k in range(j))
            factors[i][j] = (x[i,j]-cross)/diag
    return True


def spectral_bounds(x, include_minimum=True):
    numerical = np.linalg.svd(np_mid(x), compute_uv=False)
    gram = x.conjugate().transpose()*x
    ident = eye(gram.ncols())
    hi = float(numerical[0]*(1+1e-8))
    if not positive_definite(ident*(arb(hi)**2)-gram):
        hi = upper(frobenius(x))
    lo = float(numerical[-1]*(1-1e-8)) if include_minimum and x.nrows() >= x.ncols() else 0.
    if lo and not positive_definite(gram - ident*(arb(lo)**2)):
        lo = 0.
    return lo, hi


def jacobian(a,e):
    return acb_mat([[a[r,j]*e[j,s] for j in range(a.ncols())]
                    for s in range(e.ncols()) for r in range(a.nrows())])


def derivative(a,e,k,q):
    return acb_mat([[a[r,q]*k[q,j]*e[j,s]+a[r,j]*k[j,q]*e[q,s]
                    for j in range(a.ncols())]
                    for s in range(e.ncols()) for r in range(a.nrows())])


@dataclass
class Center:
    a: acb_mat
    e: acb_mat
    k: acb_mat
    y: acb_mat
    j: acb_mat
    d_j: list
    norm_k: float
    norm_a: float
    norm_e: float


def centers(blocks, material):
    for block in blocks:
        g = cmatrix(block.green)
        b = cmatrix(block.incident)
        s = cmatrix(block.receiver)
        m = block.host.copy()
        m[block.active] = material
        diag = acb_mat([[float(m[i]) if i==j else 0 for j in range(len(m))]
                        for i in range(len(m))])
        q = (eye(len(m)) - diag*g).inv()
        cur = q*diag*b
        a0, e0, k0, y0 = s*q, b+g*cur, g*q, s*cur
        active = block.active
        a = acb_mat([[a0[i,j] for j in active] for i in range(a0.nrows())])
        e = acb_mat([[e0[i,j] for j in range(e0.ncols())] for i in active])
        k = acb_mat([[k0[i,j] for j in active] for i in active])
        y = acb_mat([[y0[r,c]] for c in range(y0.ncols()) for r in range(y0.nrows())])
        j = jacobian(a,e)
        d_j = [derivative(a,e,k,i) for i in range(len(active))]
        yield Center(a, e, k, y, j, d_j, spectral_bounds(k,False)[1],
                     spectral_bounds(a,False)[1], upper(frobenius(e)))


def horizontal(left, right):
    return acb_mat([[left[i,j] for j in range(left.ncols())] +
                    [right[i,j] for j in range(right.ncols())] for i in range(left.nrows())])


def gain_columns(y, gamma):
    return horizontal(y * gamma, y * (acb(0,1)*gamma))


def single_column(x, col):
    return acb_mat([[x[i,col]] for i in range(x.nrows())])


def augment(left, right, last_row):
    return cstack([horizontal(left,right), acb_mat([last_row])])


def make_joint(y,j,d_j,gain,epsilon,gamma):
    count = j.ncols()
    gain_cols = gain_columns(y,gamma)
    matrix = augment(j*gain,gain_cols,[0]*count+[epsilon,acb(0,epsilon)])
    changes = [
        augment(d_j[i]*gain,gain_columns(single_column(j,i), gamma),[0]*(count+2))
        for i in range(count)
    ]
    zero = acb_mat(gain_cols.nrows(),2)
    changes += [
        augment(j*(gamma*z),zero,[0]*(count+2)) for z in (1,acb(0,1))
    ]
    return matrix, changes


def response_remainders(c, width):
    """Entrywise bounds of first- and second-order material Jacobian changes."""
    a, k, e = abs_matrix(c.a), abs_matrix(c.k), abs_matrix(c.e)
    p = k * arb(width)
    if any(not sum(p[i,j] for j in range(p.ncols())) < 1 for i in range(p.nrows())):
        return None
    ident = arb_mat([[int(i==j) for j in range(p.ncols())] for i in range(p.nrows())])
    u = p*(ident-p).inv()
    v = p*u
    da, de = a*u, u*e
    da2, de2 = a*v, v*e
    first = arb_mat([[da[r,i]*e[i,s]+a[r,i]*de[i,s]+da[r,i]*de[i,s] for i in range(a.ncols())]
                     for s in range(e.ncols()) for r in range(a.nrows())])
    second = arb_mat([[da2[r,i]*e[i,s]+a[r,i]*de2[i,s]+da[r,i]*de[i,s] for i in range(a.ncols())]
                      for s in range(e.ncols()) for r in range(a.nrows())])
    f_second = arb_mat([[sum(a[r,i]*de[i,s] for i in range(a.ncols()))*arb(width)]
                        for s in range(e.ncols()) for r in range(a.nrows())])
    return first, second, f_second


def rigorous_witness(blocks, material, reference, g_m=1., g_n=1., epsilon=1., gamma=0.):
    a = list(centers(blocks, material))
    b = list(centers(blocks, reference))
    h = acb_mat([[arb(float(m))-arb(float(n))] for m,n in zip(material,reference,strict=True)])
    gm, gn = cball(g_m), cball(g_n)
    ya, yb = cstack([c.y for c in a]), cstack([c.y for c in b])
    ja = cstack([c.j for c in a])
    derivative = ((ja*h)*gm + ya*(gm-gn))/arb(epsilon)
    chord = (ya*gm-yb*gn)/arb(epsilon)
    direction = (derivative.conjugate().transpose()*chord)[0,0].real
    h2 = (h.transpose()*h)[0,0].real
    if gamma:
        term = (abs(gm-gn)/arb(gamma))**2
        direction += term
        h2 += term
    return dict(direction_lo=lower(direction/h2), direction_hi=upper(direction/h2),
                chord_over_epsilon_lo=lower(frobenius(chord)))


def directions(blocks, material, gain, epsilon, gamma, observed):
    c = list(centers(blocks,material))
    d = len(material)
    y, j = (cstack([getattr(b,k) for b in c]) for k in ("y","j"))
    h = [cstack([b.d_j[i] for b in c]) for i in range(d)]
    g = cball(gain)
    obs = acb_mat([[cball(value)] for value in observed])
    center_residual = g*y-obs
    if gamma:
        center_residual = cstack([center_residual, acb_mat([[(g-1)*arb(epsilon)/arb(gamma)]])])
    center_score = upper(frobenius(center_residual)/arb(epsilon))
    if gamma:
        gain_extent = upper(1+max(abs(g.real-1),abs(g.imag))/arb(gamma))
        jc, hf = make_joint(y,j,h,g,epsilon,gamma)
    else:
        gain_extent = 0.
        jc, hf = j,h
    ref = real_jacobian(jc)
    sigma, top = spectral_bounds(ref)
    _, s, vt = np.linalg.svd(np_mid(ref).real,full_matrices=False)
    out = dict(r13_width=0., relative_width=0., relative_positive_width=0.,
               gain_half_width= gain_extent, sigma=sigma, L0=top,
               center_score=center_score, certified_arithmetic="Arb 96-bit; exact stored binary64 discrete matrices")
    if ref.nrows() < ref.ncols() or not sigma:
        return out

    w = vt.T / s
    tmat = cmatrix(w)
    abs_t = abs_matrix(tmat)
    reference_floor, _ = spectral_bounds(ref*tmat)
    if not reference_floor:
        return out
    _, map_material = spectral_bounds(cmatrix(w[:d]),False)
    projected_first = np.asarray([upper(frobenius(m*tmat)) for m in hf])
    def evaluate(width):
        radius = arb(width)
        r13, slope = None, None
        if not gamma:
            block_bounds = []
            for b in c:
                p = radius * arb(b.norm_k)
                if not p < 1:
                    break
                block_bounds.append(arb(b.norm_a)*arb(b.norm_e)*p*(2-p)/(1-p)**2)
            else:
                r13 = upper(sum(x**2 for x in block_bounds).sqrt()/arb(sigma))
        bounds = [response_remainders(b,width) for b in c]
        if any(x is None for x in bounds):
            return r13,slope
        se = astack([x[1] for x in bounds])
        terms = radius * sum(arb(v) for v in projected_first[:d])
        if gamma:
            fi, fs = (astack([x[i] for x in bounds]) for i in (0,2))
            d_g = arb(gain_extent)*arb(gamma)*arb(2).sqrt()
            se = se*abs(g) + fi*d_g
            remainder = arb_mat([[se[i,j] for j in range(d)] +
                                  [fs[i,0]*arb(gamma)]*2 for i in range(se.nrows())]+[[0]*(d+2)])
            terms += arb(gain_extent)*sum(arb(v) for v in projected_first[d:])
        else:
            remainder = se
        slope = upper((terms + frobenius(remainder*abs_t))/arb(reference_floor))
        return r13,slope

    def width_at(which, target):
        start = evaluate(0.)[which]
        if start is None or start >= target:
            return 0.
        lo, hi = 0., (30.-.5)/2
        for _ in range(21):
            trial = (lo+hi)/2
            value = evaluate(trial)[which]
            if value is not None and value < target:
                lo = trial
            else:
                hi = trial
        return lo

    out["r13_width"] = width_at(0, .7071067811865475)
    out["relative_positive_width"] = width_at(1, .7071067811865475)
    out["relative_width"] = width_at(1, .499999999)
    out["relative_eta_at_zero_material_width"] = evaluate(0.)[1]
    out["r13_available"] = not bool(gamma)
    out["r13_eta_at_0.1"] = evaluate(.1)[0]
    out["relative_eta_at_0.1"] = evaluate(.1)[1]
    out["domain_checks"] = [
        dict(material_half_width=r, r13_eta=x, relative_eta=v)
        for r in (.001,.01,.03,.1,.3,1.,3.,14.75) for x,v in [evaluate(r)]
    ]
    out["reference_floor"] = reference_floor
    out["material_transform_norm"] = map_material
    if not out["relative_width"]:
        return out
    width = out["relative_width"]
    eta = arb(evaluate(width)[1])
    base = arb(epsilon)*arb(map_material) / (arb(reference_floor)*arb(d).sqrt())
    factor = (1+eta)/(1-eta)**2
    out.update(relative_eta=upper(eta),
               joint_stationary_truth_rms_bound=upper(base*factor),
               relative_material_diameter_rms_bound=upper(2*base/(1-eta)),
               stationary_fit_rms_bound=upper(base*factor*arb(center_score)),
               can_exclude_boundary_from_feasible_component=bool(
                   (1-eta)*arb(reference_floor)*arb(width)/arb(map_material)
                   > arb(epsilon)*(1+arb(center_score))),
               normalization="raw field plus epsilon*(gain-1)/gamma when shared gain is fitted")
    return out


def local_minimum(blocks, material, observed, epsilon, truth, box):
    """Prove existence of an interior strict local minimum on a small Euclidean ball."""
    c = list(centers(blocks,material))
    f = cstack([b.y for b in c])
    jac = cstack([b.j for b in c])
    residual = f-acb_mat([[cball(y)] for y in observed])
    adj = jac.conjugate().transpose()
    grad = (adj*residual)
    grad_norm = frobenius(acb_mat([[grad[i,0].real] for i in range(grad.nrows())]))
    gram = adj*jac
    d = jac.ncols()
    mixed = cstack([cstack([b.d_j[i] for b in c]).transpose()*residual.conjugate()
                    for i in range(d)])
    hh = acb_mat([[gram[i,j].real+mixed[i*d+j,0].real for j in range(d)] for i in range(d)])
    eig = float(np.linalg.eigvalsh(np_mid(hh).real)[0])
    lam = arb(eig*.999999)
    if not lam > 0 or not positive_definite(hh-eye(d)*lam):
        return None
    _, j_norm = spectral_bounds(real_jacobian(jac),False)
    true_dist = frobenius(acb_mat([[arb(float(a))-arb(float(b))]
                                  for a,b in zip(material,truth,strict=True)]))
    residual_norm = frobenius(residual)
    for width in (1e-6, 3e-7, 1e-7, 3e-8, 1e-8, 3e-9, 1e-9):
        if np.min(np.minimum(material-box[0],box[1]-material)) <= width:
            continue
        radius = arb(width)
        q = [radius*arb(b.norm_k) for b in c]
        if any(not x < 1 for x in q):
            continue
        a_e = [arb(b.norm_a)*arb(b.norm_e) for b in c]
        beta = [2*v*arb(b.norm_k) for v,b in zip(a_e,c,strict=True)]
        j_diff = sum((v*((1-x)**-2-1))**2 for v,x in zip(a_e,q,strict=True)).sqrt()
        beta_max = sum((v/(1-x)**3)**2 for v,x in zip(beta,q,strict=True)).sqrt()
        beta_diff = sum((v*((1-x)**-3-1))**2 for v,x in zip(beta,q,strict=True)).sqrt()
        residual_diff = (arb(j_norm)+j_diff)*radius
        hess_diff = j_diff*(2*arb(j_norm)+j_diff) + residual_diff*beta_max + residual_norm*beta_diff
        floor = lam-hess_diff
        if grad_norm < floor*radius:
            return dict(verified=True, ball_radius=width, hessian_lower=lower(floor),
                        center_hessian_lower=lower(lam), gradient_upper=upper(grad_norm),
                        material_RMS_lower=lower((true_dist-radius)/arb(d).sqrt()),
                        residual_over_epsilon_lower=lower((residual_norm-residual_diff)/arb(epsilon)),
                        center=material.tolist())
    return None
