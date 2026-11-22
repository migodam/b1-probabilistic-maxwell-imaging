"""Verify simultaneous global single-illumination constraints for all unit source mixtures.

Compressed currents are embedded isometrically in the native cell space. A right
inverse corrects their mean to reproduce the stored ports exactly. The dark ray
and trace guards concern the full uncompressed 96-by-96 operator.
"""
import json
import sys
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from flint import acb, arb, ctx
from scipy.linalg import eigh, eigvalsh, qr, svd
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
sys.path.insert(0, str(STUDY / 'fine'))
sys.path.insert(0, str(STUDY / 'core'))
from all_scalar import compress
from certificate import ball, frobenius_squared, gershgorin, safe_float, star
from fourier_bound import h


def lower(a):
    return safe_float(a, upper=False)


def upper(a):
    return safe_float(a, upper=True)


def eigen_lower(matrix):
    midpoint = np.array([[complex(x) for x in row] for row in matrix.tolist()])
    basis = ball(eigh(h(midpoint), driver='evd')[1])
    gram = star(basis) * basis
    if not gershgorin(gram, upper=False) > 0:
        raise ArithmeticError('Projected sign basis is unresolved')
    return gershgorin(star(basis) * matrix * basis, upper=False) / gershgorin(gram, upper=True)


def main():
    started, cpu = perf_counter(), process_time()
    ctx.prec, ctx.threads = 144, 1
    U, T, eps = .75, 10_000, arb(2)**-40
    dark = json.loads((STUDY/'results/theory/fine_long_ratio.json').read_text())
    assert dark['cap'] == U and dark['ratio'] == T
    assert dark['verified_form_upper'] < 0 and dark['radiation_lower'] > 0
    with np.load(STUDY/'synthetic/public/fine_geometry.npz') as z:
        spec=z['spectrum']; sn=z['PW360_10x20__S']; bn=z['PW360_10x20__B']
    n=sn.shape[1]
    trace_g=sum((acb(float(z.real), float(z.imag)) for z in spec.flat), acb(0))*n/spec.size
    trace_r=trace_g.imag
    trace_combined=n-arb(U)*trace_g.real+T*trace_r
    s_norm=frobenius_squared(ball(sn))
    G,S,B,g,s,b=compress(spec,sn,bn)
    right=star(S)*(S*star(S)).inv()
    # All subsequent online calls see precisely one measured field. The geometry is shared.
    k=s.shape[0]; us,sv,vs=svd(s,full_matrices=False)
    right_np=vs.conj().T @ (us.conj().T/sv[:,None])
    null=qr(s.conj().T,mode='full')[0][:,k:]
    radiation=h(g/1j)
    ev,rotate=eigh(h(null.conj().T@radiation@null),driver='evd')
    null=null@rotate
    o_b=null.conj().T@(-.5j*b); o_g=null.conj().T@radiation
    setup=perf_counter()-started
    directory=HERE.parent/'curve/observations'
    items=json.loads((directory/'manifest.json').read_text())['observations']
    out=[]
    for item in items:
        one=perf_counter()
        with np.load(directory/item['file']) as f:y=f['Y']
        x0=right_np@y; rhs=o_g@x0+o_b
        best=None
        for beta in (1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,3e-5,3e-4):
            x=x0-null@(rhs/(ev[:,None]+beta)); je=x.conj().T@(b+g@x)
            neg_opt=-h(je/1j); neg_total=T*neg_opt-h(x.conj().T@x-U*je)
            margins=[eigvalsh(m,check_finite=False) for m in (neg_opt,neg_total)]
            score=min(float(a[0]/max(1,abs(a[-1]))) for a in margins)
            if best is None or score>best[0]:best=(score,beta,x)
        score,beta,x=best
        xa=ball(x); ya=ball(y); xa=xa+right*(ya-S*xa)
        je=star(xa)*(B+G*xa)
        neg_opt=(star(je)-je)/acb(0,2)
        neg_total=T*neg_opt-star(xa)*xa+arb(U)*(je+star(je))/2
        margin_o=eigen_lower(neg_opt); margin_c=eigen_lower(neg_total)
        # Algebraically SX=Y; this independently encloses the correction and product.
        fit=sum((abs(z).upper()**2 for z in (S*xa-ya).entries()), arb(0)).upper().sqrt(); radius=arb(item['padded_radius'])
        sector_gap=margin_c-eps*abs(trace_combined)
        optical_gap=margin_o-eps*abs(trace_r)
        residual_upper=(fit**2+eps*s_norm).sqrt()
        admissible=bool(optical_gap>0 and sector_gap>0 and residual_upper<radius)
        out.append(dict(id=item['id'],noise_pct=item['noise_pct'],cap=U,beta=beta,
                        strict_neg_optical_eigen_floor=lower(margin_o),strict_neg_combined_eigen_floor=lower(margin_c),
                        interval_mean_residual_upper=upper(fit),exact_fit_by_affine_right_inverse=True,
                        slater_epsilon=float(eps),slater_optical_numerator_lower=lower(optical_gap),
                        slater_sector_slack_lower=lower(sector_gap),slater_residual_upper=upper(residual_upper),
                        padded_radius=item['padded_radius'],verified=admissible,online_wall_seconds=perf_counter()-one))
        print(item['id'],admissible,round(lower(margin_o),5),round(lower(margin_c),3),flush=True)
    result=dict(native_dimension=n,scope='exact saved 96-grid and coherent Euclidean unit source coefficients',
                construction='Affine right-inverse correction in 3-by-3 constant-current subspace, followed by common native dark ray; lifted covariance eps I plus a dark rank-one term',
                cap=U,ratio=T,dark_witness='../../results/theory/fine_long_ratio.json',
                trace_R_upper=upper(trace_r),abs_trace_combined_upper=upper(abs(trace_combined)),
                S_norm_Frobenius_squared_upper=upper(s_norm),precision_bits=ctx.prec,
                all_verified=all(x['verified'] for x in out),observations=out,
                geometry_setup_wall_seconds=setup,wall_seconds=perf_counter()-started,cpu_seconds=process_time()-cpu)
    (HERE/'certificates.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    with threadpool_limits(1):main()
