"""Exact algebraic ambiguity in a 3-cell full-wave discrete TM model.

Richmond disk-pulse entries are evaluated in binary64 and then treated as exact
stored rational complex coefficients, matching the manuscript's stored-operator
contract. The ambiguity proof uses exact rational resultants and sign brackets;
it is not a claim of interval enclosure of the transcendental Hankel entries.
"""
import json
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.special import hankel1, j1

out = Path(__file__).parent
k, h = 2*np.pi, .10
zz = k*h/np.sqrt(np.pi)
pos = np.array([[0.,0.], [.25,0.], [0.,.25]])
c = 1j*np.pi*zz*j1(zz)/2
Gf = np.array([[1j*np.pi*zz*hankel1(1,zz)/2-1 if i == j
                else c*hankel1(0,k*np.linalg.norm(pos[i]-pos[j]))
                for j in range(3)] for i in range(3)])
Sf = c*hankel1(0,k*np.linalg.norm(pos-np.array([0.,-1.5]),axis=1))[None,:]
# Exact upward plane-wave phases at the stated cell coordinates.
Bf = np.exp(1j*k*pos[:,1])[:,None]

def cq(z):
    return sp.Rational(float(np.real(z))) + sp.I*sp.Rational(float(np.imag(z)))
G = sp.Matrix(3, 3, [cq(v) for v in Gf.ravel()])
S = sp.Matrix(1, 3, [cq(v) for v in Sf.ravel()])
B = sp.Matrix(3, 1, [cq(v) for v in Bf.ravel()])
x0 = [sp.Rational(7,10),sp.Rational(1,16),sp.Rational(1,16)]
X0 = sp.diag(*x0)
Y0 = sp.cancel((S*(sp.eye(3)-X0*G).inv()*X0*B)[0])

z,w = sp.symbols('z w', real=True)
X = sp.diag(sp.Rational(37,50),z,w)
A = sp.eye(3)-X*G
D = sp.expand(A.det())
N = 0
rhs = X*B
for j in range(3):
    Aj = A.copy()
    Aj[:,j] = rhs
    N += S[j]*Aj.det()
H = sp.expand(N-Y0*D)
hr,hi = map(sp.expand,sp.expand_complex(H).as_real_imag())
# Both equations are affine in w, so solve the real one and substitute.
ar,br = sp.expand(hr.subs(w,0)),sp.diff(hr,w)
ai,bi = sp.expand(hi.subs(w,0)),sp.diff(hi,w)
poly = sp.Poly(sp.expand(ai*br-bi*ar),z)
assert poly.degree() == 2
intervals = sp.polys.polytools.intervals(poly, eps=sp.Rational(1,10**24))

def affine_interval(expr, lo, hi):
    vals = [sp.cancel(expr.subs(z,lo)),sp.cancel(expr.subs(z,hi))]
    return min(vals),max(vals)

proofs=[]
for (lo,hi), multiplicity in intervals:
    if not (0 < lo < hi < sp.Rational(1,8)):
        continue
    aa = affine_interval(-ar,lo,hi)
    bb = affine_interval(br,lo,hi)
    if bb[0] <= 0 <= bb[1]:
        continue
    quotients = [u/v for u in aa for v in bb]
    wl,wh = min(quotients),max(quotients)
    if not (0 < wl < wh < sp.Rational(1,8)):
        continue
    # Exact sign change of the rational polynomial certifies a root.
    assert poly.eval(lo)*poly.eval(hi) < 0
    # Rational majorant on the infinity norm: |a+ib| <= |a|+|b|.
    upper_x = [sp.Rational(37,50),hi,wh]
    bound = max(upper_x[i]*sum(abs(sp.re(G[i,j]))+abs(sp.im(G[i,j]))
                              for j in range(3)) for i in range(3))
    assert bound < 1
    proofs.append({
        'background1_interval':[str(lo),str(hi)],
        'background2_interval':[str(wl),str(wh)],
        'background1_numeric':[float(lo),float(hi)],
        'background2_numeric':[float(wl),float(wh)],
        'polynomial_sign_change':True,
        'forward_infinity_norm_upper':float(bound),
        'forward_infinity_norm_upper_rational':str(bound),
        'all_background_cells_in':[0,0.125],
        'multiplicity':multiplicity,
    })
assert proofs
result = {
    'status':'EXACT_STORED_OPERATOR_ALGEBRAIC_WITNESS',
    'operator_contract':'Exact rational interpretation of stored binary64 G, B, S; no transcendental-entry enclosure claimed.',
    'cell_side_wavelengths':h,
    'positions_wavelengths':pos.tolist(),
    'receiver_position_wavelengths':[0,-1.5],
    'incident_direction':[0,1],
    'baseline_material':[float(v) for v in x0],
    'alternative_local_contrast':0.74,
    'baseline_observation':[float(sp.re(Y0)),float(sp.im(Y0))],
    'proofs':proofs,
    'proof_logic':'Real resultant has a rational sign-changing bracket. Solving the real residual equation gives a real w in (0,1/8); resultant makes imaginary residual zero. Uniform ||XG||_infinity<1 proves denominators nonzero. Therefore the two materials have exactly the same full-wave complex datum for the stored discrete operators.',
}
(out/'background_compensation_exact.json').write_text(json.dumps(result,indent=2))
np.savez(out/'background_compensation_operators.npz',G=Gf,S=Sf,B=Bf)
print(json.dumps({key:value for key,value in result.items() if key != 'proofs'},indent=2))
for p in proofs:
    print({k:v for k,v in p.items() if 'rational' not in k and not k.endswith('_interval')})
