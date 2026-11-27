"""Exact-rational validation of a nonzero same-data material tube.
The proof conditions and all comparisons are rational arithmetic.
Matrix norms are bounded by exact squared Frobenius/row norms.
"""
import sympy as s
import json,time
from pathlib import Path
Q=s.Rational
n=6
H=s.Matrix(n,n,lambda i,j:[Q(8,100),Q(6,100),Q(25,1000),Q(1,100)][min(abs(i-j),n-abs(i-j))])
B=Q(19,1000)*s.Matrix([[10,2,8],[8,10,1],[2,7,10],[-5,9,7],[6,-4,10],[-7,-5,3]])
Gamma=B*B.T+Q(24,100)*s.eye(n)
G=H+s.I*Gamma
x=s.Matrix([Q(17,10),Q(15,10),Q(20,10),Q(13,10),Q(18,10),Q(16,10)])
X=s.diag(*x)
R=(s.eye(n)-X*G).inv(method='DM')
U=2*s.I*B[:,2:3].T*R
K=G*R
E=(s.eye(n)-G*X).inv(method='DM')*B[:,:2]
# Expand just enough to reduce exact Gaussian rationals.
U=U.applyfunc(s.expand_complex)
K=K.applyfunc(s.expand_complex)
E=E.applyfunc(s.expand_complex)
cols=[]
for i in range(n):
    zz=U[:,i:i+1]*E[i:i+1,:]
    zz=zz.applyfunc(s.expand_complex)
    cols.append(s.Matrix([s.re(zz[0,j]) for j in range(2)]+[s.im(zz[0,j]) for j in range(2)]))
J=s.Matrix.hstack(*cols)
A=J[:,1:]
C=s.zeros(n,4)
C[1:,:]=A.T*(A*A.T).inv(method='DM')
h=s.zeros(n,1);h[0]=1
v=h-C*J*h

def sqnorm(M):
    return sum(s.re(t)**2+s.im(t)**2 for t in M)

kbar=Q(83,100);Dbar=Q(33,100);cbar=Q(86);vbar=Q(61,50)
checks={
    'J_C_identity':J*C==s.eye(4),
    'J_v_zero':J*v==s.zeros(4,1),
    'K_frobenius_squared_le_kbar_squared':bool(sqnorm(K)<=kbar**2),
    'U_E_product_squared_le_Dbar_squared':bool(sqnorm(U)*sqnorm(E)<=Dbar**2),
    'C_row_norms_le_cbar':all(sum(C[i,j]**2 for j in range(4))<=cbar**2 for i in range(n)),
    'v_infinity_le_vbar':all(abs(t)<=vbar for t in v)
}
a=Q(1,200);b=Q(1,50000)
d=a*vbar+b*cbar
q=kbar*d
kap=cbar*Dbar*((1-q)**(-2)-1)
ga=Dbar*kbar*(a*vbar)**2/(1-kbar*a*vbar)
checks.update(
    inverse_neumann_margin=bool(q<1),
    strict_contraction=bool(kap<1),
    ball_maps_strictly_inside=bool(ga+kap*b<b),
    entire_tube_in_material_box=bool(d<min(min(t,Q(9,2)-t) for t in x))
)
if not all(checks.values()): raise AssertionError(checks)
result={
    'rigor_status':'all listed inequalities verified by exact rational arithmetic, not interval floating point',
    'model':'six-cell passive reciprocal dissipative algebraic VIE; not a free-space Maxwell discretization',
    'all_checks_passed':all(checks.values()),'checks':checks,
    'conservative_bounds':{'K_spectral_upper':str(kbar),'D0_upper':str(Dbar),'C_2_to_infinity_upper':str(cbar),'v_infinity_upper':str(vbar)},
    'tube':{'a':str(a),'b':str(b),'d':str(d),'q':str(q),'kappa_float':float(kap),'forcing_bound_float':float(ga),'mapping_upper_float':float(ga+kap*b),'b_float':float(b)},
    'proved_task_width':str(2*a),
    'passivity':'Im G = B B^T + (6/25) I; hence full scattering I+2i B^T(X^-1-G)^-1 B is contractive for x>0',
    'data':'one receiver, two transmitters, two complex entries; four independent real data',
    'interpretation':'all t in [-1/200,1/200] admit unique correction in the chosen four-dimensional exterior chart, exactly preserving both complex readings; x_1 changes by t'
}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
