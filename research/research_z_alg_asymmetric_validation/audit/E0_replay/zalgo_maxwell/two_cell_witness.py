"""Two-cell scalar TM VIE: exact reflection symmetry, with actual Hankel kernels.

This freezes a two-basis collocation discretization, not a claim that a continuum
Maxwell field has only two inverse solutions. Symmetry itself holds in the
continuum for reflected cells and individually invariant sources/receivers.
"""
import json
from pathlib import Path
import numpy as np
from scipy.special import hankel1, jv
from scipy.integrate import quad

k = 2*np.pi
a = 0.04
centers = np.array([[-0.2, 0.0], [0.2, 0.0]])
sources = np.array([[0.0, -0.7], [0.0, -1.1]])
receivers = np.array([[0.0, 0.8], [0.0, 1.2]])

def external_cell_integral(distance):
    # k^2 integral over a circular constant basis. Valid outside the disk.
    return 0.5j*np.pi*k*a*jv(1,k*a)*hankel1(0,k*distance)

g0 = 0.5j*np.pi*k*a*hankel1(1,k*a)-1
g = external_cell_integral(0.4)
G = np.array([[g0,g],[g,g0]])

def complex_quad(f, lower, upper):
    return quad(lambda t: f(t).real,lower,upper,epsabs=1e-12,epsrel=1e-12)[0]+1j*quad(lambda t: f(t).imag,lower,upper,epsabs=1e-12,epsrel=1e-12)[0]

# Independent integration of the original disk kernels, including the self
# logarithmic singularity (integrable after the polar Jacobian factor rho).
self_quad=complex_quad(lambda rho: 0.5j*np.pi*k*k*rho*hankel1(0,k*rho),0,a)
cross_quad=complex_quad(lambda rho: complex_quad(lambda theta: 0.25j*k*k*rho*hankel1(0,k*np.sqrt(0.4**2+rho**2-2*0.4*rho*np.cos(theta))),0,2*np.pi),0,a)
E0 = 0.25j*hankel1(0,k*np.linalg.norm(centers[:,None,:]-sources[None,:,:],axis=2))
S = external_cell_integral(np.linalg.norm(receivers[:,None,:]-centers[None,:,:],axis=2))

def forward(x, incident=E0):
    A = np.eye(2)-np.diag(x)@G
    J = np.linalg.solve(A,np.diag(x)@incident)
    return S@J, J

x = np.array([1.0,3.0])
y, J = forward(x)
y_swap, J_swap = forward(x[::-1])
Y = y[0,0]/(S[0,0]*E0[0,0])
coeff = np.array([1+Y*g0, 2*(g-g0)-Y*(g0*g0-g*g)])
M = np.stack([coeff.real,coeff.imag])
pq = np.linalg.solve(M,np.array([Y.real,Y.imag]))
roots = np.roots([1,-pq[0],pq[1]])

step = 1e-5
jac=np.column_stack([(forward(x+step*np.eye(2)[i])[0][0,0]-forward(x-step*np.eye(2)[i])[0][0,0])/(2*step) for i in range(2)]).ravel()
jreal=np.stack([jac.real,jac.imag])
analytic_jac=S[0]@np.linalg.solve(np.eye(2)-np.diag(x)@G,np.diag((E0+G@J)[:,0]))
analytic_jreal=np.stack([analytic_jac.real,analytic_jac.imag])

# Break the symmetry by moving only one receiver; fixed source/receiver labels.
break_receiver=np.array([0.12,0.8])
Sbreak=external_cell_integral(np.linalg.norm(break_receiver[None,:]-centers,axis=1))[None,:]
ybreak=Sbreak@J
ybreakswap=Sbreak@J_swap

# Data-preserving dark current perturbation, but shared-material inconsistency.
E0two=E0.copy()
offaxis_source=np.array([0.1,-0.9])
E0two[:,1]=0.25j*hankel1(0,k*np.linalg.norm(centers-offaxis_source[None,:],axis=1))
ydark,Jdark=forward(x,E0two)
Jcandidate=Jdark.copy()
dark_coefficients=np.array([0.01,0.015j])
Jcandidate+=np.array([1,-1])[:,None]*dark_coefficients[None,:]
Ecandidate=E0two+G@Jcandidate
ratios=Jcandidate/Ecandidate
xhat=np.sum(np.conj(Ecandidate)*Jcandidate,axis=1)/np.sum(abs(Ecandidate)**2,axis=1)
R=Jcandidate-xhat[:,None]*Ecandidate
minors=Jcandidate[:,0]*Ecandidate[:,1]-Jcandidate[:,1]*Ecandidate[:,0]
lagrange_lhs=np.sum(abs(R)**2,axis=1)
lagrange_rhs=abs(minors)**2/np.sum(abs(Ecandidate)**2,axis=1)
yhat,Jhat=forward(xhat,E0two)
resolvent=np.linalg.inv(np.eye(2)-np.diag(xhat)@G)
bound=np.linalg.norm(S@resolvent,2)*np.linalg.norm(R,'fro')

def encode(z):
    z=np.asarray(z)
    if np.iscomplexobj(z):
        return {'real':z.real.tolist(),'imag':z.imag.tolist()}
    return z.tolist()

out={
 'contract':'2D scalar TM exp(-i omega t), k=2*pi, two circular pulse bases, exact disk-integrated Green at collocation centers; real contrast domain [0.5,4]^2',
 'radius':a,'centers':centers.tolist(),'sources':sources.tolist(),'receivers':receivers.tolist(),
 'g0':encode(g0),'g':encode(g),'material':x.tolist(),'swapped_material':x[::-1].tolist(),
 'self_disk_quadrature_absolute_error':float(abs(self_quad-g0)),
 'cross_disk_quadrature_absolute_error':float(abs(cross_quad-g)),
 'material_distance':float(np.linalg.norm(x-x[::-1])),
 'symmetry_prediction_relative_error':float(np.linalg.norm(y-y_swap)/np.linalg.norm(y)),
 'normalized_response_Y':encode(Y),'symmetric_coordinates_pq':pq.tolist(),
 'quotient_linear_matrix_det':float(np.linalg.det(M)),'recovered_roots':roots.tolist(),
 'real_jacobian_singular_values':np.linalg.svd(analytic_jreal,compute_uv=False).tolist(),
 'analytic_jacobian_finite_difference_relative_error':float(np.linalg.norm(analytic_jreal-jreal)/np.linalg.norm(analytic_jreal)),
 'resolvent_sigma_min':float(np.linalg.svd(np.eye(2)-np.diag(x)@G,compute_uv=False)[-1]),
 'offaxis_receiver_relative_difference':float(np.linalg.norm(ybreak-ybreakswap)/np.linalg.norm(ybreak)),
 'dark_current_data_relative_error':float(np.linalg.norm(S@Jcandidate-ydark)/np.linalg.norm(ydark)),
 'dark_current_source_positions':[[0.,-0.7],offaxis_source.tolist()],
 'dark_current_nullspace_complex_dimension_each_source':1,
 'dark_current_null_coefficients':encode(dark_coefficients),
 'dark_current_per_source_material_ratios':encode(ratios),
 'dark_current_material_ratio_difference_by_cell':encode(ratios[:,0]-ratios[:,1]),
 'dark_current_shared_material_residual':float(np.linalg.norm(R)),
 'dark_current_normalized_minors_by_cell':np.sqrt(lagrange_rhs).tolist(),
 'lagrange_identity_error':float(np.max(abs(lagrange_lhs-lagrange_rhs))),
 'reconstructed_material_xhat':encode(xhat),
 'projected_material_in_original_real_domain':bool(np.all(abs(xhat.imag)<1e-12) and np.all(xhat.real>=0.5) and np.all(xhat.real<=4)),
 'physical_data_defect_after_xhat':float(np.linalg.norm(yhat-ydark)),
 'physical_data_defect_upper_bound':float(bound),
 'cost_ledger':{'forward_solves_two_complex_rhs':8,'analytic_sensitivity_solves_two_complex_rhs':1,'resolvent_inverse_for_error_certificate':1,'quotient_real_2x2_solves':1,'independent_green_disk_quadratures':2,'optimization_steps':0,'training_steps':0},
 'scope':'Certified symbolically by swap symmetry and the invertible 2x2 real quotient system for this frozen two-cell material family. Numeric values alone do not certify the continuum or wider materials.'
}
Path(__file__).with_name('two_cell_witness.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
