"""Synthetic finite-dimensional algebra checks, NOT electromagnetic experiments.
Run: python checks/check_finite_source_transport.py
Requires numpy and scipy. All dense solves below are synthetic matrix operations.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from numpy.linalg import norm, solve, svd, slogdet
from scipy.linalg import null_space, eigh, subspace_angles

RNG = np.random.default_rng(20261017)

def cm(shape: tuple[int, ...]) -> np.ndarray:
    return RNG.normal(size=shape) + 1j * RNG.normal(size=shape)

def rel(a: np.ndarray, b: np.ndarray) -> float:
    return float(norm(a-b) / max(norm(a), norm(b), 1e-14))

def op(a: np.ndarray) -> float:
    return float(norm(a, 2))

def qop(d: np.ndarray, h: np.ndarray) -> np.ndarray:
    return solve(np.eye(len(d)) - np.diag(d) @ h, np.diag(d))

def galerkin(l: np.ndarray, v: np.ndarray) -> np.ndarray:
    return v @ solve(v.conj().T @ l @ v, v.conj().T)

def real_matrix(a: np.ndarray) -> np.ndarray:
    return np.block([[a.real, -a.imag], [a.imag, a.real]])

n, s, p, m = 18, 7, 3, 5
x0 = .4 + .4j + .04 * cm((n,))
g = cm((n,n)); g *= .45 / op(g)
l0 = np.eye(n) - np.diag(x0) @ g
r0 = solve(l0, np.eye(n))
e = cm((n,p)); j0 = solve(l0, x0[:,None]*e)
e0 = e + g @ j0
u = np.eye(n)[:, :s]
sensor = cm((m,n)); sensor /= op(sensor)
pfull = r0 @ u
hfull, ofull, b = u.conj().T @ g @ pfull, sensor @ pfull, u.conj().T @ e0
basis, _ = np.linalg.qr(cm((n,12)))
v, c = basis[:,:7], basis[:,7:11]
rv = galerkin(l0,v); rvnew = galerkin(l0,np.column_stack([v,c]))
pv, pn = rv @ u, rvnew @ u
barc = (np.eye(n)-rv@l0)@c
gamma = c.conj().T @ l0 @ barc
nu = c.conj().T @ (u-l0@pv)
ep = barc @ solve(gamma,nu)
hv, hn = u.conj().T@g@pv, u.conj().T@g@pn
ov, on = sensor@pv, sensor@pn
results = {}
results['background_port_exchange_relative_error'] = rel(pn-pv,ep)

d = .008 * cm((s,))
dx = u @ np.diag(d) @ u.conj().T
lh = l0-dx@g
j1 = solve(lh,(np.diag(x0)+dx)@e)
full_delta = sensor@(j1-j0)
qf, qv, qn = qop(d,hfull), qop(d,hv), qop(d,hn)
results['full_finite_port_response_relative_error'] = rel(full_delta, ofull@qf@b)
fold, fnew = ov@qv@b, on@qn@b
finite_exchange = (sensor + ov@qv@u.conj().T@g)@barc@solve(gamma,nu@qn@b)
results['doubly_dressed_finite_exchange_relative_error'] = rel(fnew-fold,finite_exchange)
kvh = galerkin(lh,v)@dx@e0
knh = galerkin(lh,np.column_stack([v,c]))@dx@e0
results['reduced_port_vs_updated_state_galerkin_relative_error'] = rel(sensor@kvh,fold)
barh = (np.eye(n)-galerkin(lh,v)@lh)@c
gh = c.conj().T@lh@barh
nh = c.conj().T@(dx@e0-lh@kvh)
results['direct_finite_state_schur_relative_error'] = rel(knh-kvh,barh@solve(gh,nh))
rho_u = u-l0@pv
rho_h = dx@e0-lh@(pv@qv@b)
results['finite_state_residual_factorization_relative_error'] = rel(rho_h,rho_u@qv@b)
alpha0 = float(svd(l0,compute_uv=False)[-1])
ep_bound = op(rho_u)/alpha0
eps_o, eps_h = op(sensor)*ep_bound, op(u.conj().T@g)*ep_bound
theta = op(qv)*eps_h
if not theta < 1:
    raise AssertionError('Chosen synthetic certificate domain is invalid.')
cert = (eps_o+op(ov@qv)*eps_h)*norm(qv@b)/(1-theta)
actual = float(norm(full_delta-fold))
results['finite_certificate'] = {'theta':theta,'true_error':actual,'bound':float(cert),'covered':bool(actual<=cert*(1+1e-10))}
if not actual <= cert*(1+1e-10):
    raise AssertionError('The finite-response bound was violated.')

# Exact change of anchor on the same physical ports.
p1 = solve(lh,u)
p1_formula = pfull @ solve(np.eye(s)-np.diag(d)@hfull,np.eye(s))
b1 = u.conj().T@(e+g@j1)
b1_formula = solve(np.eye(s)-hfull@np.diag(d),b)
results['anchor_port_transition_relative_error'] = rel(p1,p1_formula)
results['anchor_field_transition_relative_error'] = rel(b1,b1_formula)
rback = p1_formula@solve(np.eye(s)+np.diag(d)@(u.conj().T@g@p1_formula),np.eye(s))
results['anchor_roundtrip_relative_error'] = rel(rback,pfull)

# Single scalar source: the dark equivalent source is an exact rational chart.
z = null_space(ofull)
t = z@(.01*cm((z.shape[1],)))
bref = b[:,0]; et = bref+hfull@t
dchart = t/et
actual_t = qop(dchart,hfull)@bref
chart_dx = u@np.diag(dchart)@u.conj().T
jchart = solve(l0-chart_dx@g,(np.diag(x0)+chart_dx)@e)
results['scalar_source_chart_inverse_relative_error'] = rel(actual_t,t)
results['single_source_dark_output_absolute_error'] = float(norm(sensor@(jchart[:,0]-j0[:,0])))
results['other_sources_output_change'] = float(norm(sensor@(jchart[:,1:]-j0[:,1:])))
results['single_source_chart_min_total_imaginary_material'] = float(np.min((x0+u@dchart).imag))
dhdq = np.diag(1/et)@(np.eye(s)-np.diag(dchart)@hfull)
eps = 1e-6
fd = np.zeros((2*s,2*s))
for k in range(2*s):
    dt = np.zeros(s,dtype=complex); dt[k%s] = eps*(1 if k<s else 1j)
    dp = (t+dt)/(bref+hfull@(t+dt))
    dm = (t-dt)/(bref+hfull@(t-dt))
    diff = (dp-dm)/(2*eps)
    fd[:,k] = np.r_[diff.real,diff.imag]
results['source_chart_real_jacobian_relative_error'] = rel(fd,real_matrix(dhdq))
_, ld = slogdet(real_matrix(dhdq))
_, lc = slogdet(np.eye(s)-np.diag(dchart)@hfull)
ld_formula = 2*lc-2*np.log(np.abs(et)).sum()
results['source_chart_log_volume_error'] = float(abs(ld-ld_formula))

# Several independent sources: a deliberately constructed row-aligned algebraic control.
# H and O are abstract port matrices, not a realizable antenna acquisition claim.
bs = cm((s,2)); a = .25+.02*cm((s,)); gam = .04*cm((s,))
t1 = a[:,None]*bs
pseud = np.linalg.pinv(t1)
h_align = (gam[:,None]*bs)@pseud + .015*cm((s,s))@(np.eye(s)-t1@pseud)
o_align = null_space(t1.conj().T).conj().T
steps = [.05,.2,.5]
errs=[]; darks=[]
for tau in steps:
    dt = tau*a/(1+tau*gam)
    tt = qop(dt,h_align)@bs
    errs.append(rel(tt,tau*t1)); darks.append(float(norm(o_align@tt)))
results['multisource_aligned_path_relative_error_max'] = max(errs)
results['multisource_aligned_path_dark_error_max'] = max(darks)
results['multisource_control_source_rank'] = int(np.linalg.matrix_rank(bs))
# Row-misalignment gives an exact, not asymptotic, output identity.
h_pert = h_align + .04*cm((s,s))
ht = h_pert@t1
row_energy = np.sum(np.abs(bs)**2,axis=1)
gam_pert = np.sum(ht*bs.conj(),axis=1)/row_energy
res = ht-gam_pert[:,None]*bs
tau=.35; dt=tau*a/(1+tau*gam_pert)
true_response = o_align@qop(dt,h_pert)@bs
pred_response = tau*o_align@t1+tau*o_align@qop(dt,h_pert)@res
results['multisource_feedback_misalignment_identity_relative_error'] = rel(true_response,pred_response)
results['multisource_feedback_misalignment_norm'] = float(norm(res))

# The feedback-compatible injection space is a joint kernel of linear maps.
obs_cols = []
mis_cols = []
gamma_cols = []
for index in range(s):
    aa = np.eye(s,dtype=complex)[:,index]
    tt = aa[:,None]*bs
    hh = h_align@tt
    gg = np.sum(hh*bs.conj(),axis=1)/row_energy
    rr = hh-gg[:,None]*bs
    obs_cols.append((o_align@tt).reshape(-1))
    mis_cols.append(rr.reshape(-1))
    gamma_cols.append(gg)
joint = np.vstack([np.column_stack(obs_cols),np.column_stack(mis_cols)])
active = null_space(joint,rcond=1e-11)
gamma_map = np.column_stack(gamma_cols)
results['feedback_closed_dark_injection_dimension_complex'] = int(active.shape[1])
results['aligned_direction_joint_kernel_relative_error'] = float(norm(joint@a)/max(norm(joint)*norm(a),1e-14))
if active.shape[1] == 0:
    raise AssertionError('Constructed multi-source control must have a nonzero joint kernel.')
az = active@(.03*cm((active.shape[1],)))
dz = az/(1+gamma_map@az)
tz = qop(dz,h_align)@bs
results['rational_joint_kernel_chart_relative_error'] = rel(tz,az[:,None]*bs)
results['rational_joint_kernel_chart_dark_absolute_error'] = float(norm(o_align@tz))
# Fixed-phase material priors are not automatically preserved by the rational map.
phase_test = .2/(1+.2j)
results['fixed_phase_prior_counterexample_imaginary_increment'] = float(phase_test.imag)

# General multi-source realizability and its full-state residual certificate.
tp = .01*cm((s,p)); field = b+hfull@tp
rec_d = np.sum(tp*field.conj(),axis=1)/np.sum(np.abs(field)**2,axis=1)
rho = tp-rec_d[:,None]*field
true_tp = qop(rec_d,hfull)@b
results['realizability_residual_to_full_response_relative_error'] = rel(ofull@(true_tp-tp),-ofull@solve(np.eye(s)-np.diag(rec_d)@hfull,rho))

# Recover the old generalized-eigenvalue / principal-angle diagnostic.
k = cm((10,3)); ss=cm((2,10)); dark=null_space(ss)
prow=np.eye(10)-dark@dark.conj().T
lam=eigh(k.conj().T@prow@k,k.conj().T@k,eigvals_only=True)
angles=subspace_angles(k,dark)
results['generalized_sin_squared_angle_error'] = float(norm(np.sort(lam)-np.sort(np.sin(angles)**2)))

# Counterexamples: no gap, no unseen-defect certificate, and no Neumann necessity.
e1=np.eye(3)[:,:1]; bb=np.array([[0.],[1.],[0.]])
sense=np.array([[0.,1.,0.]])
jrom=sense@galerkin(np.eye(3),e1)@bb
jfull=sense@bb
results['false_weak_counterexample']={'rom_response':float(jrom.item()),'full_response':float(jfull.item()),'both_zero_residual_gn_steps':0.0}
results['stable_full_singular_galerkin_counterexample']={'full_sigma_min':1.0,'galerkin_scalar':0.0}
results['neumann_not_necessary_counterexample']={'feedback_K':2.0,'resolvent':-1.0,'neumann_converges':False}

# Acceptance criteria apply only to these algebra checks.
for key,val in results.items():
    if 'relative_error' in key and isinstance(val,float):
        threshold=1e-7 if 'jacobian' in key else 1e-10
        if val>threshold:
            raise AssertionError(f'{key}: {val} > {threshold}')
for key in ['source_chart_log_volume_error','generalized_sin_squared_angle_error','single_source_dark_output_absolute_error','multisource_aligned_path_dark_error_max','rational_joint_kernel_chart_dark_absolute_error']:
    if results[key]>1e-9:
        raise AssertionError(f'{key}: {results[key]}')
results['scope']={'seed':20261017,'physical_Maxwell_experiments':0,'GPU_training':0,'P1_read_or_assumed':False,'claim':'synthetic finite-dimensional algebra only'}
out=Path(__file__).with_name('algebra_results.json')
out.write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
