"""Original-prior decoder and real-Hilbert/GN QA; zero Maxwell, zero GT."""
import importlib.util
import json
import sys
from pathlib import Path
import numpy as np
from scipy.special import ndtr

P = Path(__file__).resolve().parents[2]
REPO = P.parents[1]
R4 = REPO.parent/'b1-p3l-r4-conditional-transport'
sys.dont_write_bytecode = True
sys.path.insert(0, str(R4/'src'))
from p1_material.prior import MaterialPrior
spec = importlib.util.spec_from_file_location('r7_math_qa_jets', P/'src/material_jets.py')
jets = importlib.util.module_from_spec(spec)
spec.loader.exec_module(jets)
cfg = json.loads((R4/'research/p3l_r4_conditional_transport/configs/TARGET.json').read_text())
points = json.loads((R4/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
prior = MaterialPrior(cfg, points)
rng = np.random.default_rng(2026102817)
report = {'scope': 'original prior only; no Maxwell import/solve; no GT', 'checks': {}, 'issues': []}

# Controlled regular points: active mixed material, all-high and all-low plateau.
coefficient_cases = []
for label, offset in [('regular', .1), ('upper_plateau', 1.7), ('lower_plateau', -1.7)]:
    x = np.zeros(66); x[0] = offset; x[64:] = [.4, -.7]
    if label == 'regular': x[1:64] = .008*rng.normal(size=63)
    v = .12*rng.normal(size=66)
    x2 = .05*rng.normal(size=66)
    x3 = .03*rng.normal(size=66)
    coefficients = jets.decoder_coefficients(prior, x, v, x2, x3)
    errors = []
    for h in [.08, .04, .02, .01]:
        actual = jets.decode(prior, jets.polynomial(x, v, h, x2, x3))[0]
        approximation = sum(h**k*c for k,c in enumerate(coefficients))
        errors.append(float(np.linalg.norm(actual-approximation)/24))
    orders = [float(np.log2(max(errors[k],1e-30)/max(errors[k+1],1e-30))) for k in range(3)]
    # Error scaling is the direct exact-decoder check, not a cancellation-prone
    # third finite difference divided by h^3 at machine precision.
    if max(errors)>1e-14:
        assert min(orders[:2]) > 3.7, (label, errors, orders)
    h = 1e-5
    first_fd = (jets.decode(prior, x+h*v)[0]-jets.decode(prior, x-h*v)[0])/(2*h)
    first_error = float(np.linalg.norm(first_fd-coefficients[1])/max(np.linalg.norm(coefficients[1]), 1e-30))
    assert first_error < 1e-7
    coefficient_cases.append(dict(case=label, steps=[.08,.04,.02,.01], RMS_errors=errors,
                                  convergence_orders=orders, first_derivative_relative=first_error))
report['checks']['decoder_coefficients'] = coefficient_cases

# A tiny-step independent degree-six complex polynomial fit recovers all Taylor
# coefficients 0..3. Scale t so the fit itself is well-conditioned.
x=np.zeros(66); x[0]=.08; x[64:]=[.9,-1.1]
v=.2*rng.normal(size=66); x2=.1*rng.normal(size=66); x3=.07*rng.normal(size=66)
exact=jets.decoder_coefficients(prior,x,v,x2,x3)
nodes=np.linspace(-1,1,17); h=.015
values=np.stack([jets.decode(prior,jets.polynomial(x,v,h*t,x2,x3))[0] for t in nodes])
fit=np.polynomial.polynomial.polyfit(nodes,values,6)
recovered=[fit[k]/h**k for k in range(4)]
errors=[float(np.linalg.norm(recovered[k]-exact[k])/max(np.linalg.norm(exact[k]),1e-30)) for k in range(4)]
assert max(errors) < 5e-6, errors
report['checks']['independent_tiny_polynomial_fit']=dict(step=h,relative_errors_orders0to3=errors)

# The material norm is a REAL Hilbert norm. No conjugation is omitted:
# [Re(a);Im(a)]^T[Re(b);Im(b)]/576 = Re(vdot(a,b))/576.
a=rng.normal(size=576)+1j*rng.normal(size=576)
b=rng.normal(size=576)+1j*rng.normal(size=576)
real_dot=float(jets.real_material(a)@jets.real_material(b))
complex_dot=float(np.vdot(a,b).real/576)
assert abs(real_dot-complex_dot) < 1e-14
R=rng.normal(size=(576,9))+1j*rng.normal(size=(576,9))
u=rng.normal(size=9)
mapping_error=float(np.linalg.norm(jets.real_material(R@u)-jets.real_material(R)@u))
assert mapping_error < 1e-13
report['checks']['real_hilbert']=dict(real_dot=real_dot,complex_real_dot=complex_dot,mapping_error=mapping_error)

# Shared real LS enforces the standard normal equation; group conflict identity
# uses the SAME lower-order b and common real columns.
A=rng.normal(size=(512,8)); b=rng.normal(size=512)
z,c=jets.common_compensation(A,b)
normal=float(np.linalg.norm(A.T@c))
within=0.; conflict=0.
for sl in (slice(0,256),slice(256,512)):
    zi,ci=jets.common_compensation(A[sl],b[sl])
    within+=float(ci@ci)
    conflict+=float((zi-z)@(A[sl].T@A[sl])@(zi-z))
defect=abs(float(c@c)-within-conflict)
assert normal < 1e-10 and defect < 1e-10
report['checks']['common_real_LS']=dict(normal_residual=normal,group_identity_absolute_defect=defect)

# GN feasibility and actual quadratic-model reduction with a real material
# pullback. This synthetic J is algebra only, not Maxwell evidence.
x0=np.zeros(66); x0[0]=.1; x=x0+.01*rng.normal(size=66)
M=jets.real_material(jets.material_derivative(prior,x))
M0=jets.real_material(jets.material_derivative(prior,x0))
v=rng.normal(size=66); v/=np.linalg.norm(M0@v)
J=rng.normal(size=(512,66)); residual=rng.normal(size=512)
step,predicted=jets.gn_step(prior,x0,x,v,J,residual)
unit=M0@v;unit/=np.linalg.norm(unit)
tangent_defect=abs(float(unit@M@step))
r=np.r_[residual,x]
independent=.5*(r@r-np.linalg.norm(r+np.r_[J@step,step])**2)
assert tangent_defect < 1e-12
assert np.linalg.norm(step) <= .5+1e-12 and np.linalg.norm(M@step)<=.01+1e-12
assert abs(predicted-independent)<1e-12 and predicted>=-1e-12
report['checks']['GN_real_constraints']=dict(signed_progress_tangent_defect=tangent_defect,
    latent_step=float(np.linalg.norm(step)),linear_material_step=float(np.linalg.norm(M@step)),
    predicted_reduction=float(predicted),independent_quadratic_reduction=float(independent))

# Concrete frozen-protocol issue: rank-two original saturated material has 64
# exact decoder gauges. Current GN retains them despite effective-only arm text.
x0=np.zeros(66);x0[0]=1.5
x=x0.copy();x[64]=.4
M=jets.real_material(jets.material_derivative(prior,x))
v=np.zeros(66);v[64]=1
J=M[:512].copy()  # Valid real material pullback; deliberately no physics.
step,predicted=jets.gn_step(prior,x0,x,v,J,np.zeros(512))
material_change=float(np.linalg.norm(jets.real_material(jets.decode(prior,x+step)[0]-jets.decode(prior,x)[0])))
_,s,Vh=np.linalg.svd(M,full_matrices=False)
good=s>max(s[0]*1e-8,1e-10)
projection=Vh[good].T@(Vh[good]@step)
kernel_step=float(np.linalg.norm(step-projection))
# Preserve the pre-fix formula as a local algebra fixture; never edit the
# implementation under review. Root has patched the material-effective frame.
from scipy.linalg import null_space
M0=jets.real_material(jets.material_derivative(prior,x0))
unit=M0@v;unit/=np.linalg.norm(unit)
oldZ=null_space((unit@M)[None,:],rcond=1e-10)
oldstep=-oldZ@np.linalg.lstsq(np.vstack((J@oldZ,oldZ)),np.r_[np.zeros(512),x],rcond=1e-10)[0]
oldstep*=min(1.,.5/max(np.linalg.norm(oldstep),1e-300),.01/max(np.linalg.norm(M@oldstep),1e-300))
oldkernel=float(np.linalg.norm(oldstep-Vh[good].T@(Vh[good]@oldstep)))
oldmaterial=float(np.linalg.norm(jets.real_material(jets.decode(prior,x+oldstep)[0]-jets.decode(prior,x)[0])))
issue=dict(priority='P1',status='RESOLVED_BY_ROOT' if kernel_step<1e-10 else 'OPEN',function='gn_step',issue='Pre-fix GN retained exact decoder gauges, violating effective-material-only frozen correction arm',
    original_prior_rank=int(good.sum()),latent_step=float(np.linalg.norm(step)),kernel_step=kernel_step,
    actual_material_RMS_change=material_change,predicted_prior_reduction=float(predicted),
    pre_fix_reproduced_kernel_step=oldkernel,pre_fix_reproduced_actual_material_RMS_change=oldmaterial,
    required_fix='Z=Veff@null_space((unit@M@Veff)[None,:]); Veff uses frozen SVD thresholds')
assert oldkernel>.49 and oldmaterial<1e-15
report['issues'].append(issue)
report['checks']['saturated_GN_gauge_counterexample']=issue

# C2 clipping does not admit a unique cubic coefficient exactly at a seam.
# For h<0 at s=1, saturation = 1 + (5/4) h^3 + O(h^4);
# for h>0 it remains exactly 1. Caller must apply regular-stratum gates.
report['issues'].append(dict(priority='GATE',function='decoder_coefficients',
    issue='Exact seam points have no single third Taylor coefficient; finite paths crossing seams require regular-stratum handling',
    module_scope='module docstring correctly disclaims global C3; caller must enforce frozen regular-stratum gate'))
assert 'a30.physics' not in sys.modules
report['new_Maxwell_RHS']=0
report['new_Maxwell_states']=0
report['status']='DECODER_AND_REAL_HILBERT_PASS_GN_PROTOCOL_FIX_REQUIRED' if kernel_step>1e-10 else 'PASS_WITH_REGULAR_STRATUM_GATE'
(P/'audit/jets_math/RESULT.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(dict(status=report['status'],issues=report['issues'],new_Maxwell_RHS=0),indent=2))
