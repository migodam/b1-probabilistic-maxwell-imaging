"""Zero-Maxwell checks of nontrivial affine-slice and gauge obligations."""
from pathlib import Path
import json
import numpy as np
from engine import Config, run_arm, effective_spectrum, PaidEvaluator, transverse_correct


class CurvedFiber:
    """An explicitly NON-Maxwell circle with one exact decoder gauge."""
    def __init__(self): self.calls = 0
    def decode(self, x): return dict(material=np.asarray(x[:2], complex), prior_legal=True)
    def evaluate(self, x, need_jac=False):
        self.calls += 1
        x = np.asarray(x, float).copy()
        r = np.array([x[0]*x[0]+x[1]*x[1]-1.])
        out = dict(x=x, residual=r, phi=.5*float(r@r), material=x[:2].astype(complex), prior_legal=True)
        if need_jac:
            out.update(jac=np.array([[2*x[0], 2*x[1], 0.]]), dmaterial=np.array([[1.,0.,0.],[0.,1.,0.]], complex))
        return out


class NonGaussianRegularized:
    """Quadratic surrogate used only to verify the native regularizer hook."""
    def evaluate(self, x, need_jac=False):
        x=np.asarray(x,float).copy()
        r=np.array([x[0]-1.])
        # The third parameter is an exact decoder gauge. The supplied native
        # regularizer is allowed to move it under full_coordinate_corrector.
        reg=.5*(x[0]**2+3*(x[1]-2.)**2+2*(x[2]-3.)**2)
        q=dict(x=x,residual=r,phi=.5*float(r@r),material=x[:2].astype(complex),
            prior_legal=True,reg_value=float(reg),
            reg_gradient=np.array([x[0],3*(x[1]-2.),2*(x[2]-3.)]),
            reg_hessian=np.diag([1.,3.,2.]))
        if need_jac:q.update(jac=np.array([[1.,0.,0.]]),dmaterial=np.array([[1.,0.,0.],[0.,1.,0.]],complex))
        return q


def main():
    problem = CurvedFiber()
    anchor = problem.evaluate(np.array([1., 0., 11.]), True)
    config = Config(arm='H1', weak_rank=1, layers=1, proposals_per_layer=2,
        layer_radii=(.1,), predictor_rms=.1, corrector_states=5,
        prior_weight=0., rms_denominator=2, material_trust=.1,
        noise_phi_max=1e-8, local_rhs_cap=1000)
    result = run_arm(problem.evaluate, anchor['x'], config,
        initial_state=anchor, decoder=problem.decode,
        random_draws=np.array([[1.],[-1.]]))
    max_fixed = max((h['weak_coordinate_defect'] for row in result['trace']
        if 'correction' in row for h in row['correction']['history']), default=0.)
    end_residual = max(float(np.linalg.norm(s['residual'])) for s in result['states'][1:])
    gauge_motion = max(float(abs(s['x'][2]-11.)) for s in result['states'])
    spectrum = effective_spectrum(anchor['jac'], anchor['dmaterial'], denominator=2)
    # A positive-residual stationary point of f(x)=x²+1 has GN-map derivative
    # -1 for T=x-2x(x²+1)/(4x²+1), while dropping DJ^T r incorrectly gives +1.
    def nonroot_map(x): return x-2*x*(x*x+1)/(4*x*x+1)
    step = 1e-6
    fd = (nonroot_map(step)-nonroot_map(-step))/(2*step)
    regularized = NonGaussianRegularized()
    regularized_anchor=regularized.evaluate(np.zeros(3),True)
    native_config=Config(corrector_states=5,latent_trust=10.,material_trust=10.,
        max_latent_radius=10.,prior_weight=999.,rms_denominator=2,
        full_coordinate_corrector=True)
    paid=PaidEvaluator(regularized.evaluate,1000,initial=regularized_anchor)
    native_endpoint,native_trace,_=transverse_correct(paid,np.zeros(3),None,native_config,
        reference=regularized_anchor,prior_center=np.full(3,100.))
    native_error=float(np.linalg.norm(native_endpoint['x']-np.array([.5,2.,3.])))
    checks = dict(status='PASS', Maxwell_RHS=0,
        curved_fiber=dict(available=len(result['states'])-1,
                          max_fixed_coordinate_defect=max_fixed,
                          max_endpoint_residual=end_residual,
                          exact_decoder_gauge_motion=gauge_motion,
                          decoded_gauge_dimension=spectrum['decoder_gauge_dimension']),
        positive_residual=dict(analytic_DT=-1., finite_difference_DT=float(fd),
                               omitting_residual_term_DT=1.),
        native_regularizer_hook=dict(known_minimizer=[.5,2.,3.],
            returned=native_endpoint['x'].tolist(),minimizer_error=native_error,
            original_gaussian_prior_ignored=True,material_gauge_move_allowed=True),
        semantics='explicit nonlinear toy only; no physical Maxwell claim')
    if len(result['states']) != 3 or max_fixed > 1e-8 or end_residual > 1e-6 or gauge_motion != 0. or abs(fd+1) > 1e-8 or native_error>1e-6:
        checks['status'] = 'FAIL'
    p = Path(__file__).resolve().parent/'ALGEBRAIC_CONTRACT_CHECKS.json'
    p.write_text(json.dumps(checks, indent=2)+'\n')
    print(json.dumps(checks))


if __name__ == '__main__': main()
