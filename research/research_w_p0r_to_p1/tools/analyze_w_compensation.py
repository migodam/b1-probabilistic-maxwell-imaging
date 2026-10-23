"""Offline finite-secant decomposition of saved full-Maxwell endpoints only."""
import json,time
from pathlib import Path
import numpy as np
from scipy.special import ndtr

root=Path(__file__).resolve().parents[1]


def main():
    start=time.process_time();source=root/'results/B2/validation/joint_seed2_narrow_fullj'
    cal=np.load(source/'cal_PAIR.npz');joint=np.load(source/'joint8_PAIR.npz')
    delta_wall=(cal['fine1']-cal['fine0'])/.01
    delta_joint=(joint['fine1']-joint['fine0'])/.01
    delta_target=delta_joint-delta_wall
    d0,d1,dj=(float(np.linalg.norm(x)) for x in (delta_wall,delta_target,delta_joint))
    cross=float(2*np.real(np.vdot(delta_wall,delta_target)))
    normsum=d0*d0+d1*d1
    out=root/'results/B2/compensation_mechanism';out.mkdir(exist_ok=True)
    np.savez_compressed(out/'FINITE_SECANTS.npz',delta_wall=delta_wall,delta_target_induced_increment=delta_target,
                        delta_joint=delta_joint,noise_quadrature_std=.01)
    answer=dict(identity='ORACLE_PHYSICAL_EXISTENCE_SAVED_ENDPOINT_DIAGNOSTIC_NOT_ONLINE',
        formula='DeltaJoint=DeltaWall+[(Joint1-Cal1)-(Joint0-Cal0)]',
        d_wall=d0,d_target_induced_increment=d1,d_joint=dj,
        realified_response_cosine=cross/(2*d0*d1),
        sum_component_squared_norms=normsum,quadratic_interference_term=cross,
        squared_identity_error=abs(dj*dj-(normsum+cross)),
        destructive_amplitude_cancellation=1-dj/(d0+d1),
        fixed_simple_hypotheses_equal_prior_testing_error=float(ndtr(-dj/2)),
        testing_error_scope='Optimal test for ONLY these two fixed means, identical Gaussian covariance, equal assigned testing priors. Not posterior odds or composite-hypothesis calibration.',
        gamma=1,noise_unchanged=True,wall_subtraction_role='offline decomposition only; inference uses full total joint data',
        posterior_probabilities=None,source_paths=['results/B2/validation/joint_seed2_narrow_fullj/cal_PAIR.npz','results/B2/validation/joint_seed2_narrow_fullj/joint8_PAIR.npz'],
        new_rhs=0,source_full_Maxwell_cost='already charged in B2; cached algebra is only incrementally zero physics',
        pure_analysis_cpu_seconds=time.process_time()-start)
    (out/'RESULTS.json').write_text(json.dumps(answer,indent=2,allow_nan=False)+'\n')
    print(json.dumps(answer,allow_nan=False))


if __name__=='__main__':main()
