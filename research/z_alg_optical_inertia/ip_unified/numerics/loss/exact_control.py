"""Rational passive reciprocal control for the partial negative split."""
from fractions import Fraction as Q
from math import sqrt

from run import HERE, write_json


def main():
    # G=diag(i/10,i/10,beta+2,beta+1/100), B=[I;0], S=[I,0].
    # The declared material is chi=(1+i/5)x, cap one on its real part.
    # Generate the data at x=(2,2,0,0).
    beta, alpha = Q(25,26), Q(5,26)
    r = Q(1,2)
    radiation = Q(1,10)
    mr = radiation + alpha
    # This product is verified below directly as A*(y_clean)=chi.
    a,b = Q(26,25), -Q(1,5)
    y_re = (2*a + Q(2,5)*b)/(a*a+b*b)
    y_im_clean = (Q(2,5)*a - 2*b)/(a*a+b*b)
    assert a*y_re-b*y_im_clean==2 and a*y_im_clean+b*y_re==Q(2,5)
    y_im = y_im_clean+Q(1,100)
    cb = Q(1,4)/mr
    hy = y_re*y_re+y_im*y_im - ((y_re+Q(1,2))**2+y_im*y_im)/(1+beta)
    tau = Q(3,50)
    gamma_controlled = Q(1,100)/alpha
    gamma_full = 2/alpha
    assert gamma_controlled < tau < gamma_full
    # Both negative centres vanish, so M_tau=tau CB in the controlled space.
    assert hy > 0 and gamma_full*cb > hy
    mixed = hy-tau*cb
    assert mixed > r*r
    # Exact in-ellipsoid port witness computed with a rational contraction.
    # Its Frobenius distance is below the same radius in both receiver channels.
    ambient = y_re*y_re + (y_im-Q(1,2)/mr)**2
    shrink = Q(21,25)
    assert shrink**2*ambient < (Q(1,2)/mr)**2
    projected_upper_sq = 2*(1-shrink)**2*ambient
    assert projected_upper_sq < r*r
    # After the corresponding complete loss shift, the smaller dark mode
    # is positive. One free negative mode remains, and both sources agree.
    k_obs = 1+beta+tau*mr
    l_re, l_im = y_re+Q(1,2), y_im+tau/2
    shifted = y_re*y_re+y_im*y_im-(l_re*l_re+l_im*l_im)/k_obs
    assert shifted >= mixed
    assert 2*Q(1,100)**2 < r*r
    vals = dict(tangent="1/5", cap="1", radius=float(r), true_re_peak=2,
                source_dimensions=2, original_inertia=2, controlled_negative=1, free_negative=1,
                pure_receiver_upper=sqrt(float(projected_upper_sq)), unshifted_distance=0,
                full_negative_common_matrix_cap_distance=0,
                mixed_distance_lower=sqrt(float(mixed)), full_shift_distance=sqrt(float(shifted)),
                parameter_exact=str(tau), physical_generating_residual_sq_exact=str(2*Q(1,100)**2),
                pure_receiver_upper_sq_exact=str(projected_upper_sq),
                mixed_distance_sq_exact=str(mixed), full_shift_distance_sq_exact=str(shifted))
    write_json(HERE/"exact_control.json",vals)
    print(vals,flush=True)


if __name__ == "__main__":
    main()
