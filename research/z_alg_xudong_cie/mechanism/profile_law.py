"""Real-scalar least-squares profile, with the frozen PixelLaw normalizations."""

import numpy as np

from cross_engine import PixelLaw


class ProfileLaw(PixelLaw):
    """Remove the positive phase penalty from the auxiliary pixel law."""

    def gradient_parts(self, J, E):
        if self.passive:
            raise ValueError("The frozen profile ablation uses synthetic real media.")
        j, e = J*self.weight, E*self.weight
        a = np.sum(abs(j)**2, axis=1, keepdims=True)
        b = np.sum(abs(e)**2, axis=1, keepdims=True)
        real_c = np.sum(j.conj()*e, axis=1, keepdims=True).real
        safe = b + self.delta
        power = np.maximum(a*b - real_c*real_c, 0.)
        direct = (b*j - real_c*e) / safe
        ge = (a*e - real_c*j) / safe - power*e/(safe*safe)
        return direct*self.weight, ge*self.weight, float(.5*np.sum(power/safe)), 0.

    def line_value(self, terms, alpha):
        if self.passive:
            raise ValueError("The frozen profile ablation uses synthetic real media.")
        a, b, c = [t0 + alpha*(t1+alpha*t2) for t0, t1, t2 in terms]
        power = np.maximum(a*b - c.real*c.real, 0.)
        return float(.5*np.sum(power/(np.maximum(b, 0.)+self.delta))), 0.
