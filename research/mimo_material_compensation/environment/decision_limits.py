"""Outward task intervals and sensitivity to a declared amplitude-dependent error."""
import json
import math
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name):
    return json.loads((HERE / name).read_text())


def up(x):
    return math.nextafter(float(x), math.inf)


def down(x):
    return math.nextafter(float(x), -math.inf)


def specifications():
    for side in (32, 96):
        proposal = load(f"low_energy_{side}/proposal.json")
        yield (side, proposal["core"], proposal["radius"], proposal["upper400"],
               proposal["lower460"], proposal["operating_energy"], proposal["beta"])
    fine = load("matched_amplitude.json")
    row = next(row for row in fine["outputs"]["460"]["budgets"] if row["beta"] == .01)
    upper = load("kyp_verified_responses_96_2.0_wave_96_2.0_1_400.json")
    yield 96, 2., fine["radius_sqrtW"], upper["upper_over_2r"], row["open_lower_over_two_r"], None, .01


def calculate(side, t, r, upper, lower, energy, beta):
    # Binary64 proposal values and the already enclosed bound endpoints are exact inputs.
    U, L = map(Fraction.from_float, (upper, lower))
    T, B = map(Fraction.from_float, (t, beta))
    result = dict(side=side, contrast_increment=t, receiver_radius=r,
                  same_frequency_pair_upper_over_2r=upper,
                  new_frequency_full_classes_lower_over_2r=lower,
                  necessary_energy_old_lower=down(1/U**2),
                  sufficient_energy_open_threshold_upper=up(1/L**2),
                  energy_ratio_lower=down((L/U)**2),
                  roi_outer_low=[0., beta],
                  roi_outer_high=[down(T-B), up(T+B)],
                  source="Open-circuit whole-class lower endpoint. Two-state verified upper endpoint.")
    if energy:
        # Below result, both class score intervals stay apart under a symmetric extra error.
        result["operating_energy"] = energy
        sqrt_lower = down(math.sqrt(energy))
        gap = Fraction.from_float(sqrt_lower)*L - 1
        result["extra_error_as_fraction_of_r_at_operating_point_below"] = down(gap)
    curves = []
    for eta_over_r in (0., .1, .25, .5, 1., 2.):
        # At energy E, total error <= r*(1 + eta_over_r*sqrt(E)).
        factor = L - Fraction.from_float(eta_over_r)
        curves.append(dict(eta_over_r=eta_over_r,
                           sufficient_E_strictly_above=up(1/factor**2) if factor>0 else None))
    result["calibration_budget_curve"] = curves
    return result


if __name__ == "__main__":
    out = dict(contract="Known nondispersive specification union. Calibration is a conditional bound. Old necessary bound only uses zero calibration error and an attainable absolute receiver noise ball.",
               records=[calculate(*item) for item in specifications()])
    (HERE/"decision_limits.json").write_text(json.dumps(out,indent=2)+"\n")
