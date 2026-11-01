#!/usr/bin/env python3
"""Exact rational audit of the existing two-cell algebraic VIE example.

Standard-library only; no SymPy, NumPy, network, or existing JSON is needed.
Run: python exact_two_cell.py [--output /path/to/exact_two_cell.json]

This checks the model already used in branches_checks.json. It does not change
that file, generate another example, or claim a Maxwell imaging benchmark.
The data are DEFINED as F(2, 7), not as a rounded decimal observation.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from functools import reduce
import json
from math import gcd, isqrt, lcm
from pathlib import Path


@dataclass(frozen=True)
class QC:
    """Complex number whose real and imaginary parts are exact rationals."""

    re: F = F(0)
    im: F = F(0)

    @staticmethod
    def coerce(value):
        return value if isinstance(value, QC) else QC(F(value))

    def __add__(self, other):
        other = self.coerce(other)
        return QC(self.re + other.re, self.im + other.im)

    __radd__ = __add__

    def __neg__(self):
        return QC(-self.re, -self.im)

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) + (-self)

    def __mul__(self, other):
        other = self.coerce(other)
        return QC(self.re * other.re - self.im * other.im,
                  self.re * other.im + self.im * other.re)

    __rmul__ = __mul__

    def conj(self):
        return QC(self.re, -self.im)

    def __truediv__(self, other):
        other = self.coerce(other)
        denominator = other.re * other.re + other.im * other.im
        if not denominator:
            raise ZeroDivisionError("Exact complex denominator is zero")
        numerator = self * other.conj()
        return QC(numerator.re / denominator, numerator.im / denominator)

    def __rtruediv__(self, other):
        return self.coerce(other) / self

    def nonzero(self):
        return bool(self.re or self.im)

    def as_json(self):
        return {"real": str(self.re), "imag": str(self.im)}


ZERO = QC()
A_G = QC(F(-3, 50), F(1, 50))
B_G = QC(F(13, 200), F(7, 1000))
D_G = QC(F(-3, 20), F(1, 50))
E = (F(1), F(7, 10))
S = (F(1), F(13, 10))
BOX = (F(1, 50), F(20))


def physical_forward(material):
    """Return S (I - diag(m) G)^(-1) diag(m) E and det(I-XG)."""
    x, z = material
    l11, l12 = 1 - x * A_G, -x * B_G
    l21, l22 = -z * B_G, 1 - z * D_G
    determinant = l11 * l22 - l12 * l21
    rhs1, rhs2 = x * E[0], z * E[1]
    current1 = (l22 * rhs1 - l12 * rhs2) / determinant
    current2 = (l11 * rhs2 - l21 * rhs1) / determinant
    return S[0] * current1 + S[1] * current2, determinant


def primitive_integer_polynomial(coefficients):
    """Clear rational denominators and common factors; leading sign positive."""
    denominator = lcm(*(v.denominator for v in coefficients))
    integers = [int(v * denominator) for v in coefficients]
    common = reduce(gcd, (abs(v) for v in integers))
    if common == 0:
        raise ValueError("The elimination polynomial vanishes identically")
    integers = [v // common for v in integers]
    if integers[0] < 0:
        integers = [-v for v in integers]
    return tuple(integers)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_suffix(".json"))
    args = parser.parse_args()

    first_material = (F(2), F(7))
    y, _ = physical_forward(first_material)
    assert y == QC(F(1826362941900, 315300667481),
                   F(157970241200, 315300667481))

    # The response numerator is A0*x + B0*z + C0*x*z; moving y*det
    # to the left gives A*x + B*z + C*x*z - y = 0.
    a0 = S[0] * E[0]
    b0 = S[1] * E[1]
    c0 = B_G * (S[0] * E[1] + S[1] * E[0]) - D_G * a0 - A_G * b0
    a = a0 + y * A_G
    b = b0 + y * D_G
    c = c0 - y * (A_G * D_G - B_G * B_G)

    # For real x, z=(y-A*x)/(B+C*x) is real iff the following numerator
    # has zero imaginary part, provided B+C*x is nonzero.
    rational_polynomial = (
        (-a * c.conj()).im,
        (y * c.conj() - a * b.conj()).im,
        (y * b.conj()).im,
    )
    p2, p1, p0 = primitive_integer_polynomial(rational_polynomial)
    assert (p2, p1, p0) == (764347498747, -4999308072894, 6941226150800)
    factor_scale, second_numerator = 764347498747, 3470613075400
    assert (p2, p1, p0) == (
        factor_scale, -second_numerator - 2 * factor_scale,
        2 * second_numerator,
    )

    discriminant = p1 * p1 - 4 * p2 * p0
    sqrt_discriminant = isqrt(discriminant)
    assert discriminant > 0 and sqrt_discriminant**2 == discriminant
    assert sqrt_discriminant == 1941918077906

    first_components = sorted((F(-p1 - sqrt_discriminant, 2 * p2),
                               F(-p1 + sqrt_discriminant, 2 * p2)))
    root_records = []
    exact_roots = []
    for x in first_components:
        recovery_denominator = b + c * x
        assert recovery_denominator.nonzero()
        z = (y - a * x) / recovery_denominator
        assert z.im == 0
        material = (x, z.re)
        response, determinant = physical_forward(material)
        derivative = 2 * p2 * x + p1
        checks = {
            "polynomial_exactly_zero": p2 * x * x + p1 * x + p0 == 0,
            "simple_polynomial_root": derivative != 0,
            "recovered_material_exactly_real": z.im == 0,
            "strictly_inside_positive_box": all(BOX[0] < v < BOX[1] for v in material),
            "recovery_denominator_nonzero": recovery_denominator.nonzero(),
            "physical_determinant_nonzero": determinant.nonzero(),
            "physical_response_exactly_equals_y": response == y,
            "bilinear_equation_exactly_zero":
                a * x + b * z + c * x * z - y == ZERO,
        }
        assert all(checks.values())
        exact_roots.append(material)
        root_records.append({
            "material_exact": [str(v) for v in material],
            "material_approx": [float(v) for v in material],
            "quadratic_derivative_exact": str(derivative),
            "checks": checks,
        })
    assert exact_roots == [first_material,
                           (F(3470613075400, 764347498747),
                            F(394249570000, 195029801279))]

    # Sylvester's criterion for H = (G+G*)/2 in this complex-symmetric model.
    h11, h12, h22 = A_G.re, B_G.re, D_G.re
    h_determinant = h11 * h22 - h12 * h12
    assert h11 < 0 and h22 < 0 and h_determinant == F(191, 40000) > 0
    pole_certificate = {
        "hermitian_part_exact": [[str(h11), str(h12)], [str(h12), str(h22)]],
        "hermitian_part_determinant_exact": str(h_determinant),
        "hermitian_part_negative_definite": True,
        "domain": "Every strictly positive diagonal real X, including [1/50,20]^2.",
        "identity": "I-XG = sqrt(X) [I-sqrt(X) G sqrt(X)] inv(sqrt(X))",
        "proof": "H<0 implies Hermitian(I-sqrt(X)Gsqrt(X))=I-sqrt(X)Hsqrt(X)>I; hence that matrix and I-XG are invertible.",
        "singular_value_caveat": "Similarity is nonunitary: sigma_min(I-XG)>=1 does not follow. A valid coarse bound is sqrt(min(m)/max(m)).",
    }

    result = {
        "scope": "Exact rational audit of the existing two-cell algebraic VIE; not an imaging benchmark or a claimed Hankel acquisition grid.",
        "arithmetic": "Python standard-library Fraction; all pass/fail assertions are exact.",
        "model": {
            "G": [[A_G.as_json(), B_G.as_json()], [B_G.as_json(), D_G.as_json()]],
            "E": [str(v) for v in E],
            "S": [str(v) for v in S],
            "box": [str(v) for v in BOX],
            "data_definition": "y := F(2,7), F(m)=S (I-diag(m)G)^(-1) diag(m) E",
            "y_exact": y.as_json(),
            "y_approx": [float(y.re), float(y.im)],
        },
        "elimination": {
            "coefficients_descending_primitive_integer": [p2, p1, p0],
            "exact_factorization": "(m1-2)*(764347498747*m1-3470613075400)",
            "discriminant_exact": str(discriminant),
            "sqrt_discriminant_exact": str(sqrt_discriminant),
            "discriminant_strictly_positive": True,
            "derivation": "Im[(y-A*m1)*conj(B+C*m1)]=0; then m2=(y-A*m1)/(B+C*m1).",
            "local_isolation": "Both roots are simple and both recovery denominators are nonzero; the original real response equations are locally equivalent to the eliminated quadratic plus recovery.",
        },
        "roots": root_records,
        "no_physical_poles_certificate": pole_certificate,
        "all_exact_checks_passed": True,
        "limitations": [
            "The two endpoints are known; no unknown endpoint generation or global optimizer success is established.",
            "The example has one existing source and one receiver; it is not a Q6/Q20 run.",
            "No floating-point Jacobian rank or grid scan is used to prove isolation or absence of positive-material poles.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"all_exact_checks_passed": True,
                      "number_of_real_legal_roots": len(root_records),
                      "output": str(args.output.resolve())}, indent=2))


if __name__ == "__main__":
    main()
