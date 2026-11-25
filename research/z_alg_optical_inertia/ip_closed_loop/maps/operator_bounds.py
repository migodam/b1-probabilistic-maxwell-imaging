"""Certify the exact binary-symbol inverse transforms by independent 128-bit Arb DFTs."""
from flint import acb, arb, ctx
import numpy as np

from common import HERE, Geometry, clock, elapsed, geom_path, save_json


def ball(z):
    return acb(float(z.real), float(z.imag))


def inverse_dft_table(rows):
    transformed = [acb.dft(row, inverse=True) for row in rows]
    return zip(*(acb.dft(list(col), inverse=True) for col in zip(*transformed)))


def discrepancy(a, b):
    return float(np.nextafter(float(abs(a - ball(b)).upper()), np.inf))


def main():
    ctx.prec = 128
    result = {}
    for side in (32, 96):
        start = clock()
        g = Geometry(geom_path(side))
        z = [[ball(a) for a in line] for line in g.symbol]
        z2 = [[acb(a.real**2 + a.imag**2) for a in line] for line in z]
        bound = {}
        for name, rows, reference in [
            ("kernel", z, g.lag), ("gram_kernel", z2, g.gram_lag)
        ]:
            err = max(discrepancy(a, b) for line, ordinary in zip(inverse_dft_table(rows), reference)
                      for a, b in zip(line, ordinary))
            bound[f"{name}_entry_error"] = err
        high = max(abs(a).upper() for line in z for a in line)
        bound["full_operator_norm_upper"] = float(np.nextafter(float(high), np.inf))
        result[str(side)] = dict(**bound, precision_bits=ctx.prec, **elapsed(start))
        print(side, bound, flush=True)
    save_json(HERE / "results" / "operator_bounds.json", result)


if __name__ == "__main__":
    main()
