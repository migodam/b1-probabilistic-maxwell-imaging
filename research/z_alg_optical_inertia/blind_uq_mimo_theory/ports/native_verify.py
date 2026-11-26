"""Native 96² forward enclosure via Arb kernel and verified inverse residual."""
import argparse
import ctypes
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from time import perf_counter, process_time

import numpy as np
from flint import arb
from scipy.linalg import lu_factor, lu_solve

from verify_witness import as_acb, exact_kernel, outward
from flint import ctx

HERE = Path(__file__).resolve().parent
unit = arb(1) / 2**53
tiny = arb(1) / 2**1022
pointer = ctypes.POINTER(ctypes.c_double)


def gamma(steps):
    u = steps * unit
    return u / (1-u)


def realparts(x):
    x = np.asarray(x, dtype=np.complex128)
    return np.ascontiguousarray(x.real), np.ascontiguousarray(x.imag)


def addr(x):
    return x.ctypes.data_as(pointer)


class ExplicitDots:
    """Check floating products with known binary64 arithmetic and order."""
    def __init__(self, dylib, threads):
        library = ctypes.CDLL(str(dylib))
        self.product_cb = library.complex_product
        self.product_cb.argtypes = 3 * [ctypes.c_size_t] + 6 * [pointer] + [ctypes.c_size_t]
        self.product_cb.restype = ctypes.c_int
        self.norm_cb = library.binary64_sum_squares
        self.norm_cb.argtypes = [ctypes.c_size_t, pointer, pointer]
        self.norm_cb.restype = ctypes.c_double
        self.threads = threads

    def norm(self, X):
        re, im = realparts(X)
        computed = self.norm_cb(re.size, addr(re), addr(im))
        if not np.isfinite(computed):
            raise ArithmeticError("Non-finite computed norm")
        # Permit subnormal flushing as error per arithmetic operation.
        return ((arb(computed) + 8*re.size*tiny) / (1-gamma(3*re.size+8))).sqrt().upper()

    def product(self, A, B):
        ar, ai = realparts(A)
        br, bi = realparts(B)
        shape = (len(A), B.shape[1])
        cr, ci = np.zeros(shape), np.zeros(shape)
        if self.product_cb(len(A), A.shape[1], B.shape[1], *map(addr, [ar, ai, br, bi, cr, ci]), self.threads):
            raise ArithmeticError("Explicit product did not receive its full worker count")
        C = np.empty(shape, complex)
        C.real, C.imag = cr, ci
        if not np.all(np.isfinite(C)):
            raise ArithmeticError("Non-finite rounded product")
        error = gamma(8*A.shape[1]+8) * self.norm(A) * self.norm(B) + 16*A.shape[1]*max(shape)*tiny
        return C, error.upper()


def exact_bare_representative(spectrum):
    enclosure = exact_kernel(spectrum)
    approx = np.array([complex(z) for z in enclosure])
    errors = [abs(z-as_acb(a)).upper() for a,z in zip(approx,enclosure,strict=True)]
    return approx, max(errors)


def certify(case_id, data, bare_lags, bare_bound, dots):
    start = perf_counter(), process_time()
    # Stored labels are read only for this forward validation. Exact zero
    # contrast eliminates exactly zero current; the inverse used no support.
    with np.load(HERE/"private"/f"{case_id}.npz") as archive:
        m = archive["material"].ravel()
        j = archive["current"]
    with np.load(HERE/"public"/"cases"/f"{case_id}.npz") as measurement:
        Y, hold = measurement["Y"], measurement["Y_hold"]
        radius, hold_radius = arb(float(measurement["radius"])), arb(float(measurement["radius_hold"]))
    ids = np.flatnonzero(m)
    assert np.all(j[m == 0] == 0), "Exact support elimination does not accept nonzero inactive current"
    assert np.min(m) >= 0 and np.max(m) <= 4.5
    d = len(ids)
    stride_y, stride_x = data["pad"]
    y, x = np.unravel_index(ids, tuple(data["grid_shape"]))
    lags = (y[:, None]-y[None, :]) % stride_y * stride_x + (x[:, None]-x[None, :]) % stride_x
    bare = bare_lags[lags]
    del lags
    lr, lr_bound = dots.product(data["low_rank_left"][ids], data["low_rank_right"][:, ids])
    G = bare + lr
    Gerror = d*bare_bound + lr_bound + gamma(5)*(dots.norm(bare)+dots.norm(lr)) + 4*d*tiny
    del lr, bare
    A = np.eye(d)-m[ids, None]*G
    max_contrast = arb(float(np.max(m)))
    # Identity diagonal and elementwise multiplication/subtraction only.
    Aerror = max_contrast*Gerror + gamma(6)*(arb(d).sqrt()+max_contrast*dots.norm(G)) + 6*d*tiny
    del G
    lstart = perf_counter()
    lu = lu_factor(np.asfortranarray(A), overwrite_a=True, check_finite=False)
    inverse = lu_solve(lu, np.eye(d, dtype=np.complex128), check_finite=False)
    del lu
    lu_seconds = perf_counter()-lstart
    inverse_bound_input = dots.norm(inverse)
    cstart = perf_counter()
    reconstructed_identity, product_error = dots.product(A, inverse)
    del inverse
    eigerr = (
        dots.norm(np.eye(d)-reconstructed_identity)
        + gamma(5)*(arb(d).sqrt()+dots.norm(reconstructed_identity)) + 4*d*tiny
        + product_error + Aerror*inverse_bound_input
    )
    dot_seconds = perf_counter()-cstart
    del reconstructed_identity
    if not eigerr < 1:
        raise ArithmeticError(f"No strict inverse defect enclosure ({eigerr})")
    inverse_bound = inverse_bound_input / (1-eigerr)
    B = m[ids, None]*data["E0"][ids]
    Berror = gamma(4)*max_contrast*dots.norm(data["E0"][ids]) + 4*len(ids)*tiny
    j = j[ids]
    rhs_prediction, state_rounding = dots.product(A, j)
    state_error = (dots.norm(rhs_prediction-B) + gamma(5)*(dots.norm(rhs_prediction)+dots.norm(B))
                   + state_rounding + Aerror*dots.norm(j) + Berror + 4*d*tiny)
    current_error = state_error*inverse_bound
    results = {}
    for name, S, Y0, observed_radius in (("train", data["S"], Y, radius), ("hold", data["S_hold"], hold, hold_radius)):
        receive = S[:, ids]
        predicted, port_rounding = dots.product(receive, j)
        display_error = dots.norm(predicted-Y0)
        arithmetic_error = gamma(5)*(dots.norm(predicted)+dots.norm(Y0)) + 4*Y0.size*tiny
        port_solution_error = dots.norm(receive)*current_error + port_rounding
        normalized_upper = (port_solution_error+display_error+arithmetic_error)/observed_radius
        results[name] = dict(distance_over_radius_upper=outward(normalized_upper)[1],
                             forward_error_upper=outward(port_solution_error)[1],
                             feasible=bool(normalized_upper < 1))
    summary = dict(case_id=case_id, native_grid_size=9216, exact_nonzero_current_state_dimension=d,
                   inverse_bound_2_from_Frobenius=outward(inverse_bound)[1], matrix_enclosure_Fro=outward(Aerror)[1],
                   approximate_right_inverse_defect_Fro_upper=outward(eigerr)[1],
                   state_residual_Fro_upper=outward(state_error)[1], **results,
                   seconds_dense_LU_and_all_inverse_columns=lu_seconds, seconds_verified_product=dot_seconds,
                   seconds_wall=perf_counter()-start[0], seconds_CPU=process_time()-start[1],
                   dense_LU=1, paid_inverse_columns=d)
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--items", nargs="+", type=int, default=(1,2))
    p.add_argument("--threads", type=int, default=4)
    args = p.parse_args()
    ctx.prec, ctx.threads = 128, 1
    with np.load(HERE/"public"/"operators"/"n96_f400_a08.npz") as archived:
        data = dict(archived)
    stage = perf_counter()
    kernel, bound = exact_bare_representative(data["spectrum"])
    seconds_kernel = perf_counter()-stage
    with TemporaryDirectory() as tmp:
        dylib = Path(tmp)/"explicit_dots.dylib"
        cmd = ["clang", "-O3", "-fno-fast-math", "-ffp-contract=off", "-dynamiclib",
               "-o", str(dylib), str(HERE/"certify_dot.c")]
        subprocess.run(cmd, check=True)
        dots = ExplicitDots(dylib, args.threads)
        summaries = []
        for index in args.items:
            row = certify(f"n96_c{index:02}", data, kernel, bound, dots)
            summaries.append(row)
            print(row, flush=True)
    path = HERE/"results"/"native_forward_arb.json"
    merged = {item["case_id"]: item for item in json.loads(path.read_text())["cases"]} if path.exists() else {}
    merged.update({item["case_id"]: item for item in summaries})
    path.write_text(json.dumps(dict(precision_Arb=ctx.prec, exact_operator="binary64 stored spectrum + network factors",
                                    bare_kernel_wall_seconds=seconds_kernel,
                                    bare_kernel_complex_error_bound=outward(bound)[1],
                                    explicit_binary64_product_threads=args.threads,
                                    cases=list(merged.values())), indent=2))


if __name__ == "__main__":
    main()
