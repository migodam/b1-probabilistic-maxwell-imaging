"""Upper-bound the native certificate by restricting its bordered Hermitian pencil.

The trial space consists of normalized indicators of uniform tiles on the 96-grid.
All accepted signs come from midpoint-radius arithmetic on the original symbol.
The optimization grid is never used to prove an upper or a lower bound here.
"""

import sys
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

PAPER = Path(__file__).resolve().parents[2] / "z_alg_optical_inertia"
sys.path.insert(0, str(PAPER))
from core.certificate import ball, gershgorin, safe_float, star


def midpoint(matrix):
    return np.fromiter((complex(x) for x in matrix.entries()), complex).reshape(
        matrix.nrows(), matrix.ncols()
    )


def block_sums(x, norm):
    side = int(np.sqrt(x.shape[0]))
    tiling = side // norm
    a = ball(x)
    cols = a.ncols()
    vals = [
        sum(
            (
                a[(ty * norm + y) * side + tx * norm + j, c]
                for y in range(norm)
                for j in range(norm)
            ),
            acb(0),
        )
        / norm
        for ty in range(tiling)
        for tx in range(tiling)
        for c in range(cols)
    ]
    return acb_mat(tiling**2, cols, vals)


class NativeRestriction:
    def __init__(self, tiling=16):
        started, cpu = perf_counter(), process_time()
        self.tiling = tiling
        ctx.prec = 112
        with np.load(
            PAPER / "synthetic/public/fine_geometry.npz", allow_pickle=False
        ) as f:
            symbol, sensor, incident = (
                f[k] for k in ("spectrum", "PW360_10x20__S", "PW360_10x20__B")
            )
        n = int(np.sqrt(incident.shape[0]))
        p = len(symbol)
        q = n // tiling
        if n % tiling or tiling % 8:
            raise ValueError("Trial grid must align with the two partitions.")
        self.s = incident.shape[1]
        self.sensor = block_sums(sensor.T, q).transpose()
        self.incident = block_sums(incident, q)
        horizontal = [
            acb.dft([acb(float(z.real), float(z.imag)) for z in row], inverse=True)
            for row in symbol
        ]
        spatial = [acb.dft([r[x] for r in horizontal], inverse=True) for x in range(p)]
        table = {}
        for dx in range(1 - tiling, tiling):
            for dy in range(1 - tiling, tiling):
                table[dy, dx] = sum(
                    (
                        (q - abs(i))
                        * (q - abs(j))
                        * spatial[(q * dx + j) % p][(q * dy + i) % p]
                        for i in range(1 - q, q)
                        for j in range(1 - q, q)
                    ),
                    acb(0),
                ) / (q * q)
        positions = list(np.ndindex(tiling, tiling))
        k = tiling**2
        self.green = acb_mat(
            k, k, [table[i - v, j - w] for i, j in positions for v, w in positions]
        )
        self.mid_g = midpoint(self.green)
        self.mid_s, self.mid_b = midpoint(self.sensor), midpoint(self.incident)
        self.ss = self.mid_s.conj().T @ self.mid_s
        self.setup = dict(
            wall_s=perf_counter() - started,
            cpu_s=process_time() - cpu,
            trial_dim=k,
            native_dim=n * n,
            precision=ctx.prec,
        )

    def check(self, y, u, mu, lam, radius):
        started, cpu = perf_counter(), process_time()
        mu = np.asarray(mu).flatten()
        lam = np.asarray(lam).flatten()
        w = int(np.sqrt(len(mu)))
        k = self.tiling**2
        if len(mu) != w * w or self.tiling % w or len(lam) != len(mu) or np.min(mu) < 0:
            raise ValueError("Invalid block profile.")
        step = self.tiling // w
        ix = np.repeat(
            np.repeat(np.arange(w * w).reshape(w, w), step, axis=0), step, axis=1
        ).ravel()
        mm = mu[ix]
        ll = lam[ix]
        z = (u * mm + 1j * ll) / 2
        scale = max(float(np.mean(mu)), 1e-16)
        root = np.sqrt(scale)
        with threadpool_limits(limits=1, user_api="blas"):
            zg = z[:, None] * self.mid_g
            K = self.ss - zg - zg.conj().T
            K.flat[:: k + 1] += mm
            L = self.mid_s.conj().T @ y + z[:, None] * self.mid_b
            block = np.block(
                [[K / scale, -L / root], [-L.conj().T / root, y.conj().T @ y]]
            )
            block = (block + block.conj().T) / 2
            yb = ball(y)
            ab_mu = [arb(float(x)) for x in mu]
            ab_z = [
                acb(arb(float(u)) * x, arb(float(v))) / 2 for x, v in zip(ab_mu, lam)
            ]
            checks = []

            def negative_at(a, state_only=False):
                matrix = block[:k, :k].copy() if state_only else block.copy()
                if not state_only:
                    matrix[k:, k:] -= a * np.eye(self.s)
                vals, vectors = eigh(
                    matrix,
                    subset_by_index=(0, self.s - 1),
                    driver="evr",
                    check_finite=False,
                )
                count = int(np.count_nonzero(vals < -1.0e-13))
                if not count:
                    return 0
                x = ball(vectors[:k, :count] / root)
                b = acb_mat(self.s, count) if state_only else ball(vectors[k:, :count])
                res = self.sensor * x - yb * b
                total = self.green * x + self.incident * b
                mx = acb_mat(
                    k,
                    count,
                    [
                        ab_mu[idx] * x[i, j]
                        for i, idx in enumerate(ix)
                        for j in range(count)
                    ],
                )
                ze = acb_mat(
                    k,
                    count,
                    [
                        ab_z[idx] * total[i, j]
                        for i, idx in enumerate(ix)
                        for j in range(count)
                    ],
                )
                cross = star(x) * ze
                gram = (
                    star(res) * res
                    + star(x) * mx
                    - cross
                    - star(cross)
                    - arb(float(a)) * star(b) * b
                )
                while count:
                    part = acb_mat(
                        count,
                        count,
                        [gram[i, j] for i in range(count) for j in range(count)],
                    )
                    edge = gershgorin(part, upper=True)
                    if edge < 0:
                        checks.append(
                            dict(
                                a=a,
                                state_only=state_only,
                                negative_at_least=count,
                                congruence_upper=safe_float(edge, upper=True),
                            )
                        )
                        return count
                    count -= 1
                return 0

            q = self.s - negative_at(0.0)
            upper = None
            negative_state_lower = None
            if q == 0:
                upper = 0.0
                negative_state_lower = negative_at(0.0, state_only=True)
            else:
                # q bounds the count above zero, independently of K having a spectral gap.
                a = float(np.nextafter(radius**2 / q * (1 - 1.0e-11), 0.0))
                if a > 0.0 and negative_at(a) == self.s:
                    candidate = safe_float((arb(q) * arb(a)).sqrt(), upper=True)
                    if candidate <= radius:
                        upper = candidate
        return dict(
            status="PROVED_ZERO"
            if upper == 0
            else "PROVED_BELOW_RADIUS"
            if upper is not None
            else "INCONCLUSIVE",
            lower=0.0,
            upper=upper,
            positive_terms_at_most=q,
            nu_lower=negative_state_lower,
            tiling=self.tiling,
            checks=checks,
            wall_s=perf_counter() - started,
            cpu_s=process_time() - cpu,
        )
