"""Native outward certificates for arbitrary row-major 4x4 and 8x8 multiplier profiles.

Uses the archived Paper 1 native actions, stationarity defect and Fourier envelope.
Only fixed geometry is loaded here. The caller owns observations and their radii.
"""

from collections import OrderedDict
from contextlib import contextmanager
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import eigh, LinAlgError
from scipy.sparse.linalg import ArpackError
from threadpoolctl import threadpool_limits


PAPER = Path(__file__).resolve().parents[2] / "z_alg_optical_inertia"
CURVE = PAPER / "submission/curve"
GEOMETRY = PAPER / "synthetic/public/fine_geometry.npz"
CAPS = (.5, .75, 1.)


def _paper_modules():
    # The archived scripts use local imports with generic names. Load a private
    # copy of each module, then restore the caller's import paths and bindings.
    bindings = {name: sys.modules.get(name) for name in ("curve_inputs", "search", "native")}
    path = sys.path[:]
    old_bytecode = sys.dont_write_bytecode

    def load(name, source):
        spec = spec_from_file_location(name, CURVE / (source + ".py"))
        mod = module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
        return mod

    try:
        sys.dont_write_bytecode = True
        load("curve_inputs", "curve_inputs")
        load("search", "search")
        native = load("_zalg_paper_native", "native")
        sys.modules["native"] = native
        return native, load("_zalg_paper_envelope", "envelope")
    finally:
        sys.dont_write_bytecode = old_bytecode
        sys.path[:] = path
        for name, value in bindings.items():
            if value is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = value


_paper, _envelope = _paper_modules()
Weights, ball, star = _paper.Weights, _paper.ball, _paper.star
gershgorin, high, low = _paper.gershgorin, _paper.high, _paper.low


def hermitian(a):
    return (a + a.conj().T) / 2


@contextmanager
def timed(receipt, label):
    wall, cpu = perf_counter(), process_time()
    try:
        yield
    finally:
        receipt[label] = {
            "wall_seconds": perf_counter() - wall,
            "cpu_seconds": process_time() - cpu,
        }


class BoundedEnvelope(_envelope.FourierEnvelope):
    """Use the same exact lower envelope, with a bounded low-memory trial search."""

    def __init__(self, geometry, max_cutoff, max_eigensolves):
        self.g = geometry
        self.roots = [acb(arb(2) * k / geometry.pad).exp_pi_i() for k in range(geometry.pad)]
        self.soph, self.dirichlet, self.numerical_grams = {}, {}, {}
        self.ss = geometry.sensor_ball * star(geometry.sensor_ball)
        self.precise_sensor_blocks = None
        width = geometry.side // geometry.block_side
        phase = 2j * np.pi * np.arange(geometry.pad) / geometry.pad
        base = np.exp(phase[:, None] * np.arange(width)[None, :]).sum(axis=1)
        self.numerical_dirichlet = np.stack([
            base * np.exp(phase * j * width) for j in range(geometry.block_side)
        ])
        self.numerical_ss_blocks = np.stack([
            geometry.sensor[:, geometry.blocks == j] @ geometry.sensor_star[geometry.blocks == j]
            for j in range(geometry.block_count)
        ])
        self.max_cutoff = max_cutoff
        self.max_eigensolves = max_eigensolves
        self.trial_eigensolves = 0

    def _port_modes(self, modes):
        g = self.g
        width = g.side // g.block_side
        phase = np.exp(-2j * np.pi * modes[:, None] * np.arange(g.side)[None, :] / g.pad)
        fields = g.sensor_star.reshape(g.side, g.side, -1)
        projected = []
        for block in range(g.block_count):
            row, col = divmod(block, g.block_side)
            ys, xs = slice(row * width, (row + 1) * width), slice(col * width, (col + 1) * width)
            local = np.einsum("ay,yzr,bz->abr", phase[:, ys], fields[ys, xs],
                              phase[:, xs], optimize=True)
            projected.append(local.reshape(len(modes)**2, len(g.sensor)) / g.pad)
        return np.stack(projected)

    def _mode_data(self, cut):
        if cut not in self.numerical_grams:
            g = self.g
            modes = np.arange(-cut, cut + 1) % g.pad
            y, x = np.meshgrid(modes, modes, indexing="ij")
            iy, ix = y.ravel(), x.ravel()
            omitted = np.ones(g.spec.shape, dtype=bool)
            omitted[iy, ix] = False
            deltas = np.arange(-2 * cut, 2 * cut + 1) % g.pad
            table = self.numerical_dirichlet[:, deltas]
            # Build each Dirichlet coefficient once for each mode difference,
            # never as a separate mode-by-mode dense matrix for every block.
            block_diffs = np.einsum("yi,xj->yxij", table, table).reshape(
                g.block_count, len(deltas), len(deltas)
            ) / g.pad**2
            dy = np.tile(np.repeat(np.arange(2 * cut + 1), 2 * cut + 1), (len(iy), 1))
            dx = np.tile(np.tile(np.arange(2 * cut + 1), 2 * cut + 1), (len(iy), 1))
            dy = dy - dy.T + 2 * cut
            dx = dx - dx.T + 2 * cut
            self.numerical_grams[cut] = (
                iy, ix, omitted, dy, dx, block_diffs, self._port_modes(modes)
            )
        return self.numerical_grams[cut]

    def _proposal(self, w, cut, beta, shift, gap):
        g = self.g
        iy, ix, omitted, dy, dx, block_diffs, ports = self._mode_data(cut)
        t, z, diagonal = self.numeric_coeffs(w, 0., beta, shift)
        iso = w.is_global
        symbol = 1 - w.cap * g.spec.real + t * g.spec.imag
        cross = g.spec - shift
        score = symbol if iso else symbol - abs(cross)**2 / beta
        floor = float(np.floor(256 * score[omitted].min()) / 256)
        alpha = floor + diagonal - gap
        kept = score[iy, ix] - floor
        if min(alpha) <= 1.e-9 or min(abs(kept)) < 1.e-11:
            return None
        return alpha, kept, symbol[iy, ix] - floor, cross[iy, ix], z, floor

    def _matrix(self, cut, proposal, w, beta):
        _, _, _, dy, dx, block_diffs, ports = self._mode_data(cut)
        alpha, kept, symbol, cross, z, _ = proposal
        inv = 1 / alpha

        def gram(coefficients):
            table = np.einsum("a,axy->xy", coefficients, block_diffs)
            return table[dy, dx]

        r = len(kept)
        idx = np.arange(r)
        receiver = np.einsum("a,aij->ij", inv, ports)
        ss = np.einsum("a,aij->ij", inv, self.numerical_ss_blocks)
        ss += w.scale * np.eye(len(ss))
        ff = gram(inv)
        if w.is_global:
            ff[idx, idx] += 1 / kept
            return np.block([[ff, receiver], [receiver.conj().T, ss]])
        fz = gram(inv * z)
        zz = gram(inv * abs(z)**2)
        pz = np.einsum("a,aij->ij", inv * z.conj(), ports)
        det = beta * kept
        ff[idx, idx] += beta / det
        zz[idx, idx] += symbol / det
        fz[idx, idx] += cross.conj() / det
        return np.block([[ff, fz, receiver], [fz.conj().T, zz, pz],
                         [receiver.conj().T, pz.conj().T, ss]])

    def select(self, w, nu):
        self.trial_eigensolves = 0
        shifts = (0., .25, .375, .5, .625)
        betas = (1., 2., 4., 8., 16., 32., 64., 128., 256.)
        settings = ((1., 0.),) if w.is_global else tuple(
            (beta, shift) for shift in shifts for beta in betas
        )
        cuts = tuple(c for c in (6, 8, 10, 12, 13) if c <= self.max_cutoff)
        for gap in (.03125, .0078125, .001953125, .0009765625):
            for cut in cuts:
                candidates = []
                for beta, shift in settings:
                    proposal = self._proposal(w, cut, beta, shift, gap)
                    if proposal is not None:
                        need = max(0, int(np.count_nonzero(proposal[1] < 0)) - nu)
                        strength = float(min(proposal[0]))
                        candidates.append((need, -strength, beta, shift, proposal))
                candidates.sort(key=lambda entry: entry[:4])
                for need, _, beta, shift, proposal in candidates[:2]:
                    if need:
                        if self.trial_eigensolves >= self.max_eigensolves:
                            raise ArithmeticError("Fourier proposal eigensolve budget exhausted.")
                        matrix = self._matrix(cut, proposal, w, beta)
                        self.trial_eigensolves += 1
                        val = eigh(hermitian(matrix), eigvals_only=True, check_finite=False,
                                   subset_by_index=(need - 1, need - 1))[0]
                        if val >= -1.e-8:
                            continue
                    iy, ix, *_ = self._mode_data(cut)
                    retained = tuple(zip(iy.tolist(), ix.tolist()))
                    return proposal[-1], beta, shift, gap, retained, w.is_global, need
        raise ArithmeticError("Fourier lower envelope did not propose the native spectral gap.")


class RegionalWitness:
    """An interval compression of each regional weight on a shared native basis."""

    def __init__(self, geometry, vectors, name):
        self.name = name
        basis = ball(vectors)
        self.gram = star(basis) * basis
        sensor = geometry.sensor_ball * basis
        self.receiver = star(sensor) * sensor
        green = geometry.exact_green(basis)
        self.regional = []
        rank = basis.ncols()
        for block in range(geometry.block_count):
            pixels = np.flatnonzero(geometry.blocks == block)
            sub = acb_mat(len(pixels), rank, [basis[int(i), j] for i in pixels for j in range(rank)])
            prop = acb_mat(len(pixels), rank, [green[int(i), j] for i in pixels for j in range(rank)])
            product = star(sub) * prop
            sector = (product + star(product)) / 2
            optical = (product - star(product)) / acb(0, 2)
            self.regional.append((star(sub) * sub, sector, optical))

    def negative_gap(self, w, required):
        projected = self.receiver
        cap, scale = arb(w.cap), arb(w.scale)
        for mu, lam, (gram, sector, optical) in zip(w.real, w.optical, self.regional, strict=True):
            projected = projected + arb(mu) * (gram - cap * sector) + arb(lam) * optical
        projected = projected / scale
        _, vectors = eigh(hermitian(_paper.mid(projected)), subset_by_index=(0, required - 1),
                          check_finite=False)
        basis = ball(vectors)
        negative = gershgorin(star(basis) * (projected * basis), upper=True)
        if negative < 0:
            gram_upper = gershgorin(star(basis) * (self.gram * basis), upper=True)
            return low(-negative / gram_upper)
        return None


class Native:
    """Score a fixed physical profile on the full 96² operator using interval arithmetic."""

    def __init__(self, *, precision=112, threads=2, max_cutoff=13, max_eigensolves=10,
                 max_profiles=192, max_witnesses=3):
        if not (6 <= max_cutoff <= 13 and max_eigensolves >= 1 and max_profiles >= 1
                and max_witnesses >= 1 and precision >= 80 and threads >= 1):
            raise ValueError("The verifier budget or interval precision is invalid.")
        self.precision, self.threads = precision, threads
        self.max_cutoff, self.max_eigensolves = max_cutoff, max_eigensolves
        self.max_profiles, self.max_witnesses = max_profiles, max_witnesses
        self.profiles = OrderedDict()
        self.geometries, self.envelopes, self.witnesses = {}, {}, {}
        self.operations = {"profile_hits": 0, "weyl_hits": 0, "witness_hits": 0}
        self.setup = {}
        with timed(self.setup, "load_fixed_geometry"):
            with np.load(GEOMETRY, allow_pickle=False) as data:
                self.spectrum = data["spectrum"]
                self.sensor = data["PW360_10x20__S"]
                self.incident = data["PW360_10x20__B"]
        assert self.sensor.shape == (20, 9216) and self.incident.shape == (9216, 10)
        self.sources = self.incident.shape[1]
        self._green_bounds = None
        self.witness_sources = {}

    @contextmanager
    def _precision(self):
        old_prec, old_threads = ctx.prec, ctx.threads
        ctx.prec, ctx.threads = self.precision, self.threads
        try:
            with threadpool_limits(limits=1):
                yield
        finally:
            ctx.prec, ctx.threads = old_prec, old_threads

    def _geometry(self, side, stages):
        if side not in self.geometries:
            with timed(stages, "partition_geometry"):
                self.geometries[side] = _paper.Geometry(
                    self.spectrum, self.sensor, self.incident, block_side=side
                )
                self.witnesses[side] = OrderedDict()
        return self.geometries[side]

    def _green_norm_bounds(self):
        if self._green_bounds is None:
            self._green_bounds = []
            for center in (0., .375, .75, 1.5):
                c = arb(center)
                norm2 = max(
                    ((arb(float(z.real)) - c)**2 + arb(float(z.imag))**2).upper()
                    for z in self.spectrum.flat
                )
                self._green_bounds.append((c, norm2.sqrt().upper()))
        return self._green_bounds

    @staticmethod
    def _expand(weights, count):
        if len(weights) == count:
            return np.asarray(weights)
        return np.asarray(weights).reshape(4, 4).repeat(2, axis=0).repeat(2, axis=1).ravel()

    def _weyl(self, w):
        anchors = [key for key, proof in reversed(self.profiles.items())
                   if proof["cap"] == w.cap
                   and (proof.get("positive_gap") is not None or proof.get("nu_lower", 0) >= self.sources)]
        if not anchors:
            return None
        approx = []
        for key in anchors:
            count = max(len(w.real), len(key.real))
            w_mu = self._expand(w.real, count)
            k_mu = self._expand(key.real, count)
            w_lam = self._expand(w.optical, count)
            k_lam = self._expand(key.optical, count)
            diff = (w.cap * (w_mu - k_mu) + 1j * (w_lam - k_lam)) / 2.
            approx.append((float(max(abs(diff))), key, w_mu, k_mu, w_lam, k_lam))
        for _, key, w_mu, k_mu, w_lam, k_lam in sorted(approx, key=lambda entry: entry[0])[:8]:
            proof = self.profiles[key]
            scale = arb(w.scale)
            ratio = arb(key.scale) / scale
            cap = arb(w.cap)
            dmu = [(arb(float(x)) - arb(float(y))) / scale for x,y in zip(w_mu,k_mu,strict=True)]
            dlam = [(arb(float(x)) - arb(float(y))) / scale for x,y in zip(w_lam,k_lam,strict=True)]
            dz = [acb(cap * x, y) / 2 for x,y in zip(dmu,dlam,strict=True)]
            cross = max(abs(x).upper() for x in dz)
            positive_error, negative_error = None, None
            for center, norm in self._green_norm_bounds():
                d = [x - 2 * center * z.real for x,z in zip(dmu,dz,strict=True)]
                interaction = 2 * cross * norm
                pos = (interaction - min(x.lower() for x in d)).upper()
                neg = (interaction + max(x.upper() for x in d)).upper()
                positive_error = pos if positive_error is None else min(pos, positive_error)
                negative_error = neg if negative_error is None else min(neg, negative_error)
            gap = proof.get("positive_gap")
            neg = proof.get("negative_gap")
            new_pos = ratio * arb(gap) - positive_error if gap is not None else None
            new_neg = ratio * arb(neg) - negative_error if neg is not None else None
            witness = new_neg is not None and new_neg > 0
            full = new_pos is not None and new_pos > 0 and (proof["nu"] == 0 or witness)
            exhausted = proof["nu_lower"] >= self.sources and witness
            if full or exhausted:
                self.operations["weyl_hits"] += 1
                return {
                    "method": "interval_weyl",
                    "cap": w.cap,
                    "nu": proof["nu"] if full else None,
                    "nu_lower": proof["nu_lower"],
                    "positive_gap": low(new_pos) if full else None,
                    "negative_gap": low(new_neg) if witness else None,
                    "anchor": key.as_dict(),
                    "lower_shift_cost_upper": high(positive_error),
                    "upper_shift_cost_upper": high(negative_error),
                }
        return None

    def _shared_zero(self, w, side):
        for name, basis in reversed(self.witnesses[side].items()):
            gap = basis.negative_gap(w, self.sources)
            if gap is not None and gap > 0:
                self.operations["witness_hits"] += 1
                return {"method": "regional_negative_compression", "cap": w.cap,
                        "nu": None, "nu_lower": self.sources, "positive_gap": None,
                        "negative_gap": gap, "witness": name}
        return None

    def _save_witness(self, geometry, w, action, stages):
        side = geometry.block_side
        pool = self.witnesses[side]
        if len(pool) >= self.max_witnesses:
            return
        name = f"{side}x{side}-basis-{len(pool) + 1}"
        with timed(stages, "build_regional_witness"):
            pool[name] = RegionalWitness(geometry, action.eig_basis, name)
        self.witness_sources[name] = w.as_dict()

    def _profile(self, w, geometry, action, stages, nu_hint):
        if w in self.profiles:
            self.operations["profile_hits"] += 1
            self.profiles.move_to_end(w)
            receipt = dict(self.profiles[w])
            receipt["cache_hit"] = True
            return receipt
        with timed(stages, "verify_cached_geometry"):
            proof = self._weyl(w)
            if proof is None:
                proof = self._shared_zero(w, geometry.block_side)
        if proof is not None:
            return proof
        with timed(stages, "native_binary64_inertia_proposal"):
            modes = action.modes(max(self.sources, nu_hint or 0))
        proposed_nu = action.nu
        if proposed_nu:
            verify_count = min(proposed_nu, self.sources)
            with timed(stages, "verify_native_negative_subspace"):
                negative_upper, negative_gap = action.ensure_negative(verify_count, return_gap=True)
        else:
            negative_upper, negative_gap = None, None
        proof = {"cap": w.cap, "nu": None, "nu_lower": min(proposed_nu, self.sources),
                 "positive_gap": None, "negative_gap": negative_gap,
                 "native_nu_proposal": proposed_nu,
                 "native_scaled_K_eigenvalues": modes["numerical_smallest_scaled_K"],
                 "native_scaled_K_defect": modes["numerical_eigen_defect"],
                 "negative_block_upper": negative_upper}
        if proposed_nu >= self.sources:
            proof["method"] = "native_negative_subspace"
            if self.witnesses[geometry.block_side] or w.cap == 1.:
                self._save_witness(geometry, w, action, stages)
            return proof
        with timed(stages, "verify_Fourier_lower_envelope"):
            side = geometry.block_side
            if side not in self.envelopes:
                self.envelopes[side] = BoundedEnvelope(
                    geometry, self.max_cutoff, self.max_eigensolves
                )
            envelope = self.envelopes[side]
            try:
                fourier = envelope.verify(w, proposed_nu)
            except (ArithmeticError, LinAlgError, ValueError) as error:
                proof.update(method="unresolved_positive_gap", explanation=str(error),
                             envelope_proposal_eigensolves=envelope.trial_eigensolves)
                return proof
            proof.update(method="native_Fourier_and_negative_subspace", nu=proposed_nu,
                         positive_gap=fourier["positive_gap_lower"],
                         envelope_proposal_eigensolves=envelope.trial_eigensolves,
                         envelope=fourier)
        return proof

    def _weights(self, u, mu, lam):
        u = float(u)
        mu, lam = np.asarray(mu, dtype=float), np.asarray(lam, dtype=float)
        if u not in CAPS:
            raise ValueError("Expected a cap of 0.5, 0.75, or 1.")
        shapes = {(16,), (64,), (4, 4), (8, 8)}
        if mu.shape not in shapes or lam.shape not in shapes or mu.size != lam.size:
            raise ValueError("Expected matching 16 or 64 row-major block weights.")
        if not np.isfinite(mu).all() or not (mu > 0).all() or not np.isfinite(lam).all():
            raise ValueError("Sector weights must be finite and positive; optical weights finite.")
        scale = float(mu.flat[0]) if np.all(mu == mu.flat[0]) else float(np.mean(mu))
        if not np.isfinite(scale) or scale <= 0.:
            raise ValueError("The positive normalization scale is not representable.")
        return Weights(u, scale, tuple(map(float, mu.flat)), tuple(map(float, lam.flat)))

    def score(self, Y, u, mu, lam, *, radius=None, nu_hint=None):
        """Return native and interval evidence; UNKNOWN never means a proved zero."""
        w = self._weights(u, mu, lam)
        y = np.asarray(Y, dtype=complex)
        if y.shape != (20, 10) or not np.isfinite(y).all():
            raise ValueError("A finite 20-receiver by 10-illumination observation is required.")
        if radius is not None and (not np.isfinite(radius) or radius < 0.):
            raise ValueError("The radius must be finite and nonnegative.")
        if nu_hint is not None and int(nu_hint) != nu_hint:
            raise ValueError("The numerical inertia hint must be an integer.")
        side = int(np.sqrt(len(w.real)))
        stages = {}
        result = {"status": "UNKNOWN", "lower": 0., "upper": None, "nu": None, "nu_lower": 0,
                  "radius": None if radius is None else float(radius), "strictly_excludes": False,
                  "native_numeric": None, "profile": w.as_dict(), "proof": None,
                  "precision_bits": self.precision, "interval_threads": self.threads, "timing": stages}
        with timed(stages, "total"), self._precision():
            geometry = self._geometry(side, stages)
            start_counters = geometry.cost.copy()
            try:
                with timed(stages, "profile_operator_setup"):
                    action = _paper.Native(geometry, w)
                proof = self._profile(w, geometry, action, stages, nu_hint)
                result.update(proof=proof, nu=proof["nu"], nu_lower=proof["nu_lower"])
                if not proof.get("cache_hit"):
                    self.profiles[w] = proof
                    if len(self.profiles) > self.max_profiles:
                        self.profiles.popitem(last=False)
                if proof["nu_lower"] >= self.sources:
                    result.update(status="PROVED_ZERO", upper=0.)
                elif proof["nu"] is not None:
                    with timed(stages, "native_binary64_stationarity_proposal"):
                        trial, vectors, vals, numeric = action.numerical_distance(y, proof["nu"])
                    result["native_numeric"] = numeric
                    with timed(stages, "native_Arb_stationarity_and_source_spectrum"):
                        interval = action.ensure_source(
                            y, proof["nu"], proof["positive_gap"], trial, vectors, vals,
                            negative_gap=proof["negative_gap"] if proof["nu"] else proof["positive_gap"]
                        )
                    lower, upper = interval["distance_lower"], interval["distance_upper"]
                    perturbation = arb(interval["upper_form_perturbation"])
                    if all(arb(float(v)) + perturbation <= 0. for v in vals[proof["nu"]:]):
                        lower, upper = 0., 0.
                        interval.update(distance_lower=lower, distance_upper=upper)
                    if lower > upper:
                        raise ArithmeticError("The verified source interval has inconsistent endpoints.")
                    result.update(lower=lower, upper=upper, source_interval=interval,
                                  strictly_excludes=bool(radius is not None and lower > radius))
                    if lower > 0.:
                        result["status"] = "CERTIFIED_POSITIVE"
                    elif upper == 0.:
                        result["status"] = "PROVED_ZERO"
                    else:
                        result["explanation"] = "No positive source lower endpoint or exact zero proved."
                else:
                    result["explanation"] = proof["explanation"]
            except (ArithmeticError, ArpackError, LinAlgError, OverflowError, ValueError) as error:
                result["explanation"] = str(error)
            result["operations"] = dict(geometry.cost - start_counters)
        return result

    def cache_summary(self):
        return {"setup": self.setup, "operation_totals": self.operations.copy(),
                "resident_proofs": len(self.profiles), "witness_sources": self.witness_sources.copy(),
                "configured_proof_budget": {"max_cutoff": self.max_cutoff,
                                            "max_proposal_eigensolves": self.max_eigensolves,
                                            "max_profiles": self.max_profiles,
                                            "max_witnesses_per_partition": self.max_witnesses}}
