"""Reciprocal scalar-TM line-port model with fixed terminated Green operator."""
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from scipy.fft import fftn, ifftn
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import LinearOperator, gmres
from scipy.special import hankel1, jv

MODEL_PATH = Path(__file__).resolve()
P = MODEL_PATH.parents[1]
ROOT = MODEL_PATH.parents[4] if len(MODEL_PATH.parents) > 4 else P
LEGACY = ROOT / "research/delegated/max_struct_r1/r8_worktree"
R8_ENTRY = Path("research/research_max_struct_r8/src/paper_physics.py")
for source in (LEGACY, P / "vendor_r8", ROOT):
    if (source / R8_ENTRY).is_file():
        sys.path.insert(0, str(source))
        break
from research.research_max_struct_r8.src.paper_physics import EqualAreaGreen

C = 3e8
R = 50.0
MAIN_STRENGTH = 48.0
COUNTS = (6, 8, 12, 16)
BANDS = (400, 460)
STRENGTHS = (16.0, MAIN_STRENGTH, 144.0)


class DiskGreen:
    """Use the existing disk-pulse operator, with the same native cell ordering."""

    def __init__(self, side, band):
        self.shape = (side, side)
        self.n = side**2
        self.band = band
        axis = (np.arange(side) + 0.5) * (2 / side) - 1
        yy, xx = np.meshgrid(axis, axis, indexing="ij")
        self.centers = np.column_stack((xx.ravel(), yy.ravel()))
        self.k = (2 * np.pi / C) * (band * 1e6)
        self.radius = (2 / side) / np.sqrt(np.pi)
        self.pulse = 2 * np.pi * self.k * self.radius * jv(1, self.k * self.radius)
        self.kernel = EqualAreaGreen(self.centers, self.shape, self.radius, self.k)

    def dense(self):
        distance = np.linalg.norm(self.centers[:, None] - self.centers[None, :], axis=-1)
        np.fill_diagonal(distance, 1)
        result = self.pulse * point(distance, self.k)
        np.fill_diagonal(result, self.kernel.diagonal)
        return result


def point(distance, k):
    return 0.25j * hankel1(0, k * distance)


def aperture(count):
    """Separate transmit and receive arcs with identical endpoints at all sizes."""

    def arc(radius, angles):
        angles = np.deg2rad(angles)
        return radius * np.column_stack((np.cos(angles), np.sin(angles)))

    tx = arc(3.75, np.linspace(95, 265, count))
    rx = arc(3.75, np.linspace(-85, 85, count))
    hold = arc(3.95, -85 + 170 * (np.arange(8) + 0.5) / 8)
    return np.concatenate((tx, rx, hold))


def port_point_matrices(disk, positions):
    """Point-to-pulse excitation, pulse-integrated readout, and point transfer."""
    distance = np.linalg.norm(disk.centers[:, None] - positions[None, :], axis=-1)
    excitation = point(distance, disk.k)
    port_distance = np.linalg.norm(positions[:, None] - positions[None, :], axis=-1)
    np.fill_diagonal(port_distance, 1)
    self_regularized = point(port_distance, disk.k)
    np.fill_diagonal(self_regularized, 0.25j)
    return excitation, disk.pulse * excitation.T, self_regularized


def network_terms(disk, count, strength, available_power=1.0, variation=None):
    positions = aperture(count)
    first_receive = count
    first_hold = 2 * count
    np_ = len(positions)
    excitation, readout, transfer = port_point_matrices(disk, positions)
    normalized_index = np.arange(np_)
    calibration = np.ones(np_)
    ohmic = np.full(np_, 8.0)
    reactance = 12 * (disk.band / 400 - 400 / disk.band) + 2 * np.cos(normalized_index)
    if variation is not None:
        calibration = calibration + variation["gain"]
        ohmic = ohmic + variation["resistance"]
        reactance = reactance + variation["reactance"]
    gain = np.sqrt(strength) * calibration
    zrad = -1j * gain[:, None] * transfer * gain[None, :]
    antenna = zrad + np.diag(ohmic + 1j * reactance)
    matrix = antenna + R * np.eye(np_)
    drive = (2 * np.sqrt(R * available_power / count)) * np.eye(np_, count)
    inv = np.linalg.solve(matrix, np.eye(np_))
    port_excitation = excitation * gain[None, :]
    port_readout = gain[:, None] * readout
    return dict(positions=positions, drive=drive, matrix=matrix, inv=inv, antenna=antenna,
                gain=gain, ohmic=ohmic, reactance=reactance,
                rx=np.arange(first_receive, first_hold), hold=np.arange(first_hold, np_),
                excitation=excitation, readout=readout, port_excitation=port_excitation,
                port_readout=port_readout, available_power=available_power)


def operator(disk, network, dense=None):
    E_bare = network["port_excitation"]
    S_bare = network["port_readout"]
    inverse = network["inv"]
    drive = network["drive"]
    left = 1j * (E_bare @ inverse)
    full_receive = -1j * np.sqrt(R) * (inverse @ S_bare)
    diagonal = disk.kernel.diagonal + np.einsum("ij,ji->i", left, S_bare)
    out = dict(E0=E_bare @ (inverse @ drive), S=full_receive[network["rx"]],
               S_hold=full_receive[network["hold"]], E0_bare=E_bare, S_bare=S_bare,
               low_rank_left=left, low_rank_right=S_bare,
               spectrum=disk.kernel.spectrum, pad=np.asarray(disk.kernel.pad),
               diagonal=disk.kernel.diagonal, diagonal_dressed=diagonal,
               grid_shape=np.asarray(disk.shape), centers=disk.centers, k=disk.k,
               disk_radius=disk.radius, pulse=disk.pulse, frequency_MHz=disk.band,
               M_nominal=network["matrix"], H=drive, Z_AA=network["antenna"],
               W=0.5*(network["antenna"] + network["antenna"].conj().T),
               receive_index=network["rx"], hold_index=network["hold"],
               source_index=np.arange(drive.shape[1]), positions=network["positions"],
               reference_ohms=R, gain=network["gain"],
               network_scale=network["gain"][0]**2, available_power=network["available_power"])
    if dense is not None:
        out["G"] = dense + left @ S_bare
    return out


class Green:
    """Exact FFT + low-rank fixed-G action and adjoint for released operators."""

    def __init__(self, operator_data):
        self.shape = tuple(operator_data["grid_shape"])
        self.pad = tuple(operator_data["pad"])
        self.spectrum = operator_data["spectrum"]
        self.diagonal = operator_data["diagonal_dressed"]
        self.left = operator_data["low_rank_left"]
        self.right = operator_data["low_rank_right"]

    def apply(self, J, adjoint=False):
        unstack = np.ndim(J) == 1
        if unstack:
            J = J[:, None]
        ny, nx = self.shape
        stack = np.zeros(self.pad + (J.shape[1],), complex)
        stack[:ny, :nx] = J.reshape(self.shape + (J.shape[1],))
        spectrum = self.spectrum.conj() if adjoint else self.spectrum
        padded = ifftn(fftn(stack, axes=(0, 1)) * spectrum[..., None], axes=(0, 1))
        result = padded[:ny, :nx].reshape(J.shape)
        result += self.right.conj().T @ (self.left.conj().T @ J) if adjoint else self.left @ (self.right @ J)
        return result[:, 0] if unstack else result

    def adjoint(self, J):
        return self.apply(J, adjoint=True)


def physical_solve(material, green, E, diagonal, dense=None):
    """Solve the unmodified full-wave current equation with independent residuals."""
    start = perf_counter(), process_time()
    m = np.asarray(material).ravel()
    B = m[:, None] * E
    if dense is not None:
        A = np.eye(len(m)) - m[:, None] * dense
        lu = lu_factor(A, overwrite_a=True, check_finite=False)
        J = lu_solve(lu, B, check_finite=False)
        steps = []
    else:
        def mv(j):
            return j - m * green(j)

        linop = LinearOperator((len(m), len(m)), mv, dtype=complex)
        precond = LinearOperator(linop.shape, lambda j: j / (1 - m * diagonal), dtype=complex)
        J = np.empty_like(B)
        steps = []
        for column, b in enumerate(B.T):
            trace = []
            J[:, column], reason = gmres(linop, b, M=precond, rtol=1e-11, atol=0, restart=80,
                                         maxiter=30, callback=trace.append, callback_type="pr_norm")
            if reason:
                raise RuntimeError(f"GMRES source {column} failed with code {reason}")
            steps.append(len(trace))
    residual = np.linalg.norm(J - m[:, None] * (E + green(J)), axis=0)
    relative = residual / np.maximum(np.linalg.norm(B, axis=0), 1e-300)
    if np.max(relative) > 6e-10:
        raise RuntimeError(f"Full-wave state residual {max(relative):.3g}")
    return J, dict(seconds_wall=perf_counter()-start[0], seconds_CPU=process_time()-start[1],
                   iterations=steps, state_residual=relative.tolist(), physical_rhs=E.shape[1],
                   LU=int(dense is not None))


def antenna_forward(bare_response, network):
    """Same physics via the independent bare VIE + object-loaded network path."""
    g = network["gain"]
    delta_z = -1j * g[:, None] * bare_response * g[None, :]
    loaded_antenna = network["antenna"] + delta_z
    loaded_network = network["matrix"] + delta_z
    current = np.linalg.solve(loaded_network, network["drive"])
    empty = network["inv"] @ network["drive"]
    y = np.sqrt(R) * (empty - current)
    return y, current, loaded_antenna


def power_report(current, antenna, drive):
    columns = np.arange(drive.shape[1])
    all_port_voltage = antenna @ current
    driven_input = np.real(np.conj(current[columns, columns]) * all_port_voltage[columns, columns])
    antenna_power = np.real(np.sum(current.conj() * all_port_voltage, axis=0))
    terminations = R * (np.sum(abs(current)**2, axis=0) - abs(current[columns, columns])**2)
    available = np.sum(abs(drive)**2, axis=0) / (4 * R)
    reflected = abs(drive[columns, columns] / (2 * np.sqrt(R)) - np.sqrt(R) * current[columns, columns])**2
    return dict(available_W=available.tolist(), driven_accepted_W=driven_input.tolist(),
                antenna_W=antenna_power.tolist(), other_termination_W=terminations.tolist(),
                reflected_W=reflected.tolist(),
                energy_balance_max_W=float(np.max(abs(driven_input - antenna_power - terminations))),
                incident_balance_max_W=float(np.max(abs(available - reflected - driven_input))),
                reciprocity_relative=float(np.linalg.norm(antenna-antenna.T)/np.linalg.norm(antenna)),
                impedance_min_hermitian_eigen_ohms=float(np.linalg.eigvalsh((antenna+antenna.conj().T)/2)[0]))
