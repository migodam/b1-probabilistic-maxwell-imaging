"""R5 fixed-prior joint-energy iteration and explicitly labelled free control.

The physical state uses source columns and either complex tensors or MPS-safe
real/imag pairs. No factorization, projection, SVD or inverse occurs here.
Autograd follows the learned prior, row minimizer, metric, and physical paths;
only Armijo branch decisions are discrete. A qualified step is an energy
descent step, not a claim of global convergence or correct material recovery.
"""
from __future__ import annotations

from pathlib import Path
import math

import numpy as np
import torch
from torch import Tensor, nn

try:
    from .architectures import PriorCNN, MetricCNN, parameter_counts
except ImportError:
    from architectures import PriorCNN, MetricCNN, parameter_counts

GRID, NC, NS, MATERIAL_SCALE = 24, 576, 12, .4


class FixedOperators(nn.Module):
    """Read-only original G/E/S buffers, with real/imag matmul and adjoints.

    Construction only loads operators.npz; it never writes an R4 cache and
    never constructs a singular basis. Move this module once with .to(device).
    Counts represent forward fixed complex operator applications, not LU RHS.
    """
    def __init__(self, path: str | Path):
        super().__init__()
        self.source_path = str(Path(path).resolve())
        with np.load(path, allow_pickle=False) as data:
            arrays = {"G": np.asarray(data["G"]),
                      "E": np.asarray(data["E0"] if "E0" in data else data["E"] if "E" in data else data["U0"]),
                      "S": np.asarray(data["S"] if "S" in data else data["H"])}
        if (arrays["G"].shape != (NC, NC) or arrays["E"].shape != (NC, NS)
                or arrays["S"].shape != (48, NC)):
            raise ValueError("expected G576x576, E576x12 and S48x576")
        for name, value in arrays.items():
            if not np.isfinite(value).all():
                raise ValueError(f"nonfinite fixed operator {name}")
            self.register_buffer(name + "_re", torch.from_numpy(np.array(value.real, dtype=np.float32)))
            self.register_buffer(name + "_im", torch.from_numpy(np.array(value.imag, dtype=np.float32)))
        self.call_counts = {"G_matmul_calls": 0, "G_adjoint_matmul_calls": 0,
                            "S_matmul_calls": 0, "S_adjoint_matmul_calls": 0}

    def matmul(self, name: str, z: Tensor, adjoint: bool = False) -> Tensor:
        if name not in ("G", "S"):
            raise ValueError("matmul name must be G or S")
        key = f"{name}_{'adjoint_' if adjoint else ''}matmul_calls"
        self.call_counts[key] += 1
        real, imag = getattr(self, name + "_re"), getattr(self, name + "_im")
        if z.is_complex():
            matrix = torch.complex(real, imag)
            return (matrix.mH if adjoint else matrix) @ z
        if z.shape[-1] != 2:
            raise ValueError("real operator inputs require a trailing real/imag axis of length 2")
        if adjoint:
            real, imag = real.T, -imag.T
        zr, zi = z[..., 0], z[..., 1]
        return torch.stack((real @ zr - imag @ zi, real @ zi + imag @ zr), dim=-1)

    def field(self, current: Tensor) -> Tensor:
        incident = (torch.complex(self.E_re, self.E_im) if current.is_complex()
                    else torch.stack((self.E_re, self.E_im), dim=-1))
        return incident + self.matmul("G", current)

    # Compatibility names do not instantiate or call the R4 cache class.
    _matmul = matmul
    total_field = field


def abs_square(value: Tensor) -> Tensor:
    return value.real.square() + value.imag.square() if value.is_complex() else value.square().sum(dim=-1)


def norm2(value: Tensor) -> Tensor:
    return abs_square(value).reshape(len(value), -1).sum(dim=1)


def row_norm2(value: Tensor) -> Tensor:
    return abs_square(value).sum(dim=-1)


def material_product(material: Tensor, value: Tensor) -> Tensor:
    return material[..., None] * value if value.is_complex() else material[..., None, None] * value


def batch_product(scale: Tensor, value: Tensor) -> Tensor:
    return scale.reshape(len(value), *([1] * (value.ndim - 1))) * value


def real_inner(left: Tensor, right: Tensor) -> Tensor:
    product = (left.conj() * right).real if left.is_complex() else left * right
    return product.reshape(len(left), -1).sum(dim=1)


def row_material(current: Tensor, field: Tensor, prior: Tensor,
                 lam: Tensor, b: Tensor, rho: float = 1.) -> Tensor:
    """Exact box-constrained real row minimizer of the specified joint energy."""
    numerator = ((field.conj() * current).real.sum(dim=-1) if current.is_complex()
                 else (field * current).sum(dim=(-1, -2)))
    numerator = rho * numerator / b[:, None] + lam * prior
    denominator = rho * row_norm2(field) / b[:, None] + lam
    # If q and Q0 both vanish, the row energy is constant in m: c is a minimizer.
    safe = torch.where(denominator > 0, denominator, torch.ones_like(denominator))
    return torch.where(denominator > 0, numerator / safe, prior).clamp(.02, 20.)


def refreshed_features(features: Tensor, material: Tensor, field: Tensor, stats: dict) -> Tensor:
    """Refresh only raw-m/.4 and TRAIN-normalized Q; gb/noise remain anchored."""
    qr, qi = (field.real, field.imag) if field.is_complex() else (field[..., 0], field[..., 1])
    raw_q = torch.cat((qr.transpose(1, 2), qi.transpose(1, 2)), dim=1).reshape(-1, 24, GRID, GRID)
    mean = torch.as_tensor(stats["mean"][2:26], dtype=features.dtype, device=features.device)[None, :, None, None]
    std = torch.as_tensor(stats["std"][2:26], dtype=features.dtype, device=features.device)[None, :, None, None]
    q_features = (raw_q - mean) / std.clamp_min(1e-20)
    return torch.cat(((material / MATERIAL_SCALE).reshape(-1, 1, GRID, GRID),
                      features[:, 1:2], q_features, features[:, 26:27]), dim=1)


def joint_energy(receiver_residual: Tensor, defect: Tensor, material: Tensor,
                 prior: Tensor, lam: Tensor, a: Tensor, b: Tensor, rho: float = 1., tolerance=None) -> Tensor:
    data = .5*norm2(receiver_residual)/a if tolerance is None else .5*torch.relu(
        (norm2(receiver_residual)/a).clamp_min(1e-30).sqrt()-tolerance).square()
    return (data + .5 * rho * norm2(defect) / b
            + .5 * (lam * (material - prior).square()).sum(dim=1))


class FeedbackCNN(nn.Module):
    def __init__(self, representation: str = "pixel", closed: bool = True, rho: float = 1., noise_ball: bool = False):
        super().__init__()
        if not math.isfinite(float(rho)) or float(rho) <= 0:
            raise ValueError("rho must be a positive fixed public scalar")
        self.representation, self.closed, self.rho = representation, bool(closed), float(rho)
        self.noise_ball=bool(noise_ball)
        self.prior_network = PriorCNN(representation)
        self.metric_network = MetricCNN(representation)

    def parameter_counts(self) -> dict[str, int]:
        return parameter_counts(self)

    def forward(self, ops: FixedOperators, features: Tensor, m0: Tensor, j0: Tensor,
                q0: Tensor, Y: Tensor, stats: dict, step_size: Tensor | float,
                steps: int = 4) -> dict:
        if not isinstance(steps, int) or isinstance(steps, bool) or steps < 1:
            raise ValueError("steps must be a positive integer")
        if m0.ndim != 2 or m0.shape[1] != NC:
            raise ValueError("m0 must have shape (B,576)")
        if features.shape != (len(m0), 27, GRID, GRID):
            raise ValueError("features must have shape (B,27,24,24)")
        if j0.shape != q0.shape or j0.shape[:3] != (len(m0), NC, NS):
            raise ValueError("j0 and q0 must have matching source-column shapes")
        a_raw, b_raw = norm2(Y), norm2(j0)
        a, b = a_raw.clamp_min(1e-20), b_raw.clamp_min(1e-20)
        lam = self.rho * row_norm2(q0) / b[:, None]
        initial_tau = torch.as_tensor(step_size, dtype=m0.dtype, device=m0.device).reshape(-1)
        initial_tau = initial_tau.expand(len(m0)) if initial_tau.numel() == 1 else initial_tau
        if initial_tau.shape != (len(m0),):
            raise ValueError("step_size must be a scalar or one public value per case")
        noise = (features[:, 26:27] * float(stats["std"][26]) + float(stats["mean"][26]))
        tolerance=torch.where(noise[:,0,0,0]>0,.003,.001) if self.noise_ball else None
        def data_energy(r):
            return .5*norm2(r)/a if tolerance is None else .5*torch.relu(
                (norm2(r)/a).clamp_min(1e-30).sqrt()-tolerance).square()
        q0_row = row_norm2(q0).clamp_min(1e-20).sqrt()
        j0_row = row_norm2(j0).clamp_min(1e-20).sqrt()
        cost = {"G_matmul_calls": 0, "G_adjoint_matmul_calls": 0,
                "S_matmul_calls": 0, "S_adjoint_matmul_calls": 0}

        def apply(name: str, value: Tensor, adjoint: bool = False) -> Tensor:
            cost[f"{name}_{'adjoint_' if adjoint else ''}matmul_calls"] += 1
            return ops.matmul(name, value, adjoint=adjoint)

        # Evaluate q from the actual fixed operator; Q0 remains the public metric
        # and lambda reference, rather than silently replacing a physical field.
        j, material = j0, m0
        q = ops.field(j)
        cost["G_matmul_calls"] += 1
        receiver_residual = apply("S", j) - Y
        fixed_prior = ((m0 + MATERIAL_SCALE * self.prior_network(features).flatten(1)).clamp(.02, 20.)
                       if self.closed else None)
        materials, currents, fields, defects, priors, row_materials, gain_history = [], [], [], [], [], [], []
        diagnostics = []
        for iteration in range(steps):
            prior = (fixed_prior if self.closed else
                     (material + MATERIAL_SCALE * self.prior_network(
                         refreshed_features(features, material, q, stats)).flatten(1)).clamp(.02, 20.))
            energy_before_row = joint_energy(
                receiver_residual, j - material_product(material, q), material, prior, lam, a, b, self.rho)
            if iteration == 0:
                initial_energy = energy_before_row
            material = row_material(j, q, prior, lam, b, self.rho)
            defect = j - material_product(material, q)
            if tolerance is None:
                data_scale=1/a
            else:
                relative=(norm2(receiver_residual)/a).clamp_min(1e-30).sqrt()
                data_scale=torch.relu(1-tolerance/relative)/a
            gradient = (batch_product(data_scale, apply("S", receiver_residual, True))
                        + self.rho * batch_product(1 / b, defect - apply(
                            "G", material_product(material, defect), True)))
            gradient_rms = (norm2(gradient) / (NC * NS)).clamp_min(1e-20).sqrt()
            metric_maps = torch.stack((material / MATERIAL_SCALE, prior / MATERIAL_SCALE,
                                       row_norm2(defect).clamp_min(1e-20).sqrt() / q0_row,
                                       row_norm2(j).clamp_min(1e-20).sqrt() / j0_row,
                                       row_norm2(q).clamp_min(1e-20).sqrt() / q0_row,
                                       row_norm2(gradient).clamp_min(1e-20).sqrt() / gradient_rms[:, None]), dim=1)
            metric_features = torch.cat((metric_maps.reshape(-1, 6, GRID, GRID), noise), dim=1)
            direction, gains = self.metric_network(metric_features, gradient)
            sd, gd = apply("S", direction), apply("G", direction)
            ld = direction - material_product(material, gd)
            directional_derivative = real_inner(gradient, direction)
            phi_before = data_energy(receiver_residual) + .5 * self.rho * norm2(defect) / b
            prior_energy_before = .5 * (lam * (material - prior).square()).sum(dim=1)
            tau = torch.zeros_like(a)
            accepted = torch.zeros_like(a, dtype=torch.bool)
            backtracks = torch.full_like(a, 6, dtype=torch.long)
            attempts = torch.zeros_like(a, dtype=torch.long)
            # Initial proposal plus at most six halvings. Candidate objectives
            # reuse Sd and Ld, so line search performs no extra operator call.
            for halvings in range(7):
                candidate_tau = initial_tau * (0.5 ** halvings)
                candidate_r = receiver_residual - batch_product(candidate_tau, sd)
                candidate_e = defect - batch_product(candidate_tau, ld)
                candidate_phi = data_energy(candidate_r) + .5 * self.rho * norm2(candidate_e) / b
                pending = ~accepted
                attempts = attempts + pending.to(attempts.dtype)
                qualifies = (pending & (candidate_tau > 0) & torch.isfinite(candidate_tau)
                             & (directional_derivative > 0) & torch.isfinite(candidate_phi)
                             & (candidate_phi <= phi_before - 1e-4 * candidate_tau * directional_derivative))
                tau = torch.where(qualifies, candidate_tau, tau)
                backtracks = torch.where(qualifies, torch.full_like(backtracks, halvings), backtracks)
                accepted = accepted | qualifies
            def accepted_or_previous(candidate: Tensor, previous: Tensor) -> Tensor:
                mask = accepted.reshape(len(accepted), *([1] * (candidate.ndim - 1)))
                return torch.where(mask, candidate, previous)

            # Explicit selection preserves the previous state on rejection,
            # including a rejected nonfinite direction (0*NaN is still NaN).
            next_j = accepted_or_previous(j - batch_product(tau, direction), j)
            next_q = accepted_or_previous(q - batch_product(tau, gd), q)
            next_r = accepted_or_previous(receiver_residual - batch_product(tau, sd), receiver_residual)
            fixed_m_defect = accepted_or_previous(defect - batch_product(tau, ld), defect)
            next_material = row_material(next_j, next_q, prior, lam, b, self.rho)
            next_defect = next_j - material_product(next_material, next_q)
            energy_before = phi_before + prior_energy_before
            fixed_m_after = joint_energy(next_r, fixed_m_defect, material, prior, lam, a, b, self.rho,tolerance)
            energy_after = joint_energy(next_r, next_defect, next_material, prior, lam, a, b, self.rho,tolerance)
            diagnostics.append({"accepted": accepted, "rejected": ~accepted,
                                "status_code": accepted.to(torch.long), "tau": tau,
                                "backtracks": backtracks, "line_search_attempts": attempts,
                                "energy_before_row": energy_before_row,
                                "energy_before": energy_before, "energy_after_fixed_m": fixed_m_after,
                                "energy_after": energy_after, "energy_change": energy_after - energy_before,
                                "gradient_norm2": norm2(gradient), "directional_derivative": directional_derivative,
                                "current_update_norm": norm2(next_j - j).sqrt(),
                                "relative_current_update": (norm2(next_j - j) / b).sqrt(),
                                "material_update_norm": (next_material - material).square().sum(dim=1).sqrt(),
                                "relative_data_residual": (norm2(next_r) / a).sqrt(),
                                "relative_constitutive_defect": (norm2(next_defect) / b).sqrt()})
            row_materials.append(material)
            material, j, q, receiver_residual = next_material, next_j, next_q, next_r
            materials.append(material); currents.append(j); fields.append(q); defects.append(next_defect)
            priors.append(prior); gain_history.append(gains)
        g_calls = cost["G_matmul_calls"] + cost["G_adjoint_matmul_calls"]
        s_calls = cost["S_matmul_calls"] + cost["S_adjoint_matmul_calls"]
        cost.update(LU=0, actual_complex_RHS=0, dense_physics_inverse_calls=0,
                    prior_network_calls=1 if self.closed else steps, metric_network_calls=steps,
                    material_row_solve_calls=2 * steps, parameter_counts=self.parameter_counts(),
                    fixed_matmul_complex_columns=len(m0) * NS * (g_calls + s_calls),
                    fixed_operator_FLOPs_batch=8 * len(m0) * NS * (g_calls * NC * NC + s_calls * 48 * NC),
                    setup_excluded=True, training_backward_cost="NOT_ESTIMATED_FORWARD_COUNTS_ONLY",
                    line_search_operator_calls=0, maximum_halvings=6)
        return {"prior": prior, "lambda": lam, "m_proposal": material, "j_aux": j,
                "initial_energy": initial_energy, "energy": energy_after,
                "q_aux": q, "e_aux": next_defect, "delta_current": j - j0,
                "stage_materials": torch.stack(materials, dim=1),
                "stage_currents": torch.stack(currents, dim=1),
                "stage_fields": torch.stack(fields, dim=1), "stage_defects": torch.stack(defects, dim=1),
                "stage_priors": torch.stack(priors, dim=1),
                "row_materials": torch.stack(row_materials, dim=1),
                "material_trajectory": torch.stack([m0, *materials], dim=1),
                "current_trajectory": torch.stack([j0, *currents], dim=1),
                "metric_gains": gain_history,
                "trace": {key: torch.stack([item[key] for item in diagnostics], dim=1)
                          for key in diagnostics[0]},
                "trace_status_labels": {0: "REJECTED", 1: "ACCEPTED"},
                "cost": cost, "preparation_counts": cost,
                "closed": self.closed,
                "noise_ball":self.noise_ball,
                "data_objective":"noise_tolerant_hinge" if self.noise_ball else "quadratic",
                "energy_contract": "FIXED_PUBLIC_ANCHOR_FIXED_PRIOR" if self.closed
                    else "FREE_CASCADE_NO_SHARED_FIXED_ENERGY_GUARANTEE",
                "zero_y_normalization": a_raw == 0, "zero_j0_normalization": b_raw == 0}
