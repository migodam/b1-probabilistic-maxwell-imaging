"""Shared-weight current/material recurrence with fixed-operator feedback.

The feedback step holds its rowLS material and rho fixed. Autograd differentiates
through the unrolled map, including the upstream material; there is no implicit
solve and no dense LU. Supports complex64 and MPS float32 real/imag pairs.
"""
from __future__ import annotations
import torch
from torch import nn
try:
    from .networks import (NC, NS, GRID, MATERIAL_SCALE, ResidualUNet,
                           abs_square, architecture_cost)
except ImportError:
    from networks import (NC, NS, GRID, MATERIAL_SCALE, ResidualUNet,
                          abs_square, architecture_cost)

MODES = ("no_feedback", "refreshed_q", "feedback1", "feedback2", "feedback_hinge")
HINGE_CLEAN_DELTA = 1e-6
HINGE_NOISY_DELTA = .0011


def _norm2(value):
    return abs_square(value).reshape(len(value), -1).sum(dim=1)


def _multiply_material(material, value):
    return material[..., None]*value if value.is_complex() else material[..., None, None]*value


def _multiply_batch(scale, value):
    return scale.reshape(len(value), *([1]*(value.ndim-1)))*value


def row_ls_material(current, field):
    """The same shared-real rowLS and [.02,20] material box as the base arms."""
    if current.is_complex():
        numerator = (field.conj()*current).sum(dim=-1).real
        denominator = field.abs().square().sum(dim=-1)
    else:
        numerator = (field*current).sum(dim=(-1, -2))
        denominator = field.square().sum(dim=(-1, -2))
    return (numerator/denominator.clamp_min(1e-20)).clamp(.02, 20.)


def feedback_step(operators, current, material, public_y, anchor_current, rho,
                  *, field=None, receiver_residual=None, y_norm2=None, j0_norm2=None, data_tolerance=0.):
    """Quadratic exact-line step or noise-ball hinge majorizer step, per case.

    e=j-m(E+Gj), r=Sj-Y, Lv=v-mGv, L*v=v-G*(mv).
    Default data_tolerance=0 preserves the original quadratic calculation.
    Positive delta uses .5 max(||r||/||Y||-delta,0)^2 and theta-scaled data
    gradient, but retains the FULL unscaled data curvature in the denominator.
    The latter is a majorizer step, not an exact hinge line minimizer.
    Cached field/residual may be passed between steps while m/rho stay fixed.
    Zero gradient/curvature produces tau=0 without ever evaluating 0/0.
    """
    rho = float(rho)
    if rho < 0:
        raise ValueError("rho must be a fixed nonnegative scalar")
    cost = {"G_matmul_calls": 0, "G_adjoint_matmul_calls": 1,
            "S_matmul_calls": 1, "S_adjoint_matmul_calls": 1}
    if field is None:
        field = operators.total_field(current)
        cost["G_matmul_calls"] += 1
    if receiver_residual is None:
        receiver_residual = operators._matmul("S", current)-public_y
        cost["S_matmul_calls"] += 1
    raw_y_norm = _norm2(public_y) if y_norm2 is None else y_norm2
    raw_j_norm = _norm2(anchor_current) if j0_norm2 is None else j0_norm2
    yn = raw_y_norm.clamp_min(1e-20)
    jn = raw_j_norm.clamp_min(1e-20)
    e = current-_multiply_material(material, field)
    lstar_e = e-operators._matmul("G", _multiply_material(material, e), adjoint=True)
    relative_before = torch.linalg.vector_norm(receiver_residual.reshape(len(current), -1), dim=1)/yn.sqrt()
    hinge_data = torch.is_tensor(data_tolerance) or (data_tolerance is not None and float(data_tolerance) != 0.)
    if hinge_data:
        delta = torch.as_tensor(data_tolerance, dtype=yn.dtype, device=yn.device).expand_as(yn)
        safe_relative = torch.where(relative_before > 0, relative_before, torch.ones_like(relative_before))
        theta = torch.where(relative_before > 0,
                            torch.relu(1-delta/safe_relative), torch.zeros_like(relative_before))
        data_scale = theta/yn
    else:
        delta = torch.zeros_like(yn)
        theta = torch.ones_like(yn)
        data_scale = yn.reciprocal()
    gradient = (_multiply_batch(data_scale, operators._matmul("S", receiver_residual, adjoint=True))
                +rho*_multiply_batch(jn.reciprocal(), lstar_e))
    sg = operators._matmul("S", gradient)
    gg = operators._matmul("G", gradient)
    cost["G_matmul_calls"] += 1
    lg = gradient-_multiply_material(material, gg)
    gradient_norm2 = _norm2(gradient)
    # The data curvature is deliberately NOT multiplied by theta.
    curvature = _norm2(sg)/yn+rho*_norm2(lg)/jn
    active = (gradient_norm2 > 0) & (curvature > 0)
    safe_curvature = torch.where(curvature > 0, curvature, torch.ones_like(curvature))
    tau = torch.where(active, gradient_norm2/safe_curvature, torch.zeros_like(curvature))
    next_current = current-_multiply_batch(tau, gradient)
    next_field = field-_multiply_batch(tau, gg)
    next_r = receiver_residual-_multiply_batch(tau, sg)
    next_e = e-_multiply_batch(tau, lg)
    relative_after = torch.linalg.vector_norm(next_r.reshape(len(current), -1), dim=1)/yn.sqrt()
    if hinge_data:
        phi_before = .5*torch.relu(relative_before-delta).square()+.5*rho*_norm2(e)/jn
        phi_after = .5*torch.relu(relative_after-delta).square()+.5*rho*_norm2(next_e)/jn
    else:
        phi_before = .5*_norm2(receiver_residual)/yn+.5*rho*_norm2(e)/jn
        phi_after = .5*_norm2(next_r)/yn+.5*rho*_norm2(next_e)/jn
    return {"current": next_current, "field": next_field, "receiver_residual": next_r,
            "tau": tau, "gradient_norm2": gradient_norm2, "curvature": curvature,
            "phi_before": phi_before, "phi_after": phi_after,
            "data_tolerance": delta, "data_theta": theta,
            "relative_data_residual_before": relative_before, "relative_data_residual_after": relative_after,
            "zero_gradient": gradient_norm2 == 0,
            "zero_curvature": curvature == 0,
            "zero_y_normalization": raw_y_norm == 0,
            "zero_j0_normalization": raw_j_norm == 0,
            "fixed_operator_cost": cost}


def refreshed_features(anchor_features, material, field, stats, refresh_q):
    """Update m/Q using the original TRAIN statistics; gb/noise stay anchored."""
    features = anchor_features.clone()
    features[:, 0] = (material/MATERIAL_SCALE).reshape(-1, GRID, GRID)
    if refresh_q:
        qr, qi = (field.real, field.imag) if field.is_complex() else (field[..., 0], field[..., 1])
        raw_q = torch.cat((qr.transpose(1, 2), qi.transpose(1, 2)), dim=1).reshape(-1, 24, GRID, GRID)
        mean = torch.as_tensor(stats["mean"][2:26], dtype=features.dtype, device=features.device)[None, :, None, None]
        std = torch.as_tensor(stats["std"][2:26], dtype=features.dtype, device=features.device)[None, :, None, None]
        features[:, 2:26] = (raw_q-mean)/std
    return features


class RecurrentFeedback(nn.Module):
    """Two passes of one j_full U-Net; matched current/feedback controls."""
    def __init__(self, mode="feedback1", rho=1.):
        super().__init__()
        if mode not in MODES:
            raise ValueError(mode)
        if float(rho) < 0:
            raise ValueError("rho must be fixed and nonnegative")
        if mode == "feedback_hinge" and float(rho) != 1.:
            raise ValueError("The frozen feedback_hinge contract fixes rho=1")
        self.mode, self.rho, self.passes = mode, float(rho), 2
        self.feedback_steps = {"no_feedback": 0, "refreshed_q": 0, "feedback1": 1, "feedback2": 2, "feedback_hinge": 1}[mode]
        self.network = ResidualUNet("j_full")

    def forward(self, operators, features, m0, j0, public_y, stats, anchor_q=None, *, public_noise=None):
        j, m = j0, m0
        data_tolerance = None
        if self.mode == "feedback_hinge":
            if public_noise is None:
                raise ValueError("feedback_hinge requires explicit public_noise from the public anchor")
            noise = torch.as_tensor(public_noise, device=m0.device, dtype=m0.dtype).reshape(len(m0))
            data_tolerance = torch.where(noise > 0, HINGE_NOISY_DELTA, HINGE_CLEAN_DELTA)
        refresh_q = self.mode != "no_feedback"
        cost = {"G_matmul_calls": 0, "G_adjoint_matmul_calls": 0,
                "S_matmul_calls": 0, "S_adjoint_matmul_calls": 0}
        q = anchor_q
        if refresh_q and q is None:
            q = operators.total_field(j)
            cost["G_matmul_calls"] += 1
        yn, jn = _norm2(public_y), _norm2(j0)
        stage_materials, stage_currents, stage_fields, raw_outputs, diagnostics = [], [], [], [], []
        for layer in range(self.passes):
            x = refreshed_features(features, m, q, stats, refresh_q)
            raw = self.network(x)
            zr = raw[:, :NS].reshape(-1, NS, NC).transpose(1, 2)
            zi = raw[:, NS:].reshape(-1, NS, NC).transpose(1, 2)
            delta = torch.complex(zr, zi) if j.is_complex() else torch.stack((zr, zi), dim=-1)
            j = j+float(stats["current_scale"])*delta
            m, q = operators.material_from_current(j)
            cost["G_matmul_calls"] += 1
            # m remains fixed for both inner steps; refresh rowLS only afterward.
            fixed_m, r = m, None
            layer_diagnostics = []
            for _ in range(self.feedback_steps):
                step = feedback_step(operators, j, fixed_m, public_y, j0, self.rho,
                                     field=q, receiver_residual=r, y_norm2=yn, j0_norm2=jn,
                                     data_tolerance=data_tolerance)
                j, q, r = step["current"], step["field"], step["receiver_residual"]
                for key, value in step["fixed_operator_cost"].items():
                    cost[key] += value
                layer_diagnostics.append({key: step[key] for key in
                    ("tau", "gradient_norm2", "curvature", "phi_before", "phi_after", "zero_gradient",
                     "zero_curvature", "zero_y_normalization", "zero_j0_normalization",
                     "data_tolerance", "data_theta", "relative_data_residual_before", "relative_data_residual_after")})
            if self.feedback_steps:
                m = row_ls_material(j, q)
            raw_outputs.append(raw)
            stage_materials.append(m)
            stage_currents.append(j)
            stage_fields.append(q)
            diagnostics.append(layer_diagnostics)
        batch = len(j0)
        g_calls = cost["G_matmul_calls"]+cost["G_adjoint_matmul_calls"]
        s_calls = cost["S_matmul_calls"]+cost["S_adjoint_matmul_calls"]
        cost.update(network_passes=self.passes, shared_parameter_count=sum(p.numel() for p in self.network.parameters()),
                    feedback_steps_per_layer=self.feedback_steps, feedback_steps_total=self.passes*self.feedback_steps,
                    material_row_ls_calls=self.passes*(1+int(self.feedback_steps > 0)),
                    LU=0, actual_complex_RHS=0, dense_physics_inverse_calls=0,
                    feedback_step_rule="hinge_majorizer" if self.mode == "feedback_hinge" else "quadratic_exact_line",
                    fixed_matmul_complex_columns=batch*NS*(g_calls+s_calls),
                    fixed_operator_FLOPs_batch=8*batch*NS*(g_calls*NC*NC+s_calls*48*NC),
                    network_conv_FLOPs_batch=batch*architecture_cost("j_full", self.passes)["network_conv_FLOPs"],
                    flop_scope="complex MAC=8 real FLOPs; scalar/rowLS/normalization/BN/ReLU excluded",
                    setup_excluded=True, training_backward_cost="NOT_ESTIMATED_FORWARD_COUNTS_ONLY")
        return {"m_proposal": m, "j_aux": j, "q_aux": q, "delta_current": j-j0,
                "raw_delta": torch.stack(raw_outputs, dim=1),
                "stage_materials": torch.stack(stage_materials, dim=1),
                "stage_currents": torch.stack(stage_currents, dim=1),
                "stage_fields": torch.stack(stage_fields, dim=1),
                "feedback_diagnostics": diagnostics, "preparation_counts": cost}
