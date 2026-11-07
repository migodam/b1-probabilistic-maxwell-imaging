"""TRAIN/VAL coupled inverse-process loss; exact R4 Maxwell VJP is reused.

No current truth is accepted. Tail differences are dynamic residual diagnostics,
not a fixed-point or convergence certificate.
"""
from __future__ import annotations

import torch
from torch import nn

from ...research_gn_nn_preservation_r4.src.train import FullVIEBridge, full_maxwell_hinge

MATERIAL_SCALE = .4
CLEAN_TOLERANCE = 1e-3
NOISY_TOLERANCE = 3e-3


def selected_stage_indices(steps: int) -> tuple[int, ...]:
    """Zero-based stages; all through four, then declared stages 2,4,last."""
    if steps < 1:
        raise ValueError("At least one inverse-process stage is required")
    return tuple(range(steps)) if steps <= 4 else tuple(sorted({1, 3, steps-1}))


def _pair(value: torch.Tensor) -> torch.Tensor:
    if value.is_complex():
        return torch.view_as_real(value.resolve_conj())
    if value.shape[-1] != 2:
        raise ValueError("Current/field must be complex or have trailing real/imag pair")
    return value


def _norm2(value: torch.Tensor) -> torch.Tensor:
    return _pair(value).square().flatten(1).sum(dim=1)


class CoupledTrainingLoss(nn.Module):
    """Loss of supervised material and measured full-forward consistency.

    ``ops`` provides field(j), or R4's total_field(j). Alternatively output must
    provide stage_fields, final_field, field, or q. Forward returns a scalar loss
    and detached tensor metrics. Use receipt() AFTER backward for adjoint costs.
    No-field controls set lambda_field=0; evaluate_field=True explicitly retains
    the validation full-forward residual while avoiding its gradient graph.
    """

    def __init__(self, operators_path, ops=None, lambda_field: float = 100.):
        super().__init__()
        if float(lambda_field) not in (0., 100.):
            raise ValueError("R5 freezes lambda_field to 100 (coupled) or 0 (no-field)")
        self.lambda_field = float(lambda_field)
        self.bridge = FullVIEBridge(operators_path)
        self.ops = ops

    def receipt(self):
        """Live actual R4 LU/source/adjoint accounting, including completed VJPs."""
        return self.bridge.receipt()

    def _field(self, output, current):
        for name in ("final_field", "field", "q"):
            if name in output:
                return output[name]
        if "stage_fields" in output:
            return output["stage_fields"][:, -1]
        if self.ops is None:
            raise ValueError("Final constitutive defect needs output field or supplied ops")
        method = getattr(self.ops, "field", None)
        if method is None:
            method = getattr(self.ops, "total_field", None)
        if method is None:
            raise TypeError("ops must provide field(j) or total_field(j)")
        return method(current)

    def forward(self, output, m_gt, public_y_pair, public_noise, j0,
                *, evaluate_field: bool = False):
        materials = output["stage_materials"]
        currents = output["stage_currents"]
        prior = output["prior"]
        if materials.ndim != 3 or materials.shape[-1] != 576:
            raise ValueError("stage_materials must have shape (B,K,576)")
        batch, steps, _ = materials.shape
        expected = (batch, steps, 576, 12)
        if (currents.shape if currents.is_complex() else currents.shape[:-1]) != expected:
            raise ValueError("stage_currents must be (B,K,576,12) complex or trailing pair")
        if prior.shape != (batch, 576) or m_gt.shape != prior.shape:
            raise ValueError("prior and m_gt must both have shape (B,576)")
        selected = selected_stage_indices(steps)
        truth = m_gt.to(device=materials.device, dtype=materials.dtype)
        selected_materials = materials[:, selected, :]
        stage_mse = (selected_materials-truth[:, None, :]).square().mean()/MATERIAL_SCALE**2
        prior_mse = (prior-truth).square().mean()/MATERIAL_SCALE**2
        material_loss = stage_mse+.5*prior_mse

        noise = torch.as_tensor(public_noise, device=materials.device).reshape(-1)
        if noise.numel() == 1:
            noise = noise.expand(batch)
        if noise.numel() != batch or not bool(torch.isfinite(noise).all()) or bool((noise < 0).any()):
            raise ValueError("public_noise must be a finite nonnegative scalar or per-case vector")
        tolerance = torch.where(noise > 0, NOISY_TOLERANCE, CLEAN_TOLERANCE)
        field_hinge = materials.new_zeros(())
        residuals = []
        if self.lambda_field:
            # The returned hinge keeps the exact full Maxwell VJP to each NN
            # stage. No detach, surrogate response, or GT-current target.
            hinges = []
            for index in selected:
                hinge, info = full_maxwell_hinge(materials[:, index], public_y_pair,
                                                  self.bridge, tolerance)
                hinges.append(hinge)
                residuals.append(info["relative_public_residual"].detach())
            field_hinge = torch.stack(hinges).mean()
        elif evaluate_field:
            # Explicit validation evaluation for the same no-field algorithm.
            with torch.no_grad():
                hinges = []
                for index in selected:
                    hinge, info = full_maxwell_hinge(materials[:, index], public_y_pair,
                                                      self.bridge, tolerance)
                    hinges.append(hinge)
                    residuals.append(info["relative_public_residual"].detach())
                field_hinge = torch.stack(hinges).mean()

        final_material, final_current = materials[:, -1], currents[:, -1]
        final_field = self._field(output, final_current)
        current_pair = _pair(final_current)
        field_pair = _pair(final_field).to(device=current_pair.device, dtype=current_pair.dtype)
        if field_pair.shape != current_pair.shape:
            raise ValueError("Final field and current shapes must match")
        anchor_norm2 = _norm2(j0).to(current_pair.device).clamp_min(1e-20)
        defect = current_pair-final_material[..., None, None]*field_pair
        constitutive = (defect.square().flatten(1).sum(1)/anchor_norm2).mean()
        if steps >= 2:
            current_delta = current_pair-_pair(currents[:, -2])
            tail_current = (current_delta.square().flatten(1).sum(1)/anchor_norm2).mean()
            tail_material = (final_material-materials[:, -2]).square().mean()/MATERIAL_SCALE**2
        else:
            tail_current = materials.new_zeros(())
            tail_material = materials.new_zeros(())
        tail = tail_current+tail_material
        loss = material_loss+self.lambda_field*field_hinge+.01*constitutive+.01*tail
        metrics = {
            "loss": loss, "material_loss": material_loss,
            "stage_material_mse_normalized": stage_mse,
            "prior_material_mse_normalized": prior_mse,
            "full_data_hinge": field_hinge,
            "constitutive_normalized": constitutive,
            "fixed_point_tail": tail, "tail_current_normalized": tail_current,
            "tail_material_mse_normalized": tail_material,
            "tail_available": materials.new_tensor(steps >= 2, dtype=torch.bool),
            "field_evaluated": materials.new_tensor(bool(residuals), dtype=torch.bool),
            "selected_stages": torch.tensor([i+1 for i in selected], device=materials.device),
            "qualification_tolerance": tolerance,
        }
        if residuals:
            relative = torch.stack(residuals, dim=1)
            metrics.update(relative_public_residual=relative.mean(),
                           relative_public_residual_per_stage=relative,
                           final_relative_public_residual=relative[:, -1],
                           final_field_qualified=relative[:, -1] <= tolerance)
        return loss, {key: value.detach() for key, value in metrics.items()}
