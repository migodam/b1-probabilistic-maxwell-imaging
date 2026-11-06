"""Residual U-Net and differentiable scalar-TM proposal heads.

The encoder/decoder structure is adapted from FAEMIS/faemis/models/unet.py
(DoubleConv, two down/up blocks), copied locally without an invasive import.
Channels are [m, gb, Re(Q[:,0:12]), Im(Q[:,0:12]), noise]. Spatial flattening
is row-major; currents have shape (..., cell, source), never source-major flat.
When flattened for a public receiver vector, use Y.T.reshape(-1).
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import torch
from torch import nn

ARMS = ("v_image", "v_physics", "j_full", "j_null", "j_null_current")
MATERIAL_SCALE = .4
GRID = 24
NS = 12
NC = GRID * GRID


def choose_device(requested="auto"):
    if requested == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        return torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    device = torch.device(requested)
    if device.type == "mps" and not torch.backends.mps.is_available():
        raise RuntimeError("MPS requested but unavailable")
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    return device


class DoubleConv(nn.Module):
    def __init__(self, cin, cout):
        super().__init__()
        self.layers = nn.Sequential(
            nn.Conv2d(cin, cout, 3, padding=1), nn.BatchNorm2d(cout), nn.ReLU(inplace=True),
            nn.Conv2d(cout, cout, 3, padding=1), nn.BatchNorm2d(cout), nn.ReLU(inplace=True))

    def forward(self, x):
        return self.layers(x)


class ResidualUNet(nn.Module):
    """27 -> (16,32,64), grid 24 -> 12 -> 6 -> 12 -> 24."""
    def __init__(self, arm):
        super().__init__()
        if arm not in ARMS:
            raise ValueError(f"Unknown arm {arm}")
        self.arm = arm
        self.stem = DoubleConv(27, 16)
        self.down1 = nn.Sequential(nn.MaxPool2d(2), DoubleConv(16, 32))
        self.down2 = nn.Sequential(nn.MaxPool2d(2), DoubleConv(32, 64))
        self.up1 = nn.ConvTranspose2d(64, 32, 2, stride=2)
        self.decode1 = DoubleConv(64, 32)
        self.up2 = nn.ConvTranspose2d(32, 16, 2, stride=2)
        self.decode2 = DoubleConv(32, 16)
        self.head = nn.Conv2d(16, 1 if arm.startswith("v_") else 24, 1)
        nn.init.zeros_(self.head.weight)
        nn.init.zeros_(self.head.bias)

    def forward(self, x):
        if x.shape[1:] != (27, 24, 24):
            raise ValueError(f"Expected (batch,27,24,24), got {tuple(x.shape)}")
        a = self.stem(x)
        b = self.down1(a)
        c = self.down2(b)
        b = self.decode1(torch.cat((b, self.up1(c)), dim=1))
        a = self.decode2(torch.cat((a, self.up2(b)), dim=1))
        return self.head(a)


def architecture_cost(arm, passes=1):
    """Analytic forward arithmetic estimate; 1 multiply+add = 2 FLOPs.

    Convolution/transpose-convolution only for the network; BN/ReLU, transfers,
    normalization and rowLS scalar work are excluded and timed in receipts.
    Complex matrix products count 8 real FLOPs per complex multiply-add.
    """
    out = 1 if arm.startswith("v_") else 24
    macs = (24**2*(27*16*9+16*16*9) + 12**2*(16*32*9+32*32*9)
            + 6**2*(32*64*9+64*64*9) + 6**2*64*32*4
            + 12**2*(64*32*9+32*32*9) + 12**2*32*16*4
            + 24**2*(32*16*9+16*16*9) + 24**2*16*out)
    is_current = not arm.startswith("v_")
    null = arm in ("j_null", "j_null_current")
    fixed_flops = (8*NC*NC*NS if is_current else 0) + (16*NC*48*NS if null else 0)
    return {"network_conv_MACs_per_pass": macs, "network_conv_FLOPs": passes*2*macs,
            "fixed_complex_matmul_FLOPs": passes*fixed_flops,
            "raw_output_real_DOF": NC*(24 if is_current else 1),
            "effective_output_real_DOF": 2*(NC-48)*NS if null else NC*(24 if is_current else 1),
            "flop_convention": "2 per real MAC; 8 per complex MAC; excludes BN/ReLU/rowLS/scalar/IO"}


def raw_features(anchor):
    m = np.asarray(anchor["m0"], dtype=np.float32).reshape(NC)
    gb = np.asarray(anchor["gb"], dtype=np.float32).reshape(NC)
    q = np.asarray(anchor["Q0"], dtype=np.complex64).reshape(NC, NS)
    noise = float(np.asarray(anchor["noise"]).reshape(-1)[0])
    return np.concatenate((m[None], gb[None], q.real.T, q.imag.T,
                           np.full((1, NC), noise, np.float32)), axis=0).reshape(27, GRID, GRID)


def fit_feature_stats(anchors):
    """Streaming, public TRAIN anchors only. Material uses the frozen .4 scale."""
    total = np.zeros(27, np.float64)
    square = np.zeros(27, np.float64)
    count = 0
    q_square = 0.
    q_count = 0
    for anchor in anchors:
        x = raw_features(anchor).reshape(27, -1).astype(np.float64)
        total += x.sum(axis=1)
        square += (x*x).sum(axis=1)
        count += x.shape[1]
        q = np.asarray(anchor["Q0"])
        q_square += float(np.sum(np.abs(q)**2))
        q_count += q.size
    if count == 0:
        raise ValueError("No training anchors for feature statistics")
    mean = total / count
    std = np.sqrt(np.maximum(square/count-mean**2, 0.))
    std[std < 1e-8] = 1.
    mean[0], std[0] = 0., MATERIAL_SCALE
    return {"mean": mean.astype(np.float32).tolist(), "std": std.astype(np.float32).tolist(),
            "current_scale": max(MATERIAL_SCALE*np.sqrt(q_square/max(q_count, 1)), 1e-8),
            "fitted_on": "train_public_anchors_only", "material_scale": MATERIAL_SCALE}


def normalized_features(anchor, stats, arm):
    x = raw_features(anchor)
    mean = np.asarray(stats["mean"], np.float32)[:, None, None]
    std = np.asarray(stats["std"], np.float32)[:, None, None]
    x = (x-mean)/std
    if arm == "v_image":
        x[1:26] = 0.
    return torch.from_numpy(x.astype(np.float32))


class OperatorCache(nn.Module):
    """Immutable G/E/S and all 48 right singular vectors of S, cached once.

    P0 Z = Z - V48 (V48^H Z). There is no weak-mode truncation. SVD is
    performed in complex128 on CPU before converting to complex64, only if the
    supplied operator NPZ does not already include V48.
    """
    def __init__(self, path: str | Path):
        super().__init__()
        with np.load(path, allow_pickle=False) as p:
            g = np.asarray(p["G"])
            e = np.asarray(p["E0"] if "E0" in p else p["E"] if "E" in p else p["U0"])
            s = np.asarray(p["S"] if "S" in p else p["H"])
            if "V48" in p:
                v = np.asarray(p["V48"])
                self.v48_source = "operators_npz"
            else:
                _, singular, vh = np.linalg.svd(s.astype(np.complex128), full_matrices=False)
                if len(singular) != 48 or singular[-1] <= np.finfo(float).eps*max(s.shape)*singular[0]:
                    raise ValueError("Expected full numerical receiver rank 48")
                v = vh.conj().T
                self.v48_source = "full_rank48_cpu_svd"
        if g.shape != (NC, NC) or e.shape != (NC, NS) or s.shape != (48, NC) or v.shape != (NC, 48):
            raise ValueError("R3 requires G576x576 E576x12 S48x576 V48=576x48")
        self._complex_cache = {}
        for name, array in (("G", g), ("E", e), ("S", s), ("V48", v)):
            if not np.isfinite(array).all():
                raise ValueError(f"Nonfinite {name}")
            self.register_buffer(name+"_re", torch.from_numpy(np.array(array.real, dtype=np.float32)))
            self.register_buffer(name+"_im", torch.from_numpy(np.array(array.imag, dtype=np.float32)))

    def _complex(self, name):
        real, imag = getattr(self, name+"_re"), getattr(self, name+"_im")
        key = (name, str(real.device))
        if key not in self._complex_cache:
            self._complex_cache[key] = torch.complex(real, imag)
        return self._complex_cache[key]

    def _matmul(self, name, z, adjoint=False):
        if z.is_complex():
            a = self._complex(name)
            return (a.mH if adjoint else a) @ z
        ar, ai = getattr(self, name+"_re"), getattr(self, name+"_im")
        if adjoint:
            ar, ai = ar.T, -ai.T
        zr, zi = z[..., 0], z[..., 1]
        return torch.stack((ar@zr-ai@zi, ar@zi+ai@zr), dim=-1)

    def null_project(self, z):
        return z - self._matmul("V48", self._matmul("V48", z, adjoint=True))

    def material_from_current(self, j):
        if j.is_complex():
            q = self._complex("E") + self._matmul("G", j)
            numerator = (q.conj()*j).sum(dim=-1).real
            denominator = q.abs().square().sum(dim=-1).clamp_min(1e-20)
        else:
            # Native float32 MPS: (...,cell,source,Re/Im), identical real algebra.
            q = torch.stack((self.E_re, self.E_im), dim=-1) + self._matmul("G", j)
            numerator = (q*j).sum(dim=(-1, -2))
            denominator = q.square().sum(dim=(-1, -2)).clamp_min(1e-20)
        return (numerator/denominator).clamp(.02, 20.), q


def proposal(model, operators, features, m0, j0, stats, passes=1):
    """No full Maxwell solve here; downstream verification is mandatory.

    Shared-weight cascade uses updated m only in channel 0, keeping gb/Q/noise
    anchored. Current stages add projected increments to accumulated j_aux.
    """
    if passes not in (1, 2):
        raise ValueError("Only one or two shared-weight passes are supported")
    if passes == 2 and model.arm not in ("v_physics", "j_null"):
        raise ValueError("Two-pass cascade enabled only for v_physics and j_null")
    m = m0
    j = j0
    q = None
    raw_outputs = []
    for _ in range(passes):
        x = features.clone()
        x[:, 0] = (m / MATERIAL_SCALE).reshape(-1, GRID, GRID)
        raw = model(x)
        raw_outputs.append(raw)
        if model.arm.startswith("v_"):
            m = (m + MATERIAL_SCALE*raw[:, 0].reshape(-1, NC)).clamp(.02, 20.)
        else:
            zr = raw[:, :NS].reshape(-1, NS, NC).transpose(1, 2)
            zi = raw[:, NS:].reshape(-1, NS, NC).transpose(1, 2)
            z = (torch.complex(zr, zi) if j.is_complex() else torch.stack((zr, zi), dim=-1))*float(stats["current_scale"])
            if model.arm in ("j_null", "j_null_current"):
                z = operators.null_project(z)
            j = j + z
            m, q = operators.material_from_current(j)
    out = {"m_proposal": m, "raw_delta": torch.stack(raw_outputs, dim=1)}
    if not model.arm.startswith("v_"):
        out.update(j_aux=j, q_aux=q, delta_current=j-j0)
    return out


def material_loss(prediction, target):
    p = prediction.reshape(-1, GRID, GRID)/MATERIAL_SCALE
    t = target.reshape(-1, GRID, GRID)/MATERIAL_SCALE
    l1 = (p-t).abs().mean()
    edge_x = ((p[:, :, 1:]-p[:, :, :-1])-(t[:, :, 1:]-t[:, :, :-1])).abs()
    edge_y = ((p[:, 1:, :]-p[:, :-1, :])-(t[:, 1:, :]-t[:, :-1, :])).abs()
    edge = (edge_x.sum()+edge_y.sum())/(edge_x.numel()+edge_y.numel())
    return l1 + .1*edge, {"material_l1_normalized": l1, "edge_l1_normalized": edge}
