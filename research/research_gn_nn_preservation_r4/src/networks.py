"""Residual U-Net and differentiable scalar-TM proposal heads.

The encoder/decoder structure is adapted from FAEMIS/faemis/models/unet.py
(DoubleConv, two down/up blocks), copied locally without an invasive import.
Channels are [m, gb, Re(Q[:,0:12]), Im(Q[:,0:12]), noise]. Spatial flattening
is row-major; currents have shape (..., cell, source), never source-major flat.
When flattened for a public receiver vector, use Y.T.reshape(-1).
"""
from __future__ import annotations

from pathlib import Path
import json
import time
import numpy as np
import torch
from torch import nn

ARMS = ("v_physics", "j_full", "j_null", "g_weak", "g_dark", "soft_dark")
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


def architecture_cost(arm, passes=2, dark_rank=67, alpha=.3):
    """Analytic convolutions and fixed complex products, excluding scalar work."""
    out = 1 if arm.startswith("v_") else 24
    macs = (24**2*(27*16*9+16*16*9) + 12**2*(16*32*9+32*32*9)
            + 6**2*(32*64*9+64*64*9) + 6**2*64*32*4
            + 12**2*(64*32*9+32*32*9) + 12**2*32*16*4
            + 24**2*(32*16*9+16*16*9) + 24**2*16*out)
    rank = {"j_null": 48, "g_weak": 19, "g_dark": dark_rank, "soft_dark": dark_rank}.get(arm, 0)
    current = arm != "v_physics"
    return {"network_conv_MACs_per_pass": macs, "network_conv_FLOPs": passes*2*macs,
            "fixed_complex_matmul_FLOPs": passes*((8*NC*NC*NS if current else 0)+16*NC*rank*NS),
            "raw_output_real_DOF": NC*out,
            "effective_output_real_DOF": (2*(NC-rank)*NS if arm in ("j_null", "g_weak", "g_dark") or (arm == "soft_dark" and alpha == 0) else NC*out),
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
    """Fixed public G/E/S, G-right V19 and S-right V48, uploaded once.

    W is an orthonormal basis of span(V19,V48), built using complex128 CPU
    SVD. All projections use Hermitian adjoints. MPS uses identical real pairs.
    """
    def __init__(self, path: str | Path):
        super().__init__()
        with np.load(path, allow_pickle=False) as p:
            g = np.asarray(p["G"])
            e = np.asarray(p["E0"] if "E0" in p else p["E"] if "E" in p else p["U0"])
            s = np.asarray(p["S"] if "S" in p else p["H"])
        if g.shape != (NC, NC) or e.shape != (NC, NS) or s.shape != (48, NC):
            raise ValueError("Expected G576x576 E576x12 S48x576")
        source = Path(path).resolve()
        cache_dir = source.parent.parent/"cache"
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path = cache_dir/"projectors.npz"
        stat = source.stat()
        signature = {"path": str(source), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
        setup_start = time.perf_counter()
        cache_hit = False
        if cache_path.exists():
            with np.load(cache_path, allow_pickle=False) as p:
                if json.loads(str(p["operator_signature"].item())) == signature:
                    v19, v48, w = (np.array(p[k]) for k in ("V19", "V48", "W"))
                    threshold = float(p["W_threshold"])
                    cache_hit = True
        if not cache_hit:
            _, _, gvh = np.linalg.svd(g.astype(np.complex128), full_matrices=False)
            _, singular, svh = np.linalg.svd(s.astype(np.complex128), full_matrices=False)
            if singular[-1] <= np.finfo(float).eps*max(s.shape)*singular[0]:
                raise ValueError("Expected full receiver rank48")
            v19, v48 = gvh[:19].conj().T, svh.conj().T
            uw, sw, _ = np.linalg.svd(np.concatenate((v19, v48), axis=1), full_matrices=False)
            threshold = np.finfo(float).eps*max(NC, 67)*sw[0]
            w = uw[:, sw > threshold]
            np.savez_compressed(cache_path, V19=v19, V48=v48, W=w,
                                W_threshold=threshold, operator_signature=json.dumps(signature))
        self.v48_source = "full_rank48_S_complex128_cpu_svd"
        self.projector_metadata = {"V19_source": "top19_G_right_singular_vectors", "V48_source": self.v48_source,
                                   "W_rank": int(w.shape[1]), "W_threshold": float(threshold),
                                   "fixed_operator_only": True, "weak_mode_truncation": "NONE",
                                   "cache_path": str(cache_path), "cache_hit": cache_hit}
        setup_receipt = {**self.projector_metadata, "setup_seconds": time.perf_counter()-setup_start,
                         "operator_signature": signature, "truth_access": "NONE"}
        (cache_dir/"projector_setup_latest.json").write_text(json.dumps(setup_receipt, indent=2)+"\n")
        if not cache_hit:
            (cache_dir/"projector_setup_creation.json").write_text(json.dumps(setup_receipt, indent=2)+"\n")
        self._complex_cache = {}
        for name, array in (("G", g), ("E", e), ("S", s), ("V19", v19), ("V48", v48), ("W", w)):
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

    def project(self, z, arm, alpha=.3):
        basis = {"j_null": "V48", "g_weak": "V19", "g_dark": "W", "soft_dark": "W"}.get(arm)
        if basis is None:
            return z
        bright = self._matmul(basis, self._matmul(basis, z, adjoint=True))
        return z-(1-float(alpha) if arm == "soft_dark" else 1)*bright

    def null_project(self, z):
        return self.project(z, "j_null")

    def total_field(self, j):
        e = self._complex("E") if j.is_complex() else torch.stack((self.E_re, self.E_im), dim=-1)
        return e+self._matmul("G", j)

    def material_from_current(self, j):
        q = self.total_field(j)
        if j.is_complex():
            numerator = (q.conj()*j).sum(dim=-1).real
            denominator = q.abs().square().sum(dim=-1).clamp_min(1e-20)
        else:
            numerator = (q*j).sum(dim=(-1, -2))
            denominator = q.square().sum(dim=(-1, -2)).clamp_min(1e-20)
        return (numerator/denominator).clamp(.02, 20.), q


def proposal(model, operators, features, m0, j0, stats, passes=2, alpha=.3):
    """Shared-weight two-pass cascade; gb/Q/noise and scales stay anchored."""
    if passes != 2:
        raise ValueError("R4 protocol requires exactly two shared-weight passes")
    if not 0 <= alpha <= 1:
        raise ValueError("Soft alpha must lie in [0,1]")
    m, j = m0, j0
    raw_outputs, stage_materials, stage_currents, stage_fields = [], [], [], []
    for _ in range(passes):
        x = features.clone()
        x[:, 0] = (m / MATERIAL_SCALE).reshape(-1, GRID, GRID)
        raw = model(x)
        raw_outputs.append(raw)
        if model.arm == "v_physics":
            m = (m+MATERIAL_SCALE*raw[:, 0].reshape(-1, NC)).clamp(.02, 20.)
        else:
            zr = raw[:, :NS].reshape(-1, NS, NC).transpose(1, 2)
            zi = raw[:, NS:].reshape(-1, NS, NC).transpose(1, 2)
            z = (torch.complex(zr, zi) if j.is_complex() else torch.stack((zr, zi), dim=-1))*float(stats["current_scale"])
            j = j+operators.project(z, model.arm, alpha)
            m, q = operators.material_from_current(j)
            stage_currents.append(j)
            stage_fields.append(q)
        stage_materials.append(m)
    out = {"m_proposal": m, "raw_delta": torch.stack(raw_outputs, dim=1),
           "stage_materials": torch.stack(stage_materials, dim=1)}
    if model.arm != "v_physics":
        out.update(j_aux=j, q_aux=q, delta_current=j-j0,
                   stage_currents=torch.stack(stage_currents, dim=1),
                   stage_fields=torch.stack(stage_fields, dim=1))
    return out


def public_regularizer(output, j0, q0_norm2, operators, loss_name):
    """Final-pass regularizers, with fixed per-observation public anchor scales."""
    if loss_name == "material":
        return torch.zeros((), device=j0.device)
    if "j_aux" not in output:
        raise ValueError("Const/internal regularizers require a current head")
    if loss_name == "const":
        q, j, m = output["q_aux"], output["j_aux"], output["m_proposal"]
        residual = j-m[..., None]*q if j.is_complex() else j-m[..., None, None]*q
        denominator = abs_square(j0).reshape(len(j0), -1).sum(dim=1)
    elif loss_name == "internal":
        residual = operators._matmul("G", output["delta_current"])
        denominator = q0_norm2
    else:
        raise ValueError(loss_name)
    numerator = abs_square(residual).reshape(len(j0), -1).sum(dim=1)
    return (numerator/denominator.clamp_min(1e-20)).mean()


def abs_square(z):
    return z.abs().square() if z.is_complex() else z.square().sum(dim=-1)


def material_loss(prediction, target):
    p = prediction.reshape(-1, GRID, GRID)/MATERIAL_SCALE
    t = target.reshape(-1, GRID, GRID)/MATERIAL_SCALE
    l1 = (p-t).abs().mean()
    edge_x = ((p[:, :, 1:]-p[:, :, :-1])-(t[:, :, 1:]-t[:, :, :-1])).abs()
    edge_y = ((p[:, 1:, :]-p[:, :-1, :])-(t[:, 1:, :]-t[:, :-1, :])).abs()
    edge = (edge_x.sum()+edge_y.sum())/(edge_x.numel()+edge_y.numel())
    return l1 + .1*edge, {"material_l1_normalized": l1, "edge_l1_normalized": edge}
