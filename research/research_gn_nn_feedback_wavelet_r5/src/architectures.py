"""Small circular CNN priors and bounded SPD wavelet/pixel metrics for R5.

All seven wavelet bands are retained. Detail CNNs share weights across the
three orientations within each scale, with three one-hot orientation maps.
There is no sharing between scales or with LL2. Macro weights are shared at
every iteration; there is no normalization, dropout, recurrent cell or resize.
"""
from __future__ import annotations

import torch
from torch import Tensor, nn

try:
    from .wavelets import WaveletPyramid2D
except ImportError:
    from wavelets import WaveletPyramid2D


class CircularCNN(nn.Module):
    def __init__(self, inputs: int, width: int):
        super().__init__()
        self.layers = nn.Sequential(
            nn.Conv2d(inputs, width, 3, padding=1, padding_mode="circular"),
            nn.SiLU(),
            nn.Conv2d(width, width, 3, padding=1, padding_mode="circular"),
            nn.SiLU(),
            nn.Conv2d(width, 1, 3, padding=1, padding_mode="circular"),
        )
        nn.init.zeros_(self.layers[-1].weight)
        nn.init.zeros_(self.layers[-1].bias)

    def forward(self, x: Tensor) -> Tensor:
        return self.layers(x)


class SubbandCNN(nn.Module):
    """One output map per retained band, using same-scale detail sharing."""
    def __init__(self, inputs: int, width: int):
        super().__init__()
        self.low = CircularCNN(inputs, width)
        self.detail = nn.ModuleDict({str(level): CircularCNN(inputs + 3, width)
                                     for level in (1, 2)})

    def forward(self, bands: dict[str, Tensor]) -> dict[str, Tensor]:
        output = {"ll2": self.low(bands["ll2"])}
        for level in (2, 1):
            for orientation, name in enumerate(("lh", "hl", "hh")):
                key = f"{name}{level}"
                x = bands[key]
                indicator = x.new_zeros((len(x), 3, *x.shape[-2:]))
                indicator[:, orientation] = 1
                output[key] = self.detail[str(level)](torch.cat((x, indicator), dim=1))
        return output


class PriorCNN(nn.Module):
    """27 fixed/public feature maps -> one full-resolution material residual."""
    def __init__(self, representation: str):
        super().__init__()
        self.representation = representation
        if representation == "pixel":
            self.wavelet = None
            self.cnn = CircularCNN(27, 38)
        elif representation in ("haar", "db2"):
            self.wavelet = WaveletPyramid2D(representation)
            self.cnn = SubbandCNN(27, 16)
        else:
            raise ValueError("representation must be pixel, haar or db2")

    def forward(self, features: Tensor) -> Tensor:
        if features.ndim != 4 or features.shape[1:] != (27, 24, 24):
            raise ValueError("prior features must have shape (B,27,24,24)")
        if self.wavelet is None:
            return self.cnn(features)
        return self.wavelet.synthesize(self.cnn(self.wavelet.analyze(features)))


def current_to_maps(current: Tensor) -> Tensor:
    """Source-as-columns (B,576,S[,2]) -> real maps (B,2*S,24,24)."""
    if current.is_complex():
        real, imag = current.real, current.imag
    else:
        if current.ndim != 4 or current.shape[-1] != 2:
            raise ValueError("real current must have shape (B,576,S,2)")
        real, imag = current[..., 0], current[..., 1]
    return torch.cat((real.transpose(1, 2), imag.transpose(1, 2)), dim=1).reshape(
        len(current), -1, 24, 24)


def maps_to_current(maps: Tensor, like: Tensor) -> Tensor:
    sources = like.shape[2]
    real, imag = maps.flatten(-2).split(sources, dim=1)
    real, imag = real.transpose(1, 2), imag.transpose(1, 2)
    return torch.complex(real, imag) if like.is_complex() else torch.stack((real, imag), dim=-1)


class MetricCNN(nn.Module):
    """P=W.T diag(p) W, p=1+.5*tanh(raw), shared across all source/RI maps.

    The pixel metric uses W=I. Every coefficient gain is strictly between .5
    and 1.5 in exact arithmetic (floating saturation may attain the endpoints).
    """
    def __init__(self, representation: str):
        super().__init__()
        self.representation = representation
        if representation == "pixel":
            self.wavelet = None
            self.cnn = CircularCNN(7, 8)
        elif representation in ("haar", "db2"):
            self.wavelet = WaveletPyramid2D(representation)
            self.cnn = SubbandCNN(7, 8)
        else:
            raise ValueError("representation must be pixel, haar or db2")

    def forward(self, features: Tensor, gradient: Tensor) -> tuple[Tensor, dict[str, Tensor]]:
        if features.shape[1:] != (7, 24, 24):
            raise ValueError("metric features must have shape (B,7,24,24)")
        maps = current_to_maps(gradient)
        if self.wavelet is None:
            gains = {"pixel": 1 + .5 * torch.tanh(self.cnn(features))}
            transformed = maps * gains["pixel"]
        else:
            raw = self.cnn(self.wavelet.analyze(features))
            gains = {name: 1 + .5 * torch.tanh(value) for name, value in raw.items()}
            coefficients = self.wavelet.analyze(maps)
            transformed = self.wavelet.synthesize(
                {name: value * gains[name] for name, value in coefficients.items()})
        return maps_to_current(transformed, gradient), gains


def parameter_counts(module: nn.Module) -> dict[str, int]:
    """Actual registered trainable parameter counts; no nominal capacity proxy."""
    return {"total": sum(p.numel() for p in module.parameters() if p.requires_grad),
            **{name: sum(p.numel() for p in child.parameters() if p.requires_grad)
               for name, child in module.named_children()}}
