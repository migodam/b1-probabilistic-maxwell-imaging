"""Lossless, periodic, critically sampled real 2D wavelet coordinates.

The final two axes are the spatial axes. Leading batch/channel axes are left
untouched. Analysis uses correlation with the explicitly documented filters;
it is not an assertion of PyWavelets' default coefficient phase convention.
No padding, resampling, coefficient selection, or learned parameter is used.
"""

from __future__ import annotations

import math
from collections.abc import Mapping

import torch
from torch import Tensor, nn


class WaveletPyramid2D(nn.Module):
    """Orthonormal Haar or db2 pyramid, defaulting to 24 -> 12 -> 6.

    ``analyze(x)`` returns ``ll2, lh2, hl2, hh2, lh1, hl1, hh1`` at
    the default two levels. In each two-letter name, the first letter refers
    to the row axis and the second to the column axis. For other valid level
    counts the same naming rule applies, ordered coarse to fine.

    Constant analysis matrices are cached per spatial size, device and dtype;
    synthesis applies their exact transpose. The cache is derived state, not
    a model parameter or checkpoint dependency. Both methods follow the input
    tensor's device/dtype, so no float64 constants are sent to MPS.
    """

    def __init__(self, kind: str = "haar", size: int = 24, levels: int = 2):
        super().__init__()
        if kind not in ("haar", "db2"):
            raise ValueError("kind must be 'haar' or 'db2'")
        if isinstance(size, bool) or not isinstance(size, int) or size <= 0:
            raise ValueError("size must be a positive integer")
        if isinstance(levels, bool) or not isinstance(levels, int) or levels <= 0:
            raise ValueError("levels must be a positive integer")
        if size % (2**levels):
            raise ValueError("size must be divisible by 2**levels (no padding)")
        self.kind, self.size, self.levels = kind, size, levels
        if kind == "haar":
            h = (1 / math.sqrt(2), 1 / math.sqrt(2))
        else:
            r3, denominator = math.sqrt(3), 4 * math.sqrt(2)
            h = tuple(v / denominator for v in (1 + r3, 3 + r3, 3 - r3, 1 - r3))
        self.low_filter = h
        self.high_filter = tuple((-1) ** t * h[-1 - t] for t in range(len(h)))
        self.band_order = (f"ll{levels}",) + tuple(
            f"{band}{level}"
            for level in range(levels, 0, -1)
            for band in ("lh", "hl", "hh")
        )
        self._matrix_cache: dict[tuple[int, torch.device, torch.dtype], Tensor] = {}

    @staticmethod
    def _validate_real(x: Tensor) -> None:
        if not isinstance(x, Tensor) or not x.is_floating_point() or x.is_complex():
            raise TypeError("wavelet inputs must be real floating-point torch tensors")

    def _matrix(self, n: int, like: Tensor) -> Tensor:
        """A[n,n]: low rows first, high rows second; periodic correlation."""
        key = (n, like.device, like.dtype)
        if key not in self._matrix_cache:
            # Keep cached constants usable by later training even if the first
            # transform call occurs inside a validation inference_mode block.
            with torch.inference_mode(False), torch.no_grad():
                # Accumulation handles periods shorter than filter support too.
                matrix = like.new_zeros((n, n))
                row = torch.arange(n // 2, device=like.device)
                for t, (low, high) in enumerate(zip(self.low_filter, self.high_filter)):
                    column = (2 * row + t) % n
                    matrix[row, column] += low
                    matrix[row + n // 2, column] += high
            self._matrix_cache[key] = matrix
        return self._matrix_cache[key]

    def analyze(self, x: Tensor) -> dict[str, Tensor]:
        """Return every wavelet coefficient without changing leading axes."""
        self._validate_real(x)
        if x.ndim < 2 or x.shape[-2:] != (self.size, self.size):
            raise ValueError(f"expected trailing shape {(self.size, self.size)}")
        details: dict[str, Tensor] = {}
        approximation = x
        for level in range(1, self.levels + 1):
            n = approximation.shape[-1]
            matrix = self._matrix(n, approximation)
            transformed = matrix @ approximation @ matrix.transpose(-2, -1)
            half = n // 2
            approximation = transformed[..., :half, :half]
            details[f"lh{level}"] = transformed[..., :half, half:]
            details[f"hl{level}"] = transformed[..., half:, :half]
            details[f"hh{level}"] = transformed[..., half:, half:]
        bands = {f"ll{self.levels}": approximation, **details}
        return {name: bands[name] for name in self.band_order}

    def synthesize(self, bands: Mapping[str, Tensor]) -> Tensor:
        """Reconstruct all pixels by the transpose of the analysis transform."""
        if set(bands) != set(self.band_order):
            raise ValueError(f"expected exactly these bands: {self.band_order}")
        approximation = bands[f"ll{self.levels}"]
        self._validate_real(approximation)
        leading = approximation.shape[:-2]
        for name in self.band_order:
            band = bands[name]
            self._validate_real(band)
            level = int(name[2:])
            side = self.size // (2**level)
            if band.shape != (*leading, side, side):
                raise ValueError(f"{name} must have shape {(*leading, side, side)}")
            if band.dtype != approximation.dtype or band.device != approximation.device:
                raise ValueError("all bands must have the same device and dtype")
        for level in range(self.levels, 0, -1):
            upper = torch.cat((approximation, bands[f"lh{level}"]), dim=-1)
            lower = torch.cat((bands[f"hl{level}"], bands[f"hh{level}"]), dim=-1)
            transformed = torch.cat((upper, lower), dim=-2)
            matrix = self._matrix(transformed.shape[-1], transformed)
            approximation = matrix.transpose(-2, -1) @ transformed @ matrix
        return approximation

    def forward(self, x: Tensor) -> dict[str, Tensor]:
        return self.analyze(x)

    def coefficient_vector(self, x: Tensor) -> Tensor:
        """Concatenate coarse-to-fine bands to shape ``(*x.shape[:-2], size**2)``."""
        return torch.cat([band.flatten(-2) for band in self.analyze(x).values()], dim=-1)

    def energy_by_band(self, error: Tensor) -> dict[str, Tensor]:
        """Squared coefficient sums over spatial axes, preserving leading axes.

        Summing the returned bands equals ``error.square().sum((-2, -1))``
        to floating-point precision. These are energies, not averages/fractions.
        """
        return {name: band.square().sum(dim=(-2, -1))
                for name, band in self.analyze(error).items()}
