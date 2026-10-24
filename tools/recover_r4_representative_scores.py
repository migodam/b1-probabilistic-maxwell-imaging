"""Additive offline repair: coerce the saved decoder matrix to an ndarray.

The original scorer and every frozen physical/inference artifact stay unchanged.
No model, candidate, physical query, threshold, or selection rule is changed.
Use a fresh output directory and a separately charged analysis receipt.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'tools'))
import numpy as np
import p3r4.probability_eval as probability_eval
import r4_fresh_representatives as original

_decoder = probability_eval.decoder_from_geometry


def ndarray_geometry(*args, **kwargs):
    geometry = _decoder(*args, **kwargs)
    geometry['A'] = np.asarray(geometry['A'], dtype=float)
    return geometry


if __name__ == '__main__':
    probability_eval.decoder_from_geometry = ndarray_geometry
    original.main()
