"""Public input boundary for the full-ring search and independent proof."""

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
PUBLIC = STUDY / 'synthetic/public'
CAPS = (0.5, 0.75, 1., 1.25, 1.5)


@dataclass(frozen=True)
class Observation:
    id: str
    object_id: str
    noise_pct: float
    sample_radius: float
    padded_radius: float
    y: np.ndarray


def full_ring_inputs():
    public = HERE / 'observations'
    with np.load(PUBLIC / 'fine_geometry.npz', allow_pickle=False) as f:
        spectrum = f['spectrum']
        sensor = f['PW360_10x20__S']
        incident = f['PW360_10x20__B']
    meta = json.loads((public / 'manifest.json').read_text())['observations']
    observations = []
    for row in meta:
        with np.load(public / row['file'], allow_pickle=False) as f:
            if set(f.files) != {'Y'}:
                raise ValueError('Certificate input must contain only Y.')
            y = f['Y']
        observations.append(Observation(row['id'], row['object_id'], row['noise_pct'],
                                        row['sample_radius'], row['padded_radius'], y))
    return spectrum, sensor, incident, observations
