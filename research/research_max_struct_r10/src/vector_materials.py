"""Design-only exact box fractions. Called only by a sealed-data generator.

No files, arrays, private targets, or random realizations are created on import.
The inverse runner must never import this module or receive design terms.
"""
import numpy as np


def box_fraction(points, h, minimum, maximum):
    lo, hi = np.asarray(minimum), np.asarray(maximum)
    if np.any(hi <= lo):
        raise ValueError('Empty box term')
    lengths = np.maximum(0., np.minimum(points+h/2, hi)-np.maximum(points-h/2, lo))
    return np.prod(lengths / h, axis=1)


def material_from_design(design, points, h):
    """Exact linear integral of disjoint/CSG inclusion-exclusion box terms.

    Negative coefficients represent declared cavity/replacement terms, not
    negative final physical material. This returns an N-real contrast vector.
    """
    material = np.zeros(len(points), dtype=float)
    for term in design['terms']:
        material += term['contrast_weight'] * box_fraction(points, h, term['minimum'], term['maximum'])
    if material.min() < -1e-12 or material.max() > 8+1e-12:
        raise ValueError('CSG design escaped shared chi box [0,8]')
    return np.clip(material, 0., 8.)
