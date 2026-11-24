"""Independent readouts of single-pixel and regional shared-material consistency."""
import numpy as np


def imag_h(value):
    return (value - value.conj().T) / (2j)


def local_readout(j, e, b, side, cap=20.0):
    """Scale each source by its incident rms. Evaluate unconstrained current iterates."""
    n, sources = e.shape
    gains = np.linalg.norm(b, axis=0) / np.sqrt(n)
    if np.any(gains <= 0) or len(j) != side**2:
        raise ValueError("Nonzero calibrated illuminations and an explicitly aligned grid are required.")
    j, e = j / gains, e / gains
    e2, j2 = abs(e)**2, abs(j)**2
    dot = (j.conj() * e).sum(axis=1)
    pixel_norm_sq = (j2.sum(axis=1) * e2.sum(axis=1) - (dot**2).real) / 2
    sum_local = np.sum(np.sqrt(np.maximum(pixel_norm_sq, 0)))
    q = imag_h(j.conj().T @ e)
    global_norm = np.linalg.norm(q)
    common = np.clip(dot.real / np.maximum(e2.sum(axis=1), 1e-300), 0., cap)
    separate = np.clip((j.conj() * e).real / np.maximum(e2, 1e-300), 0., cap)
    shared_error = np.linalg.norm(j - common[:, None] * e)**2
    single_error = np.linalg.norm(j - separate * e)**2
    yy, xx = np.indices((side, side))
    result = dict(
        source_count=sources, global_optical_norm=float(global_norm),
        sum_of_pixel_optical_norms=float(sum_local),
        shared_material_defect_squared=float(shared_error),
        separate_material_defect_squared=float(single_error),
        cost_of_sharing_fraction=float((shared_error-single_error)/shared_error) if shared_error else None,
        coherent_fraction=float(global_norm / sum_local) if sum_local else None,
    )
    for k in (2, 4):
        blocks = ((yy*k//side)*k+xx*k//side).ravel()
        total = 0.
        for region in range(k*k):
            mask = blocks == region
            total += np.linalg.norm(imag_h(j[mask].conj().T @ e[mask]))
        result[f'fraction_visible_in_{k*k}_regions'] = float(total / sum_local) if sum_local else None
    return result
