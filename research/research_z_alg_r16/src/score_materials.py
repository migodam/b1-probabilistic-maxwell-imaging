"""Offline material metrics for predictions saved by the R16 solvers."""

import numpy as np
from scipy.ndimage import uniform_filter


def structural_similarity(x, target, data_range):
    """Mean SSIM over complete 7×7 windows using sample covariance."""
    window = 7
    crop = np.s_[3:-3, 3:-3]
    if x.ndim != 2 or min(x.shape) < window or x.shape != target.shape:
        raise ValueError("SSIM requires identically shaped 2D images of at least 7×7")
    mx = uniform_filter(x, size=window)[crop]
    mt = uniform_filter(target, size=window)[crop]
    sample_correction = window**2 / (window**2 - 1)
    vx = sample_correction * (uniform_filter(x * x, size=window)[crop] - mx**2)
    vt = sample_correction * (uniform_filter(target * target, size=window)[crop] - mt**2)
    cov = sample_correction * (uniform_filter(x * target, size=window)[crop] - mx * mt)
    c1, c2 = (0.01 * data_range) ** 2, (0.03 * data_range) ** 2
    return float(np.mean(
        ((2 * mx * mt + c1) * (2 * cov + c2))
        / ((mx**2 + mt**2 + c1) * (vx + vt + c2))
    ))


def material_metrics(prediction, truth, shape=(24, 24), bounds=(0.02, 20.0)):
    """Score the historical amplitude NRMSE and fixed support IoU plus R16 SSIM."""
    prediction, truth = np.asarray(prediction), np.asarray(truth)
    if prediction.size != np.prod(shape) or truth.size != prediction.size:
        raise ValueError("Prediction and truth must match the scoring grid")
    if np.iscomplexobj(prediction) or np.iscomplexobj(truth):
        raise ValueError("These cohorts use real materials")
    p = prediction.astype(np.float64).reshape(shape)
    t = truth.astype(np.float64).reshape(shape)
    if not np.isfinite(p).all() or not np.isfinite(t).all():
        raise ValueError("Prediction and truth must be finite")
    amplitude = float(t.max() - 0.1)
    if amplitude <= 0:
        raise ValueError("The historical amplitude-NRMSE denominator must be positive")

    error = np.linalg.norm(p - t)
    rmse = float(error / np.sqrt(t.size))
    nrmse = rmse / amplitude
    p_support, t_support = p > 1.0, t > 1.0
    union = np.count_nonzero(p_support | t_support)
    iou = float(np.count_nonzero(p_support & t_support) / union) if union else 1.0
    lo, hi = bounds
    legal = bool(p.min() >= lo - 1e-12 and p.max() <= hi + 1e-12)
    return {
        "RMSE": rmse,
        "amplitude_NRMSE": nrmse,
        "relativeL2": float(error / np.linalg.norm(t)),
        "fixed_support_IoU": iou,
        "SSIM_fixed_range": structural_similarity(p, t, hi - lo),
        "material_legal": legal,
        "material_success": bool(legal and nrmse <= 0.15 and iou >= 0.7),
    }


def legacy_r3_material_metrics(prediction, truth, support):
    """MAX-STRUCT R2/R3 amplitude-12 NRMSE and prediction-only Otsu IoU."""
    p, t = np.asarray(prediction).ravel(), np.asarray(truth).ravel()
    support = np.asarray(support, dtype=bool).ravel()
    if p.shape != (576,) or p.shape != t.shape or p.shape != support.shape:
        raise ValueError("The historical R3 scorer uses the original 24×24 grid")
    if np.iscomplexobj(p) or np.iscomplexobj(t) or not np.isfinite(p).all() or not np.isfinite(t).all():
        raise ValueError("Prediction and truth must be finite and real")

    predicted_support = np.zeros(p.size, dtype=bool)
    if np.ptp(p) > 1e-12:
        counts, edges = np.histogram(p, bins=256)
        weights = np.cumsum(counts)
        moments = np.cumsum(counts * (edges[:-1] + edges[1:]) / 2)
        denominator = weights * (weights[-1] - weights)
        valid = denominator > 0
        scores = np.full(counts.size, -np.inf)
        scores[valid] = (moments[-1] * weights[valid] - moments[valid] * weights[-1])**2 / denominator[valid]
        index = int(np.argmax(scores))
        predicted_support = p > (edges[index] + edges[index + 1]) / 2

    union = np.count_nonzero(predicted_support | support)
    return {
        "legacy_NRMSE_amplitude12": float(np.linalg.norm(p - t) / (12 * 24)),
        "legacy_Otsu_IoU": float(np.count_nonzero(predicted_support & support) / union) if union else 1.0,
    }
