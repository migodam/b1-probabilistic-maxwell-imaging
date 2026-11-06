"""Offline metrics with explicit denominators; no physical/optimization calls."""
import numpy as np
from research.research_max_struct_r7.experiments.score_primary import otsu

def material_scores(m,truth,amplitude):
    m=np.asarray(m);truth=np.asarray(truth);support=truth>0;seg,threshold=otsu(m)
    def iou(s):
        union=np.count_nonzero(s|support);return float(np.count_nonzero(s&support)/union) if union else 1.
    rmse=float(np.linalg.norm(m-truth)/np.sqrt(m.size));relative=float(np.linalg.norm(m-truth)/np.linalg.norm(truth))
    return dict(relative_l2_material_error=relative,amplitude_normalized_nrmse=rmse/float(amplitude),RMSE_contrast_units=rmse,
        otsu_iou=iou(seg),otsu_threshold=threshold,IoU_fixed_half_contrast=iou(m>.5*amplitude),IoU_fixed_absolute_contrast_0p3=iou(m>.3),
        material_min=float(m.min()),material_max=float(m.max()),material_legal=bool(np.all((m>=-1e-12)&(m<=20+1e-12))),
        normalization='relativeL2=norm(error)/norm(truth); amplitudeNRMSE=norm(error)/sqrt(N)/(epsilon_r-1); Otsu uses prediction-only256bins inheritedR7')
