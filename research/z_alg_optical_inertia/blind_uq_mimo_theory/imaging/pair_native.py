"""Full-wave positive ROI counterpart to a sealed native zero witness."""
from __future__ import annotations

import argparse
import numpy as np
from threadpoolctl import threadpool_limits

from shared import Model, OUT, elapsed, now, write_arrays, write_json


def pair(case, roi_id):
    start = now()
    model = Model.read(case)
    path = OUT/case
    k = int(roi_id.rsplit("_",1)[1])-1
    with np.load(path/"baseline.npz",allow_pickle=False) as z:
        roi = z["rois"][k].ravel()
        template = z["chi"].ravel()[roi]
    with np.load(path/"b1_predictions.npz",allow_pickle=False) as z:
        zero = z["x_"+roi_id].ravel()
    if np.any(zero[roi] != 0):
        raise ValueError("Source witness violates zero-region contract")
    zero_info = model.audit(zero)
    if zero_info["residual_over_radius"] >= 1:
        raise ValueError("Zero-region source witness is outside the original training ball")
    zero_current = model.forward(zero)["j"].copy()
    proposals, pair = [], None
    for scale in (1., .5, .25, .125, .0625, .03125, .015625):
        x = zero.copy()
        x[roi] = scale*template
        result = model.audit(x, amplitude_fraction=scale,
                             roi_mean=float(x[roi].mean()), roi_peak=float(x[roi].max()))
        proposals.append(result)
        if result["residual_over_radius"] < .98:
            pair = result
            write_arrays(path/f"native_pair_{roi_id}.npz", zero=zero.reshape(model.shape),
                         positive=x.reshape(model.shape), roi=roi.reshape(model.shape),
                         zero_current=zero_current, positive_current=model.forward(x)["j"])
            break
    report = dict(case=case, roi=roi_id, zero=zero_info, positive=pair,
                  proposals=proposals, counts=dict(model.cost), **elapsed(start))
    write_json(path/f"native_pair_{roi_id}.json",report)
    print(case, roi_id, pair, flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("roi",help="Frozen ROI id such as roi_03")
    args=ap.parse_args()
    with threadpool_limits(limits=1,user_api="blas"):
        pair(args.case,args.roi)
