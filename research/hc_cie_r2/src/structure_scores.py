"""Post-seal amplitude, topology and rank-based shape metrics; never used for selection."""
import csv
from functools import lru_cache
from pathlib import Path
from time import perf_counter,process_time
import json
import numpy as np
from scipy.ndimage import uniform_filter,label
from imaging_score import get_case,otsu
from imaging_run import LAB

def local_ssim(a,b,amplitude):
    a=a.reshape(24,24);b=b.reshape(a.shape)
    ua=uniform_filter(a,size=7);ub=uniform_filter(b,size=7)
    va=48**-1*49*(uniform_filter(a*a,size=7)-ua*ua)
    vb=48**-1*49*(uniform_filter(b*b,size=7)-ub*ub)
    cov=48**-1*49*(uniform_filter(a*b,size=7)-ua*ub)
    c1=(.01*amplitude)**2;c2=(.03*amplitude)**2
    return float(np.mean((((2*ua*ub+c1)*(2*cov+c2))/((ua**2+ub**2+c1)*(va+vb+c2)))[3:-3,3:-3]))


def precision_area(a,fg):
    order=np.argsort(-a,kind="stable")
    last=np.flatnonzero(np.r_[a[order][:-1]!=a[order][1:],True])
    positives=np.cumsum(fg[order])[last]
    recall=np.r_[0.,positives/np.count_nonzero(fg)]
    precision=np.r_[1.,positives/(last+1)]
    return float(np.sum(np.diff(recall)*precision[1:])),float(np.trapezoid(precision,recall))


@lru_cache(None)
def truth(case):
    s=get_case(case);fg=s["support"].reshape(24,24)
    original=label(fg)[1]
    cuts=[]
    for i,j in np.argwhere(fg):
        reduced=fg.copy();reduced[i,j]=False
        if label(reduced)[1]>original:cuts.append(i*24+j)
    return s,np.asarray(cuts,dtype=int)


def main():
    ti,tc=perf_counter(),process_time()
    for v in ["v1","v2"]:
        base=LAB/"results/imaging"
        score=base/("score_"+v)
        with (score/"OUTCOMES.csv").open() as f:rows=list(csv.DictReader(f))
        last=None;out=[]
        for row in rows:
            folder=base/v/(row["case"]+"__"+row["endpoint"]+"__"+row["arm"])
            if folder!=last:
                with np.load(folder/"PREDICTIONS.npz") as z:mats=z["material"]
                last=folder
            a=mats[int(row["index"])];s,cuts=truth(row["case"]);b=s["material"];fg=s["support"];amp=s["amplitude"]
            ap,prauc=precision_area(a,fg)
            mask=a>1 if s["rule"]=="contrast > 1" else otsu(a)
            e=abs(a-b)
            out.append(dict(case=row["case"],endpoint=row["endpoint"],arm=row["arm"],source=row["source"],
                            index=row["index"],MAE=e.mean(),foreground_MAE=e[fg].mean(),background_MAE=e[~fg].mean(),
                            estimated_peak=a.max(),truth_peak=b.max(),estimated_roi_mean=a.mean(),truth_roi_mean=b.mean(),
                            estimated_target_mean=a[fg].mean(),truth_target_mean=b[fg].mean(),
                            SSIM7=local_ssim(a,b,amp),average_precision=ap,PR_AUC_trapezoid=prauc,
                            one_pixel_cut_vertices=len(cuts),fraction_cut_vertices_recovered=float(mask[cuts].mean()) if len(cuts) else None))
        with (score/"STRUCTURE.csv").open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
    (base/"structure_scoring.json").write_text(json.dumps(dict(offline_after_seal=True,wall_s=perf_counter()-ti,cpu_s=process_time()-tc),indent=2)+"\n")

if __name__=="__main__":main()
