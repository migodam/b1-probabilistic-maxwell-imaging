"""Aggregate and plot the alternating-current controls."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT/"results"
COHORTS = ("r10_twod","r12_confirm24")
NAMES = {
    "fixed_current": "Fixed current",
    "fft_dark": "FFT dark",
    "fft_bright": "FFT bright + dark",
    "full_current": "Full current",
    "voxel_only": "Full-wave voxel GN",
}
COLORS = ("#7a7a7a","#2d7490","#388866","#a9651f","#91435c")


def main():
    groups = {
        "metric_definitions": {
            "relativeL2": "L2 norm of the full material error divided by the L2 norm of truth",
            "amplitude_NRMSE": "RMSE divided by (maximum ground-truth contrast minus 0.1)",
            "fixed_support_IoU": "Intersection over union of material > 1.0 in prediction and truth",
            "half_peak_support_IoU": "Historical auxiliary IoU threshold max(maximum truth / 2, 0.02)",
        }
    }
    distributions = {}
    for family in COHORTS:
        data = json.loads((RESULTS/f"alternating_{family}.json").read_text())
        score = json.loads((RESULTS/f"alternating_{family}_scores.json").read_text())
        for method in NAMES:
            rows = [r for r in data["cases"] if r["method"]==method]
            metric = [score[r["array"]] for r in rows]
            error = np.array([r["relativeL2"] for r in metric])
            amplitude = np.array([r["amplitude_NRMSE"] for r in metric])
            cost = np.array([r["online_cost"]["cpu"] for r in rows])
            residual = np.array([r["residual_over_sigma_rms"] for r in rows])
            distributions[family,method] = amplitude,cost,residual
            groups.setdefault(family,{})[method] = dict(
                cases=len(rows),relativeL2_median=float(np.median(error)),
                relativeL2_q90=float(np.quantile(error,.9)),
                amplitude_NRMSE_median=float(np.median(amplitude)),
                amplitude_NRMSE_q90=float(np.quantile(amplitude,.9)),
                fixed_support_IoU_median=float(np.median([r["fixed_support_IoU"] for r in metric])),
                half_peak_support_IoU_median=float(np.median([r["half_peak_support_IoU"] for r in metric])),
                cpu=float(np.median(cost)),physical_residual=float(np.median(residual)),
                noise_compatible=sum(r["residual_over_sigma_rms"]<=r["near_noise_threshold"] for r in rows),
                relativeL2_lt_quarter=int(np.sum(error<.25)),
                current_rank=sorted({r["current_dimension"] for r in rows}),)
            if method in ("fft_dark","fft_bright","full_current"):
                groups[family][method]["largest_first_current_KKT"] = max(
                    r["first_current_projected_kkt_relative"] for r in rows)
                groups[family][method]["monotone_auxiliary_count"] = sum(
                    r["monotone_objective"] for r in rows)

    (RESULTS/"alternating_summary.json").write_text(json.dumps(groups,indent=2,allow_nan=False)+"\n")
    fig,axes = plt.subplots(1,3,figsize=(11.5,3.45),layout="constrained")
    xmax = 1.035*max(float(np.max(values[0])) for values in distributions.values())
    medians = [groups[f][m]["amplitude_NRMSE_median"] for f in COHORTS for m in NAMES]
    for idx,family in enumerate(COHORTS):
        for method,c in zip(NAMES,COLORS,strict=True):
            error,cost,_ = distributions[family,method]
            xx, yy = np.sort(error),np.arange(1,len(error)+1)/len(error)
            axes[idx].step(xx,yy,where="post",color=c,lw=1.75,label=NAMES[method])
        axes[idx].set(title=("(a) R10 two-dimensional","(b) R12 confirmation")[idx],
                      xlabel="Amplitude NRMSE",ylabel="Fraction at or below error" if idx==0 else "",
                      xlim=(0.,xmax),ylim=(0.,1.02))
    for method,c in zip(NAMES,COLORS,strict=True):
        for idx,family in enumerate(COHORTS):
            row = groups[family][method]
            axes[2].scatter(row["cpu"],row["amplitude_NRMSE_median"],color=c,marker=("o","s")[idx],s=55)
    axes[2].set(title="(c) Measured cost and error",xscale="log",xlim=(.006,1.0),
                ylim=(.84*min(medians),1.08*max(medians)),
                xlabel="Median online CPU time (s)",ylabel="Median amplitude NRMSE")
    axes[2].scatter([],[],marker="o",s=35,color="black",label="R10")
    axes[2].scatter([],[],marker="s",s=35,color="black",label="R12")
    axes[2].legend(loc="lower left",fontsize=7,frameon=False)
    line,labels = axes[0].get_legend_handles_labels()
    fig.legend(line,labels,loc="outside lower center",ncol=5,frameon=False,fontsize=8)
    for ax in axes:
        ax.spines[["top","right"]].set_visible(False)
    for suffix in (".png",".pdf"):
        fig.savefig(ROOT/"figures"/f"alternating_benchmark{suffix}",dpi=180,bbox_inches="tight")
    plt.close(fig)


if __name__=="__main__":
    main()
