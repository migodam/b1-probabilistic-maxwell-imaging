"""Summarize archived R15 comparisons and draw the equal-energy broadband figure."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT/"results"
HOSTS = {"holdout":("Lossless","o","-"),"holdout_debye":("Debye","s","--")}
DESIGNS = ("f1_t8_r8","f2_t4_r8","f4_t2_r8")
METHODS = {"fixed":("Fixed current","#777779"),"fft_joint":("Joint Fourier current","#2d7490"),
           "full_gn":("Full-wave voxel GN","#a9651f")}


def median(val):
    return float(np.median(val))


def run():
    inp=json.loads((RESULTS/"broadband_current.json").read_text())
    gt=json.loads((RESULTS/"broadband_current_scores.json").read_text())
    flat=json.loads((RESULTS/"broadband_current_flat.json").read_text())
    threshold = float(np.sqrt(chi2.ppf(.99,128)/128))
    result = dict(physical_residual_contract_noise_ratio=threshold,comparisons={})
    for tag,d in inp["designs"].items():
        cases=[(i,r) for i,r in d["cases"].items() if not gt[f"{tag}/{i}"]["is_empty"]]
        empty=next((i,r) for i,r in d["cases"].items() if gt[f"{tag}/{i}"]["is_empty"])
        assert len(cases)==8
        rows={}
        t=[r["weighted_trace"]["full_gn_fit"] for _,r in cases]
        d_retained=np.concatenate([r["random_direction_retained_fractions"] for r in t])
        for method in METHODS:
            records=[d["outputs"][method] for _,d in cases]
            scores=[gt[f"{tag}/{i}"]["methods"][method] for i,_ in cases]
            rows[method]=dict(
                RMSE=median([v["RMSE"] for v in scores]),
                signed_IoU=median([v["signed_support_IoU"] for v in scores]),
                CPU_seconds=median([v["cpu"] for v in records]),
                residual_noise_ratio=median([v["full_residual_noise_rms"] for v in records]),
                passes_99pct_contract=int(sum(v["full_residual_noise_rms"]<=threshold for v in records)),
                full_rhs=median([v["online_counts"].get("full_complex_rhs",0) for v in records]),
                reduced_rhs=median([v["online_counts"].get("reduced_complex_rhs",0) for v in records]),
                blank_support_cells=gt[f"{tag}/{empty[0]}"]["methods"][method]["false_support_pixels"])
        row=dict(nonempty=len(cases),methods=rows,
                 median_weighted_trace=median([r["trace"] for r in t]),
                 weighted_trace_range=[min(r["trace"] for r in t),max(r["trace"] for r in t)],
                 random_retention_q10_q50_q90=np.quantile(d_retained,[.1,.5,.9]).tolist(),
                 median_of_cell_min_frequency_weights=np.median(
                    [r["frequency_illumination_fraction_min"] for r in t],axis=0).tolist(),
                 median_of_cell_max_frequency_weights=np.median(
                    [r["frequency_illumination_fraction_max"] for r in t],axis=0).tolist(),
                 truth_residual_noise_ratio=median([gt[f"{tag}/{i}"]["truth_residual_noise_rms"] for i,_ in cases]))
        if len(d["freqs"])>1:
            alt=[flat[f"{tag}/{i}"] for i,_ in cases]
            row["equal_constitutive_weight_control"]=dict(
                RMSE=median([x["metrics"]["RMSE"] for x in alt]),
                signed_IoU=median([x["metrics"]["signed_support_IoU"] for x in alt]),
                residual_noise_ratio=median([x["physical_noise_ratio"] for x in alt]),
                CPU_seconds=median([x["cpu"] for x in alt]))
        result["comparisons"][tag]=row
    result["maximum_trace_diagonal_error"]=max(t["diagonal_error"] for d in inp["designs"].values()
        for r in d["cases"].values() for t in r["weighted_trace"].values())
    result["maximum_random_projection_error"]=max(t["quadratic_error"] for d in inp["designs"].values()
        for r in d["cases"].values() for t in r["weighted_trace"].values())
    (RESULTS/"broadband_current_summary.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    fig,axes=plt.subplots(1,4,figsize=(13.1,3.45),layout="constrained")
    x=np.arange(3)
    for pack,(host,marker,style) in HOSTS.items():
        for method,(name,color) in METHODS.items():
            series=[result["comparisons"][f"{pack}/{key}"]["methods"][method] for key in DESIGNS]
            for ax,metric in zip(axes,("RMSE","signed_IoU","CPU_seconds")):
                ax.plot(x,[r[metric] for r in series],marker=marker,linestyle=style,linewidth=1.5,
                        color=color,markersize=4,label=name if pack=="holdout" else None)
        rows=[result["comparisons"][f"{pack}/{key}"] for key in DESIGNS]
        spread=np.array([row["random_retention_q10_q50_q90"] for row in rows])
        dx=-.06 if pack=="holdout" else .06
        color="#2d7490" if pack=="holdout" else "#388866"
        axes[3].fill_between(x+dx,spread[:,0],spread[:,2],color=color,alpha=.13)
        axes[3].plot(x+dx,spread[:,1],marker=marker,color=color,ls=style,linewidth=1.5,
                     markersize=4,label=host+" random")
    axes[0].set(title="(a) Material RMSE",ylabel="Permittivity contrast RMSE")
    axes[1].set(title="(b) Defect support",ylabel="Signed IoU",ylim=(-.04,.9))
    axes[2].set(title="(c) Online cost",ylabel="CPU time (s)",yscale="log")
    axes[3].axhline(np.mean([r["median_weighted_trace"] for r in result["comparisons"].values()])/576,
                    color="#777779",lw=1,ls=":",label=r"Mean $\mathrm{tr}(K)/576$")
    axes[3].set(title="(d) Fitted current geometry",ylabel="Fraction of random tangent retained")
    axes[3].yaxis.set_major_formatter(PercentFormatter(1,decimals=1))
    axes[3].legend(loc="lower left",fontsize=6.8,frameon=False)
    for ax in axes:
        ax.set_xticks(x,["1 f × 8 Tx","2 f × 4 Tx","4 f × 2 Tx"])
        ax.tick_params(axis="x",labelsize=7.6)
        ax.spines[["top","right"]].set_visible(False)
    handles,labels=axes[0].get_legend_handles_labels()
    base = [plt.Line2D([0],[0],color="#555",marker=m,ls=ls,markersize=4,label=host)
            for host,m,ls in HOSTS.values()]
    fig.legend(handles+base,labels+[r.get_label() for r in base],loc="outside lower center",
               ncol=5,frameon=False,fontsize=8)
    fig.suptitle("Archived finite-plate holdout · eight geometries per host · eight shots and 64 complex samples",
                 fontsize=10.5)
    for suffix in ("png","pdf"):
        fig.savefig(ROOT/"figures"/f"broadband_current.{suffix}",dpi=205,bbox_inches="tight")
    plt.close(fig)


if __name__=="__main__":
    run()
