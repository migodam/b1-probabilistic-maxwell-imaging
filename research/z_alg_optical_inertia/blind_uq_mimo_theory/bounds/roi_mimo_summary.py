"""Report attained hulls without treating incomplete physical searches as exclusions."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
COUNTS=(6,8,12,16)
ERRORS=(0.,.005,.025,.06)


def main():
    fit=HERE/"roi_mimo_fits"
    groups=[]
    for q in COUNTS:
        directory=fit/f"n32_q{q:02}"
        cases=[]
        for path in sorted(directory.glob("n32_c*.json")):
            data=json.loads(path.read_text())
            cases.append(dict(case=data["case"],roi=data["roi"],cells=data["region_cells"],
                width=data["nominal_inner_width"],inner=data["nominal_inner_peak_range"],
                residual_ratio_upper={k:(a.get("residual_ratio_upper") if a else None)
                                      for k,a in data["nominal_endpoint_witnesses"].items()},
                nominal_wall_s=data["wall_s"],nominal_cpu_s=data["cpu_s"],power_W=data["incident_available_W"],
                paired_evidence=str(path.relative_to(HERE)),replays=[dict(
                    eta=z["calibration_level"],inner=z["inner_peak_range"],width=z["inner_width"],
                    wall_s=z["wall_s"],cpu_s=z["cpu_s"],hard_r=z["measurement_noise_radius"])
                    for z in data.get("calibration_replay",[])]))
        if cases:
            phases=[]
            for error in ERRORS:
                values=[a if error==0 else next((z for z in a["replays"] if z["eta"]==error),
                                                dict(width=None)) for a in cases]
                widths=[z["width"] for z in values if z["width"] is not None]
                phases.append(dict(eta=error,cases=len(cases),verified_full_width=sum(w==4.5 for w in widths),
                    verified_nonzero_width=sum(w>0 for w in widths),verified_any=len(widths),
                    verified_min_width=min(widths,default=None),
                    verified_median_width=float(np.median(widths)) if widths else None))
            setup=json.loads((directory/"setup.json").read_text())
            groups.append(dict(q=q,matrices=cases,phases=phases,setup_s=setup["setup_wall_s"],
                 setup_cpu_s=setup["setup_cpu_s"],nominal_wall_s=sum(a["nominal_wall_s"] for a in cases),
                 nominal_cpu_s=sum(a["nominal_cpu_s"] for a in cases),
                 replay_wall_s=sum(z["wall_s"] for a in cases for z in a["replays"]),
                 replay_cpu_s=sum(z["cpu_s"] for a in cases for z in a["replays"])))
    extras=[]
    for path in sorted((fit/"extended_budget").glob("n32_q*/n32_c*.json")):
        d=json.loads(path.read_text())
        extras.append(dict(case=d["case"],q=d["q"],path=str(path.relative_to(HERE)),
                           nominal_inner=d["nominal_inner_peak_range"],nominal_width=d["nominal_inner_width"],
                           cpu_s=d["cpu_s"],wall_s=d["wall_s"],
                           replays=[dict(eta=v["calibration_level"],inner=v["inner_peak_range"],width=v["inner_width"])
                                    for v in d["calibration_replay"]]))
    summary=dict(frequency_MHz=400,radiation_scale_ohms=48,selective_followups=extras,
                 outer_peak_interval=[0,4.5],matched_physical_replay_known_calibration=True,
                 endpoint_hull_only=True,methods=groups,
                 independent_materials=len({a["case"] for g in groups for a in g["matrices"]}),
                 nominal_rows=sum(len(g["matrices"]) for g in groups),
                 replay_rows=sum(len(a["replays"]) for g in groups for a in g["matrices"]))
    out=HERE/"roi_mimo_summary.json"
    out.write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    if summary["nominal_rows"]!=32 or summary["replay_rows"]!=96:
        return
    figure,axes=plt.subplots(2,2,figsize=(12.2,6.2),layout="constrained",sharex=True,sharey=True)
    cbar=None
    for ax,err in zip(axes.flat,ERRORS,strict=True):
        val=np.full((len(groups),8),np.nan)
        extra_used=np.zeros_like(val,dtype=bool)
        for i,g in enumerate(groups):
            for j,case in enumerate(g["matrices"]):
                v=case["width"] if err==0 else next(row["width"] for row in case["replays"] if row["eta"]==err)
                for extra in extras:
                    if extra["q"]!=g["q"] or extra["case"]!=case["case"]:
                        continue
                    extension=extra["nominal_width"] if err==0 else next(z["width"] for z in extra["replays"] if z["eta"]==err)
                    if extension is not None and (v is None or extension>v):
                        v=extension
                        extra_used[i,j]=True
                if v is not None:
                    val[i,j]=v
        cbar=ax.imshow(val,origin="upper",vmin=0,vmax=4.5,cmap="viridis",aspect="auto")
        ax.set_xticks(np.arange(8),[f"c{i+1:02}" for i in range(8)])
        ax.set_yticks(np.arange(len(groups)),[f"{g['q']}×{g['q']}" for g in groups])
        ax.set_title(f"Calibration error {err:.1%}"+(" · fitted" if err==0 else " · replayed"))
        for i,j in np.ndindex(val.shape):
            value=val[i,j]
            if np.isnan(value):
                label="—"
            else:
                label="1 point" if value==0 else f"{np.floor(value*100)/100:.2f}"
                label+="*" if extra_used[i,j] else ""
            ax.text(j,i,label,ha="center",va="center",fontsize=9,
                    color="#343434" if np.isnan(value) else ("#15263c" if value>2.4 else "white"))
    figure.suptitle("ROI 2 peak-interval width proved by physical witnesses\n* Selective extended search   ·   — No verified witness",fontsize=15)
    figure.colorbar(cbar,ax=list(axes.flat),label="Certified width floor",shrink=.85)
    figure.savefig(HERE/"figures"/"roi_mimo_interval_phase.png",dpi=180,bbox_inches="tight")
    plt.close(figure)


if __name__=="__main__":
    main()
