"""Plots from frozen topology, acquisition and post-acquisition scoring receipts."""
from pathlib import Path
import argparse
import csv
import json

from ports import A, B, PortEngine, cut_id, dump, load, project
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
PLOTS=HERE/"figures"
PLOTS.mkdir(exist_ok=True)
plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,
                     "axes.labelsize":9,"axes.titlesize":10,"figure.dpi":150})


def fixed():
    cards=json.loads((HERE/"fixed_pair"/"summary.json").read_text())
    powers=load(HERE/"fixed_pair"/"cut_scores.npz")
    d=powers["mixed_exact"]
    held=cards["heldout_medians"]
    labels=[r["case"].split("_")[-1] for r in cards["cases"]]
    x=np.arange(len(labels))
    fig,axs=plt.subplots(1,3,figsize=(12.6,3.9),layout="constrained")
    a,b=cut_id(A),cut_id(B)
    axs[0].boxplot(np.sqrt(d.T),positions=x,widths=.68,whis=(5,95),showfliers=False,
                    patch_artist=True,boxprops=dict(facecolor="#d7e5ef",edgecolor="#6c879b"),
                    medianprops=dict(color="#52758e",linewidth=1.4))
    axs[0].scatter(x,np.sqrt(d[:,a]),c="#303438",marker="o",label="A",s=13)
    axs[0].scatter(x,np.sqrt(d[:,b]),c="#bf562e",marker="D",label="B",s=17)
    axs[0].scatter(x,np.sqrt(d.max(axis=1)),c="#356b49",marker="^",label="Pair-specific best",s=21)
    axs[0].axhline(1,color="#555555",linestyle=":",lw=.8)
    axs[0].set(xticks=x,xticklabels=labels,xlabel="Independent historical material",
               ylabel=r"Fixed pair separation $d/(2r)$",title="6,435 balanced layouts per pair")
    axs[0].legend(frameon=False,loc="upper left",fontsize=7)
    parts=np.array([[r["edge_partition"][k]/r["energy_B_raw"] for k in
                           ("shared_A_B","new_B_same_A_tx_arc","new_B_same_A_rx_arc")] for r in cards["cases"]])
    axs[1].bar(x,parts[:,0],color="#959da2",width=.67,label="Shared with A")
    axs[1].bar(x,parts[:,1],bottom=parts[:,0],color="#587d98",width=.67,label="New in A transmit arc")
    axs[1].bar(x,parts[:,2],bottom=parts[:,:2].sum(axis=1),color="#bd643f",width=.67,label="New in A receive arc")
    axs[1].set(xticks=x,xticklabels=labels,xlabel="Independent historical material",
               ylim=(0,1.03),ylabel="Fraction of B separation energy",title="What B measures that A did not")
    axs[1].legend(frameon=True,fontsize=7,loc="lower right",framealpha=.95)
    selected={"Exact train median":"mixed_exact_median_rank","Born train maximin":"born_maximin",
              "Centered J maximin":"centered_J_maximin","Nuisance Fisher":"nuisance_fisher_median_rank",
              "Interleaved B":"B","Geometry":"geometry_median_rank","Random":"random","A":"A"}
    y=np.arange(len(selected))
    axs[2].barh(y,[held[k]["heldout_fraction_of_pair_ceiling"] for k in selected.values()],color="#587d98",height=.65)
    axs[2].set(yticks=y,yticklabels=list(selected),xlim=(0,1.02),
               xlabel="Median held-out fraction of pair ceiling",title="Leave one material out")
    axs[2].invert_yaxis()
    fig.savefig(PLOTS/"fixed_witness_competition.png",dpi=180)
    plt.close(fig)

    statistic=np.array([[r["metrics_relative_to_mixed"][key]["all_cut_rank_rho"] for key in
                          ("born","centered_J","hybrid")] for r in cards["cases"]])
    fig,axs=plt.subplots(1,2,figsize=(8.9,3.5),layout="constrained")
    for i,(k,color) in enumerate(zip(("Born","Centered Jacobian","One dressed field"),("#bd643f","#587d98","#629775"))):
        axs[0].bar(x+.22*(i-1),statistic[:,i],width=.2,label=k,color=color)
    axs[0].set(xlabel="Independent historical material",xticks=x,xticklabels=labels,
               ylim=(-.38,1.08),ylabel="Spearman correlation with exact rank",title="Linear proxies for finite-amplitude endpoint pairs")
    axs[0].legend(frameon=False,fontsize=7)
    raw=np.array([r["energy_C_raw"]/r["energy_B_raw"] for r in cards["cases"]])
    normalized=np.array([r["fixed_ratios"]["C"]/r["fixed_ratios"]["B"] for r in cards["cases"]])
    axs[1].plot(x,np.sqrt(raw),marker="o",color="#346b54",label=r"Raw $d_C/d_B$")
    axs[1].plot(x,normalized,marker="^",color="#bd643f",label=r"Normalized $(d_C/2r_C)/(d_B/2r_B)$")
    axs[1].axhline(1,color="#555555",ls=":",lw=.8)
    axs[1].set(xlabel="Independent historical material",xticks=x,xticklabels=labels,
                title="Frobenius contract reverses the raw order")
    axs[1].legend(frameon=False,fontsize=7,loc="best")
    fig.savefig(PLOTS/"surrogates_and_power.png",dpi=180)
    plt.close(fig)


def active():
    rows=[]
    truth_fields=[]
    engine=PortEngine()
    for p in sorted((HERE/"loop").glob("*/attack_results.json")):
        best=json.loads(p.read_text())
        if "stop" not in best:
            continue
        case=best["case"]
        folder=p.parent
        frozen=json.loads((folder/"design_frozen.json").read_text())
        strict_doc=folder/"strict_attacks.json"
        strict=json.loads(strict_doc.read_text())["arms"] if strict_doc.exists() else {}
        m=load(folder/"attacked_materials.npz")
        sealed=load(folder/"sealed_scoring.npz",("material","full_calibrated_delta"))
        ro=load(folder/"design_candidates.npz",("roi","patch"))
        xtrue=sealed["material"]
        st_truth=engine.state(xtrue)
        truth_dressed=st_truth["delta"]
        observed_old=load(folder/"initial_A.npz")
        observed_new=load(folder/"new_observed.npz")
        truth_old=float(np.linalg.norm(project(truth_dressed,A)-observed_old["Y"])/float(observed_old["radius"]))
        truth_new={name:float(np.linalg.norm(project(truth_dressed,ds["tx"],ds["rx"])-observed_new[name+"_Y"])/float(observed_new[name+"_r"]))
                   for name,ds in frozen["actions"].items()}
        assert max([truth_old,*truth_new.values()])<1
        raw_hold=project(sealed["full_calibrated_delta"],A,np.arange(16,24))
        easy={}
        for name,inrow in best["results"].items():
            arm=inrow.get("alias_of",name)
            row=best["results"][arm]
            cfg=frozen["actions"][name]
            hmin,hmax=row["numeric_inner_ROI_range"] or (None,None)
            short=dict(case=case,method=name,source_arm=arm,cut=cfg["cut"],
                  initial_feasible_pair=not frozen["candidate_pair_fallback"],
                  predicted_D_over_2r=cfg["predicted_worst_pair_D_over_2r"],
                  numerical_roi_width=(hmax-hmin if hmin is not None else 0.),
                  verified_ROI_width=strict.get(arm,{}).get("ROI_width_lower"),
                  numerical_compatible_count=sum(r["best"]["numerical_in_balls"] for r in row["families"].values()),
                  wall_s=row["wall_s"],cpu_s=row["cpu_s"],
                  state_solves=row["counts"].get("state_solves",0),dense_LU=row["counts"].get("dense_LU",0))
            score={}
            for f,rec in row["families"].items():
                key=arm+"_"+f
                mat=m[key]
                residual=max(rec["best"]["residual_over_each_radius"])
                yhold=project(m[key+"_delta"],A,np.arange(16,24))
                acc=dict(max_data_resid=residual,ROI_mean=rec["best"]["ROI_mean"],
                         exterior_L2_change=rec["best"].get("exterior_L2_area"),
                         material_relative_L2=float(np.linalg.norm(mat-xtrue)/np.linalg.norm(xtrue)),
                         target_RMSE=float(np.sqrt(np.mean((mat[xtrue>=1]-xtrue[xtrue>=1])**2))),
                         background_RMSE=float(np.sqrt(np.mean((mat[xtrue<1]-xtrue[xtrue<1])**2))),
                         offline_unread_ports_relative_error=float(np.linalg.norm(yhold-raw_hold)/np.linalg.norm(raw_hold)))
                score[f]=acc
                for key1,value in acc.items():
                    short[f+"_"+key1]=value
            rows.append(short)
            easy[name]=score
        dump(folder/"postseal_score.json",dict(case=case,online_selection="frozen before clean score",
                  ROI_truth_mean=float(xtrue[ro["roi"]].mean()),
                  sealed_truth_residual_old_over_r=truth_old, sealed_truth_residual_new_over_r=truth_new,
                  sealed_truth_original_state_defect=st_truth["state_defect"],
                  actions=easy,
                  scope="Unobserved terminal responses are simulator offline scoring channels only. Material truth was not read during design or refit.",
                  acquisition_id=frozen["case"]))
        truth_fields.append((case,m,ro,xtrue))
    if not rows:
        return
    with (HERE/"active_summary.csv").open("w") as fp:
        w=csv.DictWriter(fp,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    cases=list(dict.fromkeys(r["case"] for r in rows))
    groups=["A_repeat","A_reverse","B_interleaved","random","Born_geometry","Fisher",
            "nuisance_Fisher","signal","T_discrimination","compensation"]
    families=["zero","patch","mean_low","mean_high","unrestricted"]
    label=["A repeat","A reversal","B interleaved","Random","Born geometry","Fisher","Nuisance Fisher",
           "Signal","T-discrimination","Compensation"]
    fig,axs=plt.subplots(len(cases),2,figsize=(10,2.5*len(cases)+1.3),squeeze=False,layout="constrained")
    for i,c in enumerate(cases):
        sel={r["method"]:r for r in rows if r["case"]==c}
        mat=np.array([[sel[g][f+"_max_data_resid"] for f in families] for g in groups])
        axs[i,0].imshow(np.log10(mat),aspect="auto",cmap="RdYlGn_r",vmin=-.25,vmax=1.15)
        for x in range(len(groups)):
            for y in range(len(families)):
                axs[i,0].text(y,x,f"{mat[x,y]:.2f}",ha="center",va="center",fontsize=6)
        axs[i,0].set(yticks=np.arange(len(groups)),yticklabels=label,xticks=np.arange(len(families)),
                    xticklabels=["Zero","Patch","Low mean","High mean","Free"],
                    title=f"{c}  ·  first-budget maximum residual")
        vals=[sel[g]["numerical_roi_width"] for g in groups]
        axs[i,1].barh(np.arange(len(groups)),vals,color="#63869f")
        vals2=[sel[g]["verified_ROI_width"] or 0. for g in groups]
        axs[i,1].barh(np.arange(len(groups)),vals2,color="#3e6543",height=.4)
        if c=="fresh_bg_b":
            path=HERE/"loop"/c/"dual_cie_verified.json"
            if path.exists():
                extra=json.loads(path.read_text())["arms"]
                indices=[groups.index(name) for name in extra]
                vals3=[q["verified_inner_width"] for q in extra.values()]
                axs[i,1].barh(indices,vals3,color="#a86531",height=.45)
                for p,w in zip(indices,vals3):
                    axs[i,1].text(w+.13,p,"Longer refit",fontsize=7,va="center")
        axs[i,1].axvline(.5,color="#222222",ls=":",lw=.8)
        axs[i,1].set(yticks=np.arange(len(groups)),yticklabels=label,xlim=(0,4.5),
                    title=f"{c}  ·  realizable ROI mean span")
        axs[i,1].invert_yaxis()
    fig.savefig(PLOTS/"closed_loop_residuals_and_spans.png",dpi=180)
    plt.close(fig)
    for c,m,ro,x in truth_fields:
        folder=HERE/"loop"/c
        att=json.loads((folder/"attack_results.json").read_text())["results"]
        names=["A_repeat","B_interleaved","nuisance_Fisher","T_discrimination","compensation"]
        fams=["zero","patch","mean_high"]
        fig,axs=plt.subplots(len(fams),len(names)+1,figsize=(11,5.8),layout="constrained")
        side=int(np.sqrt(len(x)))
        roi=ro["roi"].reshape(side,side)
        for i,f in enumerate(fams):
            for j,n in enumerate([None,*names]):
                ax=axs[i,j]
                k=None if n is None else att[n].get("alias_of",n)
                v=(x if k is None else m[k+"_"+f]).reshape(side,side)
                pic=ax.imshow(v,vmin=0,vmax=4.5,cmap="viridis",origin="lower")
                ax.contour(roi,levels=[.5],colors="white",linewidths=.45)
                if i==0:
                    ax.set_title("Truth (offline)" if n is None else n.replace("_"," "))
                if n is None:
                    ax.set_ylabel(f.replace("_"," "))
                else:
                    res=att[k]["families"][f]["best"]["residual_over_each_radius"]
                    ax.set_xlabel(f"old {res[0]:.2f}  ·  new {res[1]:.2f}",fontsize=7)
                ax.set_xticks([])
                ax.set_yticks([])
        fig.colorbar(pic,ax=axs,shrink=.77,label="Material contrast")
        fig.savefig(PLOTS/f"materials_{c}.png",dpi=165)
        plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("part",choices=("all","fixed","active"),default="all",nargs="?")
    part=p.parse_args().part
    if part in ("all","fixed"):
        fixed()
    if part in ("all","active"):
        active()
