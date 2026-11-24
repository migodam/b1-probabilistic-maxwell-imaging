"""Aggregate the native array trial and render paper figures."""

import argparse
from collections import defaultdict
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from screen import CAPS, HERE, RECEIVERS, SOURCES, load_json, save_json


COLORS = {20:"#086a87", 32:"#b94b23"}
OBJECTS = tuple(f"optfresh_{i:02d}" for i in range(1,7))


def native_outputs():
    primary={}
    high={}
    for p in RECEIVERS:
        fp=HERE/f"native_p{p}.json"
        primary[p]=load_json(fp) if fp.exists() else {"profiles":[]}
        refinement=HERE/f"gap_refinement_p{p}.json"
        if refinement.exists():
            for extra in load_json(refinement)["profiles"]:
                match=next((r for r in primary[p]["profiles"] if r["name"]==extra["name"] and
                            r["cap"]==extra["cap"]),None)
                if match and "exact_nu" in extra["state"]:
                    match["state"].update(extra["state"])
                    replayed={r["id"] for r in extra["observations"]}
                    match["observations"]=[r for r in match["observations"] if r["id"] not in replayed]
                    match["observations"]+=extra["observations"]
        high[p]=[]
        for tag in ("","_recovery"):
            fn=HERE/f"high_grid_p{p}{tag}.json"
            if fn.exists():
                high[p].append(load_json(fn))
    return primary,high


def write_csv(path,rows):
    with path.open("w",newline="") as out:
        wr=csv.DictWriter(out,fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)


def cases(screen,primary,high):
    verified=defaultdict(list)
    out_spectra=[]
    for p in RECEIVERS:
        for is_high,pack in [(False,primary[p]),*((True,row) for row in high[p])]:
            for profile in pack["profiles"]:
                for r in profile["observations"]:
                    u = profile["cap"]
                    verified[r["id"],u].append((profile["name"],r))
                    if "approximate_source_eigenvalues" not in r:
                        continue
                    for j,(nominal,lo,hi) in enumerate(zip(
                            r["approximate_source_eigenvalues"],r["source_eigenvalues_lower"],
                            r["source_eigenvalues_upper"],strict=True),start=1):
                        nu=r["negative_state_count"]
                        out_spectra.append(dict(id=r["id"],sources=r["sources"],receivers=p,
                            cap=u,noise_pct=r["noise_pct"],profile=profile["name"],high_grid=is_high,
                            radius=r["radius"],nu=nu,index=j,kept=j>nu,eigenvalue=nominal,
                            eigenvalue_lower=lo,eigenvalue_upper=hi))
    feasibility={r["id"]:r for r in load_json(HERE/"forward_audit.json")["rows"]}
    rows=[]
    for q in screen["records"]:
        r,u=q["padded_radius"],q["cap"]
        objs=verified[q["id"],u]
        distances=[(name,val) for name,val in objs if "distance_lower" in val]
        strict=sorted(((n,a) for n,a in distances if a["strictly_excludes"]),
                      key=lambda x:x[1]["distance_lower"],reverse=True)
        best=max(distances,key=lambda x:x[1]["distance_lower"],default=(None,{}))
        fixed_max=max(d["distance"] for d in q["fixed_profiles"])
        grid_max=max((d["distance"] for d in q["scalar_grid_best"]),default=0.)
        w=feasibility[q["id"]]
        f=w["generating_material_peak"]
        rows.append(dict(id=q["id"],object_id=q["object_id"],sources=q["sources"],
            receivers=q["receivers"],noise_pct=q["noise_pct"],tested_cap=u,radius=r,
            material_cap_status="excluded" if strict else "feasible_witness" if f<=u else "unknown",
            feasible_peak_upper=f,feasible_residual_upper=w["feasible_residual_upper"],
            strictly_excludes=bool(strict),strict_best_profile=best[0],
            strict_profile_distance_lower=best[1].get("distance_lower"),
            strict_profile_distance_upper=best[1].get("distance_upper"),
            strict_profile_first_retained=best[1].get("first_retained_source_eigenvalue"),
            strict_profile_positive_tail_count=best[1].get("positive_retained_eigenvalue_count_lower"),
            strict_native_profiles=len(distances),native_negative_count=best[1].get("negative_state_count"),
            min_projected_scalar_inertia=q["min_projected_scalar_inertia"],
            projected_fixed_max=fixed_max,projected_grid_max=grid_max,
            projected_grid_ratio=grid_max/r,numerically_positive_projection=max(fixed_max,grid_max)>r))
    return rows,out_spectra


def strongest_rows(p,u,rows,high):
    s32=[r for r in rows if r["receivers"]==p and r["sources"]==32 and r["tested_cap"]==u]
    num=max(s32,key=lambda r:r["projected_grid_ratio"])
    chosen=[]
    for pack in high[p]:
        for profile in pack["profiles"]:
            if profile["cap"]!=u:
                continue
            winner=next((r for r in profile["observations"] if r["id"]==profile["selection_id"]),None)
            if winner and "distance_upper" in winner:
                chosen.append((profile["name"],winner))
    strongest=max(chosen,key=lambda a:a[1]["distance_lower"],default=(None,None))
    name,strict=strongest
    return dict(projected_max_ratio=num["projected_grid_ratio"],
                native_profile=name, native_best_projector_case=strict["id"] if strict else None,
                native_best_projector_bounds=[strict["distance_lower"]/strict["radius"],
                                              strict["distance_upper"]/strict["radius"]] if strict else None,
                first_retained_eigenvalue=strict.get("first_retained_source_eigenvalue") if strict else None,
                retained_positives=strict.get("positive_retained_eigenvalue_count_lower") if strict else None,
                native_negative_count=strict.get("negative_state_count") if strict else None)


def counts_by_design(rows):
    counts={}
    for p in RECEIVERS:
        counts[p]={}
        for n in (1,2):
            counts[p][n]={}
            for s in SOURCES:
                obs=[r for r in rows if r["receivers"]==p and r["noise_pct"]==n and r["sources"]==s]
                u75=[r for r in obs if r["tested_cap"]==.75]
                wins=[r["object_id"][-2:] for r in u75 if r["material_cap_status"]=="excluded"]
                counts[p][n][s]=dict(excluded_at_075=len(wins),object_ids=wins,
                    highest_excluded_by_object={ob[-2:]:max((r["tested_cap"] for r in obs
                                              if r["object_id"]==ob and r["strictly_excludes"]),default=None)
                                              for ob in OBJECTS})
    return counts


def line_panel(ax,x,y,label,**style):
    return ax.plot(x,y,label=label,solid_capstyle="round",**style)


def format_panel(ax,letter):
    ax.text(-.12,1.045,letter,transform=ax.transAxes,ha="left",va="top",
            color="#242a30",fontsize=10,fontweight="bold")
    ax.spines[["right","top"]].set_visible(False)
    ax.tick_params(axis="both",direction="out",length=3,pad=3)
    ax.xaxis.labelpad=4
    ax.yaxis.labelpad=5
    ax.grid(axis="y",color="#dce2e7",linewidth=.55)
    ax.set_axisbelow(True)


def save_fig(fig,name):
    fig.savefig(HERE/(name+".pdf"),bbox_inches="tight",pad_inches=.06)
    fig.savefig(HERE/(name+".png"),dpi=320,bbox_inches="tight",pad_inches=.06)
    plt.close(fig)


def figure(screen,bare,primary,summary):
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8.0,
                        "axes.labelsize":8.2,"legend.fontsize":6.6,"pdf.fonttype":42,
                        "axes.edgecolor":"#565e68","axes.linewidth":.65})
    fig,axs=plt.subplots(2,2,figsize=(7.35,5.22),layout="constrained")
    ax=axs[0,0]
    colors=["#a5acb4","#78878e","#0c7686","#174c77"]
    for (s,row),color in zip(screen["sources"].items(),colors,strict=True):
        vals=row["relative_singular_values"]
        line_panel(ax,np.arange(1,len(vals)+1),vals,f"$B$, $s={s}$",color=color,lw=1.35)
    for p,row in screen["receivers"].items():
        vals=row["relative_singular_values"]
        line_panel(ax,np.arange(1,len(vals)+1),vals,f"$S$, $p={p}$",color=COLORS[int(p)],
                   lw=1.1,ls="--")
    ax.set(yscale="log",ylim=(.0015,1.4),xlim=(1,33),xticks=[1,8,16,24,32],
           xlabel="Singular value index",ylabel="Relative singular amplitude")
    ax.axvline(2*screen["full_wave_number"]*screen["enclosing_radius_m"],ls=":",lw=.9,color="#545d67")
    ax.text(23.1,.0022,r"$2ka$",ha="right",va="bottom",color="#49535c")
    ax.legend(loc="lower left",ncol=2,frameon=False,handlelength=1.7,columnspacing=.85)
    format_panel(ax,"a")

    ax=axs[0,1]
    line_panel(ax,CAPS,[r["native_nu_at_least"] for r in bare["rows"]],"No receivers, real sector",
               color="#515862",lw=1.05,ls="--",marker="^",ms=4.2,mfc="white")
    for x,row in zip(CAPS,bare["rows"]):
        if "exact_nu" not in row:
            ax.annotate("",(x,row["native_nu_at_least"]+3),(x,row["native_nu_at_least"]+.3),
                        arrowprops=dict(arrowstyle="-|>",color="#515862",lw=.63))
    for p,blob in primary.items():
        d={r["cap"]:r["state"]["native_nu_at_least"] for r in blob["profiles"]
           if r["name"]=="stable_scalar_optical"}
        if d:
            line_panel(ax,list(d),list(d.values()),f"$p={p}$, optical weight",
                       color=COLORS[p],lw=1.3,marker="o",ms=4)
            for r in blob["profiles"]:
                if r["name"]=="stable_scalar_optical" and "exact_nu" not in r["state"]:
                    x=r["cap"]
                    y=r["state"]["native_nu_at_least"]
                    ax.annotate("",(x,y+2.25),(x,y+.35),
                                arrowprops=dict(arrowstyle="-|>",color=COLORS[p],lw=.6))
    ax.axhline(2*screen["full_wave_number"]*screen["enclosing_radius_m"],ls=":",lw=.85,color="#68737f")
    ax.text(.73,24.45,r"$2ka$",color="#68737f",va="bottom")
    ax.set(xlim=(.70,2.07),xticks=CAPS,xticklabels=["0.75","1","1.25","1.5","2"],
           ylim=(0,50),xlabel="Tested contrast cap $u$",ylabel="Negative state directions")
    ax.legend(frameon=False,loc="upper left",handlelength=1.5)
    format_panel(ax,"b")

    ax=axs[1,0]
    for p in RECEIVERS:
        for n,ls,mark,fill,offset in ((1,"-","o",COLORS[p],-.24),(2,"--","s","white",.24)):
            xx=np.asarray(SOURCES,dtype=float)+offset+(.13 if p==32 else -.13)
            vals=[summary["counts"][p][n][s]["excluded_at_075"] for s in SOURCES]
            line_panel(ax,xx,vals,f"$p={p}$, {n}% noise",color=COLORS[p],
                       ls=ls,lw=1.2,marker=mark,ms=4.5,mfc=fill,zorder=3)
    ax.set(xlim=(8.8,33.2),xticks=SOURCES,ylim=(0,6.35),yticks=range(7),
           xlabel="Plane-wave illuminations $s$",ylabel="Materials shown to exceed 0.75")
    ax.legend(loc="lower right",frameon=False,ncol=2,handlelength=1.7,columnspacing=.9)
    format_panel(ax,"c")

    ax=axs[1,1]
    for p,offset in ((20,-.013),(32,.013)):
        x=np.array(CAPS[1:])+offset
        q=summary["largest_high_cap_s32"][p]
        ys=[q[u]["projected_max_ratio"] for u in CAPS[1:]]
        line_panel(ax,x,ys,f"$p={p}$, projected scan",color=COLORS[p],lw=1.1,
                   ls="--",marker="o",mfc="white",ms=4)
        for xp,u in zip(x,CAPS[1:]):
            strict=q[u]["native_best_projector_bounds"]
            if strict is not None:
                ax.plot(xp,float(np.mean(strict)),color=COLORS[p],marker="o",ls="None",ms=3.35)
    ax.axhline(1.,ls=":",color="#354351",lw=.95)
    ax.text(1.98,.98,"Exclusion threshold",ha="right",va="top",color="#354351")
    ax.set(xticks=CAPS[1:],xticklabels=["1","1.25","1.5","2"],xlim=(.94,2.075),
           ylim=(-.035,1.08),xlabel="Tested contrast cap $u$",ylabel=r"Largest bound / noise radius")
    handles=[Line2D([0],[0],color=COLORS[p],ls="--",marker="o",mfc="white",ms=4)
             for p in RECEIVERS]
    handles.append(Line2D([0],[0],color="#364553",ls="None",marker="o",ms=3.4))
    ax.legend(handles,["20 receivers","32 receivers","Full-grid replay"],loc="upper left",
              bbox_to_anchor=(0,.89),frameon=False,handlelength=1.6)
    format_panel(ax,"d")
    save_fig(fig,"array_scaling")


def source_plot(primary,high):
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8.5,"pdf.fonttype":42})
    select=[("stable_scalar_optical",.75,10),("stable_scalar_optical",.75,16),
            ("stable_scalar_optical",.75,24),("stable_scalar_optical",.75,32)]
    one=primary[32]["profiles"]
    controls=[p for pack in high[32] for p in pack["profiles"] if p["cap"]==1.
              and p.get("state",{}).get("exact_nu") is not None]
    fig,axs=plt.subplots(1,2,figsize=(7.35,2.85),layout="constrained")
    for j,(ax,obj,controls0) in enumerate(((axs[0],"02",select),
                                           (axs[1],"04",[select[-1]]+[
                                               (p["name"],1.,32) for p in controls[:1]]))):
        for k,(label,u,s) in enumerate(controls0):
            prof=next((p for p in one+controls if p["cap"]==u and p["name"]==label),None)
            if not prof:
                continue
            case=next((r for r in prof["observations"] if r["object_id"].endswith(obj) and
                       r["sources"]==s and r["noise_pct"]==2),None)
            if not case or "approximate_source_eigenvalues" not in case:
                continue
            vals=np.asarray(case["approximate_source_eigenvalues"])/case["radius"]**2
            nu=case["negative_state_count"]
            x=np.arange(1,s+1)
            c=(["#98a4af","#728f96","#0f7584","#254a75"] if j==0 else ["#086a87","#b94b23"])[k]
            label0=f"$s={s}$" if j==0 else f"$u={u:g}$, $\\nu={nu}$, $d/r={case['distance_lower']/case['radius']:.2f}$"
            line_panel(ax,x,vals,label0,color=c,lw=1)
            keep=(x>nu)&(vals>0)
            ax.plot(x[keep],vals[keep],ls="None",marker="o",ms=3,mfc=c,mec=c,color=c)
            ax.plot(x[:nu],vals[:nu],ls="None",marker="o",ms=2.2,mfc="white",mec=c,color=c)
        ax.axhline(0.,color="#707983",lw=.65)
        ax.set(yscale="symlog",xlabel="Source eigenvalue index",ylabel=r"Source eigenvalue / $r^2$",
               xticks=[1,8,16,24,32],title=f"Material {obj}, 2% noise")
        ax.legend(loc="lower left",frameon=False,ncol=2 if j==0 else 1)
        format_panel(ax,"ab"[j])
    save_fig(fig,"source_spectra")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partial",action="store_true")
    args=parser.parse_args()
    screen=load_json(HERE/"screen.json")
    bare=load_json(HERE/"bare_inertia.json")
    refinement=HERE/"gap_refinement_p0.json"
    if refinement.exists():
        base={r["cap"]:r for r in bare["rows"]}
        for extra in load_json(refinement)["profiles"]:
            if "exact_nu" in extra["state"]:
                base[extra["cap"]].update(extra["state"])
    primary,high=native_outputs()
    done=all(set(r["cap"] for r in res["profiles"])==set(CAPS) for res in primary.values())
    if not args.partial and not done:
        raise RuntimeError("Primary native results still streaming.")
    rows,spectra=cases(screen,primary,high)
    write_csv(HERE/"array_scan.csv",rows)
    write_csv(HERE/"source_spectra.csv",spectra)
    compact=dict(primary_complete=done, counts=counts_by_design(rows),
        largest_high_cap_s32={p:{u:strongest_rows(p,u,rows,high) for u in CAPS[1:]} for p in RECEIVERS},
        case_rows=len(rows),number_of_objects=len(OBJECTS),
        native_state_counts={p:{u:dict(negative_count_lower=state["native_nu_at_least"],
                                       exact="exact_nu" in state) for u,state in q}
                             for p,q in [(0,[(r["cap"],r) for r in bare["rows"]]),
                                *((p,[(r["cap"],r["state"]) for r in primary[p]["profiles"]
                                      if r["name"]=="stable_scalar_optical"]) for p in RECEIVERS)]},
        unknown_by_cap={u:sum(r["tested_cap"]==u and r["material_cap_status"]=="unknown"
                             for r in rows) for u in CAPS},
        projected_crossings_without_strict_replay=[dict(id=r["id"],cap=r["tested_cap"])
          for r in rows if r["numerically_positive_projection"] and not r["strictly_excludes"]])
    save_json(HERE/"summary.json",compact)
    figure(screen,bare,primary,compact)
    source_plot(primary,high)


if __name__=="__main__":
    main()
