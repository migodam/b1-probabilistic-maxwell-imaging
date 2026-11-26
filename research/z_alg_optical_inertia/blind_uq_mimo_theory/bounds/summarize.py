"""Assemble the blind certified-bound receipts and reproduce the research figures."""
from collections import defaultdict
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from sector import HERE, UPPER

RESULTS=HERE/"results"
FIGURES=HERE/"figures"


def read(name):
    p=RESULTS/name
    return json.loads(p.read_text()) if p.exists() else None


def collect_cases(version):
    rows=[]
    follow=read("stretch_global.json") or {}
    strengthened={d["case"]:max((q["upper"] for q in d["proofs"] if q["proved"]),default=None)
                  for d in follow.get("rows",[])}
    for p in sorted(RESULTS.glob(f"n*_{version}.json")):
        d=json.loads(p.read_text())
        case=d["case"]
        source=d.get("sections",{})
        m=source.get("roi",{})
        if "roi_01" not in m:
            continue
        ver=read(f"{case}_{version}_verified.json") or {}
        audits=ver.get("section",{})
        low=audits.get("global_cap",{})
        joined={}
        for name in ("refined_cap","refined_cap_v2"):
            vv=read(f"{case}_{version}_{name}.json") or {}
            for k,a in vv.get("regions",{}).items():
                pr=[v["verification"] for v in a["strict_tests"] if v["verification"]["proved"]]
                if pr and pr[0]["roi_peak_lower"]>joined.get(k,{}).get("roi_peak_lower",-1):
                    joined[k]=pr[0]
        reported={}
        for k,r in m.items():
            ob=r["obstruction"]["profile_with_smallest_dark_count"]
            audit=audits.get(k+"_zero",{})
            cap=audit.get("strict_positive_cap",{}).get("roi_peak_lower")
            if k in joined:
                cap=max(cap or 0,joined[k]["roi_peak_lower"])
            numeric0=(r.get("zero_roi_test",{}).get("best") or {})
            dark=audits.get(k,{}).get("continuous_exterior_ratio_family",{})
            reported[k]=dict(cells=r["roi_cells"],numeric_smallest_negative_count=ob["negative_lower_numeric"],
                             numerically_best_profile=ob["method"],numeric_zero_margin=numeric0.get("margin"),
                             strict_zero_proved=audit.get("proved",False),
                             strict_zero_distance=audit.get("distance_lower"),strict_zero_nu=audit.get("nu"),
                             spectral_family_strict=dark.get("s_negatives_verified_for_every_ratio"),
                             largest_strict_positive_peak_floor=cap,
                             tested_refined_positive_details=joined.get(k))
        rows.append(dict(case=case,grid=d["side"],mask_version=version,radius=d["radius"],
                         sources=d["sources"],receivers=d["receivers"],wall_seconds=d["wall_seconds"],
                         cpu_seconds=d["cpu_seconds"],profile_windows=reported,
                         ordinary_roi_peak_interval=[0.,UPPER],
                         verified_global_peak_lower=low.get("upper") if low.get("proved") else None,
                         verified_global_peak_lower_after_dense_cap_scan=strengthened.get(case),
                         verified_global_peak_distance=low.get("distance_lower") if low.get("proved") else None,
                         finest_region_strict_zero=min((r["cells"] for r in reported.values()
                                                       if r["strict_zero_proved"]),default=None),
                         finest_region_numeric_zero=min((r["cells"] for r in reported.values()
                                                       if (r["numeric_zero_margin"] or 0)>0),default=None)))
    return rows


def signed_cost():
    out=[]
    for p in sorted(RESULTS.glob("n32*_main.json")):
        d=json.loads(p.read_text())
        up=d["sections"].get("roi",{}).get("roi_01",{}).get("signed_updates")
        if not up:
            continue
        r=up["rounds"]
        out.append(dict(case=d["case"],rounds=len(r),ranks=sorted(set(a["proof_size"] for a in r)),
                        inertia_matches=sum(a["nu_lowrank"]==a["nu_direct"] for a in r),
                        max_distance_difference=max(abs(a["distance_lowrank"]-a["distance_direct"]) for a in r),
                        base_seconds=up["preparation_seconds"],updates_seconds=sum(a["update_seconds"] for a in r),
                        direct_seconds=sum(a["direct_seconds"] for a in r)))
    return out


def miprobe(name):
    data=read(name)
    if not data:
        return None
    groups=defaultdict(list)
    for a in data["rows"]:
        key=(a["frequency"],a["channels"],a["coupling"],a["calibration_error"])
        groups[key].append(a)
    out={"grid":data["side"],"observations":len(data["rows"]),
         "independent_materials":len({z["case"] for z in data["rows"]}),
         "freq_products":len(data["frequency_products"]),
         "freq_product_positive_numeric":sum(z["joint_margin"]>0 for z in data["frequency_products"]),
         "ordinary_roi_max_surviving_columns_numeric":max(
             u["primary_roi_zero_largest_surviving_tail_numeric"] for u in data["geometry"]),
         "wall_seconds":data["wall_seconds"],"rows":[]}
    for k,v in sorted(groups.items()):
        a,b,c,d=k
        out["rows"].append(dict(frequency=a,q=b,coupling=c,calibration=d,n=len(v),
                           passes_num_absolute=sum(x["margin_absolute_envelope"]>0 for x in v),
                           passes_num_noise_only=sum(x["margin_nominal_radius_only"]>0 for x in v),
                           strict_count=sum(x.get("arb_global",{}).get("proved",False) for x in v),
                           minimum_margin=min(x["margin_absolute_envelope"] for x in v)))
    return out


def phase_picture(data):
    f,axes=plt.subplots(2,3,figsize=(11,6.9),sharex=True,sharey=True,layout="constrained")
    frequencies=(400,460)
    modes=(16,48,144)
    for row,freq in enumerate(frequencies):
        for col,c in enumerate(modes):
            a=axes[row,col]
            selected=[k for k in data["rows"] if k["frequency"]==freq and k["coupling"]==c]
            q=(6,8,12,16)
            d=(0.,0.005,0.025,0.060)
            z=np.full((len(d),len(q)),np.nan)
            for x in selected:
                z[d.index(x["calibration"]),q.index(x["q"])]=x["passes_num_absolute"]/x["n"]
            a.imshow(z,origin="lower",vmin=0,vmax=1,cmap="viridis",aspect="auto")
            a.set_title(f"{freq} MHz · radiation scale {c} Ω",fontsize=11)
            a.set_xticks(np.arange(len(q)),[f"{i}×{i}" for i in q])
            a.set_yticks(np.arange(len(d)),[f"{v:.3f}" for v in d])
            if col==0:
                a.set_ylabel("Normalized calibration error")
            if row==1:
                a.set_xlabel("MIMO array")
            for x in selected:
                j=d.index(x["calibration"]); i=q.index(x["q"])
                fraction=z[j,i]
                a.text(i,j,f"{x['passes_num_absolute']}/{x['n']}",ha="center",va="center",
                       color="white" if fraction<.52 else "#102137",fontsize=9)
    f.suptitle("Whole-field exclusion of contrast ≤ 0.25 · matrix lower-bound scan",fontsize=15)
    f.savefig(FIGURES/"mimo_global_phase.png",dpi=180,bbox_inches="tight")
    plt.close(f)


def picture_roi():
    image=HERE.parent/"imaging"/"results"/"n32_c01"/"stationary"/"baseline.npz"
    profile=RESULTS/"n32_c01_stationary.npz"
    cases=collect_cases("stationary")
    if not image.exists() or not profile.exists() or not cases:
        return
    chi=np.load(image)["chi"]
    masks=np.load(profile)
    row=next(v for v in cases if v["case"]=="n32_c01")
    keys=("roi_01","all_except_background_04","all_except_background_02","expanded_97")
    f,axs=plt.subplots(1,4,figsize=(12.5,4),layout="constrained")
    title=("Primary blind ROI","Leave background 4 free","Leave background 2 free","Expanded blind region")
    for ax,k,t in zip(axs,keys,title,strict=True):
        val=row["profile_windows"][k]
        im=ax.imshow(chi,cmap="viridis",vmin=0,vmax=UPPER,origin="lower")
        ax.contour(masks[k+"_mask"],levels=[.5],colors=["#ffbd65"],linewidths=2.0,origin="lower")
        L=val["largest_strict_positive_peak_floor"]
        cap=f"Certified peak > {np.floor(L*1e5)/1e5:.5f}" if L else "Certified interval [0, 4.5]"
        ax.set_title(t+f"\n{val['cells']}/1024 cells · {cap}",fontsize=10)
        ax.set(xticks=[],yticks=[])
    f.colorbar(im,ax=axs,location="bottom",shrink=.45,pad=.05,label="Stationary CIE contrast")
    f.savefig(FIGURES/"roi_peak_example.png",dpi=180,bbox_inches="tight")
    plt.close(f)


def main():
    data=dict(main=collect_cases("main"),stationary=collect_cases("stationary"),
              signed_update_cost=signed_cost(),mimo32=miprobe("mimo_32_stationary.json") or miprobe("mimo_32.json"),
              mimo96=miprobe("mimo_96_stationary.json"))
    (RESULTS/"summary.json").write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
    if data["mimo32"]:
        phase_picture(data["mimo32"])
    picture_roi()


if __name__=="__main__":
    main()
