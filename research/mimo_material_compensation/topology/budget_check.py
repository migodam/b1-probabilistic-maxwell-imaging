"""Outward arithmetic audit and same-scale image for the cap-matched witnesses."""
from pathlib import Path
import json

from ports import dump,load
from verify import fro, rd, up
import numpy as np
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
P=HERE/"loop"/"fresh_bg_b"


def run():
    p=load(P/"design_candidates.npz",("roi",))
    roi=p["roi"]
    rows={}
    pictures={}
    for name in ("B_interleaved","compensation"):
        results=json.loads((P/f"budget_replay_{name}.json").read_text())
        data=load(P/f"budget_replay_{name}.npz")
        centers=data["center"]
        before=results["center"]
        assert np.all(centers[roi]==0) and before["old"]["strict_feasible"] and before["new"]["strict_feasible"]
        assert centers.min()>=0 and centers.max()<=4.5
        area=float(roi.sum()*(2/32)**2)
        families={}
        pictures[name]={"ROI=0":centers}
        for family,checks in results["results"].items():
            trace=[]
            for row in checks:
                t=row["threshold"]
                x=data[f"{family}_{t}"]
                assert x.min()>=0 and x.max()<=4.5
                uncertainty=up(rd(2)*(fro(centers[~roi])+fro(x[~roi])))
                exterior=up((fro(x[~roi]-centers[~roi])+uncertainty)*(2/32))
                n=roi.sum()
                tmean=float(x[roi].mean())
                le=float(np.nextafter(tmean-rd(2*n)*tmean,-np.inf))
                passed=row["satisfies_data_balls"] and exterior<results["exterior_physical_area_L2_cap"]
                if t==.500001:
                    assert passed and le>.5
                    if family=="fixed_constant_ROI":
                        assert np.all(x[roi]==t)
                    pictures[name][("Constant ROI" if family=="fixed_constant_ROI" else "Shape-free mean")]=x
                trace.append(dict(mean_target=t,ROI_mean_lower=le,
                weighted_exterior_L2_upper=exterior,old_residual_upper_over_r=row["verified"]["old"]["upper_over_r"],
                new_residual_upper_over_r=row["verified"]["new"]["upper_over_r"],strict_nominal_witness=passed))
            families[family]=trace
        rows[name]=dict(ROI_cells=int(roi.sum()),ROI_area_m2=area,families=families,reference_upper_over_each_r=[
                         before["old"]["upper_over_r"],before["new"]["upper_over_r"]],
                         budget=results["exterior_physical_area_L2_cap"],CPU_s=results["cpu_s"])
    dump(P/"budget_witness_verified.json",dict(case="fresh_bg_b",layouts=rows,
          interpretation="Each layout's data-derived ROI-zero center plus legal alternatives, under its own A and acquired readout on fixed nominal 32² at 400 MHz. Each verified point is an upper bound on required exterior change. The finite mean samples are not a certified continuous path or a minimum cost."))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,5,figsize=(11,4),layout="constrained")
    for i,(arm,imgs) in enumerate(pictures.items()):
        bx=imgs["ROI=0"]
        for j,(label,x) in enumerate(imgs.items()):
            ax=axs[i,j]
            pic=ax.imshow(x.reshape(32,32),cmap="viridis",vmin=0,vmax=4.5,origin="lower")
            ax.contour(roi.reshape(32,32),levels=[.5],colors="white",linewidths=.5)
            if i==0:
                ax.set_title(label,fontsize=11)
        for j,key in enumerate(("Constant ROI","Shape-free mean")):
            ax=axs[i,j+3]
            d=imgs[key]-bx
            diff=ax.imshow(d.reshape(32,32),cmap="RdBu_r",vmin=-.75,vmax=.75,origin="lower")
            ax.contour(roi.reshape(32,32),levels=[.5],colors="black",linewidths=.5)
            if i==0:
                ax.set_title(("Change · constant" if j==0 else "Change · shape-free"),fontsize=11)
        axs[i,0].set_ylabel("Interleaved B" if arm=="B_interleaved" else "Compensation cut")
        for ax in axs[i]:
            ax.set_xticks([])
            ax.set_yticks([])
    fig.colorbar(pic,ax=axs[:,:3],shrink=.7,label="Material")
    fig.colorbar(diff,ax=axs[:,3:],shrink=.7,label="Material change")
    fig.savefig(HERE/"figures"/"fresh_bg_b_budget_witnesses.png",dpi=190)


if __name__=="__main__":
    with threadpool_limits(limits=1,user_api="blas"):
        run()
