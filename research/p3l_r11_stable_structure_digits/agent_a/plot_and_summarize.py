"""Zero-physics postprocessing of recorded paired timings."""
import os
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent


def rows(name):
    out=[]
    for r in csv.DictReader((HERE/name).open()):
        out.append(dict(support=int(r["support_size"]),jac=r["need_jac"]=="True",full=float(r["Full256_seconds"]),support_time=float(r["support_seconds"])))
    return out


def main():
    datasets=[("Initial, four property nodes",rows("R10_KERNEL_NODES.csv")),("Final solver, one property node",rows("R10_OPTIMIZED_KERNEL_NODES.csv"))]
    fig,axes=plt.subplots(2,2,figsize=(10,7),sharex=True,sharey=True)
    groups=[]
    for ir,(label,rr) in enumerate(datasets):
        for jc,jac in enumerate((False,True)):
            ax=axes[ir,jc];rs=[r for r in rr if r["jac"]==jac]
            ax.axhline(1,color="#777",lw=1,ls="--")
            ax.scatter([r["support"] for r in rs],[r["full"]/r["support_time"] for r in rs],s=22,c="#227c9d" if not jac else "#c05a36",alpha=.8)
            ax.set_title(label+"\n"+("Full occupancy + property Jacobian" if jac else "Forward only"),fontsize=10)
            ax.set_yscale("log")
            ax.grid(alpha=.2)
            if jc==0:ax.set_ylabel("Full256 time / support time")
            if ir==1:ax.set_xlabel("Strictly positive occupancy cells")
            for lo,hi,name in [(0,127,"s<=127"),(128,223,"128<=s<=223"),(224,255,"224<=s<=255"),(256,256,"s=256")]:
                sub=[r for r in rs if lo<=r["support"]<=hi]
                if sub:
                    full=sum(r["full"] for r in sub);ss=sum(r["support_time"] for r in sub)
                    groups.append(dict(batch=label,need_jac=jac,group=name,nodes=len(sub),full_seconds=full,support_seconds=ss,saving=1-ss/full,faster_nodes=sum(r["support_time"]<r["full"] for r in sub)))
    fig.suptitle("Exact support reduces state solves; dense shapes retain overhead",fontsize=13)
    fig.tight_layout(rect=(0,0,1,.95));fig.savefig(HERE/"R10_SUPPORT_COST.png",dpi=180);plt.close(fig)
    (HERE/"SUPPORT_COST_GROUPS.json").write_text(json.dumps(groups,indent=2)+"\n")
    if (HERE/"DIGIT_KERNEL_SUMMARY.json").exists() and (HERE/"N_PAIRED_INVERSE_NO_GT_AUDIT.json").exists():
        dig=json.loads((HERE/"DIGIT_KERNEL_SUMMARY.json").read_text())
        inv=json.loads((HERE/"N_PAIRED_INVERSE_NO_GT_AUDIT.json").read_text())["summary"]
        r10=json.loads((HERE/"R10_OPTIMIZED_KERNEL_SUMMARY.json").read_text())["costs"]["full_jacobian"]
        labels=["R10 final full-J kernel\n32 nodes","Supplied digit full-J kernel\n36 nodes","Native inverse physics kernels\n9 parents, 80 J each","Native inverse total wall\n9 parents, 80 J each"]
        ratios=[1-r10["support_saving"],1-dig["all_datasets"]["full_jacobian"]["support_saving"],1-inv["kernel_saving"],1-inv["wall_saving"]]
        fig,ax=plt.subplots(figsize=(10,4.2));x=np.arange(4)
        ax.bar(x-.18,np.ones(4),.36,color="#9a9a9a",label="Full256")
        bars=ax.bar(x+.18,ratios,.36,color=["#227c9d"]*3+["#c05a36"],label="Exact support")
        for bar,ratio in zip(bars,ratios):ax.text(bar.get_x()+bar.get_width()/2,ratio+.025,f"{ratio:.3f}",ha="center",fontsize=10)
        ax.set_xticks(x,labels,fontsize=9);ax.set_ylabel("Time / paired Full256 time");ax.axhline(1,color="#555",ls="--",lw=1);ax.set_ylim(0,1.48)
        ax.set_title("Kernel savings do not establish end-to-end inversion speedup");ax.legend(loc="upper left",frameon=False);ax.grid(axis="y",alpha=.2)
        fig.tight_layout();fig.savefig(HERE/"KERNEL_VS_NATIVE_INVERSE_COST.png",dpi=180);plt.close(fig)
    print(json.dumps(groups,indent=2))


if __name__=="__main__":main()
