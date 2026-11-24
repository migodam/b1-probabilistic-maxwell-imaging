"""Show paired differences for the independent phase-penalty ablation."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
METRICS = (
    ("inside_L2", "Interior material error"),
    ("background_L2", "Air/background error"),
    ("held", "Independent complex-port error"),
)
GROUPS = (("PW360_10x20", "Full ring", "#006e87"), ("PW180_10x20", "Half ring", "#a14445"))


def plot():
    data = json.loads((HERE/"scores.json").read_text(encoding="utf-8"))
    fig, (a,b) = plt.subplots(1,2,figsize=(9.4,3.8),layout="constrained")
    for i,(key,label) in enumerate(METRICS):
        for k,(array,title,color) in enumerate(GROUPS):
            sub = [row for row in data["paired"] if row["left"]=="profile" and row["right"]=="pixel"
                   and row["array"]==array]
            values=np.array([row["difference"][key] for row in sub])
            y = i + (.18 if k else -.18)
            a.scatter(values * 1000, y+np.linspace(-.03,.03,len(values)), color=color,
                      alpha=.65,s=20,label=title if i==0 else None)
            a.scatter(np.median(values) * 1000,y,color=color,marker="D",s=38)
    a.set_yticks(np.arange(len(METRICS)),[label for _,label in METRICS])
    a.invert_yaxis()
    a.set_title("Profile − Pixel")
    a.set_xlabel("Paired error difference (×10⁻³)")
    a.legend(loc="best",frameon=False)
    for j,(method,color,mark) in enumerate((("pixel","#625189","o"),("profile","#28775b","^"))):
        for i,(array,label,_) in enumerate(GROUPS):
            sub=[row for row in data["paired"] if row["left"]==method and row["right"]=="reset"
                 and row["array"]==array]
            y=np.array([row["difference"]["held"] for row in sub])
            place=i+(.15 if j else -.15)
            b.scatter(y,place+np.linspace(-.03,.03,len(y)),color=color,s=23,alpha=.65,
                      marker=mark,label=method.capitalize() if i==0 else None)
            b.scatter(np.median(y),place,color=color,marker="D",s=40)
    b.set_yticks(range(len(GROUPS)),[row[1] for row in GROUPS])
    b.invert_yaxis()
    b.set_title("Hidden field vs. Reset")
    b.set_xlabel("Paired relative complex-port error difference")
    b.legend(loc="best",frameon=False)
    for ax in (a,b):
        ax.axvline(0,color="#58616d",linewidth=1)
        ax.grid(axis="x",color="#dce3e8",linewidth=.75)
        ax.spines[["top","right"]].set_visible(False)
    fig.suptitle("Four independent test materials · Points are objects, diamonds are medians")
    for suffix in (".pdf",".png"):
        fig.savefig(HERE/f"paired_mechanism{suffix}",dpi=220,facecolor="white")


if __name__ == "__main__":
    plot()
