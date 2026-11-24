"""Measured-case development comparison."""

import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
def main():
    data=json.loads((HERE/"development_twin_scores.json").read_text())["rows"]
    labels = ["Original Dual-CIE","PR reset only","Diagonal optical","Full optical + sector","Pixel optical"]
    fig,ax=plt.subplots(figsize=(7.8,3.55),layout="constrained")
    centers=np.arange(len(data))
    ax.barh(centers-.18,[100*r["fit"]["physical_train_relative"] for r in data],
            .32,color="#3a6a9d",label="Training transmitters")
    ax.barh(centers+.18,[100*r["shape_scores"]["held_source_relative"] for r in data],
            .32,color="#bd6534",label="Held transmitters")
    ax.set_yticks(centers,labels)
    ax.invert_yaxis()
    ax.set_xlabel("Relative complex-field residual (%)")
    ax.set_title("Fresnel Twin, 8 GHz at the same 30 s limit",fontsize=12)
    ax.legend(loc="upper center",bbox_to_anchor=(.5,-.15),ncol=2,frameon=False,fontsize=8)
    ax.grid(axis="x",alpha=.25)
    ax.set_axisbelow(True)
    ax.spines[["top","right"]].set_visible(False)
    fig.savefig(HERE/"twin_development.png",dpi=170)
    fig.savefig(HERE/"twin_development.pdf")

if __name__=="__main__":
    main()
