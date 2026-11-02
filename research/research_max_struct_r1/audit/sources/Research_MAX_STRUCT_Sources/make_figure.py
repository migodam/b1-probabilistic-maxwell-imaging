"""Illustrate an exact local excitation formula; this is not a Maxwell simulation."""
from pathlib import Path
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/max_struct_matplotlib")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
(HERE / "figures").mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 160})
x = np.linspace(0, .7, 300)
curves = [np.ones_like(x), .7*np.ones_like(x), x,
          np.sqrt((1+.7**2)/2)*np.ones_like(x),
          np.sqrt((1+x*x)/2), np.sqrt((.7**2+x*x)/2)]
labels = [r"$\sigma_1$", r"$\sigma_2$", r"$\sigma_3$",
          r"$\sqrt{(\sigma_1^2+\sigma_2^2)/2}$",
          r"$\sqrt{(\sigma_1^2+\sigma_3^2)/2}$",
          r"$\sqrt{(\sigma_2^2+\sigma_3^2)/2}$"]
colors = ["#355070", "#6d597a", "#d45a30", "#9b7d2b", "#287d8e", "#5c8b5b"]
fig, ax = plt.subplots(1, 2, figsize=(10.4, 4.1), gridspec_kw={"width_ratios":[1.45,1]})
for y, label, color in zip(curves, labels, colors):
    ax[0].plot(x, y, lw=2.0, label=label, color=color)
ax[0].set(xlabel=r"Weak field singular value $\sigma_3$", ylabel="Tensor excitation singular values",
          title=r"Exact local spectrum ($\sigma_1=1,\ \sigma_2=0.7$)", xlim=(0,.7), ylim=(0,1.08))
ax[0].legend(loc="lower right", ncol=2, fontsize=8.5, frameon=False, handlelength=1.6)
ax[0].grid(alpha=.16)
rad=np.minimum(1., .02/np.maximum(x,1e-20))
ax[1].plot(x,rad, color="#d45a30",lw=2.3)
ax[1].fill_between(x,0,rad,color="#d45a30",alpha=.12)
ax[1].set(xlabel=r"Weak field singular value $\sigma_3$", ylabel="Guaranteed pair-close material radius",
          title=r"Conditional bound ($t_L=1,\ \delta/K=0.02$)", xlim=(0,.7),ylim=(0,1.08))
ax[1].grid(alpha=.16)
fig.tight_layout(pad=1.4)
fig.savefig(HERE / "figures" / "tensor_excitation_spectrum.png",dpi=220,bbox_inches="tight")
plt.close(fig)

