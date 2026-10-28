"""Analytical illustrations for R7-T. No Maxwell solver, data, or learned model."""
from pathlib import Path
from io import BytesIO
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/p3lr7_mpl")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import gamma, gammainc, roots_hermitenorm

OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": .2,
    "figure.constrained_layout.use": True,
    "pdf.fonttype": 42,
})
blue, orange, green = "#235789", "#D76A24", "#258773"
def save_pdf(fig, name):
    buffer = BytesIO()
    fig.savefig(buffer, format="pdf")
    (OUT/name).write_bytes(buffer.getvalue())

h = np.linspace(-.8, .8, 401)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.5))
ax[0].plot(h, -h*h, color=blue, label=r"$z_1(h)=-h^2$")
ax[0].plot(h, h*h, color=orange, label=r"$z_2(h)=h^2$")
ax[0].plot(h, np.zeros_like(h), "--", color=green, label=r"joint optimum $\bar z=0$")
ax[0].set(xlabel=r"material coordinate $h$", ylabel=r"required compensation $z$",
          title="Individually exact, jointly incompatible")
ax[0].legend(frameon=False)
ax[1].plot(h, np.sqrt(2)*h*h, color=green,
           label=r"$\min_z\|(z+h^2,z-h^2)\|_2=\sqrt{2}h^2$")
ax[1].axhline(.15, color="0.35", linestyle=":", label="example noise threshold")
ax[1].set(xlabel=r"material coordinate $h$", ylabel="joint response norm",
          title="Joint first derivative remains zero at the origin")
ax[1].legend(frameon=False, fontsize=8)
save_pdf(fig,"shared_compensation.pdf")
plt.close(fig)

sigma=np.logspace(-6,-1,180)
def branch_moment(s):
    q=(s+1)/6
    return (2*sigma*sigma)**q*gamma(q)*gammainc(q,1/(2*sigma*sigma))/6
i0=branch_moment(0)
mass=i0/(1+i0)
task_regular=branch_moment(1)/(1+i0)
task_decoder=branch_moment(3)/(1+i0)
fig,ax=plt.subplots(figsize=(8.6,3.8))
ax.loglog(sigma,mass,color=blue,label=r"omitted branch mass: $O(\sigma^{1/3})$")
ax.loglog(sigma,task_regular,color=orange,label=r"regular material task $h$: $O(\sigma^{2/3})$")
ax.loglog(sigma,task_decoder,color=green,label=r"decoder-flat task $h^3$: $O(\sigma^{4/3})$")
ax.set(xlabel=r"noise scale $\sigma$", ylabel="exact normalized omission error",
       title=r"Toy prior: $\frac{1}{2}\delta_0+\frac{1}{2}\,dh$ on $[0,1]$, response $h^3$, observation 0")
ax.legend(frameon=False,fontsize=9)
save_pdf(fig,"atom_task_scaling.pdf")
plt.close(fig)

z=np.linspace(-4,4,500)
fig,ax=plt.subplots(1,2,figsize=(10,3.6))
def tmean(z,t):
    u=np.exp(-t/2)*np.cosh(np.sqrt(t)*z)
    return u/(1+u)
for t,col in [(0.1,blue),(1,orange),(4,green)]:
    ax[0].plot(z,tmean(z,t),color=col,label=f"t = {t:g}")
ax[0].set(xlabel="new observation z",ylabel=r"$E[M^2\mid Y_t=z]$",
          title=r"$\mathrm{Cov}(M^2,M)=0$ for every $t$")
ax[0].legend(frameon=False)
nodes,weights=roots_hermitenorm(160)
weights=weights/np.sqrt(2*np.pi)
ts=np.logspace(-3,0,120)
infos=[]
for t in ts:
    v=0
    for m,p in [(-1,.25),(0,.5),(1,.25)]:
        q=tmean(nodes+np.sqrt(t)*m,t)-.5
        v+=p*np.sum(weights*q*q)
    infos.append(v)
ax[1].loglog(ts,infos,color=green,label="conditional-expectation information")
ax[1].loglog(ts,ts**2/32,"--",color="0.35",label=r"leading term $t^2/32$")
ax[1].set(xlabel="new measurement strength t",ylabel="mean task-risk reduction",
          title="Finite task information is positive")
ax[1].legend(frameon=False,fontsize=8)
save_pdf(fig,"zero_covariance_information.pdf")
plt.close(fig)
print("Generated three analytical figures; no Maxwell calculations.")

