#!/usr/bin/env python3
"""ANALYTICAL_UNIT_NOT_MAXWELL: scalar cubic saturation and posterior tasks.

X~Uniform[-1,1], M=max(X,0)^3, F=M, y=0,
likelihood exp(-M^2/(2*sigma^2)). No Maxwell or full66 computation.
Reproduce: python analytical_units.py
"""
import json
import os
from pathlib import Path
import time

import numpy as np
from scipy.integrate import quad
from scipy.special import gamma, gammainc
os.environ.setdefault("MPLCONFIGDIR",str(Path(__file__).resolve().parent/"report_build"/"mpl_cache"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent


def main():
    start = time.perf_counter()
    sigma = np.logspace(-3, -12, 37)
    # Ik=int_0^1 x^(3k) exp(-x^6/(2*sigma^2)) dx
    #   =Ck sigma^(k+1/3) P((3k+1)/6,1/(2*sigma^2)).
    # Scaled gamma formulas avoid integrating a vanishingly narrow peak in x.
    shapes = (3*np.arange(3)+1)/6
    constants = gamma(shapes)*2**shapes/6
    integrals = np.array([
        c*sigma**(k+1/3)*gammainc(shape, 1/(2*sigma**2))
        for k,(c,shape) in enumerate(zip(constants,shapes))
    ])
    i0,i1,i2 = integrals
    evidence = (1+i0)/2
    atom = 1/(1+i0)
    # Evaluate the complement directly, rather than subtract near-one numbers.
    nonhost = i0/(1+i0)
    mean = i1/(1+i0)
    second_moment = i2/(1+i0)
    mean_sq = mean**2
    variance = second_moment-mean_sq
    arrays = {
        "nonhost_posterior_mass": nonhost,
        "material_posterior_mean": mean,
        "material_posterior_variance": variance,
        "squared_mean_gap_from_zero": mean_sq,
    }
    expected = np.array([1/3,4/3,7/3,8/3])
    leading = np.array([constants[0],constants[1],constants[2],constants[1]**2])
    small_noise = sigma <= 1.000001e-9
    slopes = {}
    for (name,values),exponent in zip(arrays.items(),expected):
        fit_all = float(np.polyfit(np.log(sigma),np.log(values),1)[0])
        fit_small = float(np.polyfit(np.log(sigma[small_noise]),np.log(values[small_noise]),1)[0])
        slopes[name] = {
            "expected_asymptotic_exponent": float(exponent),
            "slope_all_sigma_1e_3_to_1e_12": fit_all,
            "slope_sigma_1e_9_to_1e_12": fit_small,
            "absolute_small_noise_slope_error": abs(fit_small-float(exponent)),
        }
        assert abs(fit_small-exponent) < 5e-4
    # Independent scaled-coordinate quadrature for the three leading constants.
    # t=x/sigma^(1/3), integrate [0,3]; exp(-3^6/2)<1e-158.
    scaled_quad=[]
    for k,c in enumerate(constants):
        value,err = quad(lambda t: t**(3*k)*np.exp(-t**6/2),0,3,
                         epsabs=1e-13,epsrel=1e-13)
        scaled_quad.append({"k":k,"gamma_constant":float(c),
                            "scaled_quadrature_0_to_3":value,
                            "quad_reported_error":err,
                            "absolute_difference":abs(value-c)})
        assert abs(value-c) < 1e-11
    assert np.all(variance>0)
    assert np.all((atom>0.5)&(atom<1))
    assert np.all(np.isfinite(integrals))

    # Scalar positive-layer material density after u=M/sigma scaling.
    # Its limit is u^(-2/3) exp(-u^2/2) / (3*C0).
    u=np.logspace(-4,np.log10(4),400)
    density=u**(-2/3)*np.exp(-u*u/2)/(3*constants[0])
    singular=u**(-2/3)/(3*constants[0])
    radial_slice=u<=0.01
    radial_slope=float(np.polyfit(np.log(u[radial_slice]),np.log(density[radial_slice]),1)[0])
    assert abs(radial_slope+2/3)<1e-4

    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,
                         "axes.titlesize":11,"axes.labelsize":10,
                         "legend.fontsize":9,"figure.facecolor":"white"})
    fig,axs=plt.subplots(1,2,figsize=(12.5,5.25),layout="constrained")
    colors=["#2563eb","#059669","#d97706","#9333ea"]
    labels=[r"Non-host mass: $\sigma^{1/3}$",
            r"Material mean: $\sigma^{4/3}$",
            r"Material variance: $\sigma^{7/3}$",
            r"Squared mean gap: $\sigma^{8/3}$"]
    for (name,values),c,beta,col,label in zip(arrays.items(),leading,expected,colors,labels):
        axs[0].loglog(sigma,values/c,color=col,marker="o",markersize=2.8,
                      linewidth=1.6,label=label)
        axs[0].loglog(sigma,sigma**beta,color=col,linestyle="--",linewidth=0.9,alpha=.5)
    axs[0].set_xlabel(r"Noise scale $\sigma$")
    axs[0].set_ylabel("Quantity / leading asymptotic constant")
    axs[0].set_title("A. Phase uncertainty and material tasks have different orders")
    axs[0].grid(True,which="major",alpha=.2)
    axs[0].legend(loc="lower right",framealpha=.95)
    axs[0].text(.03,.96,"Markers: exact gamma expression\nDashed: leading power law",
                transform=axs[0].transAxes,va="top",fontsize=8.5,
                bbox={"facecolor":"white","edgecolor":"none","alpha":.85})
    axs[1].loglog(u,density,color="#0f766e",linewidth=2,
                  label="Limiting positive-layer density")
    axs[1].loglog(u,singular,color="#374151",linestyle="--",linewidth=1.2,
                  label=r"Near-seam factor $u^{-2/3}/(3C_0)$")
    axs[1].set_xlabel(r"Scaled material amplitude $u=M/\sigma>0$")
    axs[1].set_ylabel("Conditional density in u")
    axs[1].set_title("B. Integrable material-density singularity at the seam")
    axs[1].grid(True,which="major",alpha=.2)
    axs[1].legend(loc="lower left",framealpha=.95)
    axs[1].text(.03,.96,"Scalar model; atom handled separately\nNo Maxwell data or full66 posterior",
                transform=axs[1].transAxes,va="top",fontsize=8.5,
                bbox={"facecolor":"white","edgecolor":"none","alpha":.85})
    fig.suptitle("ANALYTICAL UNIT — NOT MAXWELL\n"
                 r"$X\sim\mathrm{Unif}[-1,1],\ M=(\max\{X,0\})^3,\ F=M,\ y=0$",
                 fontsize=13,fontweight="bold")
    outdir=ROOT/"report_build"
    outdir.mkdir(parents=True,exist_ok=True)
    plot=outdir/"seam_task_scaling.png"
    fig.savefig(plot,dpi=220)
    plt.close(fig)

    result={
        "status":"ANALYTICAL_UNIT_NOT_MAXWELL",
        "model":{"prior":"X~Uniform[-1,1]","material":"M=max(X,0)^3",
                 "forward":"F=M","observation":0,"noise":"scalar N(0,sigma^2)",
                 "prior_host_atom_mass":0.5},
        "formula":{"Ik":"Ck sigma^(k+1/3) gammainc((3k+1)/6,1/(2sigma^2))",
                   "Ck":"2^((3k+1)/6) Gamma((3k+1)/6)/6",
                   "evidence":"(1+I0)/2","posterior_atom_mass":"1/(1+I0)",
                   "nonhost_mass":"I0/(1+I0)","mean":"I1/(1+I0)",
                   "variance":"I2/(1+I0) - (I1/(1+I0))^2"},
        "sigma":sigma.tolist(),"evidence":evidence.tolist(),
        "host_posterior_atom_mass":atom.tolist(),
        "quantities":{k:v.tolist() for k,v in arrays.items()},
        "leading_constants":{k:float(c) for k,c in zip(arrays,leading)},
        "loglog_slopes":slopes,
        "scaled_quadrature_checks":scaled_quad,
        "positive_layer_density_near_seam_fitted_exponent":radial_slope,
        "figure":str(plot),
        "cost":{"new_Maxwell_states":0,"new_Maxwell_RHS":0,
                "wall_seconds":time.perf_counter()-start},
        "qualification":[
            "Closed-form scalar analytical unit and float64 numerical checks only.",
            "No inference about original full66 posterior mass or high-SNR Maxwell scaling is made.",
            "Plot B is a scalar amplitude density, not a two-dimensional area density.",
            "Quad reported errors and slope fits are numerical checks, not directed-rounding certificates.",
        ],
    }
    (ROOT/"analytical_results.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"loglog_slopes":slopes,
                      "density_exponent":radial_slope,"cost":result["cost"],
                      "figure":str(plot)},indent=2))


if __name__=="__main__":
    main()
