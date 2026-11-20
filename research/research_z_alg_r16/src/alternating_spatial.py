"""Publishable spatial diagnostics for the identical public-A annulus observation."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from alternating_tsom import CurrentStep
from cohorts import historical_study, load_cohort, load_truth
from score_materials import material_metrics
from voxel_gn import Wave

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT/"results"
CASE = "r10_twod_0001"
COHORT = "r10_twod"


def plot_raster(axis, values, model, vmax, cmap):
    shape = model["shape"]
    grid = model["coords"].reshape(*shape,2)
    dx = (grid[0,1,0]-grid[0,0,0])/2
    dy = (grid[1,0,1]-grid[0,0,1])/2
    assert np.allclose(grid[:,:,0],grid[0,:,0]) and np.allclose(grid[:,:,1],grid[:,0,1,None])
    extent = (grid[0,0,0]-dx,grid[0,-1,0]+dx,grid[0,0,1]-dy,grid[-1,0,1]+dy)
    out = axis.imshow(values.reshape(shape),origin="lower",extent=extent,interpolation="none",
                      aspect="equal",cmap=cmap,vmin=0,vmax=vmax)
    axis.set_xticks((-.4,0,.4))
    axis.set_yticks((-.4,0,.4))
    axis.tick_params(labelsize=8,length=3)
    return out


def save(fig, name):
    for suffix in ("png","pdf"):
        fig.savefig(ROOT/"figures"/f"{name}.{suffix}",dpi=190,bbox_inches="tight")
    plt.close(fig)


def load_cie(model, case):
    a_path = Path("supplement_paired/inputs")/f"A__{CASE}.npz"
    mature_path = Path("RESULTS/twod")/f"{CASE}__cie_single__noise_ball_99"/"prediction.npz"
    with np.load(historical_study(12)/a_path,allow_pickle=False) as data:
        assert all(np.array_equal(data[k],model[k]) for k in ("G","S","E0","coords"))
        assert np.array_equal(data["Y"],case["Y"])
        assert float(data["sigma"])==case["sigma"]
        initial = data["m_init"].ravel()
    with np.load(historical_study(10)/mature_path,allow_pickle=False) as data:
        assert np.array_equal(data["m"].ravel(),initial)
    with (RESULTS/"baselines.csv").open() as stream:
        baseline = next(r for r in csv.DictReader(stream)
                        if r["case_id"]==CASE and r["method"]=="cie_single" and r["historical_view"]=="historical_final")
    return initial,dict(public_a=f"R12/{a_path}",historical_cie=f"R10/{mature_path}",
                        same_read_and_operators=True,same_public_a_material=True,
                        historical_complex_rhs=int(baseline["physical_complex_RHS_count"]),
                        historical_amplitude_NRMSE=float(baseline["amplitude_NRMSE"]))


def material_figure(model, case, mats, metrics):
    labels = {
        "truth": "True material",
        "fixed_current": "Fixed current",
        "fft_dark": "49 current modes",
        "voxel_only": "Full-wave voxel GN",
        "historical_cie": "Historical CIE",
    }
    fig,axes = plt.subplots(1,5,figsize=(12.6,3.25),layout="constrained")
    for i,(key,axis) in enumerate(zip(labels,axes,strict=True)):
        im = plot_raster(axis,mats[key],model,model["bounds"][1],"viridis")
        metric = metrics[key]
        suffix = ("\n\n" if key=="truth" else
                  f"\nAmp. NRMSE {metric['amplitude_NRMSE']:.3f} · IoU {metric['fixed_support_IoU']:.2f}"
                  f"\nMaxwell residual {metric['noise_rms']:.1f}×")
        axis.set_title(f"({chr(97+i)}) {labels[key]}"+suffix,fontsize=9,pad=6)
        if i:
            axis.tick_params(labelleft=False)
        else:
            axis.set_ylabel("y (model length unit)",fontsize=9)
    fig.supxlabel("x (model length unit)   ·   Residual shown in noise RMS units",fontsize=9)
    fig.suptitle("Public A · r10_twod_0001 · High-contrast annulus",fontsize=12)
    fig.colorbar(im,ax=list(axes),label="Material contrast (dimensionless)",shrink=.8,pad=.02)
    save(fig,"alternating_public_a")


def current_figure(model, before, after, source):
    delta = after-before
    vmax = float(np.max(np.abs(np.r_[before,after])))
    upper = float(np.max(np.abs(delta)))
    fig,axes=plt.subplots(1,3,figsize=(9.0,3.23),layout="constrained")
    plot_raster(axes[0],abs(before),model,vmax,"magma")
    im=plot_raster(axes[1],abs(after),model,vmax,"magma")
    changed=plot_raster(axes[2],abs(delta),model,upper,"cividis")
    for i,(axis,title) in enumerate(zip(axes,(
        r"(a) Start $|J_d|$",r"(b) After update $|J_d+Qa|$",r"(c) Actual change $|\delta j|$"),strict=True)):
        axis.set_title(title,fontsize=10)
        axis.set_xlabel("x",fontsize=9)
        if not i:
            axis.set_ylabel("y (model length unit)",fontsize=9)
        else:
            axis.tick_params(labelleft=False)
    fig.colorbar(im,ax=list(axes[:2]),label="Normalized current",shrink=.8,pad=.02)
    fig.colorbar(changed,ax=axes[2],label="Normalized change",shrink=.8,pad=.02)
    fig.suptitle(f"Public A · {CASE} · First current update · Source 1 at ({source[0]:.2f}, {source[1]:.2f})",
                 fontsize=11)
    save(fig,"alternating_public_a_current")


def run():
    model,cases = load_cohort(COHORT)
    case = next(c for c in cases if c["id"]==CASE)
    truth = load_truth(COHORT,CASE).ravel()
    with np.load(RESULTS/f"alternating_{COHORT}.npz",allow_pickle=False) as stored:
        mats = {"truth": truth, **{k: stored[f"{CASE}/{k}"] for k in ("fixed_current","fft_dark","voxel_only")}}
    mats["historical_cie"],historical = load_cie(model,case)
    metrics = {}
    for key,value in mats.items():
        prediction = Wave(model).state(value)[0]
        scores = material_metrics(value,truth,model["shape"],model["bounds"])
        metrics[key] = dict(
            **{k:scores[k] for k in ("relativeL2","RMSE","amplitude_NRMSE","fixed_support_IoU")},
            noise_rms=float(np.linalg.norm(prediction-case["Y"])/(case["sigma"]*np.sqrt(case["Y"].size))))
    assert abs(metrics["historical_cie"]["amplitude_NRMSE"]-historical["historical_amplitude_NRMSE"])<1e-12
    motion = CurrentStep(model,case,"fft_dark",halfwidth=3,tradeoff=.32)
    old = motion.anchor[:,0]
    next_current = motion.update(np.ones(len(model["G"])))[:,0]
    change = next_current-old
    material_figure(model,case,mats,metrics)
    source = model["sources"][0]
    current_figure(model,old,next_current,source)
    meta = dict(case=CASE,cohort=COHORT,archived_public_a=historical,
                spatial_unit="model length unit",material_unit="unitless real contrast",
                current_unit="normalized pulse VIE source-field unit",methods=metrics,
                metric_definitions=dict(relativeL2="L2(error) / L2(truth)",
                      amplitude_NRMSE="RMSE / (maximum truth minus 0.1)",
                      fixed_support_IoU="Intersection over union of material > 1.0"),
                one_current_source=dict(source_number=1,source_coordinate=source.tolist(),current_rank=motion.rank,
                    initial_material=1.,before="J_d",after="J_d+Qa",change_norm=float(np.linalg.norm(change)),
                    change_max=float(np.max(abs(change))),
                    change_relative_to_initial_norm=float(np.linalg.norm(change)/np.linalg.norm(old)),
                    observation_change_over_noise_rms=float(np.linalg.norm(model["S"]@change)/(case["sigma"]*np.sqrt(len(case["Y"])))) ))
    (RESULTS/"alternating_public_a_spatial.json").write_text(json.dumps(meta,indent=2,allow_nan=False)+"\n")
    np.savez_compressed(RESULTS/"alternating_public_a_current.npz",coordinate=model["coords"],
                        initial=old,updated=next_current,delta=change)
    print(json.dumps(meta,indent=2))


if __name__=="__main__":
    run()
