from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm,colors
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'RESULTS/vector_3d';FIG=ROOT/'FIGURES/vector_3d';FIG.mkdir(exist_ok=True,parents=True)
methods=['ordinary_joint','dual_like_joint','A_material_joint','direct_physical'];labels=['Truth','Vector joint','Dual-like joint','A material joint','Direct physical']
names=['sphere','two_spheres','non_spherical_bracket'];scores=json.loads((OUT/'SCORES.json').read_text());lookup={(r['target'],r['method']):r for r in scores}
fig=plt.figure(figsize=(16,9));norm=colors.Normalize(0,4);cmap=cm.viridis
for row,name in enumerate(names):
    folder=OUT/name;truth=np.load(folder/'truth_for_scoring.npz')['chi'];values=[truth]+[np.load(folder/(m+'_prediction.npz'))['prediction'] for m in methods]
    for col,(lab,x) in enumerate(zip(labels,values)):
        ax=fig.add_subplot(3,5,5*row+col+1,projection='3d');x=x.reshape(5,5,5);rgba=cmap(norm(x));rgba[...,3]=.85
        ax.voxels(x>=.5,facecolors=rgba,edgecolor=(0,0,0,.15),linewidth=.25)
        ax.set(xlim=(0,5),ylim=(0,5),zlim=(0,5),xticks=[],yticks=[],zticks=[]);ax.view_init(24,-50);ax.set_box_aspect([1,1,1])
        title=lab if row==0 else ''
        if col==0:title+='\n'+name.replace('_',' ')
        else:title+='\nNRMSE %.3f'%lookup[(name,methods[col-1])]['material_NRMSE']
        ax.set_title(title,fontsize=10,pad=0)
fig.suptitle('True 3D vector material inversions: 125 unknown real scalar voxel contrasts\nDisplay threshold chi >= 0.5; all values and fields saved; common physical corrector',fontsize=13)
fig.subplots_adjust(left=.01,right=.93,bottom=.02,top=.89,wspace=.02,hspace=.03)
cax=fig.add_axes([.95,.20,.012,.55]);fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,label='chi = epsilon_r - 1')
fig.savefig(FIG/'vector_3d_voxels.png',dpi=170);fig.savefig(FIG/'vector_3d_voxels.pdf');plt.close(fig)
fig,axes=plt.subplots(3,5,figsize=(14,8),layout='constrained')
for row,name in enumerate(names):
    folder=OUT/name;values=[np.load(folder/'truth_for_scoring.npz')['chi']]+[np.load(folder/(m+'_prediction.npz'))['prediction'] for m in methods]
    for col,x in enumerate(values):
        ax=axes[row,col];v=x.reshape(5,5,5);tri=np.concatenate([v[:,:,1],v[:,:,2],v[:,:,3]],axis=1).T
        im=ax.imshow(tri,origin='lower',vmin=0,vmax=4,cmap='viridis',aspect='equal');ax.set_xticks([]);ax.set_yticks([])
        if row==0:ax.set_title(labels[col],fontsize=11)
        if col==0:ax.set_ylabel(name.replace('_',' '))
        else:ax.set_xlabel('NRMSE %.3f | residual %.4f'%(lookup[(name,methods[col-1])]['material_NRMSE'],lookup[(name,methods[col-1])]['physical_data_residual']),fontsize=8)
fig.colorbar(im,ax=axes,label='chi',shrink=.8);fig.suptitle('Three central z slices (indices 1,2,3), shared material color range; no interpolation')
fig.savefig(FIG/'vector_3d_slices.png',dpi=170);fig.savefig(FIG/'vector_3d_slices.pdf');plt.close(fig)

rows=[]
for r in scores:
    c=r['counts'];rows.append('|%s|%s|%.4f|%.4f|%.5f|%.3f|%d|%d|%d|%d|%d|'%(r['target'],r['method'],r['raw_NRMSE'],r['material_NRMSE'],r['physical_data_residual'],r['wall_s'],c['physical_LU'],c['physical_source_RHS'],c['physical_adjoint_RHS'],c['current_normal_factorizations'],c['material_value_gradient_calls']))
text='''# R7: actual small 3D vector material imaging benchmark

**Status: COMPLETED, 12/12 reconstruction arrays; bounded finite-DDA development evidence.** A material did not show consistent recovery improvement: it improves the two-sphere score slightly and worsens the sphere and bracket relative to the ordinary vector-adapted joint control. Direct physical optimization produces better final truth scores with substantially more physical solves. This experiment does not establish equal-compute superiority, continuum imaging accuracy, blind recovery, uniqueness, or a mature vector CIE implementation.

## Physical model and independent qualification

The inherited R1 `TensorDDA` is reused without alteration. It uses the **full 3D electric dyadic Green**, including its longitudinal near-field term, for exp(−iωt), λ₀=Z₀=ε_b=μ_b=1 and k=2π. Every cell has three electric field/current components. Integrated cubic self blocks are `s I₃`, where `s=−1/3+i k³V/(6π)`. There are 125 independently unknown **real isotropic** material contrasts χ_i=ε_ri−1, not 125 scalar-TM fields. Passive lossless materials have ε_r in [1,5]; radiation reaction remains in s. μ_r=1. No anisotropic or absorbing material recovery is claimed.

With density current j=χE and dipole p=Vj, eliminating the Lorentz self gives the CM+RR polarizability `α(χ)=Vχ/(1−sχ)`. This is nonlinear in material; we optimize χ throughout and include its exact chain effect through the self-inclusive field equation. An independently assembled polarizability-coordinate solve agrees with density-coordinate p to **1.50e−15 relative error**. It uses one additional 375×375 LU and six complex RHS beyond the adapter-counted qualification calls.

Independent analytic Mie sphere verification is **reused**, not repeated: unchanged R1 radius .12λ, ε_r=3, h=.04/.03λ, N=123/257. Absolute far-vector relative errors decrease from **10.47% to 3.66%**, and scattering cross-section relative errors from **21.81% to 5.99%**. This is a bounded DDA qualification, not an H(curl) or continuum certificate. The present imaging cube has **h=.08λ**, coarser than that sphere check, so no quantitative continuum imaging claim follows. Exact source/receiver DDA data are used below; discretization error cancels in these development reconstructions. See [qualification and inherited receipt](RESULTS/vector_3d/QUALIFICATION.json).

The physical adjoint directional derivative relative error is **3.32e−9**. A later necessary off-root check has nonzero relative constitutive residual .2048 and tests the complete material-dependent A derivative (including dP/dχ): ordinary **3.95e−11**, Dual-like **2.21e−8**, A **8.68e−9**. The initial constitutive-root check is retained; it mainly tested the prior and was insufficient on its own. See [off-root receipt](RESULTS/vector_3d/NONROOT_DERIVATIVES.json). Both checks cost 5 physical LU/30 source RHS/18 adjoint RHS in the adapter, plus the separate CM-coordinate LU/6 RHS and 18 material calls. No formula sweep or error-floor study was run.

## Frozen experiment and methods

A 5×5×5 cube has center spacing .08λ, 125 material DOFs and a **375-dimensional vector Maxwell system**. There are six distinct transverse plane-wave illuminations from ±x, ±y, ±z, with declared polarizations and 30 external receiver points on radius .8λ measuring all three scattered electric components: 540 complex observations. The three nontrivial dev targets are a staircase sphere, unequal two spheres, and a non-spherical three-arm bracket. All use the same box χ∈[0,4], uniform χ=.25 initialization, quadratic nearest-neighbor smoothing λ=2e−4 and l₂=1e−6, source normalization, and common physical corrector (L-BFGS-B maximum 40 iterations). The complete config was written before execution. No target-specific tuning, noise, extra starts, GT-based selection, or support information enters a solver.

The three joint controls share a current update that exactly minimizes the same ordinary data/state quadratic by one complex Hermitian Cholesky factorization with all six source RHS. At each of 15 fixed cycles the material block uses maximum 12 L-BFGS-B iterations. For `Q=Einc+Gj`, `r=j−χQ`, `b=β−s`, `d=1+bχ`, `a=χ/d`, β=1:

- **ordinary_joint:** material state residual r; labeled *vector-adapted joint baseline*, not mature 2D global CIE.
- **dual_like_joint:** material residual `u=βr/d`; this is a local resolvent control, not a faithful implementation of a named historical Dual algorithm.
- **A_material_joint:** material residual `(I+D_a H)u`. H retains actual full 3×3 dyadic off-diagonal blocks inside disjoint 2×2×2 voxel blocks and has zero 3×3 self blocks. Each scalar isotropic a_i multiplies its three-vector block. This correctly retains self via d and cross-polarization via H; it does not extend the scalar arc solver to tensor-valued material.
- **direct_physical:** exact self-inclusive Maxwell adjoint L-BFGS-B, maximum 180 iterations before the same 40-iteration corrector.

Every candidate selection uses the original physical data objective plus the identical prior. The direct-control iteration cap is comparable to the maximum 15×12 material iterations but **not equal physical LU, wall time, or linear algebra work**. Candidate predictions are saved before the scoring file is loaded; GT is used only after all four methods for a target have sealed predictions. The three targets are disclosed development objects, not blind holdouts. The common corrector may substantially improve the weak raw joint candidates; raw scores are explicitly retained.

## Actual reconstruction results

Material NRMSE is `||χhat−χGT||₂/||χGT||₂`; residual is `||S j(χhat)−Y||F/||Y||F`. Lower is better. Raw means public-physical-selected candidate before the shared corrector. The arrays include all 125 material values plus full vector fields/current, not just a chosen voxel or illustrative plot.

|Target|Method|Raw NRMSE|Final NRMSE|Physical residual|Wall s|Physical LU|Source RHS|Adjoint RHS|Aux Cholesky|Material calls|
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
'''+ '\n'.join(rows)+'''

![Actual 3D voxel reconstructions](FIGURES/vector_3d/vector_3d_voxels.png)

![Three central slices, actual material arrays](FIGURES/vector_3d/vector_3d_slices.png)

## Costs, failures, and scope

The complete sequential benchmark took **3.607 s wall / 3.556 s CPU**, with **102,449,152 bytes peak process RSS** (macOS ru_maxrss; benchmark process including setup/qualification; per-method peaks NOT_MEASURED). Setup cost was .0150 s wall / .0142 s CPU, one dyadic Green assembly and one external-receiver assembly. Each joint method/target uses 15 auxiliary Cholesky factorizations and 90 current RHS, in addition to its listed physical costs. The three public-data simulation calls cost 3 physical LU and18 source RHS. The necessary later nonroot check costs .0175 s wall and is separate from the 3.607-s benchmark timing. Figure/report generation is not included in solver timing.

No failed physical or optimization calls occurred. Optimizer cap exits are preserved verbatim in per-method receipts and are **not convergence certificates**. No GPU work, GT-based online choice, repeated Mie runs, mesh sweeps, external training, or model-specific parameter scans occurred. Peak hardware telemetry beyond process RSS was NOT_MEASURED. No 3D finite-element/H(curl) continuum acceptance or independent discretization inverse test was performed.

The decisive evidence is modest: full vector 3D inversion runs, valid exact derivatives and honest physical forward validation provenance now exist; **A material recovery superiority does not survive this three-object test consistently**. Fine-grained baseline maturity and unequal physical compute restrict any ranking. These results cannot be rebranded as scalar-TM theorem validation, generic material uniqueness, or a broadly superior 3D method.

## Reproduction and artifact inventory

Run from repository root with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` and the project Python:

```sh
python research/research_max_struct_r7/experiments/vector_3d_reconstruct.py
python research/research_max_struct_r7/experiments/vector_derivative_nonroot.py
python research/research_max_struct_r7/experiments/vector_plot_report.py
```

The reconstruction runner refuses to overwrite a completed benchmark. Preserve/archive the original result folder before a deliberate reproduction. Each target directory contains `public_observations.npz`, scorer-only `truth_for_scoring.npz`, four `*_prediction.npz` files (initial/raw/final materials, j,E), and four `*_receipt.json` files. Shared config, qualification, simulation costs, scores, completion/memory receipt and off-root derivative receipt are under [RESULTS/vector_3d](RESULTS/vector_3d/). Source is [vector_inverse.py](src/vector_inverse.py); execution log is [delegated/vector/run.log](delegated/vector/run.log).
'''
(ROOT/'BENCHMARK_3D.md').write_text(text)
