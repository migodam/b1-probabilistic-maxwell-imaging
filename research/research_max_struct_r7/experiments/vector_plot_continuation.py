from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'RESULTS/vector_continuation';FIG=ROOT/'FIGURES/vector_continuation';FIG.mkdir(exist_ok=True,parents=True)
methods=['ordinary_joint','dual_like_joint','A_material_joint','direct_physical'];labels=['Vector joint','Dual-like joint','A material joint','Direct physical']
cases=['grid125_sphere','grid125_two_spheres','grid125_non_spherical_bracket','grid216_sphere_continued'];scores=json.loads((OUT/'SCORES.json').read_text());lookup={(r['case'],r['method']):r for r in scores}
matched=[]
for metric,ykey,ylab in [('material','incumbent_NRMSE','Material NRMSE (offline GT score)'),('physical','incumbent_data_residual','Original physical data residual')]:
    fig,axes=plt.subplots(4,2,figsize=(12,14),layout='constrained')
    for row,case in enumerate(cases):
        curves={kind:np.load(OUT/case/kind/'offline_scored_curve.npz') for kind in methods}
        for col,(cost,xlab) in enumerate([('total_physical_LU','Cumulative physical LU (aux Cholesky excluded)'),('total_cpu_s','Cumulative solver CPU seconds (aux work included)')]):
            ax=axes[row,col];lower=max(c[cost][0] for c in curves.values())
            caps=[curves[k][cost][-1] for k in methods if 'LIMIT' in lookup[(case,k)]['status']]
            upper=min(caps) if caps else max(c[cost][-1] for c in curves.values())
            for kind,label in zip(methods,labels):
                c=curves[kind];line=ax.plot(c[cost],c[ykey],label=label,linewidth=1.3)[0]
                if 'CONVERGENCE' in lookup[(case,kind)]['status'] and c[cost][-1]<upper:ax.plot([c[cost][-1],upper],[c[ykey][-1]]*2,'--',color=line.get_color(),linewidth=1)
            ax.axvspan(lower,upper,color='grey',alpha=.08);ax.set_xlabel(xlab);ax.set_ylabel(ylab);ax.set_title(case.replace('_',' '),fontsize=10);ax.grid(alpha=.2)
            if metric=='material':
                budgets=np.linspace(lower,upper,3)
                for budget in budgets:
                    vals={}
                    for kind in methods:
                        c=curves[kind];idx=max(0,np.searchsorted(c[cost],budget,side='right')-1);vals[kind]=dict(NRMSE=float(c['incumbent_NRMSE'][idx]),physical_residual=float(c['incumbent_data_residual'][idx]),actual_physical_LU=float(c['total_physical_LU'][idx]),actual_cpu_s=float(c['total_cpu_s'][idx]))
                    matched.append(dict(case=case,cost=cost,budget=float(budget),common_range=[float(lower),float(upper)],scores=vals))
    axes[0,0].legend(fontsize=8);fig.suptitle('Shared300-iteration physical continuation from sealed endpoints\nSolid: observed public-Phi incumbents; dashed: converged endpoint held without extra compute; shaded: common comparison range',fontsize=11)
    fig.savefig(FIG/('vector_'+metric+'_cost_curves.png'),dpi=160);fig.savefig(FIG/('vector_'+metric+'_cost_curves.pdf'));plt.close(fig)
(OUT/'MATCHED_COST_SCORES.json').write_text(json.dumps(matched,indent=2))
# Latest actual material slices, all endpoints.
fig,axes=plt.subplots(4,5,figsize=(14,12),layout='constrained')
for row,case in enumerate(cases):
    ns=6 if '216' in case else 5
    source=OUT/'grid216_sphere' if ns==6 else ROOT/'RESULTS/vector_3d'/case.removeprefix('grid125_')
    truth=np.load(source/'truth_for_scoring.npz')['chi'];values=[truth]+[np.load(OUT/case/k/'prediction_trace.npz')['prediction'] for k in methods]
    for col,x in enumerate(values):
        ax=axes[row,col];v=x.reshape(ns,ns,ns);tri=np.concatenate([v[:,:,ns//2-1],v[:,:,ns//2],v[:,:,ns//2+1]],axis=1).T
        im=ax.imshow(tri,origin='lower',vmin=0,vmax=4,cmap='viridis');ax.set_xticks([]);ax.set_yticks([])
        if row==0:ax.set_title((['Truth']+labels)[col],fontsize=10)
        if col==0:ax.set_ylabel(case.replace('_',' '),fontsize=9)
        else:ax.set_xlabel('NRMSE %.4f'%lookup[(case,methods[col-1])]['final_NRMSE'],fontsize=9)
fig.colorbar(im,ax=axes,label='chi = epsilon_r - 1',shrink=.8);fig.suptitle('Actual 3D material slices after common longer correction; no interpolation')
fig.savefig(FIG/'vector_continued_slices.png',dpi=170);fig.savefig(FIG/'vector_continued_slices.pdf');plt.close(fig)
# Update report around final evidence while preserving the original12 table.
p=ROOT/'BENCHMARK_3D.md';old=p.read_text()
old=old.replace('**Status: COMPLETED, 12/12 reconstruction arrays; bounded finite-DDA development evidence.** A material did not show consistent recovery improvement: it improves the two-sphere score slightly and worsens the sphere and bracket relative to the ordinary vector-adapted joint control. Direct physical optimization produces better final truth scores with substantially more physical solves. This experiment does not establish equal-compute superiority, continuum imaging accuracy, blind recovery, uniqueness, or a mature vector CIE implementation.', '**Status: COMPLETED: initial12 reconstructions,216-DOF sphere4 reconstructions, and16 common long physical continuations.** The initially large direct-versus-joint recovery gaps mostly disappear after longer shared physical correction. All four methods approach the same regularized material solution in each tested problem; A has no durable 3D basin or recovery advantage here. This is bounded finite-DDA development evidence, not continuum imaging accuracy, blind recovery, uniqueness, or a mature vector CIE implementation. The initial12 results and all costs remain archived below.')
old=old.replace('## Actual reconstruction results','## Initial12 reconstruction results (preserved; longer correction below)')
rows=[]
for r in scores:rows.append('|%s|%s|%.5f|%.5f|%.6f|%d|%d|%.3f|%s|'%(r['case'],r['method'],r['initial_NRMSE'],r['final_NRMSE'],r['physical_data_residual'],r['total_physical_LU'],r['total_auxiliary_Cholesky'],r['total_cpu_s'],'cap300' if 'LIMIT' in r['status'] else 'optimizer tolerance'))
mie=json.loads((OUT/'MIE_ACTUAL_GRIDS.json').read_text());frac=json.loads((ROOT/'RESULTS/vector_fractional_validation/RECEIPT.json').read_text());complete=json.loads((OUT/'COMPLETE.json').read_text())
append='''
## Shared longer physical correction and necessary216-DOF case

The initial12 were retained intact. At the reviewer's request, **each of the four sealed final endpoints** received the same additional physical L-BFGS-B allowance: maximum300 iterations, ftol1e−13, gtol1e−9, same material box and prior. No auxiliary objective was used in this continuation; every source/adjoint solve used the original complete Maxwell model. The public physicalPhi incumbent determines every selection. Trial material arrays, physical costs, and CPU times were sealed before their entire curves were scored offline against truth.

One additional sphere uses a6×6×6 grid (**216 material DOFs,648 vector components**) with the same physical cube width .4λ, h=.4/6, same sources, receivers, initialization, priors, box, original iteration caps and common correction. No geometry or prior scan was done. It is a same-grid synthetic development inversion; it is not an independent-discretization inverse validation.

|Case|Method|Before long correction NRMSE|After long correction NRMSE|Physical residual|Total physical LU|Aux Cholesky|Total solver CPU s|Stop|
|---|---|---:|---:|---:|---:|---:|---:|---|
'''+ '\n'.join(rows)+'''

The sphere solutions agree in NRMSE to roughly5e−5 across methods at125DOF and1e−4 at216DOF. Two-sphere/bracket joint continuations reach the300 cap, while direct reaches its numerical tolerance; their remaining differences are small, and **not a proof of identical mathematical minima or full convergence**. There is no robust A recovery advantage. The initial direct-control ranking largely reflected the limited40-step corrector on joint endpoints.

![Actual materials after long correction](FIGURES/vector_continuation/vector_continued_slices.png)

### Equal-cost curves

These plots show public-physicalPhi incumbents versus **cumulative physical LU** and cumulative actual solver CPU. Truth is used only to score the sealed curves. LU excludes the15 auxiliary Cholesky factorizations of each joint baseline; CPU includes auxiliary work and measured overhead. Dashed horizontal lines only retain an already-converged endpoint without inventing additional computation. Capped endpoints are not extended. Grey regions mark the ranges where all methods have observations or a converged endpoint. The initial-stage trajectories were not stored, so the comparison ranges start only after every method's continuation has begun. This prevents a claim about earlier equal-cost performance.

![Equal-cost material curves](FIGURES/vector_continuation/vector_material_cost_curves.png)

![Equal-cost physical residual curves](FIGURES/vector_continuation/vector_physical_cost_curves.png)

Three matched budgets inside each common range, for both cost axes, are saved in [MATCHED_COST_SCORES.json](RESULTS/vector_continuation/MATCHED_COST_SCORES.json). Selecting by public physicalPhi can worsen material NRMSE slightly; the truth score is never used to choose an iterate.

### Actual-grid independent Mie validation and one geometry correction

The new check uses the actual125/216 grids, radius .145λ and ε_r=3, with exact analytic Mie far vectors. Binary center-occupancy has27 occupied cells at125DOF (a3³ staircase cube),32 at216DOF. Occupied volume changes from .013824 to .0094815 versus analytic .0127701. Far-vector relative error **6.584%→27.588%**: **FAIL for monotone refinement**. Both original errors and arrays are retained in [MIE_ACTUAL_GRIDS.json](RESULTS/vector_continuation/MIE_ACTUAL_GRIDS.json). A small linear residual does not erase this physical geometry error.

A single physically motivated correction was authorized **for forward validation only**: a fixed24³ midpoint subcell rule estimates sphere volume fraction in every cell on both grids; isotropic volume-averaged contrast is χ_i=f_i(3−1). This is an interface/pulse approximation, not exact cut-cell Maxwell homogenization. There was no geometry sweep, subcell-resolution search, fitted parameter, or inverse support reduction. Effective volumes become .0127639/.0127788, and Mie far errors become **8.396%→8.199%**. The marginal decrease confirms that geometry aliasing contributes substantially; the remaining≈8% error prevents a strong continuum qualification. Neither inversion data nor inverse unknowns were altered. See [fractional validation receipt](RESULTS/vector_fractional_validation/RECEIPT.json).

### Additional costs and retained limits

The continuation/216/Mie-actual-grid process took **36.444 s wall /34.375 s CPU**, with161,185,792bytes process peak RSS. It includes four original216-grid reconstructions,16 continuations, two binary-sphere Mie forward solves, setup and scoring; plot generation is excluded. The separate fixed-fraction geometry check took **.2230 s wall /.2185 s CPU**, with268,042,240bytes peak RSS including its subcell temporary array; it used two physical LU and12 source RHS. These are separate process peaks, not summed memory. All physical costs, trial arrays and stop messages are retained; there were no failed calls, and capped exits remain visible.

Reproduce the added work, after preserving the archived outputs:

```sh
python research/research_max_struct_r7/experiments/vector_continuation.py
python research/research_max_struct_r7/experiments/vector_fractional_sphere_validation.py
python research/research_max_struct_r7/experiments/vector_plot_continuation.py
```

**Final3D scope:** actual125- and216-DOF vector material inversions are implemented and measured. Exact material-dependent A derivatives are supported. A did not establish durable recovery or basin superiority; longer physical optimization largely removes the early differences. Physical sphere accuracy is quantitatively bounded and geometry-sensitive. No continuum-qualified3D advantage is claimed.
'''
# Idempotent generation.
old=old.split('\n## Shared longer physical correction and necessary216-DOF case')[0]
p.write_text(old+append)
