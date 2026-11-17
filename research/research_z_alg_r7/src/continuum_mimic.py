"""Prescribed material in a disk and supported compensating material in an open annulus."""
import os
for key in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[key] = "2"

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import solve
from scipy.special import hankel1

from continuum_leakage import IN_R, OUT_R, TARGET_R, K, RECEIVER, ROOT, SOURCES, green, kernels, quadrature


def split_complex(values):
    return np.concatenate((values.real, values.imag), axis=-1)


def bump(points):
    s = (np.linalg.norm(points, axis=1) - IN_R) / (OUT_R - IN_R)
    values = np.zeros_like(s)
    interior = (s > 0) & (s < 1)
    values[interior] = np.exp(4 - 1 / (s[interior] * (1-s[interior])))
    return values


def make_modes(indices):
    xy, w = quadrature(IN_R, OUT_R, 36, 257, 0.18)
    gradients = split_complex(kernels(xy)[:, indices])
    gram = gradients.T @ ((w * bump(xy))[:, None] * gradients)
    eig, basis = np.linalg.eigh(gram)
    # Dual basis: its bare material-to-I/Q derivative is the identity.
    dual = (basis / eig) @ basis.T
    return dual, eig


def values_at(points, indices, dual):
    return bump(points)[:, None] * split_complex(kernels(points)[:, indices]) @ dual


def grid(ntar, nann, nangle, offset):
    tar, wt = quadrature(0, TARGET_R, ntar, nangle, offset)
    ann, wn = quadrature(IN_R, OUT_R, nann, nangle, offset)
    xy = np.concatenate((tar, ann))
    w = np.r_[wt, wn]
    distances = np.linalg.norm(xy[:, None] - xy[None, :], axis=2)
    np.fill_diagonal(distances, 1.)
    interaction = K**2 * 0.25j * hankel1(0, K*distances) * w[None, :]
    radius = np.sqrt(w / np.pi)
    np.fill_diagonal(interaction, 0.5j * np.pi*K*radius*hankel1(1,K*radius)-1)
    return xy, w, len(tar), interaction


def raw_operator(mesh, indices, dual):
    xy, w, n_tar, interaction = mesh
    modes = values_at(xy, indices, dual)
    incident = green(xy, np.vstack((SOURCES[indices], RECEIVER[None, :])))
    target = np.arange(len(xy)) < n_tar
    def evaluate(t, coefficients, derivatives=False):
        material = t * target + modes @ coefficients
        system = -interaction * material[None, :]
        system.flat[::len(system)+1] += 1
        total = solve(system, incident, assume_a="gen", check_finite=False, overwrite_a=True)
        y = ((K**2 * w * material * incident[:, -1]) @ total[:, :-1])
        if derivatives:
            product = K**2 * w[:,None] * total[:, -1:] * total[:, :-1]
            return y, split_complex(product).T @ modes
        return y
    return evaluate


def match(channels, indices, amplitude, meshes):
    dual, eig = make_modes(indices)
    coeff = np.zeros(2*len(indices))
    histories = []
    all_validations = []
    stage_coefficients = []
    for stage in (0,1):
        design = raw_operator(meshes[stage], indices, dual)
        history = []
        for _ in range(7):
            y, jac = design(amplitude, coeff, derivatives=True)
            history.append(np.linalg.norm(y))
            if history[-1] < 2e-15:
                break
            coeff -= np.linalg.solve(jac, split_complex(y))
        histories.append(history)
        stage_coefficients.append(coeff.tolist())
        validations = []
        # Keep coarse-design behavior, then evaluate the refined design without refitting.
        eval_indices = (0,1,2) if stage == 0 else (1,2,3)
        for idx in eval_indices:
            mesh = meshes[idx]
            evaluate = raw_operator(mesh, indices, dual)
            y = evaluate(amplitude, coeff)
            y0 = evaluate(amplitude, np.zeros_like(coeff))
            validations.append(dict(grid=idx, nodes=len(mesh[0]), raw_complex_l2=np.linalg.norm(y),
                                    target_only_l2=np.linalg.norm(y0), ratio=np.linalg.norm(y)/np.linalg.norm(y0),
                                    y=np.c_[y.real,y.imag].tolist(), y_target_only=np.c_[y0.real,y0.imag].tolist()))
        all_validations.append(validations)
    xy,w=quadrature(IN_R,OUT_R,46,379,.29)
    nuisance = values_at(xy,indices,dual) @ coeff
    return dict(channels=channels,primitive_indices=indices.tolist(),material_target=amplitude,
                nuisance_min=np.min(nuisance),nuisance_max=np.max(nuisance),nuisance_l2=np.sqrt(w @ nuisance**2),
                raw_gram_eigen_range=[eig[0],eig[-1]],raw_gram_condition=eig[-1]/eig[0],
                coefficients=coeff.tolist(),coarse_coefficients=stage_coefficients[0],newton_residuals=histories,
                coarse_validation=all_validations[0],validation=all_validations[1]),dual,coeff


def main():
    geometry=[(5,8,32,.11),(8,12,48,.32),(10,16,64,.47),(13,20,80,.19)]
    meshes=[grid(*v) for v in geometry]
    channel_sets={"two":np.array([0,4]),"four":np.array([0,3,6,9]),"twelve":np.arange(12)}
    rows=[]; raw={}
    for name,indices in channel_sets.items():
        for amplitude in ((0.025,0.10) if name=="two" else (0.10,)):
            row,dual,coeff=match(name,indices,amplitude,meshes)
            rows.append(row)
            raw[name+"_dual"]=dual
            raw[name+"_"+str(amplitude)+"_coeffs"]=coeff
    results=ROOT/"results"; figures=ROOT/"figures"
    np.savez_compressed(results/"continuum_mimic_coeffs.npz",**raw)
    (results/"continuum_mimic.json").write_text(json.dumps(dict(quadratures=geometry,rows=rows),indent=2)+"\n")
    picked=[r for r in rows if r["material_target"]==0.1]
    figure,axs=plt.subplots(1,3,figsize=(9.6,3.35),layout="constrained")
    nodes=np.linspace(-.86,.86,193);xx,yy=np.meshgrid(nodes,nodes)
    xy=np.c_[xx.ravel(), yy.ravel()]; rr=np.hypot(xy[:,0],xy[:,1])
    for ax,row in zip(axs,picked):
        name=row["channels"];indices=channel_sets[name]
        h=values_at(xy,indices,raw[name+"_dual"]) @ np.array(row["coefficients"])
        z=np.where(rr<TARGET_R,0.1,h)
        im=ax.imshow(z.reshape(xx.shape),origin="lower",extent=(nodes[0],nodes[-1],nodes[0],nodes[-1]),
                     cmap="RdBu_r",vmin=-.1,vmax=.1)
        ax.set(title=f'{len(indices)} coherent shots',xlabel="x / wavelength",ylabel="y / wavelength")
        for r in [TARGET_R,IN_R,OUT_R]: ax.add_patch(plt.Circle((0,0),r,fill=False,lw=.65,ls="--",color="#4c4c4c"))
    figure.colorbar(im,ax=axs,label="Contrast change",shrink=.85,extend="both")
    for ext in ("png","pdf"):figure.savefig(figures/f"continuum_mimic.{ext}",dpi=210)
    for row in rows:print(row["channels"],row["material_target"],row["raw_gram_condition"],row["nuisance_min"],row["nuisance_max"],[v["ratio"] for v in row["validation"]])


if __name__ == "__main__":
    main()
