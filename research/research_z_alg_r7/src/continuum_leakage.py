"""Open-region sensitivity of fixed coherent Helmholtz illuminations, wavelength units."""
from pathlib import Path
import csv
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from scipy.linalg import null_space
from scipy.special import hankel1, j0

ROOT = Path(__file__).resolve().parents[1]
K = 2 * np.pi
TARGET_R = 0.24
IN_R, OUT_R = 0.32, 0.84
SOURCE_ANGLES = (np.arange(12) + 0.13) * 2 * np.pi / 12
RECEIVER = np.array([2.13, 0.71])


def ring(radius, theta):
    return radius * np.stack([np.cos(theta), np.sin(theta)], axis=-1)


SOURCES = ring(2.2, SOURCE_ANGLES)


def quadrature(inner, outer, nr, nt, offset):
    nodes, weights = np.polynomial.legendre.leggauss(nr)
    rr = (outer - inner) * nodes / 2 + (outer + inner) / 2
    theta = (np.arange(nt) + offset) * 2 * np.pi / nt
    xy = (rr[:, None, None] * ring(1, theta)[None, :, :]).reshape(-1, 2)
    w = np.repeat((outer - inner) / 2 * weights * rr * 2 * np.pi / nt, nt)
    return xy, w


def green(points, locations):
    distances = np.linalg.norm(points[:, None, :] - locations[None, :, :], axis=2)
    return 0.25j * hankel1(0, K * distances)


def kernels(points):
    return K**2 * green(points, RECEIVER[None, :]) * green(points, SOURCES)


def operators(nr, nt, offset):
    target_points, target_weights = quadrature(0, TARGET_R, nr, nt, offset)
    nuisance_points, nuisance_weights = quadrature(IN_R, OUT_R, nr, nt, offset)
    c = target_weights @ kernels(target_points)
    g = kernels(nuisance_points)
    s = g.conj().T @ (nuisance_weights[:, None] * g)
    return c, s, g, nuisance_weights


def optimal_feed(c, s):
    unscaled = np.linalg.solve(s, c.conj())
    return unscaled / (c @ unscaled)


def null_feed(c, nodes):
    basis = null_space(kernels(nodes))
    attainable = c @ basis
    return basis @ attainable.conj() / np.vdot(attainable, attainable).real


def nonlinear_check(feed):
    """Centered change in a volume-integral solve for the optimum beam's real nuisance."""
    xy, weights = quadrature(IN_R, OUT_R, 10, 48, 0.37)
    incident = green(xy, SOURCES) @ feed
    receiver = green(xy, RECEIVER[None, :]).ravel()
    g = K**2 * receiver * incident
    phase = 0.5 * np.angle(weights @ (g**2))
    h = (g * np.exp(-1j * phase)).real
    h /= np.sqrt(weights @ h**2)
    analytic = weights @ (g * h)
    distances = np.linalg.norm(xy[:, None] - xy[None, :], axis=2)
    np.fill_diagonal(distances, 1.0)
    interaction = K**2 * 0.25j * hankel1(0, K * distances) * weights[None, :]
    # Integrate the logarithmic diagonal over a disk with the same cell area.
    radius = np.sqrt(weights / np.pi)
    np.fill_diagonal(interaction, 0.5j * np.pi*K*radius*hankel1(1, K*radius) - 1)
    identity = np.eye(len(xy))
    def scattered(amplitude):
        field = np.linalg.solve(identity - amplitude * interaction * h[None, :], incident)
        return K**2 * amplitude * weights @ (receiver * field * h)
    values = []
    for step in (0.02, 0.01, 0.005):
        plus = scattered(step)
        minus = scattered(-step)
        derivative = (plus - minus)/(2*step)
        values.append(dict(step=step, plus_scattered=[plus.real, plus.imag],
                           minus_scattered=[minus.real, minus.imag],
                           centered_abs=abs(derivative), relative_complex_error=abs(derivative/analytic-1)))
    return dict(quadrature=[10,48,0.37], analytic_derivative=[analytic.real,analytic.imag],
                nuisance_L2=np.sqrt(weights @ h**2), nuisance_Linf=np.max(abs(h)), values=values)


def main():
    results = ROOT / "results"
    figures = ROOT / "figures"
    results.mkdir(exist_ok=True)
    figures.mkdir(exist_ok=True)
    c, s, _, _ = operators(12, 80, 0.1)
    c_check, s_check, g_check, w = operators(34, 257, 0.43)
    radiation = j0(K * np.linalg.norm(SOURCES[:, None] - SOURCES[None, :], axis=2)) / 4
    probes = {"minimum feed effort": c.conj() / np.vdot(c, c).real}
    node_sets = {}
    for count in (4, 8, 11):
        nodes = ring(0.54, (np.arange(count) + 0.27) * 2 * np.pi / count)
        name = f"{count} point nulls"
        probes[name] = null_feed(c, nodes)
        node_sets[name] = nodes
    probes["minimum area leakage"] = optimal_feed(c, s)
    rows = []
    for name, feed in probes.items():
        gradient = g_check @ feed
        energy = np.sum(w * abs(gradient)**2)
        square_integral = w @ (gradient**2)
        phase = 0.5 * np.angle(square_integral)
        perturbation = (np.exp(-1j * phase) * gradient).real
        perturbation /= np.sqrt(w @ perturbation**2)
        real_worst = np.sqrt((energy + abs(square_integral)) / 2)
        realized = abs(w @ (gradient * perturbation))
        nodes = node_sets.get(name)
        rows.append(dict(illumination=name, target_derivative_abs=abs(c_check @ feed),
                         target_complex_error=abs(c_check @ feed - 1),
                         max_node_derivative_abs=None if nodes is None else np.max(abs(kernels(nodes) @ feed)),
                         complex_leakage_L2=np.sqrt(energy), real_worst_L2=real_worst,
                         real_witness_abs=realized, feed_norm_squared=np.vdot(feed, feed).real,
                         radiated_flux=np.vdot(feed, radiation @ feed).real))
    with (results / "continuum_leakage.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    mineig = np.linalg.eigvalsh(s_check)
    min_complex = np.sqrt(1 / (c_check @ np.linalg.solve(s_check, c_check.conj())).real)
    summary = dict(k=K, wavelength=1, target_radius=TARGET_R, nuisance_annulus=[IN_R, OUT_R],
                   sources=SOURCES.tolist(), receiver=RECEIVER.tolist(), collocation_radius=0.54,
                   design_quadrature=[12, 80, 0.1], verification_quadrature=[34, 257, 0.43],
                   gram_eigenvalue_range=[mineig[0], mineig[-1]], gram_condition=mineig[-1]/mineig[0],
                   relative_quadrature_change_S=np.linalg.norm(s-s_check)/np.linalg.norm(s_check),
                   relative_quadrature_change_c=np.linalg.norm(c-c_check)/np.linalg.norm(c_check),
                   independent_complex_lower_bound=min_complex, real_uniform_lower_bound=min_complex / np.sqrt(2),
                   max_real_witness_error=max(abs(row["real_worst_L2"]-row["real_witness_abs"]) for row in rows))
    summary["nonlinear_check_min_area"] = nonlinear_check(probes["minimum area leakage"])
    (results / "continuum_leakage.json").write_text(json.dumps(summary, indent=2)+"\n")
    np.savez_compressed(results / "continuum_feeds.npz", names=np.array(list(probes)), feeds=np.stack(list(probes.values())),
                        c=c_check, S=s_check, radiation_gram=radiation, **{f"nodes_{len(v)}":v for v in node_sets.values()})

    names = ["8 point nulls", "minimum area leakage"]
    grid = np.linspace(-OUT_R, OUT_R, 201)
    xx, yy = np.meshgrid(grid, grid)
    radius = np.hypot(xx, yy)
    fields = kernels(np.c_[xx.ravel(),yy.ravel()])
    maps = [np.where(radius <= OUT_R, np.abs(fields @ probes[name]).reshape(xx.shape), np.nan) for name in names]
    vmax = np.nanmax(maps)
    fig, axes = plt.subplots(1, 3, figsize=(10.9, 3.5), layout="constrained", gridspec_kw={"width_ratios":[1, 1, 1.2]})
    for ax, img, name, label in zip(axes, maps, names, ("A", "B")):
        im = ax.imshow(img, extent=[grid[0], grid[-1], grid[0], grid[-1]], norm=LogNorm(vmin=0.1, vmax=vmax), cmap="magma", origin="lower")
        for rr, color, style in [(TARGET_R, "white","-"), (IN_R, "lightblue","--"), (OUT_R, "lightblue","--")]:
            ax.add_patch(plt.Circle((0,0),rr,fill=False,color=color,lw=1,ls=style))
        if name in node_sets:
            n = node_sets[name]; ax.scatter(n[:,0],n[:,1],marker="+",s=40,c="cyan",lw=1.25)
        ax.set(title=f"{label}   {name.capitalize()}", xlabel="x / wavelength", ylabel="y / wavelength")
        ax.set_aspect("equal")
    fig.colorbar(im, ax=axes[:2], shrink=.78, label="Local sensitivity |g|")
    ax=axes[2]; labels=["Feed\nminimum","4\nzeros","8\nzeros","11\nzeros","Area\nminimum"]
    indices=np.arange(len(rows))
    ax.bar(indices-.17,[r["complex_leakage_L2"] for r in rows],width=.34,label="Complex L² norm",color="#125ca0")
    ax.bar(indices+.17,[r["real_worst_L2"] for r in rows],width=.34,label="Real worst case",color="#d88928")
    ax.set(title="C   Open-annulus leakage",ylabel="Response per unit nuisance L²",xticks=indices,xticklabels=labels)
    ax.set_yscale("log")
    ax.legend(frameon=False,fontsize=8)
    for ext in ("png","pdf"):
        fig.savefig(figures / f"continuum_leakage.{ext}",dpi=210)
    print(json.dumps(summary,indent=2))
    for r in rows: print(r)


if __name__ == "__main__":
    main()
