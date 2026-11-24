#!/usr/bin/env python3
"""Render the native-grid figure panels from the saved numerical arrays."""
import json
from pathlib import Path
from time import perf_counter, process_time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, PercentFormatter

HERE = Path(__file__).resolve().parent
EXTENT = (-1, 1, -1, 1)
BLUE = "#1B5E85"
ORANGE = "#AE562B"
INK = "#242E36"
FAINT = "#DAE0E4"
MODE_COLORS = ("#1B5E85", "#AE562B", "#337B62", "#754A86")
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.labelcolor": INK, "axes.edgecolor": "#67737A",
    "text.color": INK, "xtick.color": INK, "ytick.color": INK,
    "axes.linewidth": .65, "pdf.fonttype": 42, "savefig.facecolor": "white",
    "savefig.dpi": 350, "axes.spines.top": False, "axes.spines.right": False,
    "mathtext.fontset": "dejavusans",
})


def save(fig, name):
    for kind in ("pdf", "png"):
        fig.savefig(HERE / f"{name}.{kind}")
    plt.close(fig)


def image_style(ax, show_x_ticks, show_y_ticks, *, x_tick_labels=None):
    ax.set_xticks((-1, 0, 1), labels=x_tick_labels if show_x_ticks else ("", "", ""))
    ax.set_yticks((-1, 0, 1), labels=None if show_y_ticks else ("", "", ""))
    ax.tick_params(direction="out", length=2, pad=2)
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(.5)


def contrast_maps(objects, data):
    fig, axes = plt.subplots(2, 3, figsize=(7.1, 4.7))
    fig.subplots_adjust(left=.10, right=.88, bottom=.13, top=.93, hspace=.26, wspace=.11)
    images = []
    for i, (o, ax) in enumerate(zip(objects, axes.flat, strict=True)):
        images.append(ax.imshow(data["chi_" + o["index"]], cmap="cividis", vmin=0, vmax=4,
                                origin="lower", extent=EXTENT, interpolation="nearest"))
        ax.set_title(f'{o["index"]}  {o["name"]}', pad=5, fontsize=8.6)
        image_style(ax, i >= 3, i % 3 == 0)
    fig.supxlabel("$x$ (m)", y=.025, x=.49)
    fig.supylabel("$y$ (m)", x=.025)
    colorbar = fig.colorbar(images[0], cax=fig.add_axes([.892, .171, .017, .706]),
                            ticks=np.arange(5))
    colorbar.set_label(r"Dielectric contrast $\chi$", labelpad=6)
    save(fig, "native_contrast_maps")


def current_maps(data):
    fig, axes = plt.subplots(2, 4, figsize=(7.1, 3.88))
    fig.subplots_adjust(left=.11, right=.89, bottom=.13, top=.90, hspace=.15, wspace=.05)
    modes = data["mode_eigenvectors"][:, :4].reshape(96, 96, 4)
    high = np.abs(modes).max()
    colorbars = []
    for j in range(4):
        mode = modes[:, :, j]
        amp = np.abs(mode)
        phase = np.ma.masked_where(amp < .05 * amp.max(), np.angle(mode))
        a, b = axes[:, j]
        a.set_title(f'{j+1}  ({data["mode_eigenvalues"][j]:+.4f})', fontsize=8.5, pad=5)
        ami = a.imshow(amp, cmap="magma", vmin=0, vmax=high, origin="lower",
                       extent=EXTENT, interpolation="nearest")
        cyclic = matplotlib.colormaps["twilight"].with_extremes(bad="#EAEAEA")
        phi = b.imshow(phase, cmap=cyclic, vmin=-np.pi, vmax=np.pi, origin="lower",
                       extent=EXTENT, interpolation="nearest")
        image_style(a, False, j == 0)
        image_style(b, True, j == 0, x_tick_labels=(
            "$-1$" if j == 0 else "", "0", "1" if j == 3 else ""))
        if j == 0:
            colorbars.extend((ami, phi))
    fig.supylabel("$y$ (m)", x=.02)
    fig.supxlabel("$x$ (m)", x=.5, y=.032)
    ca, cb = [fig.add_axes(rect) for rect in ([.895, .568, .014, .301], [.895, .152, .014, .301])]
    cbar = fig.colorbar(colorbars[0], cax=ca, ticks=(0, .01, .02, .03))
    cbar.set_label(r"$|v|$", labelpad=3)
    cbar = fig.colorbar(colorbars[1], cax=cb, ticks=(-np.pi, 0, np.pi))
    cbar.ax.set_yticklabels((r"$-\pi$", "0", r"$\pi$"))
    cbar.set_label(r"$\arg v$", labelpad=3)
    save(fig, "negative_state_currents")


def source_spectra(objects):
    fig, axes = plt.subplots(2, 3, figsize=(7.1, 4.73), sharex=True, sharey=True)
    fig.subplots_adjust(left=.105, right=.982, bottom=.20, top=.92, hspace=.29, wspace=.18)
    for i, (o, ax) in enumerate(zip(objects, axes.flat, strict=True)):
        for noise, color, marker in ((1, BLUE, "o"), (2, ORANGE, "s")):
            row = next(s for s in o["source"] if s["noise_pct"] == noise)
            values = np.array(row["source_eigenvalues"])[4:] / row["padded_radius"]**2
            ax.plot(np.arange(5, 11), values, marker=marker, ms=3, lw=1.1,
                    color=color, label=f"{noise}% noise")
        ax.axhline(0, c="#71797D", lw=.8, zorder=0)
        ax.set_yscale("symlog", linthresh=1., linscale=.8)
        ax.set_ylim(-300, 12)
        ax.yaxis.set_major_locator(FixedLocator((-100, -10, -1, 0, 1, 10)))
        ax.set_xticks(np.arange(5, 11))
        ax.set_title(f'{o["index"]}  {o["name"]}', pad=5, fontsize=8.6)
        ax.grid(axis="y", color=FAINT, lw=.5, zorder=0)
        if i < 3:
            ax.tick_params(axis="x", length=2)
    fig.supylabel(r"Retained source eigenvalue / $r^2$", x=.015)
    fig.supxlabel("Retained source eigenvalue index", y=.096)
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="lower center",
               bbox_to_anchor=(.51, -.002), ncol=2, frameon=False)
    save(fig, "retained_source_spectra")


def fourier_bands(metrics):
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.95), sharex=True)
    fig.subplots_adjust(left=.10, right=.98, bottom=.27, top=.88, wspace=.33)
    rows = metrics["fourier"]["bands"]
    labels = [str(row["columns"]) for row in rows]
    x = np.arange(len(rows))
    raw = np.asarray([row["raw_fourier_energy"] for row in rows])
    projected = np.asarray([row["exact_ker_fourier_of_normalized_dark_energy"] for row in rows])
    ax, bx = axes
    for i, color in enumerate(MODE_COLORS):
        ax.plot(x, raw[:, i], linestyle="--", lw=.9, color=color, alpha=.8)
        ax.plot(x, projected[:, i], marker="o", ms=2.8, lw=1.2, color=color,
                label=f"Mode {i+1}")
    ax.set_title("(a) Energy captured", loc="left", pad=7, fontsize=8.8)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    ax.set_ylim(-.04, 1.06)
    angles = np.asarray([row["angles_dark_mode_space_to_projected_fourier_deg"] for row in rows])
    bx.fill_between(x, angles[:, 0], angles[:, -1], color=BLUE, alpha=.15,
                    label="Smallest to largest")
    bx.plot(x, angles[:, 0], color=BLUE, lw=1.05)
    bx.plot(x, angles[:, -1], color=BLUE, lw=1.05)
    bx.set_ylim(-3, 97)
    bx.set_yticks((0, 30, 60, 90))
    bx.set_ylabel("Principal angle (degrees)")
    bx.set_title("(b) Angle between projected spaces", loc="left", pad=7, fontsize=8.8)
    for ax in axes:
        ax.set_xticks(x, labels=labels)
        ax.set_xlabel("Fourier columns")
        ax.grid(axis="y", color=FAINT, lw=.5)
    fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center",
               bbox_to_anchor=(.51, .008), ncol=4, frameon=False)
    save(fig, "fourier_dark_subspaces")


def main():
    wall, cpu = perf_counter(), process_time()
    metrics = json.loads((HERE / "results.json").read_text())
    with np.load(HERE / "plot_data.npz", allow_pickle=False) as data:
        contrast_maps(metrics["objects"], data)
        current_maps(data)
    source_spectra(metrics["objects"])
    fourier_bands(metrics)
    timing = dict(wall_seconds=perf_counter()-wall, cpu_seconds=process_time()-cpu)
    (HERE / "plot_timing.json").write_text(json.dumps(timing, indent=2) + "\n")
    print(f"Plots: {timing['wall_seconds']:.3f} s wall, {timing['cpu_seconds']:.3f} s CPU")


if __name__ == "__main__":
    main()
