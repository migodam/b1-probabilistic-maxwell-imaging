"""Draw paired development material and sparse port-error maps."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from training import HERE, ROOT


def draw(source, output):
    with np.load(source, allow_pickle=False) as result:
        m, nominal, centers = (result[key] for key in ("m", "nominal", "centers"))
        train, held = (result[key] for key in ("train_absolute_port_error", "held_absolute_port_error"))
        train_tx, held_tx = (result[key] for key in ("train_tx", "held_tx"))
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42})
    fig, axs = plt.subplots(2, 2, figsize=(8.0, 6.6), layout="constrained")
    side = int(np.sqrt(len(m)))
    step = centers[1, 0] - centers[0, 0]
    extent = tuple(100 * np.array([centers[:, 0].min()-step/2, centers[:, 0].max()+step/2,
                                   centers[:, 1].min()-step/2, centers[:, 1].max()+step/2]))
    peak = max(float(m.max()), float(nominal.max()))
    for ax, contrast, title in zip(axs[0], (m, nominal), ("Recovered contrast", "Nominal reference"), strict=True):
        material_image = ax.imshow(contrast.reshape(side, side), origin="lower", extent=extent,
                                  vmin=0, vmax=peak, cmap="viridis", interpolation="nearest")
        ax.set_title(title)
        ax.set_xlabel("x (cm)")
        ax.set_ylabel("y (cm)")
    fig.colorbar(material_image, ax=axs[0].tolist(), label="Contrast χ", fraction=.055)
    error_cap = float(np.quantile(np.r_[train.ravel(), held.ravel()][
        np.isfinite(np.r_[train.ravel(), held.ravel()])], .99))
    cmap = plt.colormaps["magma"].with_extremes(bad="white")
    for ax, values, sources, title in zip(axs[1], (train, held), (train_tx, held_tx),
                                         ("Training transmitter channels", "Held transmitter channels"), strict=True):
        field = ax.imshow(values, origin="lower", vmin=0, vmax=error_cap, cmap=cmap, aspect="auto",
                          extent=(-.5, len(sources)-.5, .5, values.shape[0]+.5), interpolation="nearest")
        ax.set_xticks(np.arange(len(sources)), [str(n) for n in sources])
        ax.set_xlabel("Transmitter ID")
        ax.set_ylabel("Receiver ID")
        ax.set_title(title)
    fig.colorbar(field, ax=axs[1].tolist(), extend="max", label="Absolute complex field error", fraction=.055)
    fig.suptitle("Fresnel Twin 8 GHz, sealed Dual-CIE development endpoint")
    for suffix in (".png", ".pdf"):
        fig.savefig(output.with_suffix(suffix), dpi=220, facecolor="white")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path,
                        default=ROOT / ".artifacts/z_alg_xudong_cie/chambers/dev_twin_terminal_local.npz")
    parser.add_argument("--output", type=Path, default=HERE / "dev_twin_terminal_map")
    args = parser.parse_args()
    draw(args.source, args.output)
