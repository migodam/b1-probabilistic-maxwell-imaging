# Reproduce the fixed-operator closed-loop study

Use the repository root or extract the reproduction ZIP into an empty directory and work from that directory. Both layouts preserve the original `research/` paths. The ZIP contains the synthetic acquisitions and offline scoring arrays, proposed actions, numerical and interval receipts, and manuscript sources and figures. `RESEARCH_ARCHIVE_MANIFEST.json` at the archive root gives each inherited file’s original path and role. The measured-port section includes the processing code and derived scalars. The [run manifest](RUN_MANIFEST.json) indexes each cohort and its original result arrays. Facility downloads are needed only for a fresh measured-port calculation.

## Setup and original-scale replay

The study used Python 3.13, with dependency versions in [requirements.txt](requirements.txt). A fresh environment can be prepared from either root.

```sh
P=research/z_alg_optical_inertia/ip_closed_loop
python3.13 -m venv .repro-venv
. .repro-venv/bin/activate
python -m pip install -r "$P/requirements.txt"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B "$P/maps/replay.py" 96 ring M96-R01 broad
```

The last command reconstructs the trace and source–receiver blocks in Arb from the stored native 96² Fourier symbol. It reads the original weak-background observation, one saved candidate and its weights, then computes the strict margin anew. Look in `maps/results/arb_M96-R01_broad.json` for the resulting receipt. This is one candidate from the original wide-box map. Use `--interior` to check a complete in-domain footprint instead of the clipped one. These commands write receipts to the current tree, so use an extracted copy for independent work.

The isolated archive replay returned a strict margin of 1.0433876168 at centre `(0, 0)` with 120 footprint cells and two negative state directions. Two runs on a single-threaded Apple M4 took 5.66–6.55 s wall and 5.56–6.45 s CPU. The extracted code also loaded the native 16-transmitter, 16-receiver beam geometry and both laboratory readers from paths inside the archive. These checks used the installed environment with no network requests.

## Synthetic results

- For the 32² and native 96² full-field runs, follow the commands in [maps/README.md](maps/README.md). Public observations are already present. `maps/prepare.py` regenerates both data and offline truth. The online map and interval replay operate on the public data. Scoring, witness validation, and model transfer have their own entry points and receipts.
- The frozen coherent-beam study is described in [beams/README_ZH.md](beams/README_ZH.md). Its inherited operator files are listed in [beams/dependencies.txt](beams/dependencies.txt). The finite-angle and paired-energy comparison below reads saved old and future observations. The forecast-only form independently verifies one native case using only its saved old reading and proposal.
- [adaptive/experiment.md](adaptive/experiment.md) describes the separate spatial feedback experiment. The saved updates and witness counterchecks can be regenerated with the `read` and `updated_witnesses/recheck.py` entry points. Starting new synthetic acquisitions requires a fresh project copy and preserves the design-before-acquisition order described in each experiment.
- [reduction/README.md](reduction/README.md) covers the full 100/1,000-box direct comparison, its exact-singular regression and the matched controls. The [map PSD documentation](reduction/psd_maps/README.md) covers the continuous control and the exact receiver-null trials. Dense original comparator runs take substantially longer than a single replay.

```sh
python -B "$P/beams/paired_replay.py" power_96 401 --forecast-only --suffix _archive
python -B "$P/beams/paired_replay.py" power_32
python -B "$P/beams/paired_replay.py" power_96
python -B "$P/beams/angle_bound.py" power_96
python -B "$P/beams/secant_remainder_check.py"
python -B "$P/maps/nominal_control.py"
python -B "$P/beams/sparse_ports.py" 96
python -B "$P/beams/sparse_physical.py" 96
python -B "$P/beams/summarize.py"
python -B "$P/adaptive/loop.py" read 96
python -B "$P/adaptive/updated_witnesses/recheck.py" 96 --fit
python -B "$P/reduction/pivot_regression.py"
python -B "$P/reduction/psd_maps/dark_cone.py" --replay-full32
python -B "$P/reduction/summarize.py"
```

To rebuild the main figures and paper from included evidence, use the following commands. The LaTeX builder uses `pdflatex` and `tectonic` when available. It generates a separate anonymous source ZIP on demand from the included manuscript sources.

```sh
python -B "$P/src/figures.py"
python -B "$P/src/phase_figure.py"
python -B "$P/src/package_manuscript.py"
```

## Optional laboratory replay

The reported held-receiver ratios can be read directly in [lab/held_ports.csv](lab/held_ports.csv) and [lab/results.json](lab/results.json). The original facility readings and caches containing complex laboratory ports remain local to the reader. Use the facilities’ original usage conditions and retain the filenames below.

1. Download EIGOR Tomographic Dataset 1 from the University of Manitoba. Place the ZIP, still named `Tomographic_Dataset_1.zip`, in `.artifacts/z_alg_validation/manitoba/` at the extracted or repository root. [Facility data page](https://home.cc.umanitoba.ca/~lovetrij/EMILab/mwt_data.html) · [Original Dataset 1 ZIP](https://home.cc.umanitoba.ca/~lovetrij/EMILab/research/EIGOR/Tomographic_Dataset_1.zip)
2. Open the Fresnel 2005 article’s supplementary-data page and obtain `supp.zip` under the publisher’s terms. The page provides the current archive link. Extract `FoamDielIntTM.exp`, `FoamDielExtTM.exp`, and `FoamTwinDielTM.exp` to `research/delegated/a2_literature/data/2005_iop_21_6_S09/`. [Publisher’s supplementary data](https://beta.iopscience.iop.org/article/10.1088/0266-5611/21/6/S09/data) · [Original article](https://doi.org/10.1088/0266-5611/21/6/S09)

Prepare the nine Fresnel reading caches and rerun the held-receiver comparison from the root. The original loader supplies the measured scattered and empty fields with their actual missing channels. The matching source coordinates are 0° and 180°. `lab/measured.py` reads the Manitoba ZIP directly and applies the four reported reference combinations. The companion Manitoba calibration code is also included, though the held-receiver comparison uses its loader and analytic metal response directly.

```sh
python -B - <<'PY'
import numpy as np
from research.z_alg_optical_inertia.fresnel import physical

physical.RESULTS.mkdir(parents=True, exist_ok=True)
for frequency in physical.FREQUENCIES:
    for name in physical.NAMES:
        rows, scattered, empty = physical.load_file(physical.RAW / f"{name}.exp", frequency)
        sources = np.unique(rows[:, 0]).astype(int) - 1
        count = physical.SOURCES[name]
        if not np.array_equal(sources, np.arange(count)):
            raise ValueError(f"Unexpected recorded source indices for {name}")
        field = scattered[:, sources]
        np.savez_compressed(
            physical.RESULTS / f"{name}_{frequency}GHz_ports.npz",
            raw_scattered=field,
            raw_empty=empty[:, sources],
            mask=np.isfinite(field),
            source_degrees=np.arange(count) * (360 // count),
        )
PY
python -B "$P/lab/measured.py"
```

The original Fresnel data-access page uses publisher-managed download links. The synthetic and manuscript programs require no access to either facility.

## Build the reproduction archive

Run from the original project tree or from an extracted copy. The archive contains relative paths, with the case code and evidence under `research/z_alg_optical_inertia/ip_closed_loop/` and the inherited files beside it in their original locations.

```sh
python -B "$P/src/package_research.py" --inventory
python -B "$P/src/package_research.py"
```

The additional full-field hierarchical comparison is documented in [adaptive/multiscale/README_ZH.md](adaptive/multiscale/README_ZH.md). It includes a fresh original-data leaf sweep for every adaptive branch and every flat control at both resolutions. New fields on returned centres were verified against two separate noise balls in [adaptive/updated_witnesses/README_ZH.md](adaptive/updated_witnesses/README_ZH.md).

The independent [two-parameter full-wave template control](maps/template_control/summary.json) ranks candidates using the original data, saves material choices, replays them at the original resolution and then scores the centres. Follow [its commands](maps/template_control/README.md).
