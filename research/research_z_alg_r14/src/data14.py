"""Public legacy Maxwell observations, physical actions, and noise contracts."""
from dataclasses import dataclass
from functools import lru_cache
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
R12 = ROOT.parent / "research_z_alg_r12"
R10 = ROOT.parent / "research_z_alg_r10"
sys.path.insert(0, str(R12 / "src"))
from freq_model import BOX, Fit, FrequencyTM, ri

SIGMA = 3e-7
GAIN_SD = 4e-4
R10_OLD_EPS = 2e-4
R10_OLD_GAMMA = .01


def cn_radius(sigma, n_complex):
    return float(sigma * np.sqrt(chi2.ppf(.99, 2 * n_complex) / 2)) if n_complex else 0.


def complex_reads(values):
    values = np.asarray(values, dtype=float)
    return values[..., 0] + 1j * values[..., 1]


def encode_complex(values):
    values = np.asarray(values)
    return np.stack((values.real, values.imag), axis=-1).tolist()


@dataclass(frozen=True)
class OldCase:
    case_id: str
    cohort: str
    geometry: dict
    old: np.ndarray
    sigma: float
    gain_sd: float
    additive_radius: float
    gain_radius: float
    box: tuple = BOX
    hard_gain_radius: float | None = None

    def new_noise_radius(self, action):
        return cn_radius(self.sigma, self.old.size * action.exposures)


@dataclass(frozen=True)
class Action:
    name: str
    deltas: tuple
    energies: tuple

    @property
    def energy(self):
        return sum(self.energies)

    @property
    def exposures(self):
        return len(self.deltas)

    @property
    def bandwidth(self):
        """End-to-end relative frequency span, including the existing old block."""
        return max((0., *self.deltas)) - min((0., *self.deltas))

    @property
    def max_offset(self):
        return max((0., *map(abs, self.deltas)))


OLD_ONLY = Action("old_only", (), ())
REPEAT = Action("repeat", (0., 0.), (.5, .5))
HALF_SPLIT_FIRST = Action("half_split_first", (-.02,), (.5,))
BANDS = (.00025, .0005, .001, .002, .004, .005, .008, .012, .02)
MENU = (REPEAT,) + tuple(
    action for band in BANDS
    for action in (
        Action(f"one_minus_{band:g}", (-band,), (1.,)),
        Action(f"one_plus_{band:g}", (band,), (1.,)),
        Action(f"split_{band:g}", (-band, band), (.5, .5)),
    )
) + (
    Action("uniform_0.02", (-.02, -.02 / 3, .02 / 3, .02), (.25,) * 4),
)
WIDER = (
    Action("split_0.04", (-.04, .04), (.5, .5)),
    Action("split_0.08", (-.08, .08), (.5, .5)),
    Action("uniform_0.08", (-.08, -.08 / 3, .08 / 3, .08), (.25,) * 4),
)
ACTIONS = {action.name: action for action in (OLD_ONLY, *MENU, *WIDER, HALF_SPLIT_FIRST)}


def measured_reads(archive, case_id, action):
    if action == HALF_SPLIT_FIRST:
        return archive[f"{case_id}/split_0.02"][:1]
    return archive[f"{case_id}/{action.name}"]


def joint_fit(case, model, action, new):
    return Fit(model, (0., *action.deltas), (case.old, *new),
               weights=np.sqrt((1., *action.energies)),
               gain_calibration=((case.sigma / case.gain_sd) ** 2, 1.))


@lru_cache(maxsize=1)
def load_cases():
    """Give algorithms observed data and known geometry only, under opaque case IDs."""
    r12_geometry = {row["label"]: row["geometry"] for row in json.loads(
        (R12 / "results/spectral_new_layouts.json").read_text())}
    cohort = []
    for i, row in enumerate(json.loads((R12 / "results/spectral_online.json").read_text())):
        old = complex_reads(row["old_ri"])
        cohort.append(OldCase(
            f"spectral-{i:02d}", "spectral", r12_geometry[row["label"]], old,
            SIGMA, GAIN_SD, np.sqrt(2.) * cn_radius(SIGMA, old.size),
            np.sqrt(2.) * GAIN_SD, hard_gain_radius=GAIN_SD))
    r10_geometry = json.loads((R10 / "results/branch_holdout.json").read_text())["cases"]
    aperture_geometry = json.loads((R10 / "results/confirm_acquisition.json").read_text())[
        "settings"]["geometry"]
    for group, file in (("natural", "natural_inputs.npz"), ("aperture", "confirm_inputs.npz")):
        with np.load(R10 / "results" / file) as table:
            for name in sorted(table["case_ids"]):
                old = table[f"{name}/y_old"]
                geometry = r10_geometry[name]["geometry"] if group == "natural" else aperture_geometry
                # R10 used a fixed fraction of ε. Treat ε as r99 for new Gaussian blocks.
                sigma = R10_OLD_EPS / cn_radius(1., old.size)
                cohort.append(OldCase(
                    str(name), group, geometry, old, sigma, R10_OLD_GAMMA, R10_OLD_EPS, R10_OLD_GAMMA))
    for case in cohort:
        case.old.flags.writeable = False
    return tuple(cohort)
