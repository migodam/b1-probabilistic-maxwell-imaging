"""Physical frequency and MIMO port schedules; one unit of total transmit energy."""
from dataclasses import dataclass
from itertools import product

import numpy as np


FREQUENCY_SETS = {
    1: (16.,), 2: (2., 30.), 4: (2., 12., 20., 30.),
    8: (2., 6., 10., 16., 20., 24., 28., 30.),
}
SIGMA_FRACTION = .02


@dataclass(frozen=True)
class Design:
    name: str
    freqs: tuple
    tx: tuple
    rx: tuple

    @property
    def shots(self):
        return len(self.freqs) * len(self.tx)

    @property
    def samples(self):
        return self.shots * len(self.rx)

    @property
    def shot_energy(self):
        return 1. / self.shots


def spaced(total, number):
    return tuple(range(0, total, total // number))


def make(nf, nt, nr):
    return Design(f"f{nf}_t{nt}_r{nr}", FREQUENCY_SETS[nf], spaced(8, nt), spaced(16, nr))


def all_designs():
    arms = {d.name: d for d in (make(f,t,r) for f,t,r in
            product(FREQUENCY_SETS, (1,2,4,8), (1,2,4,8,16)))}
    for tag, freqs in (("narrow", (15.68,16.32)), ("high", (30.,)), ("low", (2.,))):
        for nt,nr in ((1,1), (4,8), (8,8), (8,16)):
            name = f"{tag}_t{nt}_r{nr}"
            arms[name] = Design(name, freqs, spaced(8,nt), spaced(16,nr))
    for nf in (1, 2, 4, 8):
        tag = f"oneside_f{nf}_t4_r8"
        arms[tag] = Design(tag, FREQUENCY_SETS[nf], (6,7,0,1), (12,13,14,15,0,1,2,3))
    return arms


RECOVERY_ARMS = (
    "f1_t8_r8", "f2_t4_r8", "f4_t2_r8", "f8_t1_r8",
    "narrow_t4_r8", "high_t8_r8", "f4_t8_r16", "f4_t2_r2", "f4_t2_r16",
    "oneside_f2_t4_r8",
)


def read_arm(pack, case_index, design, *, same_grid=False, sigma_fraction=SIGMA_FRACTION):
    """Return frequency by receiver by transmitter. The solver sees unit-amplitude data."""
    indices = [int(np.flatnonzero(np.isclose(pack["frequencies_ghz"], f, atol=1.e-8))[0])
               for f in design.freqs]
    selection = np.ix_(indices, design.rx, design.tx)
    clean = pack["clean_same" if same_grid else "clean"][case_index][selection]
    noise = pack["noise"][case_index][selection]
    sigma = float(pack["reference_amplitude"]) * sigma_fraction / np.sqrt(design.shot_energy)
    return clean + sigma * noise, sigma
