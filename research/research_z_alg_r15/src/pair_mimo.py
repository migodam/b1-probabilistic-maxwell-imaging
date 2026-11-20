"""Known-gain, two-hypothesis contrast of frequency, aperture and Tx waveforms."""

import json
from functools import lru_cache
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from scipy.special import j0, log_ndtr
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
R14 = ROOT.parent / "research_z_alg_r14"
R10 = ROOT.parent / "research_z_alg_r10"
WITNESS = R14 / "results/continuous_band/fresh6-16_known_coherent.json"
GEOMETRY_SOURCE = R10 / "results/branch_holdout.json"
sys.path.insert(0, str(ROOT.parent / "research_z_alg_r2" / "src"))
from physics import Model

C0 = 299792458
F0 = 16.
ENERGY = 3.
PORT_EFFORT_CAP = 2.
TX_SLOTS_PER_TONE = 6
NARROW = np.linspace(15.68, 16.32, 8)
WIDE = np.linspace(2, 30, 8)
OLD_RADIUS = 2e-4


def radius(sigma, samples):
    return sigma * np.sqrt(chi2.ppf(.99, 2 * samples) / 2)


SIGMA = OLD_RADIUS / radius(1., 3)


def encode(x):
    x = np.asarray(x)
    return np.stack((x.real, x.imag), axis=-1).tolist()


def complex_norm(y):
    return float(np.linalg.norm(y))


def layout(original):
    old_sources = np.array(original["sources"])
    old_receiver = np.array(original["receivers"])
    theta_0 = np.angle(complex(*old_receiver[0]))
    # Offset the six transmitters by half a sector from the existing receiver.
    theta_tx = theta_0 + np.pi / 6 + np.arange(6) * np.pi / 3
    theta_rx = theta_0 + np.arange(1, 8) * np.pi / 4
    tx_ring = np.mean(np.linalg.norm(old_sources, axis=1)) * np.c_[np.cos(theta_tx), np.sin(theta_tx)]
    # Seven ports on an outer ring plus the existing port in its original position.
    rx_ring = 1.6 * np.c_[np.cos(theta_rx), np.sin(theta_rx)]
    return {"old": old_sources, "uniform6": tx_ring}, np.vstack((old_receiver, rx_ring))


class Comparison:
    def __init__(self):
        self.witness = json.loads(WITNESS.read_text())
        self.materials = np.array(self.witness["materials"])
        self.original = json.loads(GEOMETRY_SOURCE.read_text())["cases"][self.witness["case_id"]]["geometry"]
        self.sources, self.receivers = layout(self.original)
        self.k0_si = 2 * np.pi * F0 * 1e9 / C0
        self.length_scale_m = self.original["k"] / self.k0_si

    @lru_cache(maxsize=None)
    def at(self, ghz, tx):
        sources = self.sources[tx] * self.length_scale_m
        k = 2 * np.pi * ghz * 1e9 / C0
        model = Model(
            centers=np.array(self.original["centers"]) * self.length_scale_m,
            sources=sources,
            receivers=self.receivers * self.length_scale_m,
            radius=self.original["radius"] * self.length_scale_m,
            k=k,
        )
        states = [model.forward(m, independent=True) for m in self.materials]
        fields = np.array([state["y"] for state in states])
        power_gram = j0(k * np.linalg.norm(sources[:, None] - sources[None, :], axis=2))
        condition = [state["resolvent_sigma_min"] for state in states]
        return fields * np.sqrt(ghz / F0), power_gram, condition

    def measure(self, tones, tx, rx, waveform="fixed", switch_phase=0):
        ns = len(self.sources[tx])
        num_shots = len(tones) * TX_SLOTS_PER_TONE
        pulse_energy = ENERGY / num_shots
        fields, rx_indices, frequencies, weights = [], [], [], []
        for tone_index, ghz in enumerate(tones):
            y, q, _ = self.at(float(ghz), tx)
            for slot in range(TX_SLOTS_PER_TONE):
                receivers = selected_receivers(rx, tone_index, slot, switch_phase)
                if waveform == "fixed":
                    coeff = np.eye(ns, dtype=complex)[slot % ns]
                else:
                    cap = PORT_EFFORT_CAP if waveform == "capped" else None
                    coeff, _ = best_beam(y[:, receivers], q, cap)
                coeff = np.sqrt(pulse_energy) * coeff
                fields.append(y[:, receivers] @ coeff)
                rx_indices.extend(receivers)
                weights.append(coeff)
                frequencies.append(float(ghz))
        fields = np.concatenate(fields, axis=1)
        return measurement_result(fields, rx_indices, frequencies, np.array(weights), self, tx)


def selected_receivers(rx, tone_index, slot, phase=0):
    if rx == "simultaneous8":
        return np.arange(8)
    if rx == "switched8":
        return np.array([(tone_index + slot + phase) % 8])
    return np.array([0])


def best_beam(y, q, cap=None):
    delta = y[0] - y[1]
    objective = delta.conj().T @ delta
    ns = len(q)
    identity = np.eye(ns)

    def top(penalty):
        _, v = eigh(objective - penalty * identity, q, subset_by_index=(ns - 1, ns - 1))
        return v[:, 0]

    coeff = top(0.)
    if cap is not None and np.vdot(coeff, coeff).real > cap:
        def slack(penalty):
            v = top(penalty)
            return np.vdot(v, v).real - cap

        hi = np.linalg.norm(objective, ord=2)
        while slack(hi) > 0:
            hi *= 2
        coeff = top(brentq(slack, 0., hi, xtol=hi * 1e-12))
    pivot = np.argmax(np.abs(coeff))
    coeff *= np.exp(-1j * np.angle(coeff[pivot]))
    return coeff, complex_norm(delta @ coeff)**2


def log_pe(gap, sigma=SIGMA):
    return float(log_ndtr(-gap / (np.sqrt(2) * sigma)) / np.log(10))


def sensitivity(mu, receiver_index, gap):
    # Static per-receiver multiplicative mismatch, constant across the tones and shots.
    # Relax the rectangular magnitude/phase set to its enclosing complex disk.
    delta = mu[0] - mu[1]
    groups = np.unique(receiver_index)
    projected_fields = np.array([
        [np.vdot(delta[receiver_index == r], m[receiver_index == r]) for r in groups]
        for m in mu
    ])
    projection_budget = np.sum(np.abs(projected_fields), axis=1)
    projection_budget_varying = np.sum(np.abs(delta.conj() * mu), axis=1)
    noise = [
        dict(sigma_multiplier=multiplier, log10_pe=log_pe(gap, multiplier * SIGMA))
        for multiplier in (1, 3, 10, 30, 100)
    ]
    calibration = []
    for amplitude, phase_deg in ((.001, .1), (.01, 1.), (.05, 5.)):
        eta = abs((1 + amplitude) * np.exp(1j * np.deg2rad(phase_deg)) - 1)
        row = dict(
            amplitude_fraction=amplitude, phase_degrees=phase_deg, enclosing_complex_disk=eta,
        )
        for mode, budget in (("static", projection_budget), ("varying", projection_budget_varying)):
            remaining = .5 * gap**2 - eta * budget
            row[mode] = dict(remaining_projection_fraction=(2 * remaining / gap**2).tolist())
            for multiple in (1, 10):
                conditional_log10 = log_ndtr(-remaining / (multiple * SIGMA * gap / np.sqrt(2))) / np.log(10)
                row[mode][f"sigma{multiple}_equal_prior_log10_pe_upper"] = float(
                    np.logaddexp(*(conditional_log10 * np.log(10))) / np.log(10) - np.log10(2)
                )
        calibration.append(row)
    return dict(noise=noise, calibration=calibration, projection_budget=projection_budget.tolist(),
                projection_budget_varying=projection_budget_varying.tolist())


def measurement_result(mu, receiver_index, frequencies, weights, experiment, tx):
    receiver_index = np.asarray(receiver_index)
    gap = complex_norm(mu[0] - mu[1])
    n = len(receiver_index)
    per_shot_energy = [float(np.vdot(w, experiment.at(f, tx)[1] @ w).real)
                       for f, w in zip(frequencies, weights, strict=True)]
    active = abs(weights) > 1e-12
    effort = np.sum(abs(weights)**2, axis=1)
    result = dict(
        tx=tx, transmitter_count=len(experiment.sources[tx]),
        distinct_receiver_count=int(np.unique(receiver_index).size),
        rx_samples=n, shots=len(weights), powered_tx_slots=int(active.sum()),
        transmit_count_per_port=active.sum(axis=0).tolist(),
        shots_by_receiver=np.bincount(receiver_index, minlength=8).tolist(),
        energy=sum(per_shot_energy), min_shot_energy=min(per_shot_energy), max_shot_energy=max(per_shot_energy),
        current_effort_per_radiated_energy=float(sum(effort) / sum(per_shot_energy)),
        field_norms=[complex_norm(m) for m in mu], max_sample_amplitude=np.max(abs(mu), axis=1).tolist(),
        gap=gap, rho99=radius(SIGMA, n), gap_over_2rho99=gap / (2 * radius(SIGMA, n)),
        log10_equal_prior_pe=log_pe(gap), sensitivity=sensitivity(mu, receiver_index, gap),
        pulse_tones_ghz=frequencies, coefficients_ri=encode(weights),
        sample_receiver_zero_based=receiver_index.tolist(), responses_ri=encode(mu),
    )
    return result


def single_tone(experiment, ghz, tx, rx):
    y, q, condition = experiment.at(float(ghz), tx)
    ids = selected_receivers(rx, 0, 0)
    y = y[:, ids]
    beam, value = best_beam(y, q)
    single_beam = np.sqrt(ENERGY) * beam
    capped_beam, capped_value = best_beam(y, q, PORT_EFFORT_CAP)
    capped_beam *= np.sqrt(ENERGY)
    fixed_gap = np.sqrt(ENERGY / y.shape[-1]) * complex_norm(y[0] - y[1])
    return dict(ghz=float(ghz), tx=tx, rx=rx, fixed_gap=fixed_gap, optimal_gap=np.sqrt(ENERGY * value),
                capped_gap=np.sqrt(ENERGY * capped_value),
                optimal_coefficients_ri=encode(single_beam),
                capped_coefficients_ri=encode(capped_beam),
                optimal_amplitudes=np.abs(single_beam).tolist(),
                optimal_phase_deg=np.angle(single_beam, deg=True).tolist(),
                optimal_responses_ri=encode(y @ single_beam),
                capped_responses_ri=encode(y @ capped_beam),
                single_tx_element_fields_ri=encode(y),
                power_gram=q.tolist(), min_power_gram_eigenvalue=float(np.linalg.eigvalsh(q)[0]),
                material_resolvent_min_sv=condition)


def arb_bound(experiment, ghz, rx):
    """Bound gap of the 48 fixed pulses at the specified single frequency."""
    from flint import acb, acb_mat, arb, ctx

    ctx.dps = 55
    original = experiment.original
    # Interpret the saved IEEE input numbers as exact, representing the frozen model.
    k = arb(original["k"]) * arb(ghz) / arb(F0)
    a = arb(original["radius"])
    cells = np.array(original["centers"])
    tx = experiment.sources["old"]
    receivers = experiment.receivers if rx == "switched8" else experiment.receivers[selected_receivers(rx, 0, 0)]
    j = acb(0, 1)

    def arg(x, y):
        return k * sum((arb(float(i)) - arb(float(l)))**2 for i, l in zip(x, y, strict=True)).sqrt()

    def hankel(x, order):
        return acb(x).bessel_j(order) + j * acb(x).bessel_y(order)

    prefactor = j * arb.pi() * k * a / 2
    disk = prefactor * acb(k * a).bessel_j(1)
    self_term = prefactor * hankel(k * a, 1) - 1
    g = acb_mat([[self_term if i == l else disk * hankel(arg(x, y), 0)
                  for l, y in enumerate(cells)] for i, x in enumerate(cells)])
    b = acb_mat([[j / 4 * hankel(arg(x, y), 0) * (arb(ghz) / arb(F0)).sqrt() for y in tx] for x in cells])
    s = acb_mat([[disk * hankel(arg(x, y), 0) for y in cells] for x in receivers])
    answers = []
    for m in experiment.materials:
        d = [arb(float(v)) for v in m]
        operator = acb_mat([[int(i == l) - d[i] * g[i, l] for l in range(len(m))] for i in range(len(m))])
        rhs = acb_mat([[d[i] * b[i, l] for l in range(len(tx))] for i in range(len(m))])
        answers.append(s * operator.solve(rhs))
    # Each old Tx occurs 16 times with energy 3/48, i.e. energy 1 per Tx.
    norm = sum(abs(value)**2 for value in (answers[0] - answers[1]).entries()).sqrt()
    n = 48 if rx == "switched8" else 48 * len(receivers)
    if rx == "switched8":
        # Each receiver sees 1/8 of the occurrences of each old Tx at a fixed tone.
        norm /= arb(8).sqrt()

    def gamma_integer_survival(x, shape):
        term = arb(1)
        total = term
        for i in range(1, shape):
            term *= x / i
            total += term
        return (-x).exp() * total

    sigma_hi = float(SIGMA * (1 + 1e-11))
    prob = arb(1) / 100
    assert gamma_integer_survival((arb(OLD_RADIUS) / arb(sigma_hi))**2, 3) > prob
    noise_hi = float(radius(sigma_hi, n) * (1 + 1e-11))
    assert gamma_integer_survival((arb(noise_hi) / arb(sigma_hi))**2, n) < prob
    sigma_lo = float(SIGMA * (1 - 1e-11))
    assert gamma_integer_survival((arb(OLD_RADIUS) / arb(sigma_lo))**2, 3) < prob
    noise_lo = float(radius(sigma_lo, n) * (1 - 1e-11))
    assert gamma_integer_survival((arb(noise_lo) / arb(sigma_lo))**2, n) > prob
    strictly_separates = bool(norm > 2 * arb(noise_hi))
    return dict(
        ghz=ghz, rx=rx, gap_lower=float(np.nextafter(float(norm.lower()), -np.inf)),
        gap_upper=float(np.nextafter(float(norm.upper()), np.inf)), sigma_upper=sigma_hi,
        two_rho99_upper=2 * noise_hi, two_rho99_lower=2 * noise_lo,
        strictly_separates=strictly_separates, certified_overlap=bool(norm < 2 * arb(noise_lo)),
        arb_dps=ctx.dps,
        rx_samples=n, shots=48, total_energy=ENERGY,
    )


def summarize(experiment):
    factorial = {}
    for name, tones in (("center", np.full(8, F0)), ("narrow", NARROW), ("wide", WIDE)):
        factorial[name] = {}
        for tx in ("old", "uniform6"):
            for rx in ("single", "switched8", "simultaneous8"):
                result = experiment.measure(tones, tx, rx)
                if rx == "switched8":
                    # Rotation is fixed before seeing either material; report sensitivity to that phase.
                    result["all_switch_phase_gaps"] = [experiment.measure(tones, tx, rx, switch_phase=i)["gap"]
                                                        for i in range(8)]
                factorial[name][f"{tx}_{rx}"] = result

    oracle, capped_oracle = {}, {}
    for name, tones in (("narrow", NARROW), ("wide", WIDE)):
        oracle[name] = {
            f"{tx}_{rx}": experiment.measure(tones, tx, rx, waveform="coherent")
            for tx in ("old", "uniform6")
            for rx in ("single", "switched8", "simultaneous8")
        }
        capped_oracle[name] = {
            f"{tx}_{rx}": experiment.measure(tones, tx, rx, waveform="capped")
            for tx in ("old", "uniform6")
            for rx in ("single", "switched8", "simultaneous8")
        }

    sweep_grid = np.unique(np.r_[np.linspace(2, 30, 1121), 15.68, 16., 16.32, NARROW, WIDE])
    sweeps, optimal_tones = {}, {}
    for tx in ("old", "uniform6"):
        for rx in ("single", "simultaneous8"):
            name = f"{tx}_{rx}"
            all_rows = [single_tone(experiment, f, tx, rx) for f in sweep_grid]
            sweeps[name] = [[r["ghz"], r["fixed_gap"], r["optimal_gap"], r["capped_gap"]] for r in all_rows]
            optimal_tones[name] = {}
            for band, lower, upper in (("narrow", 15.68, 16.32), ("wide", 2., 30.)):
                optimal_tones[name][band] = {}
                for waveform in ("optimal", "capped"):
                    winning = dict(max((r for r in all_rows if lower <= r["ghz"] <= upper),
                                       key=lambda r: r[f"{waveform}_gap"]))
                    # Repeat the winning tone/beam at equal per-shot energy for all 48 shots.
                    y, _, _ = experiment.at(winning["ghz"], tx)
                    ids = selected_receivers(rx, 0, 0)
                    coeff_ri = np.array(winning[f"{waveform}_coefficients_ri"])
                    coeff = (coeff_ri[:, 0] + 1j * coeff_ri[:, 1]) / np.sqrt(48)
                    winning["equal_48_shot_replay"] = measurement_result(
                        np.tile(y[:, ids] @ coeff, (1, 48)), np.tile(ids, 48),
                        [winning["ghz"]] * 48, np.tile(coeff, (48, 1)), experiment, tx,
                    )
                    detuning = []
                    for offset_mhz in (-100, -50, -10, -1, 0, 1, 10, 50, 100):
                        f = winning["ghz"] + offset_mhz / 1000
                        drifted_fields, drifted_gram, _ = experiment.at(f, tx)
                        delta = (drifted_fields[0, ids] - drifted_fields[1, ids]) @ coeff
                        detuning.append(dict(
                            offset_mhz=offset_mhz, ratio_to_reference_gap=float(
                                complex_norm(delta) * np.sqrt(48) / winning[f"{waveform}_gap"]
                            ),
                            drifted_radiated_energy=float(np.vdot(coeff, drifted_gram @ coeff).real * 48),
                        ))
                    winning["fixed_coefficients_frequency_detuning"] = detuning
                    optimal_tones[name][band][waveform] = winning

    # Fixed candidate tones chosen from the disclosed uniform grid; no waveform optimization here.
    rigorous = [arb_bound(experiment, 16., "single"), arb_bound(experiment, 16., "switched8"),
                arb_bound(experiment, 16., "simultaneous8"),
                arb_bound(experiment, 22., "single"), arb_bound(experiment, 22., "simultaneous8")]
    geometry_m = dict(
        centers=(np.array(experiment.original["centers"]) * experiment.length_scale_m).tolist(),
        radius=experiment.original["radius"] * experiment.length_scale_m,
        sources={k: (v * experiment.length_scale_m).tolist() for k, v in experiment.sources.items()},
        receivers=(experiment.receivers * experiment.length_scale_m).tolist(),
    )
    source_and_receiver_spacing_mm = {
        tx: float(np.min(np.linalg.norm(v[:, None, :] - experiment.receivers[None, :, :], axis=2))
                  * experiment.length_scale_m * 1000)
        for tx, v in experiment.sources.items()
    }
    return dict(
        case=experiment.witness["case_id"], materials=experiment.materials.tolist(), known_gain=[1, 0],
        inherited_narrow_continuum_bound_ratio=experiment.witness["upper_ratio"],
        material_rms=float(np.sqrt(np.mean((experiment.materials[0] - experiment.materials[1])**2))),
        frequencies_ghz=dict(center=F0, narrow=NARROW.tolist(), wide=WIDE.tolist()),
        sigma=SIGMA, total_energy=ENERGY, k0_per_m=experiment.k0_si,
        length_scale_m=experiment.length_scale_m, geometry_m=geometry_m,
        nearest_tx_rx_spacing_mm=source_and_receiver_spacing_mm,
        factorial=factorial, fixed_tones_optimal_waveform=oracle,
        fixed_tones_capped_waveform=capped_oracle, best_single_tone_on_grid=optimal_tones,
        cap_per_unit_radiated_energy=PORT_EFFORT_CAP,
        rigorous_points=rigorous,
        sweep_columns=["GHz", "fixed_gap_equally_distributed", "optimized_gap", "capped_optimized_gap"],
        sweeps=sweeps,
    )


def plot(result):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "font.size": 10})
    fig, axes = plt.subplots(1, 3, figsize=(14.1, 4.1), width_ratios=(.87, 1.42, 1.28), layout="constrained")
    geometry = result["geometry_m"]
    xy = np.array(geometry["centers"]) * 1e3
    for p in xy:
        axes[0].add_patch(Circle(p, geometry["radius"] * 1e3, ec="#393f45", fc="#bac4c8"))
    for s, lab, marker, color in (
        (geometry["sources"]["old"], "Original 3 Tx", "^", "#bf731b"),
        (geometry["sources"]["uniform6"], "Uniform 6 Tx", "x", "#8f601c"),
        (geometry["receivers"][:1], "Original Rx", "s", "#006c76"),
        (geometry["receivers"][1:], "7 added Rx", "o", "#3aa3ac"),
    ):
        axes[0].scatter(*(np.array(s).T * 1e3), label=lab, marker=marker, c=color, s=26)
    axes[0].set(aspect="equal", xlabel="x (mm)", ylabel="y (mm)", title="Fixed physical geometry")
    axes[0].legend(loc="lower left", bbox_to_anchor=(-.06, -.32), ncol=2, fontsize=8, frameon=False)

    for key, lab, color in (
        ("old_single", "Original Rx", "#bf731b"),
        ("old_simultaneous8", "8 Rx simultaneously", "#006c76"),
    ):
        sweep = np.array(result["sweeps"][key])
        axes[1].plot(sweep[:, 0], sweep[:, 1], color=color, label=lab)
        axes[1].plot(sweep[:, 0], sweep[:, 2], color=color, ls="--", alpha=.75)
    axes[1].axvspan(15.68, 16.32, color="#d4d4d4", alpha=.6)
    axes[1].set(yscale="log", xlabel="Frequency (GHz)", ylabel="Two-material signal gap",
                title="Original three transmitters")
    axes[1].legend(loc="lower right", fontsize=8, frameon=False)
    axes[1].text(.99, .98, "Solid  Equal energy per Tx\nDashed  Pair-optimized coherent beam",
                 transform=axes[1].transAxes, fontsize=8, va="top", ha="right")

    width = .23
    for i, (rx, lab, color) in enumerate((
        ("single", "Original Rx, 48 reads", "#bf731b"),
        ("switched8", "Switched 8 Rx, 48 reads", "#73b2b4"),
        ("simultaneous8", "Simultaneous 8 Rx, 384 reads", "#006c76"),
    )):
        vals = [result["factorial"][band][f"old_{rx}"]["gap"] for band in ("center", "narrow", "wide")]
        axes[2].bar(np.arange(3) + (i-1)*width, vals, width, label=lab, color=color)
    axes[2].set(yscale="log", title="Fixed waveforms and equal Tx energy", ylabel="Two-material signal gap",
                xticks=np.arange(3), xticklabels=("16 GHz", "Narrow band", "2–30 GHz"))
    axes[2].legend(loc="upper center", bbox_to_anchor=(.5, -.14), fontsize=8, frameon=False)
    for ax in axes[1:]:
        ax.grid(axis="y", alpha=.18)
        ax.set_axisbelow(True)
    stem = ROOT / "figures/pair_mimo"
    fig.savefig(stem.with_suffix(".png"), dpi=195)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout="constrained", sharey=True)
    for ax, key, label in zip(axes, ("uniform6_single", "uniform6_simultaneous8"),
                              ("Original receiver", "Eight receivers in parallel"), strict=True):
        freq, fixed, ideal, capped = np.array(result["sweeps"][key]).T
        ax.plot(freq, ideal, color="#a05451", ls="--", label="Radiated energy only")
        ax.plot(freq, capped, color="#006c76", label="Radiated energy and drive cap")
        ax.plot(freq, fixed, color="#bf731b", label="Equal energy per transmitter")
        ax.axvline(5.075, linewidth=1, color="#505050", alpha=.5)
        ax.axvspan(15.68, 16.32, color="#d4d4d4", alpha=.6)
        ax.set(yscale="log", xlabel="Frequency (GHz)", title=label)
        ax.grid(axis="y", alpha=.18)
    axes[0].set_ylabel("Two-material signal gap")
    axes[1].legend(loc="upper right", fontsize=8, frameon=False)
    stem = ROOT / "figures/pair_source_effort"
    fig.savefig(stem.with_suffix(".png"), dpi=195)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


if __name__ == "__main__":
    comparison = Comparison()
    summary = summarize(comparison)
    (ROOT / "results/pair_mimo.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    plot(summary)
    for band, runs in summary["factorial"].items():
        for key, r in runs.items():
            print(f"{band:>6} {key:>22}  N={r['rx_samples']:3d}  gap={r['gap']:.8g} "
                  f"gap/(2r99)={r['gap_over_2rho99']:.5g}  log10Pe={r['log10_equal_prior_pe']:.6g}")
    for witness in summary["rigorous_points"]:
        print(witness)
