"""Plot paired Stage 2 outcomes from the published per-case CSV."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import LogLocator, NullFormatter
import numpy as np


PILOT = Path(__file__).resolve().parents[1]
INPUT = PILOT / "results/cases.csv"
OUTPUT = PILOT / "figures"
INK = "#1e2634"
MUTED = "#53606e"
TEAL = "#137369"
BORDER = "#d6dde2"
COVERAGE = LinearSegmentedColormap.from_list(
    "certificate_coverage", ("#f1f5f3", "#a5d6c9", "#14766a")
)
DECISIONS = {
    "PASS": ("Strict exclusion", TEAL, ""),
    "BELOW": ("Upper bound at/below radius", "#dabc83", ""),
    "ZERO": ("Proved zero", "#d8dee6", ""),
    "UNKNOWN": ("Unresolved", "#a9a0bd", "///"),
}
GROUPS = (
    *(("primary", "Fresh6", noise, cap) for noise in (1, 2, 10) for cap in (0.5, 0.75, 1.0)),
    ("primary", "SRAP22", 10, 0.5),
    ("primary", "SRAP22", 10, 0.75),
    ("supplement", "SRAP22", 1, 0.75),
    ("supplement", "SRAP22", 2, 0.75),
)
CAPS = (0.5, 0.75, 1.0)
COHORTS = ("Fresh6", "SRAP22")
COHORT_NAMES = {"Fresh6": "Fresh shapes (6)", "SRAP22": "Archived SRAP shapes (22)"}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "text.color": INK,
    "axes.labelcolor": INK,
    "axes.edgecolor": BORDER,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "xtick.color": MUTED,
    "ytick.color": INK,
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})


@dataclass(frozen=True)
class Method:
    arm: str
    side: int
    name: str
    mark: str
    color: str


FIXED = Method("fixed_portfolio", 0, "Fixed portfolio", "D", "#656b74")
FIXED_COMPONENTS = (
    Method("global_real", 0, "Fixed global", "v", "#687382"),
    Method("fixed_16_0", 4, "Fixed 16 #0", "^", "#687382"),
    Method("fixed_16_2", 4, "Fixed 16 #2", "s", "#687382"),
)
DYNAMIC = (
    Method("no_feature", 8, "Geometry-only", "X", "#73834e"),
    Method("neural", 8, "Neural", "o", TEAL),
    Method("hybrid", 8, "Hybrid", "s", "#296b9e"),
    Method("gradient1", 8, "Gradient 1 s", "^", "#a76538"),
    Method("gradient5", 8, "Gradient 5 s", "v", "#835462"),
)


@dataclass(frozen=True)
class Observation:
    study: str
    cohort: str
    object_id: str
    noise: int
    cap: float

    @property
    def stratum(self):
        return self.study, self.cohort, self.noise, self.cap


@dataclass(frozen=True)
class Row:
    observation: Observation
    radius: float
    decision: str
    charged_wall: float | None
    coarse_d: float | None

    @property
    def passed(self):
        return self.decision == "PASS"


def optional_number(s):
    if s in (None, ""):
        return None
    number = float(s)
    return number if np.isfinite(number) else None


def read_cases(path, preview=False):
    trials = defaultdict(dict)
    with path.open(newline="") as stream:
        for record in csv.DictReader(stream):
            obs = Observation(record["study"], record["cohort"], record["oid"],
                              int(float(record["noise"])), float(record["u"]))
            key = record["arm"], int(float(record["side"]))
            verdict = record["decision"].upper()
            if verdict not in DECISIONS or key in trials[obs]:
                raise ValueError(f"Unexpected decision or repeated method: {obs}, {key}, {verdict}")
            data = Row(obs, float(record["radius"]), verdict,
                       optional_number(record.get("charged_wall")), optional_number(record.get("coarse_d")))
            trials[obs][key] = data
    if not trials:
        raise ValueError("No case rows are available.")
    if not preview:
        totals = Counter(obs.stratum for obs in trials)
        expected = {stratum: (6 if stratum[1] == "Fresh6" else 22) for stratum in GROUPS}
        if totals != expected:
            raise ValueError(f"{path} does not contain the frozen 142 observation/cap conditions.")
        for obs, methods in trials.items():
            required = {(m.arm, m.side) for m in (FIXED, *FIXED_COMPONENTS, *DYNAMIC)}
            if not required.issubset(methods):
                raise ValueError(f"{path} is incomplete for {obs}. Use --preview for a partial table.")
    return trials


def row_for(methods, method):
    return methods.get((method.arm, method.side))


def matched(rows, method):
    result = []
    for condition in rows:
        arm = row_for(condition, method)
        reference = row_for(condition, FIXED)
        if arm is not None and reference is not None:
            result.append((arm, reference))
    return result


def heading(fig, title, detail, preview):
    fig.text(.065, .958, title, fontsize=16, weight="bold", ha="left")
    fig.text(.065, .910, detail, fontsize=9, color=MUTED, ha="left")
    if preview:
        fig.text(.955, .960, "INCOMPLETE PREVIEW", color="#9f5f30", fontsize=8,
                 ha="right", weight="bold")


def certified_outcomes(trials, preview):
    methods = (FIXED, *DYNAMIC)
    strata = defaultdict(list)
    for obs, rows in trials.items():
        strata[obs.stratum].append(rows)
    fig, ax = plt.subplots(figsize=(12.1, 8.05))
    fig.subplots_adjust(left=.218, right=.955, top=.815, bottom=.130)
    heading(fig, "Strict exclusions on the native 96² grid",
            "Same material, noise and contrast cap paired with the frozen portfolio on every line.",
            preview)
    for index, group in enumerate(GROUPS):
        for pane, method in enumerate(methods):
            pairs = matched(strata[group], method)
            n = len(pairs)
            p = sum(a.passed for a, _ in pairs)
            unresolved = sum(a.decision == "UNKNOWN" for a, _ in pairs)
            unmatched = sum((a.passed and b.decision == "UNKNOWN") or
                            (b.passed and a.decision == "UNKNOWN") for a, b in pairs)
            gains = sum(a.passed and b.decision in ("ZERO", "BELOW") for a, b in pairs)
            losses = sum(b.passed and a.decision in ("ZERO", "BELOW") for a, b in pairs)
            fraction = p / n if n else 0
            face = COVERAGE(fraction) if n else "#fbfbfc"
            text_color = "white" if fraction >= .81 else INK
            ax.add_patch(Rectangle((pane-.475, index-.45), .95, .90, facecolor=face,
                                   ec="white", lw=.6))
            if n:
                parts = []
                if gains or losses:
                    parts.append(f"+{gains}  −{losses}")
                if unresolved:
                    parts.append(f"U{unresolved}")
                if unmatched:
                    parts.append(f"A{unmatched}")
                ax.text(pane, index-.13 if parts else index, f"{p}/{n}", ha="center", va="center",
                        fontsize=9.5, weight="bold", color=text_color)
                if parts:
                    ax.text(pane, index+.19, "   ".join(parts), ha="center", va="center", fontsize=7.2,
                            color=text_color if fraction >= .81 else MUTED)
            else:
                ax.text(pane, index, "—", ha="center", va="center", color=BORDER, fontsize=12)
    ax.set(xlim=(-.5, len(methods)-.5), ylim=(len(GROUPS)-.5, -.5))
    ax.set_xticks(range(len(methods)), [m.name + ("\nCommon reference" if m is FIXED else "\n64 regions")
                                                for m in methods])
    labels = [
        f"{cohort}    {noise:>2}%    $u$={cap:.2f}" + ("  supp." if study == "supplement" else "")
        for study, cohort, noise, cap in GROUPS
    ]
    ax.set_yticks(range(len(GROUPS)), labels)
    ax.xaxis.tick_top()
    ax.tick_params(axis="both", length=0, pad=7)
    for edge in (8.5, 10.5):
        ax.axhline(edge, color=BORDER, lw=1.1)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.text(.065, .078, "Cell shading tracks exclusions / matched cases.  + / − = confirmed gains / losses against the fixed portfolio.",
             fontsize=8.3, color=MUTED)
    fig.text(.065, .049, "U = unresolved candidate outcomes.  A = pairs with one strict exclusion and one unresolved outcome.  Each row counts distinct objects.",
             fontsize=8.3, color=MUTED)
    return fig


def timing_summary(rows, method):
    pairs = matched(rows, method)
    by_object = defaultdict(list)
    for a, b in pairs:
        if a.passed and b.passed and a.charged_wall and b.charged_wall:
            ratio = a.charged_wall / b.charged_wall
            if ratio > 0 and np.isfinite(ratio):
                by_object[a.observation.object_id].append(ratio)
    ratios = np.array([np.median(values) for values in by_object.values()])
    coverage = sum(a.passed for a, _ in pairs)
    unresolved = sum(a.decision == "UNKNOWN" for a, _ in pairs)
    quartiles = np.quantile(ratios, (.25, .50, .75)) if len(ratios) else None
    return len(pairs), coverage, unresolved, len(ratios), quartiles


def compute_reliability(trials, preview):
    methods = (FIXED, *FIXED_COMPONENTS, *DYNAMIC)
    panes = {}
    estimates = [1.]
    for cohort in COHORTS:
        for cap in (0.5, 0.75):
            records = [v for o, v in trials.items() if (o.cohort, o.cap) == (cohort, cap)]
            pane = [timing_summary(records, method) for method in methods]
            panes[cohort, cap] = pane
            for _, _, _, _, quartiles in pane:
                if quartiles is not None:
                    estimates.extend(quartiles)
    lower_exp = np.floor(np.log10(min(estimates)*.8))
    upper_exp = np.ceil(np.log10(max(estimates)*1.25))
    xmin, xmax = 10**lower_exp, 10**upper_exp

    fig, axes = plt.subplots(2, 2, figsize=(13, 8.4))
    fig.subplots_adjust(left=.16, right=.88, bottom=.145, top=.845, hspace=.47, wspace=.87)
    heading(fig, "Online cost on matched successful certificates",
            "Charged wall time relative to the fixed portfolio on the same material, cap and noise level.",
            preview)
    for i, cap in enumerate((0.5, 0.75)):
        for j, cohort in enumerate(COHORTS):
            ax = axes[i, j]
            ax.set_title(f"{COHORT_NAMES[cohort]}   $u$={cap:.2f}", loc="left",
                         fontsize=10, weight="bold", pad=11)
            ax.set_xscale("log")
            ax.set(xlim=(xmin, xmax), ylim=(len(methods)-.55, -.55))
            ax.set_yticks(range(len(methods)), [m.name for m in methods])
            ax.xaxis.set_major_locator(LogLocator(base=10., numticks=8))
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.grid(axis="x", which="major", color="#e6eaee", lw=.75)
            ax.axvline(1., color="#88929d", lw=1., ls=(0, (3, 3)))
            ax.yaxis.set_tick_params(length=0)
            ax.spines["left"].set_visible(False)
            ax.axhline(3.5, color=BORDER, lw=.75)
            right_columns = (1.035, 1.245, 1.365)
            for x, label in zip(right_columns, ("Pass", "?", "$n$"), strict=True):
                ax.text(x, 1.026, label, ha="left", va="bottom", fontsize=8,
                        color=MUTED, weight="bold", transform=ax.transAxes)
            for y, (method, summary) in enumerate(zip(methods, panes[cohort, cap], strict=True)):
                nrows, passes, unk, nobjects, quartiles = summary
                if quartiles is not None:
                    q1, med, q3 = quartiles
                    ax.plot((q1, q3), (y, y), color=method.color, lw=1.7, zorder=3)
                    ax.plot(med, y, linestyle="none", marker=method.mark, color=method.color,
                            ms=5.4, zorder=4)
                values = (f"{passes}/{nrows}" if nrows else "—", str(unk) if nrows else "—",
                          str(nobjects) if nrows else "—")
                for x, value in zip(right_columns, values, strict=True):
                    ax.text(x, y, value, ha="left", va="center", fontsize=8.1, color=MUTED,
                            transform=ax.get_yaxis_transform())
    fig.text(.065, .095, "Points mark the median object-level time ratio. Lines show the middle 50%. Ratios below 1 favor the method shown.",
             fontsize=8.3, color=MUTED)
    fig.text(.065, .065, "Repeated successful noise pairs are first collapsed within each material. $n$ counts paired materials. Pass and ? count all observed conditions.",
             fontsize=8.3, color=MUTED)
    fig.text(.065, .038, "Fixed profiles can reuse state proofs prepared offline. Their plotted figures use the recorded online charges. Missing points indicate no timed matched success.",
             fontsize=8.3, color=MUTED)
    return fig


def coarse_native_outcomes(trials, preview):
    totals = defaultdict(Counter)
    for observation, rows in trials.items():
        for method in DYNAMIC:
            result = row_for(rows, method)
            if result is not None and result.coarse_d is not None:
                screened = result.coarse_d > result.radius
                totals[observation.cap, method, screened][result.decision] += 1
    fig, axes = plt.subplots(1, 3, figsize=(12.7, 5.1), sharey=True)
    fig.subplots_adjust(left=.135, right=.965, top=.77, bottom=.285, wspace=.20)
    heading(fig, "Coarse screening and the native verdict",
            "64-region proposals. Each pair of bars uses the same arm and contrast cap.", preview)
    limit = 1.10 * max(1, max((sum(counts.values()) for counts in totals.values()), default=0))
    for cap, ax in zip(CAPS, axes, strict=True):
        ax.set_title(f"Contrast cap $u$={cap:.2f}", fontsize=10, loc="left", weight="bold", pad=9)
        for j, method in enumerate(DYNAMIC):
            for flag, offset in ((True, -.20), (False, .20)):
                counts = totals[cap, method, flag]
                start = 0
                for status, (_, color, hatch) in DECISIONS.items():
                    count = counts[status]
                    if not count:
                        continue
                    ax.barh(j+offset, count, left=start, height=.32, color=color, hatch=hatch,
                            edgecolor="#6e627d" if hatch else "white", linewidth=.4, zorder=3)
                    start += count
        ax.set(xlim=(0, limit), ylim=(len(DYNAMIC)-.5, -.5))
        ax.set_yticks(range(len(DYNAMIC)), [d.name for d in DYNAMIC])
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        ax.grid(axis="x", color="#e5e9ee", lw=.75)
    fig.text(.55, .213, "Observation / cap cases", ha="center", fontsize=9)
    handles = [Patch(label=label, facecolor=color, edgecolor="#6e627d" if hatch else color, hatch=hatch)
               for label, color, hatch in DECISIONS.values()]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.55, .147),
               frameon=False, ncol=4, columnspacing=2, fontsize=8.7)
    fig.text(.065, .103, "The upper bar in each pair counts coarse distances above the radius. The lower bar counts distances at or below it.",
             fontsize=8.3, color=MUTED)
    fig.text(.065, .064, "Counts include paired noise conditions. Unresolved native checks are shown separately from both kinds of proved non-pass.",
             fontsize=8.3, color=MUTED)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--preview", action="store_true", help="Draw a labelled preview from incomplete case rows.")
    args = parser.parse_args()
    if args.preview and (args.out is None or args.out.resolve() == OUTPUT):
        parser.error("--preview requires a separate --out directory.")
    out = args.out or OUTPUT
    trials = read_cases(args.input, preview=args.preview)
    out.mkdir(parents=True, exist_ok=True)
    for name, render in (
        ("certified_outcomes", certified_outcomes),
        ("compute_reliability", compute_reliability),
        ("coarse_native_outcomes", coarse_native_outcomes),
    ):
        fig = render(trials, args.preview)
        for extension in ("png", "pdf"):
            target = out / f"{name}.{extension}"
            fig.savefig(target, dpi=240)
            print(target, flush=True)
        plt.close(fig)


if __name__ == "__main__":
    main()
