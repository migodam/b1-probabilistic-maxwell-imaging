"""Run the independent eight-material, finite-aperture six-cell acquisition."""

from importlib.metadata import version
import json

from natural_acquisition import (
    ARMS, EPS, GAIN, HI, LO, ROOT, Experiment, acquire_extra, acquire_old,
    as_complex_pair, find_central_roots, noisy_vector, report_feed, rms, select_feeds,
)
import numpy as np
from scipy.special import j0

MATERIAL_SEED = 2026111901
NOISE_SEED = 705190
SOURCE_SEED = 865032
OLD_ANGLES = [-158., -77., 19.]
ADDED_ANGLES = np.linspace(-179., 41., 20)
SOURCE_RADIUS = 1.2


def circle(angles):
    radians = np.deg2rad(angles)
    return SOURCE_RADIUS * np.column_stack((np.cos(radians), np.sin(radians)))


def material(i):
    return np.random.default_rng(MATERIAL_SEED+i).uniform(10., 30., size=6)


def ordered_roots(rows, n):
    roots = np.asarray(rows).reshape(-1, n)
    distance = np.linalg.norm(roots-(LO+HI)/2, axis=1)
    return roots[np.argsort(distance, kind="stable")]


def main():
    geometry = dict(
        centers=[[-.42,-.23],[-.06,-.35],[.37,-.17],[-.31,.19],[.16,.05],[.29,.37]],
        radius=.053, k=7.05, receivers=[[-.20,1.03]], sources=circle(OLD_ANGLES).tolist(),
    )
    old = Experiment(geometry)
    source_positions = np.vstack((old.model.sources,circle(ADDED_ANGLES)))
    physical = Experiment(dict(geometry,sources=source_positions.tolist()))
    source_separations = np.linalg.norm(source_positions[:,None]-source_positions[None,:], axis=-1)
    gram = j0(old.model.k*source_separations)
    names = [f"confirm6-{i:02d}" for i in range(8)]
    settings = dict(
        geometry=geometry, old_source_angles_deg=OLD_ANGLES, old_aperture_deg=177.,
        added_source_angles_deg=ADDED_ANGLES.tolist(), added_aperture_deg=220.,
        external_source_radius=SOURCE_RADIUS, finite_source_positions=source_positions.tolist(),
        material_seed_base=MATERIAL_SEED, material_range=[10.,30.],
        noise_seed_base=NOISE_SEED, source_seed_base=SOURCE_SEED,
        additive_radius=float(EPS), gain_radius=float(GAIN), additive_radius_fraction=.55,
        gain_radius_fraction=.8, joint_radius_used=float(np.hypot(.55,.8)),
        pivot_index_origin=0, homotopy_pypolsys_version=version("pypolsys"),
        truth_use="Independent full VIE acquisitions and post-design scoring.",
    )
    report = {"settings":settings, "cases":{}}
    data = dict(case_ids=names, arms=ARMS, eps=EPS, gamma=GAIN, box=[LO,HI])
    choices = []
    for i, name in enumerate(names):
        rng = np.random.default_rng(NOISE_SEED+i)
        epsilon = noisy_vector(rng, 4) * EPS * .55
        gain_error = GAIN * .8 * np.exp(2j*np.pi*rng.uniform())
        y_old = acquire_old(geometry, material(i), epsilon, gain_error)
        hits = find_central_roots(geometry, y_old)
        candidates = ordered_roots(hits["validated_roots"], old.n)
        feeds, rule = select_feeds(physical, gram, old.model.ns, candidates, SOURCE_SEED+i)
        pivot = np.argsort(-np.ptp(candidates,axis=0),kind="stable")[:2] if len(candidates)>=2 else np.array([0,1])
        choices.append((y_old,epsilon,gain_error,candidates,feeds,rule,pivot,hits))
        print(name, len(candidates), rule, pivot.tolist(), flush=True)

    # The exact old synthetic read and its alternative roots are opened after choosing all feeds.
    for i, (name, selected) in enumerate(zip(names,choices)):
        y_old,epsilon,gain_error,candidates,feeds,rule,pivot,hits = selected
        truth = material(i)
        new = acquire_extra(physical, truth, feeds, epsilon, gain_error)
        clean_y = old.model.forward(truth, independent=True)["y"][0]
        clean_hits = find_central_roots(geometry, clean_y)
        clean_roots = ordered_roots(clean_hits["validated_roots"], old.n)
        record = dict(
            rule=rule, old_shared_gain_error=as_complex_pair(gain_error),
            additive_noise=[as_complex_pair(z) for z in epsilon], pivots=pivot.tolist(),
            central_roots_count=len(candidates), central_roots=candidates.tolist(),
            central_truth_RMS=[float(rms(c-truth)) for c in candidates],
            clean_witness_count=len(clean_roots), clean_witness_methods=[["POLSYS_PLP"] for _ in clean_roots],
            clean_witness_truth_RMS=[float(rms(c-truth)) for c in clean_roots],
            central_path_classifications=hits["path_classifications"],
            central_paths_attempted=hits["paths_attempted"], central_solve_status=hits["solve_status"],
            central_total_cpu_s=hits["total_cpu_seconds"], clean_path_classifications=clean_hits["path_classifications"],
            old_read=[as_complex_pair(z) for z in y_old], arms={},
        )
        payload = dict(
            B0=old.b, H=old.h, R=old.r, host=np.zeros(old.n), y_old=y_old, y_old_clean=clean_y,
            central_roots=candidates, clean_witnesses=clean_roots, truth=truth, pivots=pivot,
            source_positions=source_positions, power_gram=gram, additive_noise=epsilon, gain_error=gain_error,
        )
        data.update({f"{name}/{k}":v for k,v in payload.items()})
        predictions = [physical.response(c)[0] for c in candidates]
        for arm,feed in feeds.items():
            b = physical.b @ feed
            B = np.column_stack((old.b,b))
            y = np.r_[y_old,new[arm]]
            record["arms"][arm] = dict(report_feed(feed,gram,[v@feed for v in predictions]),
                                        reading=as_complex_pair(new[arm]),
                                        role="random fallback" if not len(candidates) else arm,
                                        rank=int(np.linalg.matrix_rank(B)))
            for k,v in dict(b=b,B=B,y=y,y_host=np.zeros(4,dtype=complex),feed=feed,pivots=pivot).items():
                data[f"{name}/{arm}/{k}"] = v
        report["cases"][name] = record
    np.savez_compressed(ROOT/"results"/"confirm_inputs.npz", **data)
    (ROOT/"results"/"confirm_acquisition.json").write_text(json.dumps(report,indent=2,allow_nan=False))


if __name__ == "__main__":
    main()
