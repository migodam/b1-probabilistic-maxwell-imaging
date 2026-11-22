"""Generate the frozen transfer16 observations and separate offline labels."""

import json
from random import SystemRandom
from time import perf_counter, process_time

from generate_fields import HERE, PILOT, arb, ctx, material, native_geometry, np, observe, threadpool_limits
from generate_srap_bands import high, low


DESIGN = HERE / "TRANSFER16_DESIGN.json"
OFFLINE = HERE / "transfer16_offline.json"
RECEIPT = HERE / "transfer16_generation.json"


def object_ids(count):
    ids = [f"transfer16_{n + 1:02d}" for n in range(count)]
    if OFFLINE.is_file():
        labels = json.loads(OFFLINE.read_text())["rows"]
        indexed = sorted(labels, key=lambda row: row["material_index"])
        if [row["material_index"] for row in indexed] != list(range(count)):
            raise ValueError("The offline index-to-ID mapping does not cover the frozen cohort.")
        saved_ids = [row["object_id"] for row in indexed]
        if sorted(saved_ids) != ids:
            raise ValueError("The saved opaque IDs are different from the frozen cohort.")
        return saved_ids
    SystemRandom().shuffle(ids)
    return ids


def write_json(path, contents):
    staged = path.with_suffix(".json.tmp")
    staged.write_text(json.dumps(contents, indent=2, allow_nan=False) + "\n")
    staged.replace(path)


def main():
    begin, cpu = perf_counter(), process_time()
    design = json.loads(DESIGN.read_text())
    ctx.prec, ctx.threads = design["precision_bits"], design["threads"]
    ids = object_ids(design["n_materials"])
    geometry = native_geometry()
    setup = {"wall_s": perf_counter() - begin, "cpu_s": process_time() - cpu}
    noise_pct = design["noise_pct"]
    observations, labels, audit = {}, [], []
    for i, oid in enumerate(ids):
        start, item_cpu = perf_counter(), process_time()
        m, family, target_peak = material(i, seed=design["geometry_seed"], high_count=design["high_count"])
        peak = float(m.max())
        limits = design["high_peak_uniform"] if i < design["high_count"] else design["low_peak_uniform"]
        if m.size != np.prod(design["native_grid"]) or not limits[0] <= peak <= limits[1]:
            raise ArithmeticError("The generated material differs from its frozen grid or peak band.")
        ys, radii, slacks = observe(
            geometry, m, index=i, noise_seed=design["noise_seed"], noise_pct=noise_pct,
        )
        margins = [low(slack) for slack in slacks]
        distance_uppers = [high(arb(r) - slack) for r, slack in zip(radii, slacks, strict=True)]
        if any(margin <= 0 or upper >= r for margin, upper, r in zip(margins, distance_uppers, radii, strict=True)):
            raise ArithmeticError(f"The outward-rounded truth-in-ball certificate failed for {oid}.")
        observations[oid] = (ys, radii)
        labels.append({
            "material_index": i, "object_id": oid,
            "geometry": design["family_order_by_index_mod_4"][family],
            "geometry_seed_sequence_entropy": [design["geometry_seed"], i],
            "noise_seed_sequence_entropy": [design["noise_seed"], i],
            "sampled_target_peak": float(target_peak), "peak": peak,
            "under_cap": peak <= design["peak_cap"],
        })
        audit.append({
            "object_id": oid, "radius_slack_lower": margins,
            "truth_distance_upper": distance_uppers,
            "wall_s": perf_counter() - start, "cpu_s": process_time() - item_cpu,
        })
        print(oid, round(audit[-1]["wall_s"], 2), flush=True)

    public_ids = sorted(observations)
    public = PILOT / design["data_output"]
    staged = public.with_suffix(".npz.tmp")
    with staged.open("wb") as stream:
        np.savez_compressed(
            stream,
            object_id=np.array(public_ids),
            Y=np.array([observations[oid][0] for oid in public_ids]),
            radius=np.array([observations[oid][1] for oid in public_ids]),
            noise_pct=np.array(noise_pct),
        )
    write_json(OFFLINE, {"cohort": design["cohort"], "design": DESIGN.name, "rows": labels})
    write_json(RECEIPT, {
        "cohort": design["cohort"], "design": DESIGN.name, "public_data": design["data_output"],
        "objects": len(ids), "observations": len(ids) * len(noise_pct), "noise_pct": noise_pct,
        "physical_grid": design["native_grid"][0], "frequency_MHz": design["frequency_MHz"],
        "array": design["array"], "precision_bits": ctx.prec,
        "native_field_rhs": len(ids) * design["sources"],
        "native_receiver_adjoint_rhs": len(ids) * design["receivers"],
        "geometry_setup": setup, "wall_s": perf_counter() - begin, "cpu_s": process_time() - cpu,
        "rows": sorted(audit, key=lambda row: row["object_id"]),
    })
    staged.replace(public)


if __name__ == "__main__":
    with threadpool_limits(limits=1):
        main()
