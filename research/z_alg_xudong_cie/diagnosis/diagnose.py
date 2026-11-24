#!/usr/bin/env python3
"""Read-only snapshots of the mature per-source Dual-CIE trajectory."""

import argparse
from collections import Counter
from contextlib import contextmanager
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import sys
from time import perf_counter, process_time

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"

import numpy as np
from scipy.linalg import eigvalsh
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
SRAP = RESEARCH / "srap_cie_array_study_20261121"
LA = RESEARCH / "limited_aperture_routed_dr_20261120/.runs"
OPTICS = RESEARCH / "z_alg_optical_inertia"
DATA = OPTICS / "synthetic"
sys.path.insert(0, str(SRAP))
import engine

Austria = "la180_austria_rot90_eps3p5"
Slotted = "optfresh_01__pw180_10x20__n02"
Shell = "s2d_004__pw180_10x20"
Twin = "fresnel_twin_8GHz_half180_tx"
CASES = (Austria, Slotted, Shell, Twin)


def read_npz(path):
    with np.load(path, allow_pickle=False) as arrays:
        return {key: arrays[key] for key in arrays.files}


def write(path, data):
    tmp = path.with_suffix(".partial")
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def get_case(name):
    """Return one of four historical settings and public training arrays."""
    if name not in CASES:
        raise ValueError(f"Case not in the four-case diagnostic: {name}")
    source = (LA / "public" if name in (Austria, Twin) else SRAP / "data/public") / f"{name}.npz"
    if name != Slotted:
        public = read_npz(source)
        case = json.loads(str(public.pop("case_json")))
        if name in (Austria, Twin):
            modes, betas, rank = engine.legacy.stage_table(case)
            case["public_schedule"] = dict(modes=modes, betas=betas, rank=rank)
        case["measurement_files"] = [str(source.relative_to(RESEARCH))]
        return case, public

    source = DATA / f"public/{name}.npz"
    geom_path = DATA / "public/pw180_10x20.npz"
    reference = SRAP / f"data/public/{Shell}.npz"
    base_case, base = get_case(Shell)
    geom = read_npz(geom_path)
    if not (np.array_equal(geom["B"], base["E"]) and np.array_equal(geom["S"], base["S"])):
        raise ValueError("The optical inertia archive and SRAP inverse geometry do not match.")
    case = {**base_case, "case": name, "object_id": "optfresh_01", "noise_relative": .02}
    radius_path = OPTICS / "results/theory/fine_feasibility.json"
    radius_record = next(item for item in json.loads(radius_path.read_text())["observations"]
                         if item["id"] == name)
    # Delta is unused by the Dual-CIE updates, but the archived value is recorded.
    case["measurement_files"] = [str(path.relative_to(RESEARCH))
                                 for path in (source, geom_path, reference, radius_path)]
    public = dict(E=geom["B"], S=geom["S"], Y=read_npz(source)["Y"],
                  centers=base["centers"], k=base["k"], radius=base["radius"],
                  delta=float(radius_record["inflated_radius"]))
    return case, public


class EngineClock(engine.legacy.Clock):
    """Stop on solver wall time and book diagnostic overhead separately."""

    def __init__(self):
        self.paused_wall = 0.
        self.paused_cpu = 0.
        super().__init__(float("inf"))

    def wall(self):
        return super().wall() - self.paused_wall

    def cpu(self):
        return super().cpu() - self.paused_cpu

    @contextmanager
    def pause(self):
        wall, cpu = perf_counter(), process_time()
        try:
            yield
        finally:
            self.paused_wall += perf_counter() - wall
            self.paused_cpu += process_time() - cpu


def quotient(numerator, denominator):
    return float(numerator / denominator) if denominator > 0 else None


def col_energy(values):
    return np.sum(abs(values)**2, axis=0)


def helm(a):
    return (a + a.conj().T) / 2


def optical_terms(J, Q):
    jj = J.conj().T @ J
    qq = Q.conj().T @ Q
    jq = J.conj().T @ Q
    O = (jq - jq.conj().T) / (2j)
    diagonal = np.diag(O).real.copy()
    off = O - np.diag(diagonal)
    n = len(diagonal)
    sj, sq = col_energy(J), col_energy(Q)
    denom = float(np.sqrt(sj.sum() * sq.sum()))
    d, f = float(np.linalg.norm(diagonal)), float(np.linalg.norm(off))
    strength = np.sqrt(sj * sq)
    pair_den = np.sqrt(np.outer(strength, strength))
    entry = np.divide(O, pair_den, out=np.zeros_like(O), where=pair_den > 0)
    entry_d = np.diag(entry).real
    entry_o = entry - np.diag(entry_d)
    optical_eig = eigvalsh(O, check_finite=False)
    off_eig = eigvalsh(off, check_finite=False)
    raw = dict(J_frob2=float(sj.sum()), E_frob2=float(sq.sum()), J_source_energy=sj.tolist(),
               E_source_energy=sq.tolist(), J_Gram_fro=float(np.linalg.norm(jj)),
               E_Gram_fro=float(np.linalg.norm(qq)), JE_Gram_fro=float(np.linalg.norm(jq)),
               optical_scale=denom, optical_diag=d, optical_offdiag=f,
               optical_diag_rms=d/np.sqrt(n), optical_offdiag_rms=f/np.sqrt(n*(n-1)),
               optical_diag_norm=quotient(d, denom), optical_offdiag_norm=quotient(f, denom),
               optical_diag_per_source=diagonal.tolist(), optical_diag_per_source_normalized=entry_d.tolist(),
               optical_diag_pair_normalized_rms=float(np.linalg.norm(entry_d)/np.sqrt(n)),
               optical_offdiag_pair_normalized_rms=float(np.linalg.norm(entry_o)/np.sqrt(n*(n-1))),
               optical_diag_energy_fraction=quotient(d*d, d*d+f*f),
               optical_eigenvalues=optical_eig.tolist(), optical_offdiag_eigenvalues=off_eig.tolist(),
               optical_operator_norm=float(np.max(abs(optical_eig))),
               optical_rank_relative_1e_9=int(np.count_nonzero(abs(optical_eig) > 1.e-9 * denom)))
    return raw, jj, jq


def pixel_arrays(J, Q, upper):
    a = (abs(J)**2).sum(axis=1)
    b = (abs(Q)**2).sum(axis=1)
    c = (J.conj()*Q).sum(axis=1)
    d = a - upper*c.real
    root = np.sqrt(d*d + upper*upper*np.maximum(a*b-abs(c)**2, 0.))
    sector = np.maximum(.5*(d + root), 0.)
    return a, b, c, sector


def audit_record(model, m, Y, step, stage, kind):
    start = (perf_counter(), process_time())
    before = model.clock.counts.copy()
    state = model.state(m)
    f = model.ports.receive(state["j"])
    d = masked(model, f-Y)
    den = float(np.linalg.norm(masked(model, Y)))
    physical_raw, _, _ = optical_terms(state["j"], state["Q"])
    front = float(np.linalg.norm(state["j"] - m[:, None]*state["Q"]))
    return dict(step=step, stage=stage, point=kind, fit=float(np.linalg.norm(d)),
                fit_relative=quotient(np.linalg.norm(d), den),
                physical_optical_diag=physical_raw["optical_diag"],
                physical_optical_offdiag=physical_raw["optical_offdiag"],
                physical_optical_diag_normalized=physical_raw["optical_diag_norm"],
                physical_optical_offdiag_normalized=physical_raw["optical_offdiag_norm"],
                physical_constitutive=front,
                physical_constitutive_relative=quotient(front, np.linalg.norm(state["j"])),
                wall_seconds=perf_counter()-start[0], cpu_seconds=process_time()-start[1],
                operations=dict(model.clock.counts-before))


def masked(model, y):
    mask = model.ports.valid_mask
    return np.where(mask, y, 0.) if mask is not None else y


def metric_record(model, J, Q, m, Y, count, stage, beta_material, beta_current, kind):
    beta = 1/beta_material
    upper = 20.
    raw, jj, jq = optical_terms(J, Q)
    sector = eigvalsh(helm(jj) - upper*helm(jq), check_finite=False)
    pos = np.maximum(sector, 0.)
    current_contract = 1/(m + 1/beta_current)
    c = J - m[:, None]*Q
    data_fit = masked(model, model.ports.receive(J)-Y)
    px_j, px_q, interaction, local_sector = pixel_arrays(J, Q, upper)
    a = upper * raw["optical_scale"] + raw["J_frob2"]
    sector_scale = beta**3/(beta+upper)
    return dict(step=count, stage=stage, point=kind, beta_engine_material=beta_material,
                beta_engine_current=beta_current, beta_user_material=beta,
                beta_user_current=1/beta_current, contrast_upper=upper,
                sector_eigenvalues=sector.tolist(), sector_positive_Fro=float(np.linalg.norm(pos)),
                sector_positive_trace=float(pos.sum()), sector_positive_rank=int(np.sum(sector > 1.e-9*a)),
                sector_positive_normalized=quotient(np.linalg.norm(pos), a),
                sector_pixel_positive_trace=float(local_sector.sum()),
                sector_pixel_positive_max=float(local_sector.max()),
                contraction_sector_positive_Fro=float(sector_scale*np.linalg.norm(pos)),
                contraction_sector_factor=sector_scale,
                contraction_optical_diag=beta*raw["optical_diag"],
                contraction_optical_offdiag=beta*raw["optical_offdiag"],
                constitutive_frob=float(np.linalg.norm(c)),
                constitutive_relative=quotient(np.linalg.norm(c), np.linalg.norm(J)+np.linalg.norm(m[:, None]*Q)),
                contracted_current_constitutive_frob=float(np.linalg.norm(current_contract[:,None]*c)),
                aux_data_frob=float(np.linalg.norm(data_fit)),
                aux_data_relative=quotient(np.linalg.norm(data_fit), np.linalg.norm(masked(model, Y))),
                m_min=float(m.min()), m_peak=float(m.max()), m_rms=float(np.sqrt(np.mean(m*m))),
                **raw)


def snapshot(path, model, m, J, Q, count, stage):
    a, b, c, sector = pixel_arrays(J, Q, 20.)
    residual = J - m[:, None] * Q
    diagonal_absolute = np.sum(abs((J.conj()*Q).imag), axis=1)
    total_defect = (Q-model.E)-model.green.apply(J)
    cross = J.conj()[:, :, None]*Q[:, None, :]
    optics = (cross - cross.conj().transpose(0, 2, 1))/(2j)
    diag = optics[:, np.arange(J.shape[1]), np.arange(J.shape[1])]
    offdiag = np.maximum(np.sum(abs(optics)**2, axis=(1,2)) - np.sum(abs(diag)**2, axis=1), 0.)
    np.savez_compressed(path, m=m, J=J, E=Q, step=count, stage=stage, incident_B=model.E,
                        current_energy=a, field_energy=b, optical_pixel_trace=c.imag,
                        optical_pixel_diagonal_abs=diagonal_absolute,
                        optical_pixel_offdiag_fro=np.sqrt(offdiag), sector_pixel_positive=sector,
                        constitutive_pixel=np.linalg.norm(residual, axis=1),
                        graph_relative_defect=np.linalg.norm(total_defect)/max(np.linalg.norm(Q), 1.e-300))


def run_case(name, limit, destination):
    case, public = get_case(name)
    Y = public["Y"]
    fold = destination / name
    fold.mkdir(parents=True, exist_ok=True)
    modes, betas, rank = engine.schedule(case)
    clock = EngineClock()
    all_start = (perf_counter(), process_time())
    model = engine.legacy.PaperModel(public, case, clock)
    audit = engine.legacy.PaperModel(public, case, engine.legacy.Clock(float("inf")))
    blank = np.zeros(model.n)
    history = []
    saved = []
    audited = []
    ref_materials = []
    mark_stage = -1
    final_step = 0
    ref = SRAP / f"results/predictions/{name}__dual.npz"
    with clock.pause():
        geom = {key: value for key, value in public.items() if key != "case_json"}
        geom["spectrum"] = model.green.spectrum
        np.savez_compressed(fold/"geometry.npz", **geom)
        audited.append(audit_record(audit, blank, Y, 0, -1, "air_physical_initial"))

    archived = {}
    if ref.is_file():
        with clock.pause():
            reference = read_npz(ref)
            if "traj_step" in reference:
                archived = {int(step): reference["traj_material"][i]
                            for i, step in enumerate(reference["traj_step"])}
    archive_matches = []
    with (fold / "updates.jsonl").open("w") as log:
        gen = engine.run_cie_arm(model, Y, case, "dual", blank.copy(), float(public["delta"]))
        try:
            for m, step, stage, event in gen:
                frame = gen.gi_frame.f_locals
                J, Q = frame["J"], frame["Q"]       # read-only access to actual yielded iterate
                bm, bj = float(frame["bm"]), float(frame["bj"])
                current_engine_wall = clock.wall()
                current_engine_cpu = clock.cpu()
                final_step = int(step)
                tick = clock.counts["PR_CG_steps"]
                if tick != step:
                    raise ValueError("Tracker missed an actual per-source update.")
                start_stage = stage != mark_stage
                with clock.pause():
                    if step == 1:
                        dpcs = metric_record(model, frame["JD"], frame["QD"], blank,
                                             Y, 0, stage, bm, bj, "dpcs_start")
                        log.write(json.dumps(dpcs, allow_nan=False) + "\n")
                        snapshot(fold/"step_0000_dpcs.npz", audit, blank, frame["JD"], frame["QD"], 0, stage)
                        print(name, "DPCS", {k: dpcs[k] for k in ("optical_diag","optical_offdiag",
                              "optical_diag_rms","optical_offdiag_rms","optical_diag_norm",
                              "optical_offdiag_norm","aux_data_relative")}, flush=True)
                    record = metric_record(model, J, Q, m, Y, step, stage, bm, bj, "update")
                    record["engine_wall_seconds"] = current_engine_wall
                    record["engine_cpu_seconds"] = current_engine_cpu
                    record["stage_last_update"] = bool(event)
                    record["operations"] = dict(clock.counts)
                    log.write(json.dumps(record, allow_nan=False) + "\n")
                    ref_materials.append(m.copy())
                    if step in archived:
                        difference = float(np.max(abs(m - archived[step])))
                        archive_matches.append(dict(step=step, linf=difference))
                    do_snapshot = (step in (1, 4, 16, 64, 128, 256, 512, 1024)
                                   or start_stage or event or step % 160 == 0 or current_engine_wall >= limit)
                    if do_snapshot:
                        filename = f"step_{step:04d}.npz"
                        snapshot(fold/filename, audit, m, J, Q, step, stage)
                        saved.append(dict(step=step, stage=stage, path=filename))
                    do_audit = (step in (1, 16, 64, 128, 256, 512, 1024) or start_stage
                                or event or step % 160 == 0 or current_engine_wall >= limit)
                    if do_audit:
                        audited.append(audit_record(audit, m, Y, step, stage, "material_update"))
                    history.append((int(step), int(stage), float(current_engine_wall)))
                    if step in (1, 16, 64) or start_stage:
                        print(name, "step", step, "stage", stage, "diag/off",
                              record["optical_diag"], record["optical_offdiag"],
                              "peak", record["m_peak"], "fit", record["aux_data_relative"],
                              flush=True)
                mark_stage = stage
                if current_engine_wall >= limit:
                    status = "SOLVER_WALL_CUT_AFTER_COMPLETE_UPDATE"
                    break
            else:
                status = "SCHEDULE_COMPLETE"
            if final_step == 0:
                raise ValueError("Engine returned no actual material update.")
            with clock.pause():
                if not saved or saved[-1]["step"] != final_step:
                    filename = f"step_{final_step:04d}.npz"
                    snapshot(fold/filename, audit, m, J, Q, final_step, stage)
                    saved.append(dict(step=final_step, stage=stage, path=filename))
                if audited[-1]["step"] != final_step:
                    audited.append(audit_record(audit, m, Y, final_step, stage, "material_update"))
        finally:
            gen.close()
    with clock.pause():
        np.savez_compressed(fold/"materials.npz", materials=np.asarray(ref_materials),
                            step=np.array([i for i, _, t in history]),
                            stage=np.array([i for _, i, t in history]),
                            engine_wall=np.array([t for _, _, t in history]))
    result = dict(case=case, engine_source=str((SRAP/"engine.py").relative_to(RESEARCH)),
                  method="dual", measurement_shape=list(Y.shape), raw_noise_radius=float(public["delta"]),
                  max_engine_wall_seconds=limit, schedule_modes=modes, schedule_betas=betas,
                  schedule_requested_rank=rank, solver_events=model.trace_route,
                  stop=status, final_step=final_step, optimizer_wall_seconds=current_engine_wall,
                  optimizer_cpu_seconds=current_engine_cpu, clock_split_wall_seconds=clock.paused_wall,
                  clock_split_cpu_seconds=clock.paused_cpu,
                  measured_full_wall_seconds=perf_counter()-all_start[0],
                  measured_full_cpu_seconds=process_time()-all_start[1],
                  operations=dict(clock.counts), audit_counts=dict(audit.clock.counts),
                  accepted_state_audits=audited, material_snapshots=saved,
                  dpcs_snapshot="step_0000_dpcs.npz", archived_step_matches=archive_matches,
                  scope="Existing exact coarse inverse VIE and original training channels",
                  truth_access="NONE", sealed_utc=datetime.now(timezone.utc).isoformat())
    write(fold/"receipt.json", result)
    print(name, status, final_step, "physical final residual", audited[-1]["fit_relative"],
          "wall", result["optimizer_wall_seconds"], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", choices=CASES)
    parser.add_argument("--solver-wall", default=30., type=float)
    parser.add_argument("--output", default=HERE/"trajectories", type=Path)
    args = parser.parse_args()
    with threadpool_limits(1):
        for name in args.case or CASES:
            run_case(name, args.solver_wall, args.output)


if __name__ == "__main__":
    main()
