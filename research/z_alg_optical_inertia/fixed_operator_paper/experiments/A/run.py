"""Decision-only Phase A entrypoint. This module cannot read the sealed labels."""
import argparse
import itertools
import json
import resource
import sys
from time import perf_counter

from common import BGS, PAPER, PUBLIC, RATIOS, REGION, RESULTS, SEAL
from common import check_preregistration, digest, load, require_p0, save, write_json
import numpy as np
from threadpoolctl import threadpool_limits
from bounds import (Acquisition, BASE_WEIGHTS, LOCAL, OPTICAL, assemble, tail_proposal,
                    threshold_options, verify_inertia)
from baseline_neumann import neumann
from baseline_psd import psd_qcqp
from compressed import TraceBank, interval_search, local_update

sys.path.append(str(PAPER.parent / "ip_guided_v4/pool"))
from pool_common import FullWave

PARAMS = list(itertools.product(map(float, BASE_WEIGHTS), LOCAL, OPTICAL))
BUDGET = 10.0
METHODS = "ABCDE"


class ArrayStore:
    def __init__(self):
        self.arrays = {}
        self.number = 0

    def add(self, matrix):
        key = f"array_{self.number:06d}"
        self.arrays[key] = np.asarray(matrix)
        self.number += 1
        return key


def make_witness(acq, ell, upp, save_frame):
    t0 = perf_counter()
    full = FullWave(acq.G, acq.S, acq.B, acq.Y)
    center = (ell+upp)/2
    try:
        misfit = full.evaluate(center)[0]
        if np.linalg.norm(misfit-acq.Y) <= acq.radius:
            material, search = center, dict(stop_reason="center_in_noise_proposal")
        else:
            material, _, search = full.refine(center, ell, upp, 45)
        proposal_seconds = perf_counter()-t0
        strict_start = perf_counter()
        audit = acq.strict_witness(material, ell, upp)
        return dict(search=search, audit=audit,
                    material=save_frame(material) if audit["proved"] else None,
                    proposal_seconds=proposal_seconds, strict_seconds=perf_counter()-strict_start)
    except (ArithmeticError, ValueError, np.linalg.LinAlgError) as exc:
        return dict(audit=dict(proved=False, error=str(exc)), proposal_seconds=perf_counter()-t0,
                    strict_seconds=0.)


def numerical(acq, ell, upp, bank, reduced_bases=None):
    start = perf_counter()
    if reduced_bases is None:
        states = acq.n
        matrices = [(params, assemble(acq, ell, upp, *params)) for params in PARAMS]
    else:
        states = len(REGION)*2
        matrices = [(params, reduced_bases[params[0]]+local_update(acq, ell, upp, *params))
                    for params in PARAMS if params[0] in reduced_bases]
    values = [(params, border, tail_proposal(border, states, acq.radius))
              for params, border in matrices]
    return values, perf_counter()-start


def inertia(acq, ell, upp, values, save_frame, *, reduced=None, exact_bases=None,
            deadline=np.inf):
    start = perf_counter()
    errors, proof = [], None
    states = acq.n if reduced is None else len(REGION)*2
    ranked = sorted(values, key=lambda result: result[2]["strength"], reverse=True)
    for params, border, candidate in ranked[:4]:
        options = list(threshold_options(candidate, acq.radius))
        if not options:
            errors.append("tail_below_noise")
            continue
        if perf_counter() >= deadline:
            errors.append("wall_budget")
            break
        try:
            if reduced is None:
                exact = assemble(acq, ell, upp, *params, strict=True)
            else:
                if params[0] not in exact_bases:
                    errors.append("trace_base_could_not_be_certified")
                    continue
                exact = exact_bases[params[0]]+local_update(acq, ell, upp, *params, strict=True)
            for terms in options:
                proof = verify_inertia(exact, border, states, terms, acq.radius, save_frame)
                if proof is not None:
                    break
            if proof is not None:
                if reduced is not None:
                    proof["trace_bank_mu"] = params[0]
                    proof["positive_eliminated_states"] = acq.n - states
                proof["weights"] = list(params)
                proof["numeric_K_negative_count"] = candidate["nu"]
                break
            errors.append("one_sided_frame_check_failed")
        except (ArithmeticError, ZeroDivisionError, ValueError) as exc:
            errors.append(str(exc))
    best = ranked[0][2] if ranked else dict(strength=0., nu=None, K_gap=None)
    return dict(excluded=proof is not None, proof=proof, best_numeric_tail_over_r=best["strength"],
                best_numeric_K_negative=best["nu"], best_numeric_K_absolute_gap=best["K_gap"],
                failures=errors, strict_seconds=perf_counter()-start,
                proposed_weight_trials=len(values)), ranked


def chart_stage(acq, ell, upp, ranked, reduced, exact, save_frame, other_exclusion,
                witnessed, deadline):
    if other_exclusion:
        return dict(excluded=True, stage="inherited_D", strict_seconds=0., proposal_seconds=0.)
    if witnessed:
        return dict(excluded=False, stage="audited_physical_witness", strict_seconds=0., proposal_seconds=0.)
    start = perf_counter()
    attempts = []
    stages = {}
    proof_cost = 0.
    if not ranked:
        return dict(excluded=False, stage="no_trace_bank", strict_seconds=0., proposal_seconds=0.)
    # Fixed legal ranking of rows. Zero-width laws are least costly. Next take
    # the largest box-relative contribution to local constitutive uncertainty.
    l = ell[REGION]
    u = upp[REGION]
    order = np.argsort((u-l) / np.maximum(1., .5*(l+u)), kind="stable")
    for q in (1, 2, 4, 8):
        if perf_counter() >= deadline:
            stages[q] = "wall_budget"
            continue
        selected = tuple(map(int, order[:q]))
        # Best relaxed multiplier first plus one different-strength anchor.
        for params, _, _ in ranked[:2]:
            mu = params[0]
            if mu not in exact:
                stages[q] = "trace_bank_not_certified"
                continue
            elapsed = perf_counter()
            result = interval_search(acq, reduced[mu], exact[mu], ell, upp, selected,
                                     params, save_frame, deadline)
            result["total_wall_seconds"] = perf_counter()-elapsed
            proof_cost += result["elapsed_components"]["proof"]
            attempts.append(result)
            if result["excluded"]:
                stages[q] = "excluded"
                total = perf_counter()-start
                return dict(excluded=True, stage="added_physical_rows", strict_seconds=proof_cost,
                            proposal_seconds=max(0., total-proof_cost), levels=stages, attempts=attempts,
                            strict_distance_over_radius_lower=min(
                                leaf["verified"]["strict_distance_over_radius_lower"] for leaf in result["leaves"]),
                            rows_at_first_exclusion=q)
        stages[q] = stages.get(q, "not_excluded_at_cap")
    total = perf_counter()-start
    return dict(excluded=False, stage="no_complete_interval_cover", strict_seconds=proof_cost,
                proposal_seconds=max(0., total-proof_cost), levels=stages, attempts=attempts)


def run_pair(a, pair, global_bank):
    rows = []
    arrays = ArrayStore()
    pairdir = RESULTS / "raw"
    pairdir.mkdir(exist_ok=True)
    out = pairdir / f"pair_{pair+1:02d}.json"
    if out.exists():
        return
    start = perf_counter()
    preprocessed = []
    for bg_idx, cap in enumerate(BGS):
        for target in range(2):
            acq_start = perf_counter()
            acq = Acquisition(a["G"], a["S"], a["B"], a["Y"][pair, target, bg_idx],
                              float(a["radii"][pair, target, bg_idx]))
            setup = perf_counter()-acq_start
            if cap not in global_bank:
                global_bank[cap] = {mu: TraceBank(acq, cap, mu) for mu in BASE_WEIGHTS}
                preprocessed.extend(dict(cap=cap, mu=float(mu), seconds=bank.preprocess,
                                         verified=bank.available, detail=bank.verification)
                                    for mu, bank in global_bank[cap].items())
            banks = global_bank[cap]
            t0 = perf_counter()
            reduced = {mu: bank.base(acq, False) for mu, bank in banks.items() if bank.numeric_ready}
            reduce_base_numeric = perf_counter()-t0
            t0 = perf_counter()
            exact = {mu: bank.base(acq, True) for mu, bank in banks.items() if bank.available}
            reduce_base_ball = perf_counter()-t0
            for w, ratio in enumerate(RATIOS):
                for slot in range(2):
                    ell, upp = a["boxes"][pair, bg_idx, w, slot]
                    witness = make_witness(acq, ell, upp, arrays.add)
                    operation_zero = perf_counter()
                    C0, proposal_C = numerical(acq, ell, upp, banks)
                    C, crank = inertia(acq, ell, upp, C0, arrays.add, deadline=operation_zero+BUDGET)
                    C["proposal_seconds"] = proposal_C
                    operation_zero = perf_counter()
                    D0, proposal_D = numerical(acq, ell, upp, banks, reduced)
                    D, drank = inertia(acq, ell, upp, D0, arrays.add, reduced=banks,
                                       exact_bases=exact, deadline=operation_zero+BUDGET)
                    D["proposal_seconds"] = proposal_D
                    A = neumann(acq, ell, upp, arrays.add, perf_counter()+BUDGET)
                    B = psd_qcqp(acq, ell, upp, [(x, y, z["nu"]) for x,y,z in C0],
                                 arrays.add, perf_counter()+max(0., BUDGET-proposal_C))
                    B["reused_full_proposal_seconds"] = proposal_C
                    E = chart_stage(acq, ell, upp, drank, reduced, exact, arrays.add,
                                    D["excluded"], witness["audit"]["proved"], perf_counter()+BUDGET)
                    statuses = {}
                    rec = dict(zip(METHODS, (A,B,C,D,E), strict=True))
                    for method, value in rec.items():
                        if value["excluded"] and witness["audit"]["proved"]:
                            raise ArithmeticError(f"certificate conflicts with strict witness at {pair,target,bg_idx,w,slot}, {method}")
                        if value["excluded"]:
                            statuses[method] = "EXCLUDED"
                        elif witness["audit"]["proved"]:
                            statuses[method] = "WITNESS_FOUND"
                        else:
                            statuses[method] = "UNRESOLVED"
                    t_mark = [target, bg_idx, w, slot]
                    rows.append(dict(target_idx=target, bg_idx=bg_idx, width_idx=w, box_idx=slot,
                                     pair_id=str(a["pair_ids"][pair]),
                                     target_id=str(a["target_ids"][pair,target]),
                                     condition_id=f"a.{float(cap):.5f}.{float(ratio):.2f}.{slot}",
                                     background_cap=float(cap), width_ratio=float(ratio),
                                     box_positions=t_mark, candidate_peak_lower=float(max(ell[REGION])),
                                     candidate_peak_upper=float(max(upp[REGION])),
                                     status=statuses, witness=witness, methods=rec,
                                     eligible_C_D_max_numeric_difference=max(
                                         (abs(cval[2]["strength"]-dval[2]["strength"])
                                          for cval,dval in zip(C0,D0,strict=True)), default=None)
                                     if len(C0)==len(D0) else None,
                                     acquisition_setup_seconds=setup,
                                     reduced_base_numeric_seconds=reduce_base_numeric,
                                     reduced_base_ball_seconds=reduce_base_ball))
            print(f"pair {pair+1:02d}/24, background {cap}, acquisition {target}, "
                  f"{len(rows)} boxes", flush=True)
    save(pairdir / f"pair_{pair+1:02d}_arrays.npz", **arrays.arrays)
    write_json(out, dict(pair_index=pair, pair_id=str(a["pair_ids"][pair]),
                         elapsed_wall_seconds=perf_counter()-start, shared_trace_preprocess=preprocessed,
                         peak_resident_memory_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                         rows=rows))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--begin", type=int, default=1)
    parser.add_argument("--end", type=int, default=24)
    args = parser.parse_args()
    check_preregistration()
    require_p0()
    if not 1 <= args.begin <= args.end <= 24:
        raise ValueError("Invalid pair shard.")
    seal = json.loads(SEAL.read_text())
    if digest(PUBLIC) != seal["files"][PUBLIC.name]:
        raise ArithmeticError("The public decision input differs from the frozen acquisition.")
    a = load(PUBLIC)
    bank = {}
    with threadpool_limits(1):
        for pair in range(args.begin-1, args.end):
            run_pair(a, pair, bank)


if __name__ == "__main__":
    main()
