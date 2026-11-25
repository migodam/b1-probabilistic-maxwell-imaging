"""Paired hard-box sourcewise Dual-CIE on existing 24² scalar-TM observations."""

import argparse
from collections import Counter
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
from flint import ctx
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
POOL = HERE.parent
GUIDED = POOL.parent
PROJECT = GUIDED.parents[2]
SOURCE = PROJECT / "research"
sys.path[:0] = [str(POOL), str(SOURCE / "srap_cie_array_study_20261121"),
                str(SOURCE / "delegated/max_struct_r1/r13_worktree")]
import engine
from research.research_max_struct_r2.src.cie_global import _baseline_readout

from box_certificates import Bounds
from pool_common import FullWave, clock, elapsed, write_arrays, write_json
from structured_pool import BRANCHES, N_STARTS, RX, TX, filenames, proposals

spec = importlib.util.spec_from_file_location("guided_cie_threshold", GUIDED/"threshold/run.py")
threshold = importlib.util.module_from_spec(spec)
spec.loader.exec_module(threshold)
PROTOCOL = json.loads((HERE/"PROTOCOL.json").read_text())
DEVELOPMENT = POOL/"structured_results"
CONFIRMATION = GUIDED/"threshold/confirmation"
WEIGHTS = (1/32, 1/128)


class HardBoxCIE:
    """Reuse inherited DPCS, projected native Fourier maps, and material readout."""

    def __init__(self, g, s, b, y, radius, physical):
        tic = clock()
        self.g, self.s, self.b, self.y = g, s, b, y
        self.radius, self.physical = radius, physical
        self.history_clock = engine.legacy.Clock(np.inf)
        self.cost = self.history_clock.counts
        adapter = SimpleNamespace(E=b, ports=engine.legacy.Ports(s), ns=b.shape[1],
                                  clock=self.history_clock)
        self.jd, groups = engine.dpcs(adapter, y, PROTOCOL["dpcs_rank"], cap=False)
        self.qd = b + g@self.jd
        self.cost["matrix_green_forward_rhs"] += b.shape[1]
        self.dpcs_rank = [v.shape[1] for _, v in groups]
        self.bases = [engine.FourierCurrents((7,7), mf, groups, self.cost)
                      for mf in PROTOCOL["fourier_orders"]]
        self.weights_data = 1/np.maximum(np.sum(abs(y)**2, axis=0), 1e-300)
        self.weights_state = 1/np.maximum(np.sum(abs(b)**2, axis=0), 1e-300)
        self.auditor = Bounds(g, s, b, y)
        self.setup = dict(dpcs_rank=self.dpcs_rank, fourier_columns=[u.k for u in self.bases], **elapsed(tic))

    def refine(self, initial, lo, hi):
        begin = clock()
        base_counts = self.cost.copy()
        phys_counts = self.physical.counts.copy()
        start_audit = self.physical.evaluate(initial)
        j_physical = start_audit[1]
        m = initial.copy()
        best_residual = float(np.linalg.norm(start_audit[0]-self.y))
        best_material, best_data = m.copy(), start_audit[0]
        best_verified = None
        all_audits = [dict(stage=-1, step=0, residual=best_residual)]
        stages = []
        b_current = PROTOCOL["beta_current"]
        cfg = dict(cache_material_value_gradient=True)
        no_penalty = lambda z: (0., np.zeros_like(z))
        iteration = 0
        reason = "CIE stage schedule complete"
        current_start = None
        audit_time = dict(wall_seconds=0., cpu_seconds=0.)

        def verify_point(material):
            t = clock()
            result = self.auditor.point_upper(material)
            for key, value in elapsed(t).items():
                audit_time[key] += value
            self.cost["strict_physical_ball_LU"] += 1
            return result

        if best_residual <= self.radius:
            best_verified = verify_point(best_material)
        feasible = best_verified is not None and best_verified <= self.radius
        if feasible:
            reason = "start already strictly feasible"

        prior_coefficients = {}
        for stage, basis in enumerate(self.bases):
            if feasible:
                break
            target_beta = PROTOCOL["beta_material"][stage]
            tolerance = PROTOCOL["relative_coefficients_squared_tolerance"][stage]
            budget = PROTOCOL["stage_cycles"][stage]
            frequencies = [tuple(p) for p in basis.frequencies]
            if stage == 0:
                coefficient = basis.project_coefficients(j_physical-self.jd)
            else:
                coefficient = np.array([prior_coefficients.get(f, np.zeros(self.b.shape[1], complex))
                                        for f in frequencies])
            fluctuation = basis.apply(coefficient)
            j, q = self.jd+fluctuation, self.qd+self.g@fluctuation
            self.cost["matrix_green_forward_rhs"] += self.b.shape[1]
            if stage == 0:
                current_start = float(np.linalg.norm(j-j_physical)/max(np.linalg.norm(j_physical),1e-300))
            old_gradient, old_direction = None, None
            stage_reason = "cycle budget"
            for st in range(1, budget+1):
                contraction = (b_current/(1+b_current*m))[:,None]
                port_err = self.s@j-self.y
                constitutive = contraction*(j-m[:,None]*q)
                v = contraction*(constitutive*self.weights_state)
                gradient = basis.adjoint(self.s.conj().T@(port_err*self.weights_data)+
                                         v-self.g.conj().T@(m[:,None]*v))
                direction = -gradient
                if old_gradient is not None:
                    pr = np.maximum(0.,np.sum((gradient.conj()*(gradient-old_gradient)).real,axis=0) /
                                    np.maximum(np.sum(abs(old_gradient)**2, axis=0),1e-300))
                    direction += old_direction*pr
                    bad = np.sum((gradient.conj()*direction).real,axis=0)>=0.
                    direction[:,bad] = -gradient[:,bad]
                jdir = basis.apply(direction)
                qdir = self.g@jdir
                dd, dr = self.s@jdir, contraction*(jdir-m[:,None]*qdir)
                energy = self.weights_data*np.sum(abs(dd)**2,axis=0)+self.weights_state*np.sum(abs(dr)**2,axis=0)
                alpha = -np.sum((gradient.conj()*direction).real,axis=0)/np.maximum(energy,1e-300)
                change = direction*alpha
                delta2 = float(np.mean(np.sum(abs(change)**2,axis=0) /
                                       np.maximum(np.sum(abs(coefficient)**2,axis=0),1e-300)))
                coefficient += change
                j += jdir*alpha
                q += qdir*alpha
                m, _ = _baseline_readout(j,q,m,target_beta,self.weights_state,1.,(lo,hi),
                                         no_penalty,False,cfg,self.cost)
                if not (np.all(np.isfinite(m)) and np.all((lo<=m)&(m<=hi))):
                    raise ArithmeticError("The CIE readout escaped its declared branch.")
                old_gradient, old_direction = gradient,direction
                iteration += 1
                self.cost["PR_CG_steps"] += 1
                self.cost["matrix_green_forward_rhs"] += self.b.shape[1]
                self.cost["matrix_green_adjoint_rhs"] += self.b.shape[1]
                tolerance_end = st>=PROTOCOL["minimum_stage_cycles"] and delta2<tolerance
                if st % PROTOCOL["audit_cadence"] == 0 or st == budget or tolerance_end:
                    prediction, _, _ = self.physical.evaluate(m)
                    residue = float(np.linalg.norm(prediction-self.y))
                    all_audits.append(dict(stage=stage, step=iteration, residual=residue))
                    if residue<best_residual:
                        best_material,best_data,best_residual = m.copy(), prediction, residue
                        best_verified = None
                        if residue<=self.radius:
                            best_verified = verify_point(best_material)
                            if best_verified <= self.radius:
                                feasible = True
                                reason = stage_reason = "strict physical noise ball"
                                break
                if tolerance_end:
                    stage_reason = "relative squared coefficient tolerance"
                    break
            prior_coefficients = dict(zip(frequencies, coefficient, strict=True))
            stages.append(dict(stage=stage, cycles=st, stop=stage_reason, final_relative_coefficients_squared=delta2,
                               fourier_columns=basis.k, beta_material=target_beta,
                               stopping_squared_tolerance=tolerance))
        return best_material, best_data, dict(
            residual=best_residual, strict_upper=best_verified, feasible=feasible, stop=reason, stages=stages,
            current_initialization_relative_difference=current_start, original_audits=all_audits,
            physical_counts=dict(self.physical.counts-phys_counts),
            cie_counts=dict(self.cost-base_counts), strict_audit_time=audit_time, **elapsed(begin))


def get_inputs(cohort):
    groups = []
    if cohort in ("all","development"):
        groups.extend(("development", seed, name, DEVELOPMENT/"public"/f"{name}.npz",
                       DEVELOPMENT/"offline_input"/f"{name}.npz") for seed,name in filenames())
    if cohort in ("all","confirmation"):
        metadata = json.loads((CONFIRMATION/"public/cases.json").read_text())
        groups.extend(("confirmation", row["seed"], row["id"], CONFIRMATION/"public"/f"{row['id']}.npz",
                       CONFIRMATION/"scoring"/f"{row['id']}.npz") for row in metadata)
    return groups


def one_arm(method, g,s,b,y,radius, seed, physical_correction=False):
    start = clock()
    solver = FullWave(g,s,b,y)
    samples, lo, hi, prep = proposals(g,s,b,y,solver,seed)
    order = list(range(len(BRANCHES)))
    screening_rows = []
    tic = clock()
    if method == "point_to_fit":
        field = np.stack([solver.evaluate(m)[0] for m in samples])
        residues = np.linalg.norm((field-y).reshape(len(samples),-1),axis=1)
        best = residues.reshape(-1,N_STARTS).min(axis=1)
        order = np.argsort(best,kind="stable").tolist()
        screening_rows = best.tolist()
    elif method == "neumann_directional":
        cert = Bounds(g,s,b,y)
        screening_rows = [cert.neumann(l,h,radius) for l,h in zip(lo,hi,strict=True)]
        order = [i for i,row in enumerate(screening_rows) if row["directional_lower"]<=radius]
    elif method == "shifted_threshold":
        cert = threshold.Border(g,s,b)
        cert.measurement(y)
        screening_rows = [cert.check(l,h,y,radius,weights=WEIGHTS,jvalues=(1,))
                          for l,h in zip(lo,hi,strict=True)]
        order = [i for i,row in enumerate(screening_rows) if not row["rejected"]]
    screening = elapsed(tic)
    estimator = HardBoxCIE(g,s,b,y,radius,solver)
    rows, x_values, ports = {}, {}, {}
    winner = None
    for i in order:
        x, f, receipt = estimator.refine(samples[i*N_STARTS],lo[i],hi[i])
        if physical_correction and not receipt["feasible"]:
            x, f, correction = solver.refine(x,lo[i],hi[i],max_calls=42)
            receipt["pure_cie_residual"] = receipt["residual"]
            receipt["full_wave_correction"] = correction
            receipt["residual"] = correction["residual"]
            if correction["residual"]<=radius:
                t = clock()
                strict = estimator.auditor.point_upper(x)
                estimator.cost["strict_physical_ball_LU"] += 1
                receipt["strict_upper"] = strict
                receipt["feasible"] = bool(strict<=radius)
                receipt["corrected_audit_time"] = elapsed(t)
        rows[i],x_values[i],ports[i] = receipt,x,f
        if winner is None or receipt["residual"]<rows[winner]["residual"]:
            winner = i
        if method=="point_to_fit" and receipt["feasible"]:
            winner = i
            break
    endpoint = clock()
    upper = None
    if winner is not None:
        upper = rows[winner]["strict_upper"]
        if upper is None:
            upper = estimator.auditor.point_upper(x_values[winner])
            estimator.cost["strict_physical_ball_LU"] += 1
    endpoint_cost = elapsed(endpoint)
    serialized = dict(winner=winner, visited=list(rows), allowed=order, preparation=prep,
                      screening=screening, screening_rows=screening_rows, branch_records=rows,
                      setup=estimator.setup, estimated_winner_residual=rows[winner]["residual"] if winner is not None else None,
                      verified_winner_upper=upper, endpoint_extra_audit=endpoint_cost, radius=radius,
                      strict_data_feasible=bool(upper is not None and upper<=radius),
                      physical_counts=dict(solver.counts), cie_counts=dict(estimator.cost), **elapsed(start))
    data = dict(lower=lo, upper=hi, sample_center=samples[::N_STARTS],
                **{f"branch_{i}_material":x for i,x in x_values.items()},
                **{f"branch_{i}_train":f for i,f in ports.items()})
    if winner is not None:
        data.update(material=x_values[winner], train_prediction=ports[winner])
    return serialized, data


def online(out, cohort, physical_correction=False):
    ctx.prec,ctx.threads=112,1
    t=clock()
    if (out/"online.json").exists():
        raise FileExistsError(out/"online.json")
    with np.load(DEVELOPMENT/"public/geometry.npz",allow_pickle=False) as a:
        g,s,b = a["green"],a["sensor"][RX],a["incident"][:,TX]
    geometry_setup = elapsed(t)
    out.mkdir(parents=True,exist_ok=True)
    results, arrays = {}, {}
    stream = np.random.default_rng(610419)
    for tag, seed, key, published, _ in get_inputs(cohort):
        with np.load(published,allow_pickle=False) as case:
            y = case["train"] if tag=="development" else case["Y"]
            r = float(case["radius"])
        order = stream.permutation(PROTOCOL["arms"]).tolist()
        rows = {}
        for name in order:
            result, vals = one_arm(name,g,s,b,y,r,seed+1007,physical_correction)
            rows[name] = result
            arrays.update({f"{key}__{name}__{k}": v for k,v in vals.items()})
            print(f"{key} {name}: {len(result['visited'])}/16; winner/r "
                  f"{(result['verified_winner_upper'] or float('nan'))/r:.4f}; "
                  f"{result['wall_seconds']:.4f}s LU {result['physical_counts']['lu']}",
                  flush=True)
        results[key] = dict(cohort=tag, method_execution_order=order, methods=rows)
        write_json(out/"in_progress.json",dict(cases=results, geometry_read=geometry_setup))
    write_arrays(out/"predictions.npz",**arrays)
    write_json(out/"online.json",dict(cases=results, geometry_read=geometry_setup,
                                      corrected_after_cie=physical_correction,total=elapsed(t)))
    (out/"in_progress.json").unlink()


def score(out, cohort):
    ctx.prec,ctx.threads=112,1
    online_result = json.loads((out/"online.json").read_text())
    with np.load(DEVELOPMENT/"public/geometry.npz",allow_pickle=False) as data:
        g,s,b = data["green"],data["sensor"],data["incident"]
    truth_boxes = {}
    actual = {}
    recv = np.ix_(np.setdiff1d(np.arange(s.shape[0]),RX),TX)
    source = np.ix_(RX,np.setdiff1d(np.arange(b.shape[1]),TX))
    with np.load(out/"predictions.npz",allow_pickle=False) as arr, \
            np.load(CONFIRMATION/"predictions.npz",allow_pickle=False) as frozen_confirmed:
        for tag, _, key, public_case, private_case in get_inputs(cohort):
            with np.load(private_case,allow_pickle=False) as labels, np.load(public_case,allow_pickle=False) as public:
                material = labels["material"] if tag=="development" else labels["truth"]
                clean = labels["noiseless_ports"] if tag=="development" else labels["ports"]
                y = public["train"] if tag=="development" else public["Y"]
                radius = float(public["radius"])
            methods = online_result["cases"][key]["methods"]
            ref = FullWave(g,s,b,clean)
            cert = threshold.Border(g,s[RX],b[:,TX])
            cert.measurement(y)
            base = Bounds(g,s[RX],b[:,TX],y).point_upper(material)
            if base>radius:
                raise ArithmeticError("Known reference does not fit the finite contract.")
            entries, checks = {}, []
            for arm, rec in methods.items():
                base_key = f"{key}__{arm}"
                lo,hi = arr[base_key+"__lower"],arr[base_key+"__upper"]
                if tag=="development":
                    with np.load(DEVELOPMENT/f"{key}.npz",allow_pickle=False) as frozen:
                        assert np.array_equal(lo,frozen["lower"]) and np.array_equal(hi,frozen["upper"])
                else:
                    assert np.array_equal(lo,frozen_confirmed[key+"__lower"]) and np.array_equal(hi,frozen_confirmed[key+"__upper"])
                protected = np.flatnonzero(np.all((lo<=material)&(material<=hi),axis=1)).tolist()
                truth_boxes[key] = protected
                if arm=="shifted_threshold":
                    for i,row in enumerate(rec["screening_rows"]):
                        if row["rejected"]:
                            assert i not in protected
                            verify = cert.check(lo[i],hi[i],y,radius,verify=True,weights=WEIGHTS,jvalues=(1,))
                            assert verify["arb_congruence_lower"]>0
                            checks.append(dict(branch=i,**verify))
                if rec["winner"] is None:
                    entries[arm] = None
                    continue
                restored = arr[base_key+"__material"]
                predicted, _, _ = ref.evaluate(restored)
                actual[base_key+"__all_ports"] = predicted
                truth_error = float(np.linalg.norm(restored-material)/np.linalg.norm(material))
                entries[arm] = dict(full_relative_l2=truth_error, inside_relative_l2=truth_error,
                    background_rms=0., heldout_receivers=float(np.linalg.norm((predicted-clean)[recv])/np.linalg.norm(clean[recv])),
                    heldout_sources=float(np.linalg.norm((predicted-clean)[source])/np.linalg.norm(clean[source])),
                    strict_data_feasible=rec["strict_data_feasible"],
                    online_verified_winner_upper=rec["verified_winner_upper"],
                    re_solved_train_residual=float(np.linalg.norm(predicted[np.ix_(RX,TX)]-y)))
            ref_physical = dict(cie_counts=ref.counts)
            truth_boxes[key] = dict(contains_truth=protected, methods=entries, truth_verified_upper=base,
                                    radius=radius, arb_checked_shifted_exclusions=checks, scoring_counts=ref_physical)
    write_arrays(out/"offline_ports.npz",**actual)
    write_json(out/"offline.json",truth_boxes)


if __name__=="__main__":
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("phase",choices=("online","score"))
    args.add_argument("--cohort",default="all",choices=("development","confirmation","all"))
    args.add_argument("--out",default=HERE/"results",type=Path)
    args.add_argument("--physical-correction",action="store_true")
    a = args.parse_args()
    with threadpool_limits(limits=1):
        if a.phase=="online":
            online(a.out,a.cohort,a.physical_correction)
        else:
            score(a.out,a.cohort)
