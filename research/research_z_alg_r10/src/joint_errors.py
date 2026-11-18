"""At complete material vectors, compare separate and shared-noise cofactor tests."""

from collections import Counter
from itertools import combinations
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq

from bounded_certificate import feasible_witness

ROOT = Path(__file__).resolve().parents[1] / "results"
EPS, GAMMA = 2e-4, .01
BOX = (.5, 30.)
TARGET, MAX_STEP = 1.2, .5
PAIRS = list(combinations(range(4),2))


def scaled_spectrum(M):
    s = np.linalg.svd(M,compute_uv=False)
    return float(s[-1]/s[0])


class Case:
    def __init__(self, data, name):
        self.H, self.R = (data[f"{name}/{key}"] for key in ("H","R"))
        self.B, self.host, self.y = (data[f"{name}/active/{key}"] for key in ("B","y_host","y"))
        n,s = self.B.shape
        U, values, Vh = np.linalg.svd(self.B,full_matrices=True)
        self.N = U[:, s:].conj().T
        self.B_inv = (Vh.conj().T/values) @ U[:, :s].conj().T
        self.B_gap = float(values[-1]/values[0])
        self.pivots = data[f"{name}/pivots"].tolist()
        self.nonpivots = [i for i in range(n) if i not in self.pivots]

    def forward(self, q):
        L = np.diag(1/q)-self.H
        return self.host+self.R @ np.linalg.solve(L,self.B), L

    def profile(self, q):
        f,_ = self.forward(q)
        residual = self.y-f
        g = np.vdot(f,residual)/(np.vdot(f,f).real+(EPS/GAMMA)**2)
        return float(np.linalg.norm((residual-g*f)/EPS)**2+abs(g/GAMMA)**2)

    def evaluate(self, q, label):
        U = np.eye(len(q))-self.H*q
        W = self.N @ U
        Ks, material_predictions = [], []
        for i in self.nonpivots:
            indices = [*self.pivots,i]
            cols = W[:,indices]
            weights = np.array([(-1)**j*np.linalg.det(np.delete(cols,j,axis=1)) for j in range(len(indices))])
            Ks.append(U[:,indices] @ weights)
            material_predictions.append((q[indices]*self.R[indices]) @ weights)
        K = np.column_stack(Ks)
        A = self.B_inv @ K
        c = self.host @ A + np.asarray(material_predictions)
        r = self.y @ A-c
        design = np.vstack((GAMMA*c, EPS*A))
        col_scale = np.linalg.norm(design,axis=0)
        design /= col_scale
        r /= col_scale

        def qr_score(index):
            # (row normalized noise) @ design = residual, solve in the shared noise coordinates.
            T = np.linalg.qr(design[:,index],mode="reduced")[1]
            return float(np.linalg.norm(np.linalg.solve(T.conj().T,r[index].conj()))**2)

        f,L = self.forward(q)
        direct_system = np.column_stack((GAMMA*f,EPS*np.eye(len(f))))
        forward_score = float(np.linalg.norm(np.linalg.lstsq(direct_system,self.y-f,rcond=None)[0])**2)
        singles = (abs(r)**2).tolist()
        pair_scores = [qr_score(list(pair)) for pair in PAIRS]
        full = qr_score(list(range(len(r))))
        maximum_single,maximum_pair = max(singles),max(pair_scores)
        gap_tol = 1e-10*max(1.,full)
        assert maximum_single <= maximum_pair+gap_tol and maximum_pair <= full+gap_tol
        discrepancy = abs(full-forward_score)/max(1.,forward_score)
        assert discrepancy < 1e-8, (q,discrepancy)
        return dict(label=label,q=q.tolist(), forward_score=forward_score, scalar_scores=singles,
                    pair_scores=pair_scores, maximum_scalar=maximum_single,maximum_pair=maximum_pair,
                    full_cofactor_score=full, full_vs_forward_relative=discrepancy,
                    score_category=("forward_feasible" if full<=1 else "scalar_reject" if maximum_single>1
                                    else "pair_reject_after_scalars_pass" if maximum_pair>1
                                    else "full_reject_after_pairs_pass"),
                    singular_ratios=dict(forward=scaled_spectrum(L),B=self.B_gap,
                                         pivot_system=scaled_spectrum(W[:,self.pivots]),
                                         A=scaled_spectrum(A),scaled_noise_design=scaled_spectrum(design),
                                         worst_scaled_pair_noise_design=min(scaled_spectrum(design[:,pair]) for pair in PAIRS)),
                    smallest_raw_noise_column=float(min(col_scale)),
                    synthesis_relative=float(np.linalg.norm(self.B@A-K)/np.linalg.norm(K)),
                    physical_prediction_relative=float(np.linalg.norm(f@A-c)/np.linalg.norm(c)))


def scan_case(data, name):
    case = Case(data,name)
    truth = data[f"{name}/truth"]
    others = data[f"{name}/clean_witnesses"]
    original = [case.evaluate(truth,"saved_truth"), case.evaluate(np.full(len(truth),sum(BOX)/2),"box_center")]
    original.extend(case.evaluate(q,f"old_witness_{j}") for j,q in enumerate(others))
    scans, misses = [], []
    for j in range(len(truth)):
        for sign in (-1,1):
            limit = min(MAX_STEP,truth[j]-BOX[0] if sign<0 else BOX[1]-truth[j])
            direction = np.eye(len(truth))[j]*sign
            start, end = case.profile(truth), case.profile(truth+limit*direction)
            if limit<=0 or not start<TARGET<=end:
                misses.append(dict(axis=j,sign=sign,limit=limit,start_score=start,end_score=end))
                continue
            step = brentq(lambda t:case.profile(truth+t*direction)-TARGET,0,limit,xtol=1e-12)
            rec = case.evaluate(truth+step*direction,"truth_coordinate_path")
            rec.update(axis=j,sign=sign,step=step,score_at_end=end)
            scans.append(rec)
    return dict(case=name,pivots=case.pivots,nonpivots=case.nonpivots, baseline_points=original,
                diagnostic_points=scans,unbracketed_paths=misses),case


def main():
    start = time.process_time()
    reports,examples = {},{}
    for cohort in ("natural","confirm"):
        with np.load(ROOT/f"{cohort}_inputs.npz") as data:
            for name in data["case_ids"].tolist():
                result,case = scan_case(data,name)
                result["cohort"] = cohort
                reports[name] = result
                for p in result["diagnostic_points"]:
                    category = p["score_category"]
                    if category in ("pair_reject_after_scalars_pass","full_reject_after_pairs_pass") and category not in examples:
                        witness = feasible_witness(np.array(p["q"]),case.H,case.R,case.B,case.host,case.y)
                        assert witness["outside"], witness
                        examples[category] = dict(case=name, pivots=case.pivots, nonpivots=case.nonpivots,
                                                  strict_arb_forward=witness, point=p)
    all_points = [p for c in reports.values() for key in ("diagnostic_points","baseline_points") for p in c[key]]
    summary = {}
    for cohort in ("natural","confirm"):
        cases = [c for c in reports.values() if c["cohort"]==cohort]
        scans = [p for c in cases for p in c["diagnostic_points"]]
        summary[cohort] = dict(case_count=len(cases), diagnostic_points=len(scans),
                               unbracketed_paths=sum(len(c["unbracketed_paths"]) for c in cases),
                               diagnostic_categories=dict(Counter(p["score_category"] for p in scans)),
                               cases_with_scalar_pass_and_joint_rejection=sum(any(
                                   p["score_category"] in ("pair_reject_after_scalars_pass","full_reject_after_pairs_pass")
                                   for p in c["diagnostic_points"]) for c in cases))
    out = dict(noise_scales=[EPS,GAMMA],pair_column_indices=PAIRS,diagnostic_target_forward_score=TARGET,
               design="Use the stored old-data pivots. Enumerate all 12 signed coordinate directions at each saved truth. Search a bracket within a material displacement of 0.5 for forward score 1.2. This post-hoc diagnostic uses truth and active observations, and is not a blind experiment or a regional certificate.",
               summary=summary, checked_points=len(all_points),
               maximum_full_vs_forward_relative=max(p["full_vs_forward_relative"] for p in all_points),
               minimum_singular_ratios={key:min(p["singular_ratios"][key] for p in all_points)
                                        for key in all_points[0]["singular_ratios"]},
               examples=examples,cases=reports,cpu_seconds=time.process_time()-start)
    (ROOT/"joint_errors.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:out[k] for k in ("summary","checked_points","maximum_full_vs_forward_relative",
                                         "minimum_singular_ratios","examples")},indent=2))


if __name__ == "__main__":
    main()
