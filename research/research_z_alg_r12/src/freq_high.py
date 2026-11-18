"""High-contrast shared-material, shared-gain comparison on archived TM geometry."""
from argparse import ArgumentParser
from collections import Counter
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from threadpoolctl import threadpool_limits

from high_endpoints import Model, PAIRED, R11, R12, Regularizer, bounds, historical_path, load_endpoint
from research.research_z_alg_r4_diversity.src.physics import build_operators

OFFSET = .10
PLANS = {"repeat_old": (0.,0.), "negative": (-OFFSET,-OFFSET),
         "positive": (OFFSET,OFFSET), "symmetric": (-OFFSET,OFFSET)}
GAIN_SIGMA = .002
GAIN_ACTUAL = 1. + .0013*np.exp(.63j)
HALF_POWER_AMPLITUDE = np.sqrt(.5)
OUT = R12/"results"/"freq_high"


def power_operators(public, delta):
    k0 = float(public["k"])
    geometry = {k: public[k] for k in ("centers","sources","receivers","radius")}
    operators = build_operators(**geometry,k=k0*(1+delta))
    # Each independent external line source has free-background P=|C|²/(8 Z k).
    operators["E0"] *= np.sqrt(1+delta)
    return operators


def forward_complex(model, material):
    return model.forward(material)["y"].T.reshape(-1)


def old_candidates(name):
    endpoint, public, start = load_endpoint(name)
    path = R11/f"{endpoint.case_id}.npz"
    with np.load(path) as data:
        old_path = data["m"]
        candidate = old_path[-1] if old_path.ndim == 2 else old_path[0,-1]
    return endpoint, public, np.array([start,candidate])


def predict_gap(a,b):
    def profiled(x,y):
        return np.linalg.norm(y-x*(np.vdot(x,y)/np.vdot(x,x)))
    return float(min(profiled(a,b),profiled(b,a)))


def design(name, public, candidates, models):
    sigma = float(public["sigma"])
    responses = {d: [forward_complex(m,x) for x in candidates] for d,m in models.items()}
    w = sigma/GAIN_SIGMA
    old_axes = [np.r_[y,w] for y in responses[0.]]
    old_gap = predict_gap(*old_axes)
    entries = {}
    for plan, ds in PLANS.items():
        yn = [np.concatenate([HALF_POWER_AMPLITUDE*responses[d][i] for d in ds])
              for i in range(2)]
        joined = [np.r_[responses[0.][i],yn[i],w] for i in range(2)]
        gap = predict_gap(*joined)
        entries[plan] = dict(raw_new_gap=float(np.linalg.norm(yn[1]-yn[0])),
                             joint_profile_gap=gap,
                             novelty_L2=float(np.sqrt(max(0., gap**2-old_gap**2))))
    candidate_plan = max(("negative","positive","symmetric"),key=lambda key:entries[key]["novelty_L2"])
    return dict(case=name, selected=candidate_plan,
                candidate_pair_material_RMS=float(np.linalg.norm(candidates[1]-candidates[0])/np.sqrt(candidates.shape[1])),
                old_profile_gap=old_gap, plans=entries, sigma=sigma,
                old_data_residual_L2=[float(np.linalg.norm(y-public["Y"].T.reshape(-1))) for y in responses[0.]],
                design_RHS=sum(m.counts["actual_complex_RHS"] for m in models.values()))


def acquire(name):
    endpoint, public, candidates = old_candidates(name)
    rebuilt = power_operators(public,0.)
    relative = {k:float(np.linalg.norm(rebuilt[k]-public[k])/np.linalg.norm(public[k]))
                for k in ("G","E0","S")}
    models = {d:Model(public if not d else power_operators(public,d)) for d in (-OFFSET,0.,OFFSET)}
    metadata = design(name, public, candidates, models)
    metadata.update(case_id=endpoint.case_id, k0=float(public["k"]), relative_k0_transfer_differences=relative,
                    archive=historical_path(R11/f"{endpoint.case_id}.npz"),
                    gain_calibration_sigma=GAIN_SIGMA,
                    gain_actual=[GAIN_ACTUAL.real,GAIN_ACTUAL.imag])
    (OUT/f"{name}_design.json").write_text(json.dumps(metadata,indent=2)+"\n")
    sigma = float(public["sigma"])
    # This is the only truth read in acquisition. It is absent from the inverse solver.
    with np.load(PAIRED/"offline_truth"/"A"/f"{endpoint.case_id}.npz") as source:
        truth = source["m"]
    true_response = {d:forward_complex(m,truth) for d,m in models.items()}
    rng = np.random.default_rng(20261119+int(endpoint.case_id[-4:]))
    shape = (2, public["Y"].size)
    noise = sigma/np.sqrt(2)*(rng.normal(size=shape)+1j*rng.normal(size=shape))
    # Keep the archived first exposure (including its original noise), with the same
    # simulated receiver gain as the new exposures. Paired arms reuse the draws by index.
    data = {"old":GAIN_ACTUAL*public["Y"].T.reshape(-1), "new_noise":noise}
    for plan, offsets in PLANS.items():
        data[plan] = np.array([GAIN_ACTUAL*HALF_POWER_AMPLITUDE*true_response[d]+noise[i]
                               for i,d in enumerate(offsets)])
    np.savez_compressed(OUT/f"{name}_acquisition.npz", **data)
    print(name,"acquired",metadata["selected"],{k:round(v["novelty_L2"],5) for k,v in metadata["plans"].items()},flush=True)


class HighFit:
    def __init__(self, models, blocks, reg, box, sigma, scale):
        self.models, self.blocks, self.reg, self.box = models,blocks,reg,box
        self.unique = tuple(sorted({d for d,_,_ in blocks}))
        self.observed = np.concatenate([y for _,_,y in blocks])
        self.sqrt_lambda = sigma/GAIN_SIGMA
        self.scale = scale
        self.cache = {}
        self.count = Counter()

    def cost(self, jacobian):
        model = self.models[self.unique[0]]
        return len(self.unique)*(model.ns + (model.nr if jacobian else 0))

    def evaluate(self, material, jacobian=False):
        key = material.tobytes()
        last = self.cache.get(key)
        if last and (last[1] is not None or not jacobian):
            return last
        outputs = {}
        for d in self.unique:
            model = self.models[d]
            count = model.counts.copy()
            if model.state is None or not np.array_equal(model.state["m"],material):
                model.forward(material)
            if jacobian:
                model.jacobian()
            self.count.update(model.counts-count)
            state = model.state
            outputs[d] = state["y"].T.reshape(-1), state.get("Jcomplex")
        y = np.concatenate([w*outputs[d][0] for d,w,_ in self.blocks])
        lam = self.sqrt_lambda**2
        norm_y2 = np.vdot(y,y).real+lam
        gain = (np.vdot(y,self.observed)+lam)/norm_y2
        resid = np.r_[gain*y-self.observed,self.sqrt_lambda*(gain-1)]/self.scale
        j = None
        if jacobian:
            columns = np.vstack([w*outputs[d][1] for d,w,_ in self.blocks])
            gp = (columns.conj().T@self.observed - gain*2*np.real(y.conj()@columns))/norm_y2
            dc = np.vstack([gain*columns+y[:,None]*gp,self.sqrt_lambda*gp])/self.scale
            j = np.vstack([dc.real,dc.imag])
        regularizer = self.reg.evaluate(material)
        misfit = np.vdot(resid,resid).real/2
        value = float(misfit+regularizer[0])
        answer = (value,j,np.r_[resid.real,resid.imag],regularizer,gain,
                  float(misfit),float(np.linalg.norm(resid[:-1])*self.scale))
        if len(self.cache) > 12:
            self.cache.clear()
        self.cache[key] = answer
        return answer

    def solve(self, start, rhs_budget):
        tic = perf_counter()
        count0 = self.count.copy()
        m = start.copy()
        trace = []
        damping = 1e-3
        for _ in range(100):
            if self.count["actual_complex_RHS"]-count0["actual_complex_RHS"]+self.cost(True) > rhs_budget:
                break
            loss,j,r,prior,_,_,_ = self.evaluate(m,True)
            grad = j.T@r+prior[1]
            low,high = self.box
            face = ((m<=low+1e-9)&(grad>0)) | ((m>=high-1e-9)&(grad<0))
            projected = np.where(face,0.,grad)
            trace.append([self.count["actual_complex_RHS"]-count0["actual_complex_RHS"],
                          loss,float(np.linalg.norm(projected,np.inf)),int(np.sum(face))])
            if np.linalg.norm(projected,np.inf)<1e-10:
                break
            free = ~face
            jf = j[:,free]
            D = self.reg.D[:,free]
            curvature = jf.T@jf+(D.T@D.multiply(prior[2][:,None])).toarray()
            diagonal = np.maximum(np.diag(curvature),1e-8)
            accepted = False
            for _ in range(7):
                if self.count["actual_complex_RHS"]-count0["actual_complex_RHS"]+self.cost(False) > rhs_budget:
                    break
                matrix = curvature + np.diag(damping*diagonal)
                direction = cho_solve(cho_factor(matrix, lower=True, check_finite=False),
                                      -grad[free], check_finite=False)
                trial = m.copy()
                trial[free] = np.clip(m[free]+direction,low,high)
                actual = self.evaluate(trial)[0]
                if actual < loss:
                    m = trial
                    accepted = True
                    damping = max(1e-7, damping*.35)
                    break
                damping *= 8
            if not accepted:
                break
        value, _, _, _, gain, data_fit, raw_residual = self.evaluate(m)
        setup = self.count-count0
        return m, dict(value=value,data_loss=data_fit,data_residual_L2=raw_residual,
                       gain_estimate=[float(gain.real),float(gain.imag)],
                       counts=dict(setup), RHS_budget=rhs_budget,line_trace=trace,
                       elapsed=perf_counter()-tic)


def solve(name, rhs_budget):
    endpoint, public, candidates = old_candidates(name)
    with np.load(OUT/f"{name}_acquisition.npz") as source:
        acquired = {k:source[k] for k in source.files}
    metadata = json.loads((OUT/f"{name}_design.json").read_text())
    box = bounds(public)
    sigma = float(public["sigma"])
    scale = np.linalg.norm(acquired["old"])
    starts = {"bad":candidates[0],"old_polished":candidates[1],
              "cold":np.full(candidates.shape[1],np.clip(1.,*box))}
    regularizer = Regularizer(public["grid_shape"],smooth=1e-4,tv=1e-4,epsilon=.01)
    arrays = dict(starts)
    report = {"selected":metadata["selected"],"arms":{}}

    def create(blocks):
        # Fresh states make actual RHS counts comparable across initializations.
        return HighFit({d:Model(operators[d]) for d in set([0.]+[row[0] for row in blocks])},
                       blocks,regularizer,box,sigma,scale)

    operators = {d:public if not d else power_operators(public,d) for d in (-OFFSET,0.,OFFSET)}
    for plan, offsets in PLANS.items():
        blocks = [(0.,1.,acquired["old"])] + [
            (d,HALF_POWER_AMPLITUDE,acquired[plan][i]) for i,d in enumerate(offsets)]
        arm = {}
        finals = {}
        for start in ("bad","cold","old_polished"):
            fit = create(blocks)
            finals[start],arm[start] = fit.solve(starts[start],rhs_budget)
        for start in ("bad","cold","old_polished"):
            fit = create(blocks)
            finals[f"{start}_half"],arm[f"{start}_half"] = fit.solve(starts[start],rhs_budget//2)
        for rule,pool in (("ordinary_two_starts",("bad_half","cold_half")),
                          ("old_pair_two_starts",("bad_half","old_polished_half"))):
            choice = min(pool,key=lambda start:arm[start]["value"])
            finals[rule] = finals[choice]
            counts = Counter()
            for k in pool:
                counts.update(arm[k]["counts"])
            arm[rule] = dict(value=arm[choice]["value"],data_loss=arm[choice]["data_loss"],
                             data_residual_L2=arm[choice]["data_residual_L2"],counts=dict(counts),
                             RHS_budget=rhs_budget,selected_start=choice)
        if plan == "symmetric":
            for start in ("bad","cold"):
                stage1 = create(blocks[1:2])
                phase1,cost1 = stage1.solve(starts[start],rhs_budget//3)
                stage2 = create(blocks[:2])
                phase2,cost2 = stage2.solve(phase1,rhs_budget//3)
                stage3 = create(blocks)
                result,cost3 = stage3.solve(phase2,rhs_budget-2*(rhs_budget//3))
                finals[f"cumulative_{start}"] = result
                counts = Counter(cost1["counts"])+Counter(cost2["counts"])+Counter(cost3["counts"])
                arm[f"cumulative_{start}"] = dict(value=cost3["value"],data_loss=cost3["data_loss"],
                                                 data_residual_L2=cost3["data_residual_L2"],
                                                 counts=dict(counts),RHS_budget=rhs_budget,
                                                 stages=[cost1,cost2,cost3],order=[-OFFSET,0.,OFFSET])
        report["arms"][plan] = arm
        arrays.update({f"{plan}__{key}":m for key,m in finals.items()})
        print(name,plan,{k:round(arm[k]["value"],7) for k in (
            "bad","cold","old_polished","ordinary_two_starts","old_pair_two_starts")},flush=True)
        (OUT/f"{name}_solves.json").write_text(json.dumps(report,indent=2)+"\n")
        np.savez_compressed(OUT/f"{name}_fits.npz", **arrays)


def score(name):
    endpoint, _, _ = old_candidates(name)
    with np.load(PAIRED/"offline_truth"/"A"/f"{endpoint.case_id}.npz") as data:
        truth = data["m"]
    with np.load(OUT/f"{name}_fits.npz") as fits:
        truth_norm = np.linalg.norm(truth)
        errors = {k:dict(relative_L2=float(np.linalg.norm(fits[k]-truth)/truth_norm),
                         RMS=float(np.linalg.norm(fits[k]-truth)/np.sqrt(truth.size))) for k in fits.files}
    (OUT/f"{name}_scores.json").write_text(json.dumps(errors,indent=2)+"\n")


WIDE = (-.50, -.25)


def acquire_wide(name):
    endpoint, public, _ = old_candidates(name)
    with np.load(OUT/f"{name}_acquisition.npz") as primary:
        old, noise = primary["old"], primary["new_noise"]
    with np.load(PAIRED/"offline_truth"/"A"/f"{endpoint.case_id}.npz") as source:
        truth = source["m"]
    observations = []
    for i, delta in enumerate(WIDE):
        model = Model(power_operators(public,delta))
        observations.append(
            GAIN_ACTUAL * HALF_POWER_AMPLITUDE * forward_complex(model,truth) + noise[i])
    np.savez_compressed(OUT/f"{name}_wide_acquisition.npz",old=old,wide=observations,new_noise=noise)
    print(name,"wide band acquired",flush=True)


def solve_wide(name, rhs_budget):
    _, public, candidates = old_candidates(name)
    sigma = float(public["sigma"])
    box = bounds(public)
    regularizer = Regularizer(public["grid_shape"],smooth=1e-4,tv=1e-4,epsilon=.01)
    with np.load(OUT/f"{name}_acquisition.npz") as source:
        narrow = source["symmetric"]
    with np.load(OUT/f"{name}_wide_acquisition.npz") as source:
        old, wide = source["old"], source["wide"]
    scale = np.linalg.norm(old)
    operators = {d:public if not d else power_operators(public,d) for d in (*WIDE,-OFFSET,0.,OFFSET)}

    def new_fit(blocks):
        return HighFit({d:Model(operators[d]) for d in {x[0] for x in blocks}},
                       blocks,regularizer,box,sigma,scale)

    starts = {"bad": candidates[0], "cold": np.full(candidates.shape[1],np.clip(1.,*box))}
    wide_blocks = [(d,HALF_POWER_AMPLITUDE,wide[i]) for i,d in enumerate(WIDE)] + [(0.,1.,old)]
    narrow_blocks = [(-OFFSET,HALF_POWER_AMPLITUDE,narrow[0]),(0.,1.,old),
                     (OFFSET,HALF_POWER_AMPLITUDE,narrow[1])]
    report, finals = {}, {}
    for label,start in starts.items():
        vals = new_fit(wide_blocks)
        finals[f"wide_direct_{label}"],report[f"wide_direct_{label}"] = vals.solve(start,rhs_budget)
        for band, blocks in (("wide",wide_blocks),("narrow",narrow_blocks)):
            key = f"{band}_cumulative_{label}"
            m = start
            costs = []
            for i,budget in enumerate((rhs_budget//3,rhs_budget//3,rhs_budget-2*(rhs_budget//3)),1):
                fit = new_fit(blocks[:i])
                m, stats = fit.solve(m,budget)
                costs.append(stats)
            counts = Counter()
            for stage in costs:
                counts.update(stage["counts"])
            report[key] = dict(value=costs[-1]["value"],data_loss=costs[-1]["data_loss"],
                               data_residual_L2=costs[-1]["data_residual_L2"],
                               counts=dict(counts), RHS_budget=rhs_budget, stages=costs,
                               order=[x[0] for x in blocks])
            finals[key] = m
        print(name,rhs_budget,label,
              {k:round(v["value"],7) for k,v in report.items() if k.endswith(label)},flush=True)
    (OUT/f"{name}_wide_{rhs_budget}_solves.json").write_text(json.dumps(report,indent=2)+"\n")
    np.savez_compressed(OUT/f"{name}_wide_{rhs_budget}_fits.npz",**finals)


def score_wide(name,rhs_budget):
    endpoint, _, _ = old_candidates(name)
    with np.load(PAIRED/"offline_truth"/"A"/f"{endpoint.case_id}.npz") as source:
        truth = source["m"]
    with np.load(OUT/f"{name}_wide_{rhs_budget}_fits.npz") as fits:
        norm = np.linalg.norm(truth)
        scores = {k:float(np.linalg.norm(fits[k]-truth)/norm) for k in fits.files}
    (OUT/f"{name}_wide_{rhs_budget}_scores.json").write_text(json.dumps(scores,indent=2)+"\n")
    print(name,rhs_budget,{k:round(v,3) for k,v in scores.items()},flush=True)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["acquire","solve","score","all","wide-acquire","wide-solve","wide-score","wide-all"])
    parser.add_argument("cases", nargs="+", choices=["ring","shell","bridge"])
    parser.add_argument("--rhs-budget",type=int,default=1440)
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    with threadpool_limits(limits=1):
        for name in args.cases:
            if args.stage in ("acquire","all"):
                acquire(name)
            if args.stage in ("solve","all"):
                solve(name,args.rhs_budget)
            if args.stage in ("score","all"):
                score(name)
            if args.stage in ("wide-acquire","wide-all"):
                acquire_wide(name)
            if args.stage in ("wide-solve","wide-all"):
                solve_wide(name,args.rhs_budget)
            if args.stage in ("wide-score","wide-all"):
                score_wide(name,args.rhs_budget)
