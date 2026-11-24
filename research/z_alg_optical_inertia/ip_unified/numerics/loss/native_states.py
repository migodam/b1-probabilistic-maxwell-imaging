"""Screen legal native lossy state profiles and replay selected winners outward."""
from collections import defaultdict
from types import SimpleNamespace
import argparse
import json
import sys

from flint import arb, ctx
import numpy as np
from threadpoolctl import threadpool_limits

from run import HERE, STUDY, clock, elapsed, cap_coefficients, native_geometry, write_json

sys.path[:0] = [str(STUDY / "physical_revision/array"), str(STUDY / "submission/curve")]
from envelope import FourierEnvelope
from native import Geometry, Weights
from search import Projected
from verify_native import NativeState, analyze_row, check_state, reusable_stationarity


MU = np.logspace(-5., -2., 11)
ETA = (0., .025, .05, .1, .2, .4, .8, 1.6, 3.2, 8., 16., 32., 64.)
PHASE = (-1., -.5, 0., .2, .4, .65, .8, .9, 1., 1.1, 1.2, 1.5, 2.)


def legal_coefficients(cap, tangent, mode, weight):
    beta, alpha = cap_coefficients(cap, tangent)
    contraction = 1 - 1.e-11
    if mode == "old":
        return contraction / beta, contraction * weight / beta
    if mode == "shift":
        a = beta + weight*alpha
        return contraction / a, contraction * weight / a
    lam = weight / tangent
    return contraction * (1. - weight) / beta, contraction * lam / beta


def declaration(cap, tangent, mode, weights):
    beta = 1 / (arb(str(cap)) * (1 + arb(str(tangent))**2))
    alpha = arb(str(tangent)) * beta
    ua, ea = arb(weights.cap), arb(weights.optical[0])/arb(weights.scale)
    if mode == "old":
        slack = 1 - beta*ua
        valid = ua >= 0 and ea >= 0
    elif mode == "shift":
        slack = 1 - beta*ua - alpha*ea
        valid = ua >= 0 and ea >= 0
    else:
        sector = ua + arb(str(tangent))*ea
        slack = 1 - beta*sector
        valid = sector > 0
    if not valid or not slack > 0:
        raise ArithmeticError("Profile not enclosed in the declared loss contract.")
    return float(slack.lower())


def screen(geo, database, ys, *, resume):
    start = clock()
    proj = Projected(geo.spectrum, geo.s, geo.b, 7)
    output = HERE/"native_state_screen.json"
    existing = json.loads(output.read_text()) if resume and output.exists() else {}
    profiles = existing.get("profiles",[])
    already = {(row["tangent"],row["cap"]) for row in profiles}
    for t in sorted(set(row["tangent"] for row in database["observations"])):
        rows = [o for o in database["observations"] if o["tangent"] == t]
        data = [SimpleNamespace(y=ys[o["id"]]) for o in rows]
        radius = np.array([o["rho"] for o in rows])
        for cap in sorted(set(row["cap"] for row in database["receivers"] if row["tangent"] == t)):
            if (t,cap) in already:
                continue
            best = {o["id"]: {mode: dict(ratio=0.) for mode in ("old","shift","fixed_phase")} for o in rows}
            for mode, weights in (("old",ETA),("shift",ETA),("fixed_phase",PHASE)):
                for weight in weights:
                    ua, ea = legal_coefficients(cap,t,mode,weight)
                    proposals = proj.scalar(ua, ea, data, mu_grid=MU)
                    for mu, (nu, ds, first) in proposals.items():
                        for row, ratio, h in zip(rows, ds/radius,first,strict=True):
                            entry = best[row["id"]][mode]
                            if "nu" not in entry or ratio > entry["ratio"]:
                                entry.update(ratio=float(ratio), nu=int(nu), mu=float(mu), effective_cap=ua,
                                             effective_optical=ea, input_weight=weight, first_retained=float(h))
            profiles.append(dict(cap=cap,tangent=t,best=best))
            print("projected", t,cap,[(mode,round(max(r[mode]["ratio"] for r in best.values()),3))
                                       for mode in ("old","shift","fixed_phase")], flush=True)
            cost = elapsed(start)
            counts = dict(proj.cost)
            if existing:
                cost = {key: value + existing[key] for key,value in cost.items()}
                counts = {key: value + existing["counts"].get(key,0) for key,value in counts.items()}
            write_json(output,dict(profiles=profiles, dimension=proj.dimension, counts=counts, **cost))


def replay(geo, database, ys, *, focus):
    start = clock()
    src = json.loads((HERE/"native_state_screen.json").read_text())
    g = Geometry(geo.spectrum, geo.s, geo.b)
    envelope = FourierEnvelope(g)
    out = dict(cases=[], rule="Best projected candidates. A pass requires native interval strict lower margin.", audited_forward=None)
    radii = {item["id"]:item["rho"] for item in database["observations"]}
    if focus:
        profiles = [info for info in src["profiles"] if info["tangent"]==.2 and info["cap"]==1.]
        out["rule"] = "Largest old and shifted proposals at t=.2,u=1, with same observation radius."
    else:
        profiles = src["profiles"]
    for info in profiles:
        cap, t = info["cap"], info["tangent"]
        jobs = defaultdict(list)
        if focus:
            names = [name for name in info["best"] if not name.split("__")[0].endswith("_low")]
            selections = [(max(names,key=lambda name:info["best"][name][mode]["ratio"]),mode)
                          for mode in ("shift","old")]
        else:
            selections = [(name,mode) for name,modes in info["best"].items() for mode,values in modes.items()
                          if values["ratio"] >= .985 and not name.split("__")[0].endswith("_low")]
        for name,mode in selections:
            values = info["best"][name][mode]
            jobs[(mode, values["effective_cap"], values["effective_optical"], values["mu"], values["nu"])].append(name)
        for (mode, ua, ea, mu, nu), names in jobs.items():
            w = Weights.global_(ua,mu,ea)
            slack = declaration(cap,t,mode,w)
            engine = NativeState(g,w)
            try:
                state = check_state(engine,envelope,nu,equilibrate=True)
                item = dict(cap=cap,tangent=t,arm=mode,weights=w.as_dict(),physical_slack_lower=slack,state=state)
                if "exact_nu" in state:
                    s, defect, slots, stationarity = reusable_stationarity(engine,{geo.inputs:geo.b})
                    item["stationarity"] = stationarity
                    item["rows"] = [dict(id=name, result=analyze_row(ys[name],state,mu,s,defect,slots[geo.inputs],radii[name]))
                                    for name in names]
                out["cases"].append(item)
            except (ArithmeticError, np.linalg.LinAlgError) as ex:
                out["cases"].append(dict(cap=cap,tangent=t,arm=mode,weights=w.as_dict(),cases=names,error=str(ex)))
            print("replay",t,cap,mode,names,flush=True)
            write_json(HERE/"native_state_verified.json",{**out,"total_green_columns":dict(g.cost), **elapsed(start)})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay", action="store_true")
    parser.add_argument("--focus", action="store_true")
    parser.add_argument("--resume-screen", action="store_true")
    args = parser.parse_args()
    geo = native_geometry()
    data = json.loads((HERE/"native.json").read_text())
    with np.load(HERE/"native_inputs.npz",allow_pickle=False) as ys:
        if args.replay:
            replay(geo,data,ys,focus=args.focus)
        else:
            screen(geo,data,ys,resume=args.resume_screen)


if __name__ == "__main__":
    ctx.prec, ctx.threads = 128,1
    with threadpool_limits(1):
        main()
