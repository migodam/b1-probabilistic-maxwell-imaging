"""Offline one-sided replay of the simulator material after design and refit are sealed."""
from pathlib import Path
import argparse
import json

from ports import A,PortEngine,dump,elapsed,load,now
from verify import VerifiedWave,audit_response,audit_published_a,audit_projection
import numpy as np
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent


def check(case):
    started=now()
    p=HERE/"loop"/case
    attack=json.loads((p/"attack_results.json").read_text())
    assert "stop" in attack
    spec=json.loads((p/"design_frozen.json").read_text())
    old=load(p/"initial_A.npz")
    new=load(p/"new_observed.npz")
    x=load(p/"sealed_scoring.npz",("material",))["material"]
    engine=PortEngine()
    v=VerifiedWave(dict(engine.op,E0=engine.b,S=engine.c))
    st=engine.state(x)
    q,e,receipt=audit_response(v,x,st["j"])
    if case.startswith("n32_"):
        a=audit_published_a(v,engine,x,st["j"],receipt["inverse"]["inverse_norm_upper"],old["Y"],float(old["radius"]))
    else:
        a=audit_projection(q,e,A,np.arange(8,16),old["Y"],float(old["radius"]))
    actual={}
    for name,action in spec["actions"].items():
        actual[name]=audit_projection(q,e,action["tx"],action["rx"],new[name+"_Y"],float(new[name+"_r"]))
    assert a["strict_feasible"] and all(r["strict_feasible"] for r in actual.values())
    dump(p/"postseal_truth_verified.json",dict(case=case,scope="Offline scoring only. This material and score were never read by topology selection or the joint numerical fit.",old=a,new=actual,
          strict_original_operator_receipt=receipt,**elapsed(started)))
    print(case,round(a["upper_over_r"],6),round(max(q["upper_over_r"] for q in actual.values()),6),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("cases",nargs="+")
    with threadpool_limits(limits=1,user_api="blas"):
        for c in p.parse_args().cases:
            check(c)
