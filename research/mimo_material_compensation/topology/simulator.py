"""Isolated old/new acquisition and post-design clean scoring."""
from pathlib import Path
import argparse
import json

from ports import A, PRIOR, PortEngine, WATT, dump, elapsed, load, now, project, save
import numpy as np
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
CFG=json.loads((HERE/"campaign.json").read_text())
OUT=HERE/"loop"


def simulator_material(case,engine):
    if not case.startswith("fresh_"):
        return load(PRIOR/"ports"/"private"/f"{case}.npz",("material",))["material"].ravel()
    index={"fresh_bg_a":0,"fresh_bg_b":1}[case]
    rng=np.random.default_rng((3357412997,3689502347)[index])
    x,y=engine.op["centers"].T
    background=.34+.22*np.cos(2.9*x+2.3*y)+.20*np.sin(4.2*y-1.1*x)
    # These variants include spatially extended smooth background and shape changes.
    value=background
    for (cx,cy),(a,b),contrast,turn in (
        ((-.33,.05),(.47,.39),2.90,.76),
        ((.39,.48),(.33,.25),2.42,1.54),
        ((.32,-.43),(.38,.30),2.65,-.65)):
        cx,cy=np.array([cx,cy])+rng.uniform(-.22,.22,2)
        u=(x-cx)*np.cos(turn)+(y-cy)*np.sin(turn)
        v=(y-cy)*np.cos(turn)-(x-cx)*np.sin(turn)
        footprint=((u/a)**2+(v/b)**2<1+.21*np.cos(3*np.arctan2(v/b,u/a)+rng.uniform(0,6.28)))
        value=value+contrast*footprint*(1+.13*np.cos(4*x+5*y))
    return np.clip(value,0,4.5)


def noise_draw(case,shot):
    row=[r["case"] for r in CFG["objects"]].index(case)
    rng=np.random.default_rng(np.random.SeedSequence([28462026,row,shot]))
    return rng.normal(size=(24,24))+1j*rng.normal(size=(24,24))


def noisy_block(full,tx,rx,radius,direction):
    clean=project(full,tx,rx)
    nd=direction[np.asarray(rx)[:,None],np.asarray(tx)[None,:]]
    noise=(2*radius/3)*nd/np.linalg.norm(nd)
    return clean+noise,clean,noise


def initial(case):
    cfg=next(x for x in CFG["objects"] if x["case"]==case)
    side=cfg["side"]
    path=OUT/case
    start=now()
    tx,rx=A,np.arange(8,16)
    if case.startswith("fresh_"):
        engine=PortEngine(side)
        material=simulator_material(case,engine)
        full,_,cost=engine.independent_response(material)
        # Predeclared absolute radius for new backgrounds. This does not use truth magnitude.
        radius=0.003
        y,_,_=noisy_block(full,tx,rx,radius,noise_draw(case,0))
        observation_id=f"mcomp26_{case}_f400_a08_c48_initialA"
        provenance="new independent bare VIE and loaded circuit"
    else:
        d=load(PRIOR/"ports"/"public"/"cases"/f"{case}.npz",("Y","radius","observation_id"))
        y,radius=d["Y"],float(d["radius"])
        observation_id=str(d["observation_id"].item())
        provenance="original immutable public noisy A"
        cost={}
    save(path/"initial_A.npz",Y=y,radius=radius,tx=tx,rx=rx)
    dump(path/"initial_receipt.json",dict(observation_id=observation_id,provenance=provenance,
                  calibration="exact nominal fixed instrument",generator_cost=cost,**elapsed(start)))


def acquire(case):
    path=OUT/case
    proposal_path=path/"design_frozen.json"
    proposal_bytes=proposal_path.read_bytes()
    proposal=json.loads(proposal_bytes)
    if proposal["case"]!=case or not (path/"design_candidates.npz").exists():
        raise RuntimeError("A saved case-specific design and fixed predictions are required before new acquisition")
    start=now()
    old=load(path/"initial_A.npz")
    cfg=next(x for x in CFG["objects"] if x["case"]==case)
    engine=PortEngine(cfg["side"])
    material=simulator_material(case,engine)
    full,total,cost=engine.independent_response(material)
    current=(engine.eye-total)/np.sqrt(engine.zref)
    ant=np.linalg.solve(current.T,engine.drive.T).T-engine.zref*engine.eye
    records={}
    arrays={}
    sealed_arrays={}
    for name,record in proposal["actions"].items():
        tx,rx=np.asarray(record["tx"]),np.asarray(record["rx"])
        radius=float(old["radius"])
        y,clean,noise=noisy_block(full,tx,rx,radius,noise_draw(case,1))
        j=current[:,tx]
        accepted=(np.conj(j[tx,np.arange(8)])*(ant@j)[tx,np.arange(8)]).real
        antenna=np.real(np.sum(np.conj(j)*(ant@j),axis=0))
        reflected=abs(total[tx,tx])**2
        others=np.sum(abs(total[:,tx])**2,axis=0)-reflected
        records[name]=dict(observation_id=f"mcomp26_{case}_f400_a08_c48_{name}",
            power_per_pulse_available_W=WATT,source_energy_sum_W_times_pulse_duration=1.,
            real_driven_accepted_W_sum=float(WATT*accepted.sum()),terminated_other_ports_W_sum=float(WATT*others.sum()),
            radiated_and_absorbed_W_sum=float(WATT*antenna.sum()),reflected_W_sum=float(WATT*reflected.sum()),
            radius=radius,actual_noise_over_radius=float(np.linalg.norm(noise)/radius),pulses=8,complex_reads=64,
            port_role_changes_from_old_A=int(np.count_nonzero(np.isin(np.arange(16),A)!=np.isin(np.arange(16),tx))))
        arrays[name+"_Y"]=y
        arrays[name+"_r"]=radius
        sealed_arrays[name+"_clean"]=clean
        sealed_arrays[name+"_noise"]=noise
    save(path/"new_observed.npz",**arrays)
    save(path/"sealed_scoring.npz",material=material,full_total_S=total,full_calibrated_delta=full,
         full_port_current=current,loaded_antenna=ant,**sealed_arrays)
    dump(path/"acquisition_receipt.json",dict(case=case,actions=records,independent_generator_cost=cost,
        simulation_accounting="One bare 24-RHS material solve simulates all paid interventions; pulses, reads, energy, port switches stated for each counterfactual action independently",
        frozen_design_mtime_ns=proposal_path.stat().st_mtime_ns,
        observation_mtime_ns=(path/"new_observed.npz").stat().st_mtime_ns,**elapsed(start)))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("phase",choices=("initial","acquire"))
    p.add_argument("case", choices=tuple(x["case"] for x in CFG["objects"]))
    args=p.parse_args()
    with threadpool_limits(limits=1,user_api="blas"):
        (initial if args.phase=="initial" else acquire)(args.case)
