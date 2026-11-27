"""Independent bare-wave and loaded-circuit replay of decisive packet witnesses."""
from bound import CFG, HERE, block, elapsed, load_case, now, np, up, write_json
from quick import Response
from model import power_report
import json


def indexed(stem):
    data=json.loads((HERE/f"budget_{stem}.json").read_text())
    return {(r["layout"],r["contract"],r["kind"],r["t"],r["beta"]):r for r in data["records"]}


def focused_best():
    records=indexed("focus_joint32")
    choices=[]
    for key,row in records.items():
        lay,contract,kind,t,beta=key
        if lay!="A" or row["status"]!="hideable" or "x_key" not in row:
            continue
        level=up(row["exterior_norm"]["max"])
        for k,b in records.items():
            if k[:4] == ("B",contract,kind,t) and k[4]>=level and b["status"]=="unhideable":
                choices.append((k[4]-level,level,beta,t,row,b))
    if not choices:
        return None
    return max(choices,key=lambda q:q[0])


def independently_loaded(base,side,x,br):
    j=br.full(x)["j"]
    baredelta=base["S_bare"]@j
    loaded=base["M_nominal"]-1j*baredelta
    antenna=base["Z_AA"]-1j*baredelta
    drive=2*np.sqrt(float(base["reference_ohms"])/8)*np.eye(24,16)
    empty=np.linalg.solve(base["M_nominal"],drive)
    current=np.linalg.solve(loaded,drive)
    p=np.sqrt(float(base["reference_ohms"]))*(empty-current)
    receipt=power_report(current,antenna,drive)
    result={}
    for tag,spec in CFG["layout"].items():
        ids=spec["tx"]
        result[tag]=dict(
            available_W=float(sum(receipt["available_W"][i] for i in ids)),
            driven_accepted_W=float(sum(receipt["driven_accepted_W"][i] for i in ids)),
            reflected_W=float(sum(receipt["reflected_W"][i] for i in ids)),
            dissipated_or_radiated_W=float(sum(receipt["antenna_W"][i] for i in ids)),
            fixed_terminal_dissipation_W=float(sum(receipt["other_termination_W"][i] for i in ids)),
            sources=len(ids),reads=len(ids)*len(spec["rx"]))
    green=br.model.green.apply(j)
    rhs=x[:,None]*base["E0_bare"]
    residual=np.linalg.norm(j-x[:,None]*green-rhs,axis=0)
    relative=0. if not np.any(rhs) else float(np.max(residual/np.maximum(np.linalg.norm(rhs,axis=0),1e-300)))
    return p,j,current,dict(layouts=result,relative_current_residual=relative,
            network_power_balance_max_W=receipt["energy_balance_max_W"],
            source_power_balance_max_W=receipt["incident_balance_max_W"],
            reciprocity_relative=receipt["reciprocity_relative"])


def main():
    total=now()
    requests=[
        ("joint32","A","prospective","max",0.25,0.001),
        ("joint32","B","actual","max",0.045,0.0001),
        ("air32","B","prospective","max",0.56,0.003),
        ("air96","B","prospective","max",0.58,0.003),
        ("native96","B","prospective","max",0.05,0.00001)
    ]
    focus_receipt=None
    if (HERE/"budget_focus_joint32.json").is_file():
        best=focused_best()
        if best:
            _,bmin,point,t,row,b=best
            requests.append(("focus_joint32","A","prospective","max",t,point))
            focus_receipt=dict(t=t,minimum_budget_of_saved_A_witness_upper=bmin,
                               audited_A_witness_grid_budget=point,
                               B_lower_certified_until=b["beta"],
                               A_point_upper_over_two_r=row["witness_upper_over_threshold"],
                               A_no_compensation_lower_over_two_r=indexed("focus_joint32")[("A","prospective","max",t,0.)]["lower_over_threshold"],
                               B_all_edits_lower_over_two_r=b["lower_over_threshold"],
                               B_norm_control_lower_over_two_r=b["norm_control_lower_over_threshold"],
                               B_empirical_direction_lower_over_two_r=b["empirical_direction_lower_over_threshold"],
                               exteriormaterial_norms=row["exterior_norm"])
    rows=[]
    arrays={}
    for count,(stem,lay,contract,kind,t,beta) in enumerate(requests):
        start=now()
        key=stem.removeprefix("focus_")
        source=indexed(stem)[(lay,contract,kind,t,beta)]
        assert source["status"]=="hideable"
        side,base,exc,rec,roi,c,r,obs_id,observed=load_case(key)
        with np.load(HERE/f"raw_{stem}.npz",allow_pickle=False) as raw:
            prefix=f"t{source['t_index']:02d}_"
            x_name=source.get("x_key",prefix+"anchor")
            x=raw[x_name]
            predicted=raw[x_name[:-2]+"_p"] if "x_key" in source else raw[prefix+"p_anchor"]
            reference=raw["p_reference"]
        assert np.all(x[roi]==t) and np.min(x)>=0 and np.max(x)<=4.5
        maxdiff=up(np.max(np.abs((x-c)[~roi]))) if np.any((x-c)[~roi]) else 0.
        assert maxdiff<=beta
        # Remove all circuit loading from volume Green and solve for bare port
        # columns. The object-dependent 24-port loading is then applied explicitly.
        bare=dict(base,low_rank_left=np.zeros((len(c),1),complex),
                  low_rank_right=np.zeros((1,len(c)),complex),
                  diagonal_dressed=base["diagonal"])
        bare.pop("G",None)
        bare_response=Response(side,bare,base["E0_bare"],base["S_bare"])
        p0,j0,i0,i0log=independently_loaded(base,side,c,bare_response)
        p1,j1,i1,i1log=independently_loaded(base,side,x,bare_response)
        a,b=block(p0,lay),block(p1,lay)
        allocated=2*r if contract=="prospective" else r
        residual=float(np.linalg.norm((b-a) if contract=="prospective" else (b-observed[lay]))/allocated)
        direct_prior=np.linalg.norm(block(p0,lay)-block(reference,lay))/r
        direct_new=np.linalg.norm(block(p1,lay)-block(predicted,lay))/r
        if max(direct_prior,direct_new)>3e-6 or residual>=1:
            raise RuntimeError(f"Independent replay disagreed with stored witness: {stem} {lay}")
        key_id=f"w{count:02d}_"
        arrays.update({key_id+"x_reference":c,key_id+"x_witness":x,key_id+"roi":roi,
                       key_id+"p0":p0,key_id+"p1":p1,key_id+"currents_bare0":j0,key_id+"currents_bare1":j1,
                       key_id+"current_ports0":i0,key_id+"current_ports1":i1,
                       key_id+"midpoint_if_prospective":(a+b)*.5 if contract=="prospective" else observed[lay]})
        rows.append(dict(case=stem,layout=lay,contract=contract,t=t,beta=beta,raw_prefix=key_id,
                        observation_id=obs_id,r=r,strict_witness_upper=source["witness_upper_over_threshold"],
                        independently_replayed_over_threshold=residual,
                        independent_to_stored_reference_difference_over_r=float(direct_prior),
                        independent_to_stored_witness_difference_over_r=float(direct_new),
                        exterior_norms=source["exterior_norm"],verified_max_budget=float(maxdiff),
                        reference=i0log,witness=i1log,independent_cost=dict(bare_response.model.cost),**elapsed(start)))
    shape=json.loads((HERE/"shape_witness.json").read_text())["records"][-1]
    side,base,exc,rec,roi,c,r,obs_id,observed=load_case("joint32")
    with np.load(HERE/"shape_witness_raw.npz",allow_pickle=False) as raw:
        x=raw[shape["raw_prefix"]+"_x"]
        independent_shape_nominal=raw[shape["raw_prefix"]+"_p"]
    assert np.array_equal(x[~roi],c[~roi]) and x[roi].mean()>=0.27
    bare=dict(base,low_rank_left=np.zeros((len(c),1),complex),
              low_rank_right=np.zeros((1,len(c)),complex),diagonal_dressed=base["diagonal"])
    bare.pop("G",None)
    br=Response(side,bare,base["E0_bare"],base["S_bare"])
    p0,j0,i0,i0log=independently_loaded(base,side,c,br)
    p1,j1,i1,i1log=independently_loaded(base,side,x,br)
    distance=float(np.linalg.norm(block(p1-p0,"B"))/(2*r))
    difference=float(np.linalg.norm(block(p1,"B")-block(independent_shape_nominal,"B"))/r)
    if distance>=1 or difference>3e-6:
        raise RuntimeError("Independent shape-witness replay failed")
    arrays.update(shape_x_reference=c,shape_x_witness=x,shape_roi=roi,shape_p0=p0,shape_p1=p1,
                  shape_currents_bare0=j0,shape_currents_bare1=j1,
                  shape_current_ports0=i0,shape_current_ports1=i1)
    shape_replay=dict(t=shape["t"],observation_id=obs_id,r=r,raw_prefix="shape_",
                      independent_norm_over_two_r=distance,strict_upper_over_two_r=shape["outward_witness_upper_over_two_r"],
                      independent_to_stored_witness_difference_over_r=difference,
                      reference=i0log,witness=i1log,cost=dict(br.model.cost))
    write_json(HERE/"independent_replay.json",dict(rows=rows,focus_window=focus_receipt,
               shape_replay=shape_replay,**elapsed(total)))
    with (HERE/"independent_replay.npz").open("wb") as f:
        np.savez_compressed(f,**arrays)
    print(focus_receipt,flush=True)
    print([(z["case"],z["independently_replayed_over_threshold"],z["independent_to_stored_witness_difference_over_r"]) for z in rows],flush=True)


if __name__=="__main__":
    main()
