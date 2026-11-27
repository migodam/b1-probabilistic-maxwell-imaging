"""Recompute each bounded nonuniform witness via the bare wave and loaded circuit."""
from pathlib import Path
import sys
import json

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from bound import Response, block, elapsed, load_case, now, np, write_json
from replay import independently_loaded


def main():
    start=now()
    receipt=json.loads((HERE/"receipt.json").read_text())
    side,base,_,_,roi,c,r,obs,_=load_case("joint32")
    bare=dict(base,low_rank_left=np.zeros((len(c),1),complex),
              low_rank_right=np.zeros((1,len(c)),complex),diagonal_dressed=base["diagonal"])
    bare.pop("G",None)
    physics=Response(side,bare,base["E0_bare"],base["S_bare"])
    p0,j0,i0,ref=independently_loaded(base,side,c,physics)
    arrays=dict(x_reference=c,p_reference=p0,currents_reference_bare=j0,currents_reference_ports=i0)
    rows=[]
    with np.load(HERE/"raw.npz",allow_pickle=False) as archive:
        for q in receipt["records"]:
            witness=q["witness"]
            if not witness:
                continue
            key=witness["raw_prefix"]
            x=archive[key+"_x"]
            p,j,i,audit=independently_loaded(base,side,x,physics)
            d=np.linalg.norm(block(p-p0,"B"))/(2*r)
            agreement=np.linalg.norm((block(p,"B")-block(archive[key+"_p"],"B")))/r
            if d>=1 or agreement>=3e-6:
                raise ArithmeticError("Bridge witness disagrees with independent physical replay")
            arrays.update({key+"_bare_currents":j,key+"_port_currents":i,key+"_response":p})
            rows.append(dict(beta=q["beta"],alpha=witness["alpha_upper"], B_norm_over_two_r=d,
                            independent_vs_loaded_difference_over_r=agreement, audit=audit))
    with (HERE/"replay_raw.npz").open("wb") as f:
        np.savez_compressed(f,**arrays)
    write_json(HERE/"replay.json",dict(obs_id=obs, reference=ref,rows=rows,
               cost=dict(physics.model.cost),**elapsed(start)))
    print([(v["beta"],v["alpha"],v["B_norm_over_two_r"]) for v in rows],flush=True)


if __name__=="__main__":
    main()
