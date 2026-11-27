"""Independent one-scalar simulations and later blind decisions."""
from fractions import Fraction
import argparse
import json
import secrets

from coherent import (HERE, ENV, BETA, BAND, elapsed, fro, lower_ball, now, np, packet,
                      rd, stamp, up, upper, write_json)
from probe import Circuit, at_load
from passive_replay import VerifiedWave, load_cover, stored_response


TAGS=("capture_37","capture_82")
NOISE_FRACTION=.85
ACTIONS=("four_by_four","single_pair")


def receive_u(a):
    return np.array([complex(*x) for x in a["u"]])


def incident_v(a):
    return np.array([complex(*x) for x in a["v"]])


def entries(side):
    p=json.loads((HERE/f"proposal_{side}.json").read_text())
    return p,{k:a for k,a in p["best"].items() if a["passes_fixed_unit_energy"]}


def synthesize(base,ix,side,seed):
    rng=np.random.default_rng(seed)
    x=base.copy()
    row,col=np.divmod(ix,side)
    modes=rng.uniform(.11,.43,size=(2,2))
    phases=rng.uniform(0.,2*np.pi,size=2)
    wave=np.cos(modes[:,0,None]*row+modes[:,1,None]*col+phases[:,None]).mean(axis=0)
    x[ix]=np.clip(x[ix]+.006*wave+.002*rng.uniform(-1.,1.,size=len(ix)),0.,4.5)
    left=np.setdiff1d(np.arange(len(base)),ix)
    n=max(16,side//2)
    j=left[np.linspace(12,len(left)-13,n,dtype=int)]
    x[j]=rng.uniform(.001,.0075,size=len(j))
    return x,np.unique(np.r_[ix,j])


def physical_scalar(s,u,v):
    transmit_wave=np.zeros(16,complex)
    transmit_wave[:8]=v
    ports=at_load(s,1.)@transmit_wave
    return np.vdot(u,ports[8:16]),ports


def absolute_noise(r,seed):
    rng=np.random.default_rng(seed)
    q=complex(*rng.normal(size=2))
    return (NOISE_FRACTION*r/abs(q))*q


def save(dest,**data):
    with dest.open("xb") as f:
        np.savez_compressed(f,**data)


def acquire(side):
    t=now()
    p,possible=entries(side)
    if not possible:
        print(side,"no certified unit-energy action",flush=True)
        return
    dst=HERE/f"readings_{side}.npz"
    if dst.exists():
        raise FileExistsError("A frozen acquisition exists")
    xs,ix,roi,r,obs,source,dual,receipt=packet(side)
    rng=secrets.SystemRandom()
    observed={}
    truth={}
    records={}
    raw={}
    for name,action in possible.items():
        order=[0,1]
        rng.shuffle(order)
        u,v=receive_u(action),incident_v(action)
        counts={}
        for tag,k in zip(TAGS,order):
            h=f"{name}_{tag}"
            seed=rng.getrandbits(128)
            noise_seed=rng.getrandbits(128)
            x,ids=synthesize(xs[k],ix,side,seed)
            w=Circuit(side,BAND,ids)
            s=w.solve(x[ids])
            val,port=physical_scalar(s,u,v)
            y=val+absolute_noise(r,noise_seed)
            observed[h]=y
            raw[h+"_material"]=x
            raw[h+"_full_forward_scattering_only_for_simulation"]=s
            raw[h+"_simulated_outgoing_16port_wave"]=port
            raw[h+"_numerical_clean_scalar"]=val
            truth[h]=dict(reference_index=k,label="high" if k==0 else "low",
                          sample_seed_hex=f"{seed:032x}",noise_seed_hex=f"{noise_seed:032x}")
            counts[h]=dict(w.count,complex_channels_actually_read=1,simultaneous_shots=1,
                           coherent_driven_ports=int(np.count_nonzero(v)), receive_feeds=int(np.count_nonzero(u)),
                           total_incident_wave_norm_squared=float(np.sum(abs(v)**2)),
                           numerical_receive_output_power=float(abs(val)**2))
        records[name]=counts
    save(HERE/f"offline_material_{side}.npz",roi=roi,**raw)
    write_json(HERE/f"offline_labels_{side}.json",dict(generated_utc=stamp(),samples=truth))
    save(dst,**observed)
    write_json(HERE/f"acquisition_{side}.json",dict(sealed_utc=stamp(),side=side,frequency_MHz=BAND,
                noise_fraction=NOISE_FRACTION,physical_scalar="u^* S_open v",action_counts=records,
                simulation="25-port bare VIE and circuit assembly offline; one composite transmit column, one harvested complex scalar per trial.",
                **elapsed(t)))
    print(side,"acquired",len(observed),"scalar samples",flush=True)


def decide(side):
    t=now()
    p,possible=entries(side)
    if not possible:
        return
    name=HERE/f"decisions_{side}.json"
    if name.exists():
        raise FileExistsError("A frozen prediction exists")
    ans={}
    with np.load(HERE/f"readings_{side}.npz",allow_pickle=False) as f:
        for key in f.files:
            mode=next(v for v in possible if key.startswith(v+"_"))
            chosen=possible[mode]
            h=chosen["frozen_threshold"]
            # The selected functional is Re(z), an exact stored binary64
            # coordinate of a single complex readout.
            scalar=complex(f[key])
            v=scalar.real
            label="high" if v>h else "low" if v<h else "unresolved"
            ans[key]=dict(action=mode,label=label,value_real=v,value_imag=scalar.imag,
                          threshold=h,margin_abs=abs(v-h),
                          task_interval=([float(np.nextafter(1.99,-np.inf)),float(np.nextafter(2.01,np.inf))]
                                         if label=="high"
                                         else [0.,float(np.nextafter(.01,np.inf))] if label=="low" else [0.,4.5]))
    write_json(name,dict(sealed_utc=stamp(),inputs=[f"proposal_{side}.json",f"readings_{side}.npz"],results=ans,
                         scope="Fixed scalar comparison only. The sampler labels and materials remain unread.",**elapsed(t)))
    print(side,"blind predictions sealed",flush=True)


def score(side):
    t=now()
    p,possible=entries(side)
    if not possible:
        return
    if (HERE/f"scores_{side}.json").exists():
        raise FileExistsError("Sealed scores exist")
    decision=json.loads((HERE/f"decisions_{side}.json").read_text())
    acquired=json.loads((HERE/f"acquisition_{side}.json").read_text())
    if not p["freeze_utc"]<acquired["sealed_utc"]<decision["sealed_utc"]:
        raise ArithmeticError("Freeze/acquire/decide rule broken")
    offline=json.loads((HERE/f"offline_labels_{side}.json").read_text())["samples"]
    xs,ix,roi,r,obs,source,dual,previous=packet(side)
    results={}
    proof_arrays={}
    with np.load(HERE/f"offline_material_{side}.npz",allow_pickle=False) as f, \
         np.load(HERE/f"readings_{side}.npz",allow_pickle=False) as y:
        for tag,ans in decision["results"].items():
            actor=possible[ans["action"]]
            u,v=receive_u(actor),incident_v(actor)
            key=offline[tag]
            x=f[tag+"_material"]
            base=xs[key["reference_index"]]
            seed=int(key["sample_seed_hex"],16)
            z,ids=synthesize(base,ix,side,seed)
            if not np.array_equal(z,x):
                raise ArithmeticError("Independent material generator is not reproducible")
            recorded=complex(y[tag])
            approximate=f[tag+"_numerical_clean_scalar"]
            if approximate+absolute_noise(r,int(key["noise_seed_hex"],16))!=recorded:
                raise ArithmeticError("Scalar readout generator is not reproducible")
            max_change=max(abs(Fraction(a)-Fraction(b)) for a,b in zip(x,base))
            if max_change>Fraction(1,100):
                raise ArithmeticError("Specimen outside its own class")
            w=Circuit(side,BAND,ids)
            n=side**2
            basis=dict(spectrum=w.symbol,grid_shape=(side,side),S=w.c,E0=w.b,
                       low_rank_left=np.zeros((n,1)),low_rank_right=np.zeros((1,n)))
            verifier=VerifiedWave(basis)
            s,err,slog=stored_response(w,verifier,x[ids])
            op,operr,_=load_cover(s,err,1.,0.,0.)
            exact_proposal=np.vdot(u,op[8:16,:8]@v)
            scalar_error=up(operr+rd(20)*fro(op))
            strict=up((abs(exact_proposal-recorded)+scalar_error+rd(5)*(abs(exact_proposal)+abs(recorded)))/r)
            if strict>=1:
                raise ArithmeticError("Actual measurement outside the declared scalar ball")
            mean=sum(Fraction(a) for a in x[roi])/int(np.sum(roi))
            valid=(Fraction(ans["task_interval"][0])<=mean<=Fraction(ans["task_interval"][1]))
            correct=key["label"]==ans["label"]
            if not valid or not correct:
                raise ArithmeticError("Declared class threshold did not meet the scored material")
            results[tag]=dict(action=ans["action"],sealed_prediction=ans["label"],
                              revealed_class=key["label"],correct=correct,true_region_mean=float(mean),
                              mean_in_interval=valid,distance_infinity_to_own_spec_upper=up(float(max_change)),
                              actual_noise_ratio_upper=strict,independent_scalar_error_upper=scalar_error,
                              verified_inverse_rhs=len(ids),verified_physical_rhs=25,verified_material_lu=1,
                              independent_wave_certificate=slog)
            proof_arrays[tag+"_independent_forward_scalar"]=exact_proposal
            proof_arrays[tag+"_observed_scalar"]=recorded
    save(HERE/f"verified_scalar_{side}.npz",**proof_arrays)
    write_json(HERE/f"scores_{side}.json",dict(scored_utc=stamp(),side=side,results=results,**elapsed(t)))
    print(side,[(a,b["sealed_prediction"],b["true_region_mean"],b["actual_noise_ratio_upper"])
                   for a,b in results.items()],flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage",choices=("acquire","decide","score"))
    p.add_argument("side",type=int,choices=(32,96),nargs="+")
    a=p.parse_args()
    for side in a.side:
        globals()[a.stage](side)
