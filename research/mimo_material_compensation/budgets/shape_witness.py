"""Test whether reshaping only the six ROI pixels already defeats the flat-task result."""
from bound import CFG, HERE, Response, audited_state, block, down, elapsed, fro, load_case, now, np, rd, up, write_json
from budget_proof import VerifiedWave
from scipy.optimize import minimize


def main():
    start=now()
    side,base,exc,rec,roi,c,r,obs_id,obs=load_case("joint32")
    rec=rec[:16]
    physical=Response(side,base,exc,rec)
    wave=VerifiedWave(dict(base,S=rec,E0=exc))
    reference=audited_state(wave,physical,rec,exc,c)
    target=block(reference.p,"B")
    bx=CFG["layout"]["B"]
    rows=[]
    arrays=dict(reference=c,roi=roi,p_reference=reference.p)
    n=int(roi.sum())
    factor=2**20
    for i,t in enumerate([0.21875,0.25,np.ceil(0.27*factor)/factor]):
        def objective(p):
            x=c.copy()
            x[roi]=p
            state=physical.full(x,derivatives=True)
            diff=(block(state["p"],"B")-target)/(2*r)
            u=state["u"][bx["rx"]][:,roi]
            e=state["e"][roi][:,bx["tx"]]
            grad=np.einsum("ij,ik,kj->k",diff.conj(),u,e,optimize=True).real/(2*r)
            return .5*np.linalg.norm(diff)**2, grad
        solutions=[]
        for init in [np.full(n,t),*list(n*t*np.eye(n))]:
            ans=minimize(objective,init,jac=True,method="SLSQP",bounds=[(0.,4.5)]*n,
                         constraints=[dict(type="eq",fun=lambda x:x.sum()-n*t,jac=lambda x:np.ones(n))],
                         options=dict(maxiter=100,ftol=1e-12))
            solutions.append(ans)
        best=min(solutions,key=lambda z:z.fun)
        ticks=np.clip(np.rint(best.x*factor).astype(int),0,int(4.5*factor))
        missing=int(n*t*factor)-int(ticks.sum())
        if missing>0:
            ticks[int(np.argmin(ticks))]+=missing
        if missing<0:
            ticks[int(np.argmax(ticks))]+=missing
        assert ticks.sum()==int(n*t*factor) and ticks.min()>=0 and ticks.max()<=int(4.5*factor)
        x=c.copy()
        x[roi]=ticks/factor
        audit=audited_state(wave,physical,rec,exc,x)
        d=block(audit.p,"B")-target
        err=up(audit.response_error+reference.response_error+rd(3)*(fro(target)+fro(block(audit.p,"B"))))
        upper=up((fro(d)+err)/(2*r))
        arrays[f"t{i:02d}_x"]=x
        arrays[f"t{i:02d}_p"]=audit.p
        arrays[f"t{i:02d}_j"]=audit.j
        rows.append(dict(t=t,raw_prefix=f"t{i:02d}",roi_values=(ticks/factor).tolist(),
                         exterior_max=0.,unconstrained_proposals=7,optimizer_results=[
                         dict(iterations=int(z.nit),method=z.message,normalized_norm=float(np.sqrt(2*z.fun))) for z in solutions],
                         outward_witness_upper_over_two_r=upper,
                         hideable=upper<1.,**elapsed(start)))
    write_json(HERE/"shape_witness.json",dict(obs_id=obs_id,r=r,records=rows,cost=dict(physical.model.cost),**elapsed(start)))
    with (HERE/"shape_witness_raw.npz").open("wb") as f:
        np.savez_compressed(f,**arrays)
    print([(q["t"],q["outward_witness_upper_over_two_r"],q["roi_values"]) for q in rows],flush=True)


if __name__=="__main__":
    main()
