"""Development replay to measure availability of the next routed proposal and audit."""
import json
import time
from check_interventions import CASE_NAMES, HISTORICAL, HERE,engine,model_for,np,read_public
import lazy_srap


def run():
    records=[]
    for name in CASE_NAMES:
        pub,c=read_public(name,HISTORICAL)
        methods=('srap','graph','srap_lazy') if name==CASE_NAMES[-1] else ('srap','graph')
        for method in methods:
            model=model_for(pub,c)
            zero=np.zeros(model.n)
            model.objective(zero,pub['Y'])
            next_audit=model.clock.wall()
            events=[]
            run_arm=lazy_srap.run_cie_arm if method=='srap_lazy' else engine.run_cie_arm
            for m,step,stage,force in run_arm(model,pub['Y'],c,method,zero,pub['delta']):
                if force or model.clock.wall()>=next_audit:
                    ready=model.clock.wall()
                    phi,res,_=model.objective(m,pub['Y'])
                    available=model.clock.wall()
                    next_audit=available+1.
                    if force:
                        record=dict(step=step,stage=stage,proposal_wall=ready,audited_wall=available,
                                    phi=phi,relative_physical_residual=res,
                                    affine=next(v for v in reversed(model.trace_route) if v['event']=='AP'))
                        events.append(record)
                        if len(events)>=2:
                            break
            records.append(dict(case=name,method=method,events=events,counts=dict(model.clock.counts)))
            print(name,method,[(r['step'],round(r['audited_wall'],3)) for r in events],flush=True)
            (HERE/'next_checkpoint.json').write_text(json.dumps(records,indent=2,default=float))
    return records


if __name__=='__main__':
    run()
