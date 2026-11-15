"""Post-seal public-data stationarity probes, never an online candidate search.

Fixed exposed cases and steps; no truth, optimizer, candidate replacement or
stopping-threshold calibration. Every trial is validation work and is charged.
"""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
import argparse,csv,json,sys,time
from pathlib import Path


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    root=a.root.resolve();out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    protocol=dict(case_indices=[7,16,29,37,44,3,17],pairs=['00','01','10','11'],
        material_step_Linf=[.1,.01,.001],selection='saved primary, no replacement',
        truth_access='NONE',scope='Post-seal exposed-case public-physics diagnostic; not solver continuation')
    (out/'PROTOCOL_BEFORE_PROBES.json').write_text(json.dumps(protocol,indent=2)+'\n')
    sys.path[:0]=[str(root/'src'),str(root.parents[1])]
    import score_dualmetric as scoring
    seal,tasks,_,account=scoring.authorize(root,root/'audit/PREDICTION_SEAL_AUTHORIZED.json')
    import numpy as np
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    cpu=time.process_time();wall=time.perf_counter();rows=[];trials=[];costs=[]
    wanted={f'r13_fresh_{i:04d}' for i in protocol['case_indices']}
    for task in tasks:
        if task['case_id'] not in wanted:continue
        ident=task['task_id'];folder=root/'RESULTS'/task['campaign']/ident
        with np.load(root/task['input_file'],allow_pickle=False) as z:public={k:z[k].copy() for k in z.files}
        with np.load(folder/'primary.npz',allow_pickle=False) as z:m=z['m'].copy()
        paid=json.loads((folder/'primary.json').read_text());model=MaxwellModel(public)
        Y=public['Y'];lo=float(public['lower']);hi=float(public['upper']);local_cpu=time.process_time()
        value,g=model.data_value_gradient(m,Y);state=model.state(m)
        pg=m-np.clip(m-g,lo,hi);scale=float(np.max(np.abs(pg)))
        direction=-pg/scale if scale else np.zeros_like(pg)
        dJ=model.tangent(m,direction)[0];dY=model.S@dJ
        gauss_newton_curvature=float(np.vdot(dY,dY).real/np.vdot(Y,Y).real)
        _,_,weights=model.reg.evaluate(m);dedge=model.reg.D@direction
        prior_curvature=float(np.dot(weights,dedge*dedge))
        row=dict(task_id=ident,case_id=task['case_id'],pair=task['pair'],
            saved_Phi=paid['Phi'],recomputed_Phi=value,Phi_recompute_difference=value-paid['Phi'],
            PG_Linf=scale,PG_L2=float(np.linalg.norm(pg)),gradient_L2=float(np.linalg.norm(g)),
            lower_active=int(np.count_nonzero(m<=lo+1e-8)),upper_active=int(np.count_nonzero(m>=hi-1e-8)),
            directional_derivative=float(g@direction),GN_data_curvature=gauss_newton_curvature,
            prior_curvature=prior_curvature,termination=json.loads((folder/'RECEIPT.json').read_text())['termination_reason'])
        rows.append(row)
        np.savez_compressed(out/(ident+'.npz'),m=m,gradient=g,projected_gradient=pg,
            direction=direction,current_tangent=dJ,port_tangent=dY)
        for step in protocol['material_step_Linf']:
            proposal=np.clip(m+step*direction,lo,hi);actual_step=proposal-m
            phi=model.physical_score(proposal,Y)
            trials.append(dict(task_id=ident,case_id=task['case_id'],pair=task['pair'],
                requested_Linf_step=step,actual_Linf_step=float(np.max(np.abs(actual_step))),
                Phi_before=value,Phi_after=phi,actual_Phi_change=phi-value,
                first_order_prediction=float(g@actual_step),
                second_order_GN_prediction=float(g@actual_step)+.5*step*step*(gauss_newton_curvature+prior_curvature),
                unprojected_ray_model=bool(np.array_equal(proposal,m+step*direction)),
                trial_source='OFFLINE_DIAGNOSTIC_ONLY; never inserted into paid online candidate pool'))
        costs.append(dict(task_id=ident,CPU_seconds=time.process_time()-local_cpu,
            physical_cost=dict(model.cost),source_RHS_dimension=model.ns,matrix_dimension=model.n))
    for name,data in [('ENDPOINTS',rows),('FINITE_DESCENT_PROBES',trials),('COSTS',costs)]:
        (out/(name+'.json')).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
        keys=list(dict.fromkeys(k for x in data for k in x))
        with (out/(name+'.csv')).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
            for x in data:w.writerow({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in x.items()})
    receipt=dict(status='PUBLIC_ENDPOINT_DIAGNOSTIC_COMPLETE',expected_endpoints=28,
        actual_endpoints=len(rows),trials=len(trials),truth_access='NONE',optimizer_calls=0,
        online_selections_changed=False,CPU_seconds=time.process_time()-cpu,
        wall_seconds=time.perf_counter()-wall,counts=costs,
        interpretation='A finite descent probe is not a deployment result. Nonzero PG or descent disproves stationarity at the measured scale; no successful probe does not prove a strict minimum or identifiability.')
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='counts'}))


if __name__=='__main__':main()
