"""Longer paired mechanism continuation after the prespecified 500-step probe.

All three positive/negative cases continue, retaining the first study outputs.
This starts each branch afresh from the saved AP state, not a rerun of AP.
An exposed exploratory comparison; no GT participates in this program.
"""
from normal_memory_probe import *

def regularizer(m):
    a=m.reshape(24,24)
    d=np.r_[np.diff(a,axis=0).ravel(),np.diff(a,axis=1).ravel()]/19.5
    return 1e-4*np.mean(d*d)+1e-4*np.mean(np.sqrt(d*d+.01**2)-.01)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--case',required=True);ap.add_argument('--steps',type=int,default=3000);ap.add_argument('--integer-only',action='store_true');args=ap.parse_args()
    out=ROOT/'work/normal_memory'/args.case;old=json.loads((out/'result.json').read_text());ss=dict(np.load(out/'states.npz'))
    pub=dict(np.load(BASE/'inputs'/(args.case+'.npz')));G,E,S,Y=pub['G'],pub['E0'],pub['S'],pub['Y']
    P=fc.DenseFieldCurrentProjectors(G,E,S,Y,old['delta'],scale=6.,lower=.02,upper=20.,grouped_ports=True)
    cx,cq,dx,dq,m=[ss[k]for k in ['cX','cQ','dX','dQ','m']];crit=old['first_tcrit'];first=old['first_cell']
    rows=[];arrays={};receipts={};normY2=float(np.vdot(Y,Y).real)
    def score(mm):
        f=lu_factor(np.eye(len(G))-mm[:,None]*G,check_finite=False)
        jj=lu_solve(f,mm[:,None]*E,check_finite=False);r=S@jj-Y
        D=.5*float(np.vdot(r,r).real)/normY2;R=regularizer(mm)
        return dict(Phi=float(D+R),D=D,R=float(R),residual_over_delta=float(np.linalg.norm(r)/P.delta))
    init=score(pub['m_init']);arrays['initial_material']=pub['m_init'];arrays['AP_material']=m
    arms=[('integer_event',max(1.,float(np.floor(crit)+1)))] if args.integer_only else [('ordinary',1.),('event',max(1.,1.001*crit))]
    for arm,t0 in arms:
        cpu=time.process_time();wall=time.perf_counter();counts=dict(P.costs)
        zx,zq=cx+t0*dx,cq+t0*dq;best=dict(init,iteration=-1);bestm=pub['m_init'].copy()
        first_box=None;pred_cell_box=None;first50=[];nLU=0
        for k in range(args.steps):
            px,pq,pm=P.material(zx,zq)
            # Large slope transitions are distinct from small smooth motion in
            # interior material cells caused by a nonzero AP pair defect.
            boxjump=(m<=.020000001)&(pm>=19.999999999)
            if first_box is None and boxjump.any():first_box=k
            if pred_cell_box is None and boxjump[first]:pred_cell_box=k
            ax,aq=P.affine_ball(2*px-zx,2*pq-zq)
            if k<50:first50.append(dict(iteration=k,predicted_cell_m=float(pm[first]),full_box_jumps=int(boxjump.sum()),max_material_change=float(np.max(abs(pm-m))),normal_ray_error=pnorm(zx-(cx+(t0+k)*dx),zq-(cq+(t0+k)*dq))))
            if k%50==0 or k==args.steps-1:
                sc=score(pm);nLU+=1
                if sc['Phi']<best['Phi']:best=dict(sc,iteration=k);bestm=pm.copy()
                rows.append(dict(arm=arm,iteration=k,cpu_seconds=time.process_time()-cpu,shadow_gap=pnorm(ax-px,aq-pq),**sc))
                arrays[f'{arm}_m_{k:05d}']=pm.copy()
            zx,zq=zx+ax-px,zq+aq-pq
        arrays[arm+'_best_m']=bestm;arrays[arm+'_last_m']=pm
        receipts[arm]=dict(best=best,first_box_jump=first_box,predicted_cell_box_jump=pred_cell_box,first50=first50,
            cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
            additional_projection_cost={k:P.costs[k]-counts[k]for k in counts},physical_LU=nLU,physical_RHS=nLU*E.shape[1])
    prefix='integer' if args.integer_only else 'followup'
    np.savez_compressed(out/(prefix+'_materials.npz'),**arrays)
    result=dict(case=args.case,steps=args.steps,threshold=crit,initial=init,arms=receipts,trajectory=rows,
        shared_AP_cpu_seconds_excluded_from_incremental_branch_times=old['cpu_seconds'],
        accounting_note='Original probe CPU includes AP3000 and both500 continuations; it must not be claimed to equal AP-only setup. Followup branch times exclude published projector setup and shared AP; exact projection/physical calls are reported.',
        selection='minimum original Phi over initial and scheduled full physical candidates; no corrector; no truth read',
        scope='exposed exploratory followup, not independent confirmation or full G19/DR benchmark')
    (out/(prefix+'_result.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(case=args.case,arms={k:{q:v[q]for q in ['best','first_box_jump','predicted_cell_box_jump','cpu_seconds']}for k,v in receipts.items()}),indent=2),flush=True)

if __name__=='__main__':main()
