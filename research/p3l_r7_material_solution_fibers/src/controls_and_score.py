"""Separate post-seal diagnostics; never change the eight frozen algorithms."""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from physics import OUT,Book,Target,write,save,totals
from material_jets import geometry,match_distance,real_material,tasks
from run_e2 import load_parent,correct,metrics,counters


def require_seal():
    seal=json.loads((OUT/'E2/PREDICTIONS_SEALED.json').read_text())
    for row in seal['files']:
        path=Path(row['path'])
        if hashlib.blake2b(path.read_bytes()).hexdigest()!=row['blake2b']:
            raise RuntimeError('sealed prediction altered')
    return seal


def gauge_controls():
    require_seal()
    with np.load(OUT/'audit/copula_full66/CONTROLLED_DISCRETE_FIBRE_ARRAYS.npz',allow_pickle=False) as z:
        parents=z['parent_id'].copy();correct_x=z['x_correct_curve'].copy();cross_x=z['x_fixed_off_diagonal'].copy()
    planned=4*4*8
    if totals()['RHS']+planned>50000:raise RuntimeError('no gauge-control budget')
    with Book('E2','NONUNIFORM_GAUGE_COPULA_CONTROL',role='CONTROLLED_DISCRETE_FIBRE_NOT_POSTERIOR',max_wall=300) as book:
        target=Target(book)
        rows=[]
        try:
            for i,parent in enumerate(parents):
                points=[correct_x[i,0],correct_x[i,-1],cross_x[i,0],cross_x[i,1]]
                labels=['MATCHED_FIRST','MATCHED_LAST','CROSSED_FIRST_LAST','CROSSED_LAST_FIRST']
                reference=None
                for label,x in zip(labels,points):
                    p=target.evaluate(x,str(parent)+'_'+label)
                    if reference is None:reference=(p.state.chi.copy(),p.state.data.copy())
                    delta=target.white(p.state.data-reference[1])
                    row=dict(parent=str(parent),label=label,
                       material_RMS=float(np.linalg.norm(p.state.chi-reference[0])/24),
                       response_white_L2=float(np.linalg.norm(delta)),endpoint=p.key,
                       scope='CONTROLLED_DISCRETE_FIBRE; equal four-point law, not conditional prior/posterior')
                    rows.append(row)
                    save(OUT/'E2/controls/gauge'/f'{parent}_{label}.npz',x=p.x,q=p.q,
                        material=p.state.chi,prediction=p.state.data,current=p.state.current,field=p.state.field,
                        whitened_response_change=delta)
            write(OUT/'E2/controls/gauge/RESULTS.json',dict(rows=rows,cost=counters(book),
                same_material_gauge_never_counted_as_new_material=True))
        finally:target.close()


def profiles():
    require_seal()
    if totals()['RHS']+1616>50000:raise RuntimeError('no registered development-profile budget')
    rows=[]
    for parent in ['R6_CDEV_0000','R6_CDEV_0002']:
        folder=OUT/'E2/controls/profiles'/parent
        with Book('E2','ORDINARY_PROFILE20_'+parent,parent=parent,role='EXPOSED_DEVELOPMENT_PROFILE_REFERENCE',max_wall=600) as book:
            target=Target(book)
            try:
                base,y=load_parent(target,parent,folder)
                g=geometry(target.prior,base.x,base.J,2026102815)
                xp,match=match_distance(target.prior,base.x,lambda h:base.x+h*g['v'],.02)
                if xp is None:raise RuntimeError('registered development radius infeasible')
                endpoint=correct(target,base,y,xp,g['v'],.02,folder/'FullJ20',max_states=20)
                row=metrics(target,endpoint,base,y,g['weak_material'],.02)
                row.update(parent=parent,cost=counters(book),profile_qualification='BUDGETED20_NOT_PROVEN_CONVERGED',GT_used=False)
                rows.append(row);write(folder/'RESULT.json',row)
            finally:target.close()
    write(OUT/'E2/controls/profiles/SUMMARY.json',rows)


def score():
    seal=require_seal()
    rows=[]
    # This is the only R7 evaluation process that opens hidden truth.
    for parent in seal['parents']:
        d=OUT/'E2/parents'/parent
        with np.load(OUT/'raw/offline_truth'/f'{parent}.npz',allow_pickle=False) as z:
            truth=z['material'].copy();xtrue=z['x'].copy()
        for arm in ['ANCHOR','A_FROZEN','B_FULLJ','C_JETS','D_RANDOM']:
            path=d/'ANCHOR.npz' if arm=='ANCHOR' else d/arm/'RMS_0.0200/ENDPOINT.npz'
            if not path.exists():
                rows.append(dict(parent=parent,arm=arm,status='UNAVAILABLE_NO_GT_SUBSTITUTION'));continue
            with np.load(path,allow_pickle=False) as z:
                m=z['material'].copy();x=z['x'].copy()
            if m.dtype!=np.complex128 or truth.dtype!=np.complex128 or x.dtype!=np.float64:raise ValueError('dtype contract')
            # ROI is determined by the public geometry/config, never by truth phase.
            from p3r4.core import MaterialPrior,ROOT as R4,OUT as R4OUT
            cfg=json.loads((R4OUT/'configs/TARGET.json').read_text())
            points=json.loads((R4/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
            prior=MaterialPrior(cfg,points)
            e=m-truth;roi=prior.mask
            rows.append(dict(parent=parent,arm=arm,status='SCORED_AFTER_SEAL',
                full_MSE=float(np.mean(abs(e)**2)),full_MAE=float(np.mean(abs(e))),
                ROI_MSE=float(np.mean(abs(e[roi])**2)),ROI_MAE=float(np.mean(abs(e[roi]))),
                latent_distance_to_GT=float(np.linalg.norm(x-xtrue)),
                predicted_material_tasks=tasks(prior,m).tolist(),GT_material_tasks=tasks(prior,truth).tolist()))
    write(OUT/'E2/OFFLINE_GT_SCORES.json',dict(scored_utc=time.time(),sealed_utc=seal['utc'],
       rows=rows,GT_selected_candidate=False,primary_imaging_action='ANCHOR remains default; candidate tasks are geometry diagnostics',
       qualification='8-parent offline geometry only; not Bayesian risk, no oracle best-of-K claim'))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['gauge','profiles','score'])
    a=parser.parse_args()
    {'gauge':gauge_controls,'profiles':profiles,'score':score}[a.command]()
