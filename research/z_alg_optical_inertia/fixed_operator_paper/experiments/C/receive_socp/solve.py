"""Post-study single-beam receiver-direction proposals from frozen old data."""
import sys
from pathlib import Path
from time import perf_counter, process_time

import cvxpy as cp
from flint import arb, ctx
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
C=HERE.parent
PAPER=C.parents[1]
sys.path.insert(0,str(C))

from design import envelope, form_border, source_pair, split_query
from model import (ANGLES, OLD, geometry, hash_arrays, hash_file, herm, mid, noise_radius,
                   read_json, write_json)

SOURCE=C/'diagnostic'
INPUTS=('registration.json','proposal.json','public_model.npz','full_menu_replay.json','menu_query_replay.npz')


def bundle():
    proposal=read_json(SOURCE/'proposal.json')
    registry=read_json(SOURCE/'registration.json')
    with np.load(SOURCE/'public_model.npz',allow_pickle=False) as public:
        for plan in proposal['source_rows']:
            scene,contract=plan['scene'],plan['energy_contract']
            with np.load(SOURCE/'old'/f'{scene}__{contract}.npz',allow_pickle=False) as data:
                y=data['Y']
            if hash_arrays(y)!=plan['old_data_hash']:
                raise ArithmeticError('A frozen old-data block changed.')
            lo,hi=[public[scene+'__'+k] for k in ('lo','hi')]
            if hash_arrays(lo,hi)!=plan['public_candidate_cover']:
                raise ArithmeticError('Frozen candidate boxes changed.')
            g=geometry(plan['layout'])
            if hash_arrays(g.g,g.s,g.b)!=plan['operator_model_hash']:
                raise ArithmeticError('Frozen operator recipe changed.')
            old_power=arb(1) if contract=='added' else arb(1)/2
            nr=noise_radius(registry,scene)
            z=source_pair(OLD,old_power.sqrt(),y)
            options=[]
            proofs=[]
            for i,entry in enumerate(plan['old_classes'].values()):
                chosen=entry['chosen_envelope']
                if not chosen:
                    raise ArithmeticError('The frozen condition had no bounded envelope.')
                state, proof=form_border(g,lo[i],hi[i],chosen['weight'])
                if state is None:
                    raise ArithmeticError('Cannot independently rebuild the frozen state border.')
                decision,meta=envelope(state,g,z,nr)
                if decision is None:
                    raise ArithmeticError('Cannot independently rebuild the frozen envelope.')
                options.append(decision)
                proofs.append(dict(state=proof,conditional=meta))
            yield plan,g,nr,options,proofs


def solve_direction(option,angle,sigma):
    """Solve the unit-norm SOCP on the same two fixed receiving ellipsoids."""
    centers,radii,transforms=[],[],[]
    for side in option:
        centers.append(mid(side.center)[:,angle])
        radii.append(np.sqrt(max(0.,float(side.gram[angle,angle].real.mid()))))
        vals,vec=eigh(herm(mid(side.inverse)),check_finite=False)
        transforms.append(np.sqrt(np.maximum(vals,0.))[:,None]*vec.conj().T)
    d=centers[1]-centers[0]
    x=cp.Variable(len(d),complex=True)
    caps=cp.Variable(3,nonneg=True)
    bill=radii[0]*caps[0]+radii[1]*caps[1]+2*sigma*caps[2]
    constraints=[cp.norm(transforms[0]@x)<=caps[0],cp.norm(transforms[1]@x)<=caps[1],
                 cp.norm(x)<=caps[2],caps[2]<=1.]
    p=cp.Problem(cp.Maximize(cp.real(d.conj()@x)-bill),constraints)
    value=p.solve(solver='CLARABEL',max_iter=200,tol_gap_abs=1.e-10,tol_gap_rel=1.e-10,tol_feas=1.e-10)
    if x.value is None or not np.all(np.isfinite(x.value)):
        return None,dict(status=p.status,objective_proposal=None)
    proposed=np.asarray(x.value,complex)
    # Preserve a binary64 vector strictly inside the requested norm bound.
    w=(proposed*((1.-2.**-45)/max(1.,np.linalg.norm(proposed))))[:,None]
    return w,dict(status=p.status,objective_proposal=float(value),iterations=p.solver_stats.num_iters,
                  candidate_norm_numeric=float(np.linalg.norm(w)),solver_seconds=p.solver_stats.solve_time)


def main():
    wall,cpu=perf_counter(),process_time()
    outputs={}
    data=[]
    frozen=read_json(SOURCE/'full_menu_replay.json')
    audited={(entry['scene'],entry['contract']):entry['by_batch']['1'] for entry in frozen['rows']}
    with np.load(SOURCE/'menu_query_replay.npz',allow_pickle=False) as directions:
        for plan,g,nr,options,_ in bundle():
            scene,contract=plan['scene'],plan['energy_contract']
            prior=audited[scene,contract]
            assert len(prior)==10
            for index,row in enumerate(prior):
                started=perf_counter()
                angle=row['source_indices'][0]
                if [angle]!=plan['batches']['1']['all_action_menu'][index]:
                    raise ArithmeticError('Changed frozen physical source menu.')
                _,center=split_query(options,(angle,),arb(1),nr)
                if center is None:
                    raise ArithmeticError('No center direction for a frozen menu candidate.')
                old=row['verified_support']
                if old is not None:
                    key=old['support_query_index']
                    archived=directions[key]
                    if hash_arrays(archived)!=old['support_query_hash'] or not np.array_equal(center,archived):
                        raise ArithmeticError('Center-directed W disagrees with the actual frozen archive.')
                    center=archived
                try:
                    optimal, proposed=solve_direction(options,angle,float(arb(nr)))
                except cp.error.SolverError as exc:
                    optimal,proposed=None,dict(status='SOLVER_FAILED',reason=str(exc))
                label=f'{scene}__{contract}__{angle}'
                outputs[label+'__center']=center
                if optimal is not None:
                    outputs[label+'__socp']=optimal
                data.append(dict(key=label,scene=scene,contract=contract,source_index=angle,
                   source_angle=ANGLES[angle],original_archive_status=row['status'],
                   center_provenance='archived_binary_direction' if old else 'same_frozen_center_rule_reconstructed_from_old',
                   proposal=proposed,direction_seconds=perf_counter()-started))
    np.savez_compressed(HERE/'directions.npz',**outputs)
    write_json(HERE/'proposals.json',dict(status='POST_STUDY_OLD_ONLY_COUNTERFACTUAL',source_round='diagnostic',
        future_data_read=False,source_input_identities={x:hash_file(SOURCE/x) for x in INPUTS},
        method='CVXPY / Clarabel; numerical direction proposals only',
        proposal_bits=ctx.prec,rows=data,wall_seconds=perf_counter()-wall,cpu_seconds=process_time()-cpu))
    print(f'Saved {len(data)} old-only single-beam direction pairs.')


if __name__=='__main__':
    ctx.prec=128
    ctx.threads=1
    with threadpool_limits(1):
        main()
