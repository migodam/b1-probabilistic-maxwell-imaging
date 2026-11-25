"""Compute and freeze all beam proposals with access to the old data only."""
import argparse
from datetime import datetime, timezone
import json

import numpy as np
from scipy.linalg import lu_solve
from threadpoolctl import threadpool_limits

from beam_common import (HERE, OLD, clocks, elapsed, atomic_write, require_absent,
                         save, digest, geometry, forward, normalize, herm, eigh)
from quadratic import ResidualForm, quad_search, envelope


POLICIES=('uniform_angle','random_angle','fisher_angle','direct_coherent','source_eigenbeam','source_shaped')
FLUENCE_TARGET=2.40


def robust_response(y, d, width, field, eta, action, direction):
    half=np.linalg.norm((direction.conj()@d)*width)*np.linalg.norm(field@action)/(1-eta)
    center=np.vdot(direction,y@action).real
    return float(center-half),float(center+half)


def for_class(g,s,b,center,width):
    y,e,d,j,born,lu=forward(g,s,b,center,auxiliary=True)
    t=lu_solve(lu,g.T,trans=1,check_finite=False).T
    eta=np.linalg.norm(t*width[None,:])
    return y,e,d,t,eta,j,born


def exact_best_angle(delta, gram):
    cov=delta.conj().T@delta
    ev,vec=eigh(herm(cov[1:,1:]),gram[1:,1:],check_finite=False)
    return np.r_[0.,vec[:,-1]],float(ev[-1])


def support_score(pair, widths, action, sigma):
    diff=(pair[1][0]-pair[0][0])@action
    if np.linalg.norm(diff)==0:
        return -np.inf,np.zeros(len(diff))
    direction=diff/np.linalg.norm(diff)
    ends=[robust_response(y,d,widths,e,eta,action,direction) for y,e,d,_,eta,*_ in pair]
    return float(ends[1][0]-ends[0][1]-2*sigma),direction


def optical_score(candidates, action, direction, m, sigma):
    intervals=[]
    for options in candidates:
        proposals=[]
        for _,record,q,mode,source in options[:12]:
            receive,center,extra,_=envelope(q,m,action)
            if extra<0:
                continue
            reach=np.sqrt(extra*max(0,np.vdot(direction,np.linalg.solve(receive,direction)).real))
            projection=np.vdot(direction,center).real
            proposals.append((float(projection-reach),float(projection+reach),record))
        if not proposals:
            return None
        intervals.append(proposals)
    idx1=min(intervals[0],key=lambda v:v[1])
    idx2=max(intervals[1],key=lambda v:v[0])
    return float(idx2[0]-idx1[1]-2*sigma),idx1[2],idx2[2]


def design_case(geom,gram,bank,record,Yold):
    tick=clocks()
    i=str(record['index'])
    cells,centers,widths,fg,bg=[bank[i+'__'+v] for v in ('cells','centers','widths','fg','bg')]
    g,s,b=geom.restricted(cells)
    pair=[for_class(g,s,b,ctr,widths) for ctr in centers]
    y0,y1=pair[0][0],pair[1][0]
    delta=y1-y0
    direct,eigenvalue=exact_best_angle(delta,gram)
    r_old=record['old_radius']
    sigma=record['candidate_new_noise_radius']
    fluence=record.get("fixed_fluence", (FLUENCE_TARGET*sigma/np.sqrt(eigenvalue))**2)
    sources=b.shape[1]
    standard=np.eye(sources,dtype=complex)
    uniform=(1+sources)//2
    randomized=int(np.random.default_rng(np.random.SeedSequence([202611260472,int(i)])).integers(1,sources))
    # A local Fisher score for the material-contrast direction with two uncertain
    # background amplitudes. All columns refer to one co-normalized plain beam.
    nominal=(centers[0]+centers[1])/2
    _,field,d,*_=forward(g,s,b,nominal,auxiliary=True)
    d_fish=np.einsum('pi,is,i->ps',d,field,centers[1]-centers[0],optimize=True)
    bg_ids=np.flatnonzero(bg)
    two=np.array_split(bg_ids,2)
    nuis=np.stack([np.einsum('pi,is->ps',d[:,ids],field[ids],optimize=True)*.045 for ids in two],axis=-1)
    fish=[]
    for j in range(sources):
        c=herm(sigma**2*np.eye(s.shape[0])+nuis[:,j,:]@nuis[:,j,:].conj().T)
        fish.append(float(np.vdot(d_fish[:,j],np.linalg.solve(c,d_fish[:,j])).real/gram[j,j].real))
    fisher=1+np.argmax(fish[1:])
    base=dict(uniform_angle=standard[:,uniform],random_angle=standard[:,randomized],
              fisher_angle=standard[:,fisher],direct_coherent=direct)
    option_sets=[]
    all_trials=[]
    for y,e,d,t,eta,*_ in pair:
        row_power=np.linalg.norm(d,axis=0)*np.linalg.norm(e,axis=1)
        examples=np.flatnonzero(fg)
        # Four observable cells and four weak cells detect both useful and
        # dark optical directions with a finite predeclared proposal budget.
        ordered=examples[np.argsort(row_power[examples])]
        active=np.r_[ordered[:4],ordered[-4:]]
        form=ResidualForm(widths,e,y,t,d)
        opts,trials=quad_search(form,active,Yold,r_old,sigma,gram)
        option_sets.append(opts)
        all_trials.append(trials)
    modes=[item[3] for options in option_sets for item in options[:12]]
    if modes:
        beam_mode=max(modes,key=lambda a:support_score(pair,widths,normalize(a,gram,fluence),sigma)[0])
    else:
        beam_mode=direct
    base['source_eigenbeam']=beam_mode
    pool=[direct,*modes]
    for mode in modes:
        phase=np.angle(np.vdot(direct,gram@mode))
        phase_matched=mode*np.exp(-1j*phase)
        pool.extend((direct+ .15*phase_matched,direct-.15*phase_matched,direct+.35j*phase_matched))
    candidates=[normalize(a,gram,fluence) for a in pool]
    opt_record=[]
    for a in candidates:
        gap,w=support_score(pair,widths,a,sigma)
        optical=optical_score(option_sets,a,w,sources,sigma)
        # Use old-conditioned optical form if it can certify a positive
        # two-sided threshold. Otherwise use the full-wave secant support.
        opt_record.append((gap,optical))
    score=lambda item: (item[1] is not None and item[1][0]>0,
                        item[1][0] if item[1] is not None and item[1][0]>0 else item[0])
    winner=max(range(len(candidates)),key=lambda k:score(opt_record[k]))
    base['source_shaped']=candidates[winner]
    actions=np.column_stack([normalize(base[name],gram,fluence) for name in POLICIES])
    output={}
    dirs=[]
    q_per_action=[]
    for j,name in enumerate(POLICIES):
        action=actions[:,j]
        gap,direction=support_score(pair,widths,action,sigma)
        dirs.append(direction)
        trial=optical_score(option_sets,action,direction,sources,sigma)
        q_per_action.append(np.stack([next(v[2] for v in block if v[1] is rec) for block,rec in zip(option_sets,trial[1:])]) if trial is not None else np.full((2,2*sources,2*sources),np.nan))
        output[name]=dict(nominal_incident_fluence=float(np.vdot(action,gram@action).real),
                          old_compatible_center_witness_new_distance=float(np.linalg.norm(delta@action)),
                          normalized_predicted_box_direction_gap=gap/sigma,
                          both_center_witnesses_ambiguous_at_midpoint=bool(np.linalg.norm(delta@action)<2*sigma),
                          optical_old_conditioned_gap=None if trial is None else trial[0],
                          optical_support_candidate=None if trial is None else list(trial[1:]))
    stats=dict(policies=list(POLICIES),selected_policy='source_shaped',action_policy=list(output),
               selected_source_was_optical=bool(opt_record[winner][1] is not None and opt_record[winner][1][0]>0),
               eigenbeam_fallback_to_direct=not bool(modes),
               source_shaped_equals_direct=bool(abs(np.vdot(actions[:,3],gram@actions[:,-1]))/fluence > 1-1.e-12),
               support_candidates=len(candidates),baselines_angular_index=dict(uniform=uniform,random=randomized,fisher=int(fisher)),
               same_fluence=fluence, nominal_fluence_calibration=('fixed across independent confirmation cohort' if 'fixed_fluence' in record else 'pre-acquisition old-compatible class-center contrast at 2.40 new noise radii for ideal direct one-shot coherent beam'),
               old_center_residual_over_radius=[float(np.linalg.norm(y[:,OLD]-Yold)/r_old) for y in (y0,y1)],
               prescribed_noise_balls={'old':r_old,'new':sigma},
               nominal_Born_error_bin_full_coherent_menu=[float(np.linalg.norm(p[0]-p[6])/np.linalg.norm(p[0])) for p in pair],
               class_material_difference_rmse_active=float(np.sqrt(np.mean((centers[0]-centers[1])**2))),
               old_threshold_uses_only_r_old_squared=True,
               prospective_quadratic_is_for_clean_future=True,
               stored_conditional_q_trial_forms=[[v[1] for v in opt[:8]] for opt in option_sets],
               optical_trials=all_trials, viable_optical_count=[len(o) for o in option_sets], designs=output,
               **elapsed(tick))
    return stats,dict(actions=actions,directions=np.column_stack(dirs),
                      center_fields=np.stack([p[0] for p in pair]),
                      conditional_clean_q_by_action=np.stack(q_per_action),
                      conditional_clean_q_trials=np.stack([np.stack([v[2] for v in options[:8]]) for options in option_sets]) if all(len(o)>=8 for o in option_sets) else np.empty((2,0,2*sources,2*sources)),
                      potential_fisher=np.array(fish), old=Yold,
                      optical_candidate_source_modes=np.column_stack(modes) if modes else np.empty((sources,0)),
                      all_support_candidates=np.column_stack(candidates),
                      all_candidate_secant_gaps=np.array([r[0] for r in opt_record]))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('side',type=int,choices=(32,96))
    parser.add_argument('--tag',default='main')
    args=parser.parse_args()
    folder=HERE/f'{args.tag}_{args.side}'
    require_absent(folder/'proposal.json',folder/'proposal_arrays.npz',folder/'future_readings.npz')
    registered=json.loads((folder/'registration.json').read_text())
    records=[rec for rec in registered['cases'] if rec['selected']]
    bank=np.load(folder/'public_classes.npz')
    old=np.load(folder/'old_readings.npz')
    geom=geometry(args.side)
    gram=bank['power']
    start=clocks()
    proposals={}
    arrays={}
    for rec in records:
        idx=str(rec['index'])
        props,data=design_case(geom,gram,bank,rec,old[idx+'__Y'][:,0])
        proposals[idx]=props
        arrays.update({idx+'__'+k:v for k,v in data.items()})
        print(idx,props['viable_optical_count'],list(np.round([p['normalized_predicted_box_direction_gap'] for p in props['designs'].values()],4)),flush=True)
    save(folder/'proposal_arrays.npz',**arrays)
    atomic_write(folder/'proposal.json',dict(stage='FROZEN_BEFORE_FUTURE',
                   frozen_utc=datetime.now(timezone.utc).isoformat(),
                   old_file_digest=digest(np.frombuffer((folder/'old_readings.npz').read_bytes(),dtype=np.uint8)),
                   design_binary_digest=digest(np.frombuffer((folder/'proposal_arrays.npz').read_bytes(),dtype=np.uint8)),
                   rows=proposals,**elapsed(start)))


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
