#!/usr/bin/env python3
"""Post-terminal saved-array A0 analysis; zero PDE, no online truth access.

Candidate-only noiseless pair search and optimization diagnostics. Execute only
through the accounting wrapper; this file never launches a physical solver.
"""
from pathlib import Path
import argparse,csv,json
from collections import defaultdict
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
LOWER=np.array([1.5,0.,-1.,1.6,-.01,-.01])
UPPER=np.array([6.,.25,1.,2.6,.01,.01])
SEED=2026102510


def realpack(a):
    a=np.asarray(a).ravel();return np.r_[a.real,a.imag]/.01


def complex_json(value):
    a=np.asarray(value)
    if np.iscomplexobj(a):return a.astype(np.complex128)
    if a.ndim<1 or a.shape[-1]!=2:raise ValueError('Expected complex [real,imag] JSON encoding')
    return np.asarray(a[...,0],float)+1j*np.asarray(a[...,1],float)


def field_grid(support,eps,ell):
    s=np.asarray(support,bool)
    if s.shape not in ((2,2),(3,3),(6,24)):raise ValueError('Unexpected support grid')
    e=np.broadcast_to(np.asarray(eps,float),s.shape);l=np.broadcast_to(np.asarray(ell,float),s.shape)
    ry,rx=6//s.shape[0],24//s.shape[1]
    def expand(a):return np.repeat(np.repeat(a,ry,0),rx,1)
    return expand(np.where(s,e-1,0.)),expand(np.where(s,l,0.)),expand(s)


def json_default(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--output',type=Path,default=None)
    parser.add_argument('--clean-mean-convention',choices=['includes-gamma','unmultiplied'],default='includes-gamma',
        help='Default verified from a0_prepare: mean=gamma*ev.data; means.append(mean).')
    args=parser.parse_args();root=args.root.resolve()
    fitpath=root/'results/a0/optimization/RESULTS.json';scorepath=root/'results/a0/scoring/RESULTS.json'
    fits=json.loads(fitpath.read_text());scoring=json.loads(scorepath.read_text())
    if fits.get('status') in ('RUNNING',None) or scoring.get('status') in ('RUNNING',None):
        raise RuntimeError('Both fitting and separate scoring must be terminal before offline analysis')
    records=fits['records'];scores=scoring['rows']
    if len(records)!=72 or len(scores)!=72:raise RuntimeError('Expected all72 frozen method/observation records, including failures')
    indexed={(int(r['observation_id']),r['method']):r for r in scores}
    if len(indexed)!=72:raise ValueError('Duplicate scoring identities')
    labels=json.loads((root/'inputs/offline/CONTINUOUS_TRUTHS.json').read_text())
    truth={int(r['parent_id']):r for r in labels}
    with np.load(root/'inputs/offline/CONTINUOUS_CLEAN_MEANS.npz',allow_pickle=False) as z:clean=z['data'].copy()
    if clean.shape!=(12,4,16):raise ValueError('Expected12parents×4TX×16RX clean means')
    with np.load(root/'inputs/online/A0_OBSERVATIONS.npz',allow_pickle=False) as z:
        observations=z['base'].copy();parent_ids=z['parent_id'].copy();repeat_ids=z['noise_repeat'].copy()
    with np.load(root/'inputs/online/COARSE_GRID.npz',allow_pickle=False) as z:
        legal=z['legal_supports'].copy();support_prior=z['support_prior'].copy()
    # This read occurs only after terminal results. gamma is applied exactly once.
    true_gamma={p:(1+r['instrument']['amplitude_offset'])*np.exp(1j*r['instrument']['phase']) for p,r in truth.items()}
    if args.clean_mean_convention=='unmultiplied':
        clean=np.stack([true_gamma[p]*clean[p] for p in range(12)])
    out=(args.output or root/'results/a0/analysis').resolve();out.mkdir(parents=True,exist_ok=False)
    def resolve(path):
        p=Path(path)
        if p.exists():return p
        parts=p.parts
        if 'results' in parts:
            q=root.joinpath(*parts[parts.index('results'):])
            if q.exists():return q
        raise FileNotFoundError(str(p))
    rows=[];candidates=[];arrays={};missing=[]
    for record in records:
        obs=int(record['observation_id']);parent=int(record['parent_id']);repeat=int(record['noise_repeat']);method=record['method']
        if int(parent_ids[obs])!=parent or int(repeat_ids[obs])!=repeat:raise ValueError('Observation identity mismatch')
        score=indexed[(obs,method)];selected=record.get('selection')
        row=dict(observation_id=obs,parent_id=parent,noise_repeat=repeat,method=method,
            online_status=record.get('status'),score_status=score.get('score_status'),
            fit_status=selected.get('fit_status') if selected else None,
            primary_fit_converged=bool(selected.get('success')) if selected else None,
            scipy_status=selected.get('status') if selected else None,
            selected_nfev=selected.get('nfev') if selected else None,
            actual_method_nfev=record.get('actual_nfev'),actual_method_njev=record.get('actual_njev'),
            optimizer_message=selected.get('message') if selected else None,
            scipy_optimality='NOT_SAVED_IN_FIT_SCHEMA',gradient_available=False,
            finer_status=score.get('finer_status'),finer_prediction_difference_white=score.get('finer_prediction_difference_white'),
            training_noise_envelope_compatible=score.get('training_noise_envelope_compatible'),
            holdout_noise_envelope_compatible=score.get('holdout_noise_envelope_compatible'),
            both_noise_envelopes_compatible=score.get('both_noise_envelopes_compatible'),
            training_noise_energy=score.get('training_noise_energy'),holdout_noise_energy=score.get('holdout_noise_energy'),
            noiseless_mean_scope='OFFLINE_ALL4TX_INSTRUMENT_MEAN; not noisy compatibility or algorithm acceptance',
            selected_candidate_available=selected is not None,
            noiseless_pair_status='NOT_EVALUABLE_MISSING_SAVED_PREDICTIONS')
        label=truth[parent]
        instrument=label['instrument'];prior_penalty=.5*((instrument['amplitude_offset']/.005)**2+(instrument['phase']/.005)**2)
        true_raw=realpack(clean[parent,:2]-observations[obs])
        truth_objective=float(.5*true_raw@true_raw+prior_penalty)
        matches=[i for i,s in enumerate(legal) if np.array_equal(s,np.asarray(label['support'],bool))]
        if len(matches)!=1:raise ValueError('True support does not match declared2x2 prior')
        true_hier=truth_objective-float(np.log(support_prior[matches[0]]))+float(np.log(2.))
        row.update(true_cached_data_energy=float(true_raw@true_raw),true_cached_nuisance_prior_penalty=float(prior_penalty),
                   true_cached_MAP_objective=truth_objective,true_cached_hierarchical_MAP_objective=true_hier,
                   truth_objective_scope='cached true physical mean and actual noisy base data; offline optimization diagnostic only')
        if selected:
            theta=np.asarray(selected['theta'],float)
            row.update(selected_MAP_objective=float(selected['objective']),
                selected_hierarchical_MAP_objective=float(selected['hierarchical_map_objective']),
                fit_minus_true_MAP_objective=float(selected['objective'])-truth_objective,
                fit_minus_true_hierarchical_MAP_objective=float(selected['hierarchical_map_objective'])-true_hier,
                target_pose_change=float(np.linalg.norm(theta[2:4]-np.asarray(label['target_pose']))),
                fitted_gain_offset=float(theta[4]),fitted_phase=float(theta[5]),true_gain_offset=float(instrument['amplitude_offset']),
                true_phase=float(instrument['phase']),distance_to_lower_bounds=theta-LOWER,distance_to_upper_bounds=UPPER-theta)
            contrast,loss,mask=field_grid(selected['support'],theta[0],theta[1])
            tc,tl,tm=field_grid(label['support'],label['epsinf'],label['ell'])
            row.update(material_contrast_L2=float(np.sqrt(.125**2*np.sum((contrast-tc)**2))),
                material_loss_L2=float(np.sqrt(.125**2*np.sum((loss-tl)**2))),
                support_symmetric_difference_area=float(.125**2*np.sum(mask!=tm)),
                support_iou=float(np.sum(mask&tm)/np.sum(mask|tm)) if np.any(mask|tm) else 1.)
            basepred=complex_json(selected['pred']).reshape(2,16)
            try:
                path=resolve(score['prediction_path'])
                with np.load(path,allow_pickle=False) as z:
                    outer=z['holdout_prediction'].copy();saved_base=z['train_prediction'].copy()
                if basepred.shape!=(2,16) or outer.shape!=(2,16):raise ValueError('Unexpected source layout')
                if not np.array_equal(basepred,saved_base):raise ValueError('Selected base prediction differs from primary score archive')
                mean=np.concatenate((basepred,outer),axis=0)
                delta=mean-clean[parent];distance=float(np.linalg.norm(realpack(delta)))
                candidate=bool(distance<=1. and row['target_pose_change']>=.1)
                row.update(noiseless_joint_mean_white_distance=distance,
                    noiseless_base_mean_white_distance=float(np.linalg.norm(realpack(delta[:2]))),
                    noiseless_outer_mean_white_distance=float(np.linalg.norm(realpack(delta[2:]))),
                    noiseless_pair_status='CANDIDATE_ONLY_PENDING_BOTH_WORLD_INDEPENDENT_NUMERICAL_VALIDATION' if candidate else 'NOT_FOUND_FOR_THIS_SAVED_ENDPOINT',
                    offline_pair_candidate=candidate,
                    candidate_condition='d_all4_mean<=1 and target_pose_change>=.1; screening only, not a scientific gate')
                key=f'obs{obs:02d}_{method}'
                arrays[key+'_fitted_all4_mean']=mean;arrays[key+'_true_all4_mean']=clean[parent]
                arrays[key+'_noiseless_difference']=delta
                if candidate:candidates.append(dict(row))
            except (KeyError,FileNotFoundError,ValueError) as error:
                row['noiseless_pair_error']=type(error).__name__+': '+str(error)
                missing.append(dict(observation_id=obs,method=method,kind='primary_prediction',error=str(error)))
            # Exact final-coordinate cached J only. No interpolation, new forward,
            # finite differences, or assumed optimizer terminal optimality.
            try:
                trace=resolve(selected['trace_path']);final=trace.with_name(trace.stem+'_final.npz')
                if not final.exists():raise FileNotFoundError(str(final))
                with np.load(final,allow_pickle=False) as z:
                    keys=[k for k in z.files if k.endswith('_physical_theta')]
                    match=[k for k in keys if np.array_equal(z[k],theta[:4])]
                    if len(match)!=1:raise ValueError('No unique exact final physical coordinate in cache')
                    prefix=match[0][:-len('physical_theta')]
                    if prefix+'J' not in z:raise KeyError('Final physical Jacobian was not saved')
                    physical_j=z[prefix+'J'].copy();f=z[prefix+'data'].copy()
                gamma=(1+theta[4])*np.exp(1j*theta[5]);phase=np.exp(1j*theta[5])
                cj=np.column_stack((gamma*physical_j,phase*f.ravel(),1j*gamma*f.ravel()))
                jac=np.vstack((cj.real,cj.imag))/.01
                prior_j=np.zeros((2,6));prior_j[:,4:]=np.eye(2)/.005
                jac=np.vstack((jac,prior_j))
                resid=np.r_[realpack(gamma*f.reshape(2,16)-observations[obs]),theta[4:]/.005]
                grad=jac.T@resid
                mapping=theta-np.clip(theta-grad,LOWER,UPPER)
                row.update(gradient_available=True,gradient=grad,gradient_infinity_norm=float(np.max(np.abs(grad))),
                    projected_gradient_mapping_step1=mapping,projected_gradient_mapping_infinity_norm=float(np.max(np.abs(mapping))),
                    gradient_scope='reconstructed exact final cached J; physical-coordinate step1 mapping, not saved scipy optimality or convergence gate')
            except (KeyError,ValueError,FileNotFoundError) as error:
                row['gradient_unavailable_reason']=type(error).__name__+': '+str(error)
        rows.append(row)
    def write_csv(name,data):
        if not data:return
        keys=list(dict.fromkeys(k for row in data for k in row))
        with (out/(name+'.csv')).open('w',newline='') as f:
            wr=csv.DictWriter(f,keys);wr.writeheader()
            for row in data:
                wr.writerow({k:json.dumps(v,default=json_default) if isinstance(v,(list,dict,np.ndarray)) else v for k,v in row.items()})
    write_csv('ALL72_OFFLINE_ANALYSIS',rows)
    # Explicit headers remain present when no candidates exist.
    if candidates:write_csv('CANDIDATE_ONLY_PAIRS',candidates)
    else:(out/'CANDIDATE_ONLY_PAIRS.csv').write_text('observation_id,parent_id,noise_repeat,method,noiseless_joint_mean_white_distance,target_pose_change,status\n')
    np.savez_compressed(out/'SAVED_MEAN_PAIRS.npz',**arrays)
    metrics=['noiseless_joint_mean_white_distance','target_pose_change','material_contrast_L2','material_loss_L2',
             'support_iou','fit_minus_true_hierarchical_MAP_objective','training_noise_envelope_compatible',
             'holdout_noise_envelope_compatible','both_noise_envelopes_compatible','primary_fit_converged']
    groups=defaultdict(list)
    for row in rows:groups[(row['parent_id'],row['method'])].append(row)
    parentrows=[]
    for (parent,method),items in sorted(groups.items()):
        row=dict(parent_id=parent,method=method,saved_repeats=len(items))
        for metric in metrics:
            valid=[float(r[metric]) for r in items if r.get(metric) is not None and np.isfinite(float(r[metric]))]
            row[metric]=float(np.mean(valid)) if valid else None;row[metric+'_valid_repeats']=len(valid)
        parentrows.append(row)
    bootstrap=[]
    for metric in metrics:
        values=defaultdict(dict)
        for row in parentrows:
            if row[metric+'_valid_repeats']==2:values[row['method']][row['parent_id']]=row[metric]
        for method in sorted(values):
            ids=sorted(set(values[method])&set(values.get('FULL_J_TR',{})))
            differences=np.array([values[method][p]-values['FULL_J_TR'][p] for p in ids])
            row=dict(metric=metric,method=method,baseline='FULL_J_TR',paired_parent_count=len(ids),
                     descriptive_only=True,draws=2000,seed=SEED,
                     mean_difference=float(differences.mean()) if len(ids) else None)
            if len(ids)==12:
                rng=np.random.default_rng(SEED);draws=differences[rng.integers(0,12,size=(2000,12))].mean(axis=1)
                row.update(low=float(np.quantile(draws,.025)),high=float(np.quantile(draws,.975)),status='DESCRIPTIVE_12_PARENT_PAIRED')
            else:row.update(low=None,high=None,status='NOT_RUN_INCOMPLETE_12_PARENT_REPEAT_PAIRS')
            bootstrap.append(row)
    methodstats=[]
    for method in sorted(set(r['method'] for r in rows)):
        rr=[r for r in rows if r['method']==method]
        row=dict(method=method,records=len(rr),selected=sum(r['selected_candidate_available'] for r in rr),
            converged=sum(r['primary_fit_converged'] is True for r in rr),gradient_available=sum(r['gradient_available'] for r in rr),
            candidate_only_pairs=sum(r.get('offline_pair_candidate') is True for r in rr),
            finer_checks=sum(r['finer_status']=='SCORED' for r in rr))
        for key in ('training_noise_envelope_compatible','holdout_noise_envelope_compatible','both_noise_envelopes_compatible'):
            row[key+'_count']=sum(r.get(key) is True for r in rr);row[key+'_scored']=sum(r.get(key) is not None for r in rr)
        for metric in metrics:
            pp=[r[metric] for r in parentrows if r['method']==method and r[metric+'_valid_repeats']==2]
            row[metric+'_parent_mean']=float(np.mean(pp)) if pp else None;row[metric+'_parent_count']=len(pp)
        methodstats.append(row)
    fine=[r['finer_prediction_difference_white'] for r in rows if r['finer_status']=='SCORED' and r['finer_prediction_difference_white'] is not None]
    evaluable=sum('noiseless_joint_mean_white_distance' in row for row in rows)
    summary=dict(status='OFFLINE_SAVED_ARRAY_ANALYSIS_COMPLETE' if evaluable==72 else 'INCOMPLETE_NOISELESS_PAIR_ANALYSIS',
        noiseless_pair_evaluable_count=evaluable,noiseless_pair_diagnosis_complete=evaluable==72,
        zero_new_PDE=True,new_rhs=0,new_factorizations=0,
        main_terminal_status=fits['status'],scoring_terminal_status=scoring['status'],record_count=len(rows),
        clean_mean_convention=args.clean_mean_convention,
        clean_mean_provenance='tools/a0_prepare.py: mean=gamma*ev.data; means.append(mean); default already includes true gamma',
        extra_gamma_multiplication_performed=args.clean_mean_convention=='unmultiplied',
        candidate_only_pair_count=len(candidates) if evaluable==72 else None,candidate_condition='d_true_vs_fitted_all4_noiseless_mean<=1 and pose_change>=.1',
        candidates_are_not_validated_witnesses=True,both_world_independent_numerical_validation='NOT_RUN_BY_THIS_ZERO_PDE_TOOL',
        numerical_fixed_parent_checks=dict(saved_count=len(fine),predeclared_maximum=24,
            count_exceeds_predeclared_maximum=len(fine)>24,maximum_saved_difference=max(fine) if fine else None,
            scope='candidate outer-prediction refinement only; not independent validation of both all4-world means'),
        method_stats=methodstats,missing_saved_inputs=missing,
        optimization_diagnosis='Fit versus cached true noisy objective can expose optimization limitations; neither objective gap nor failed fitting proves physical nonexistence or identifiability.',
        no_new_acceptance_gate=True,no_algorithm_selected_by_GT=True,online_use_forbidden=True,
        costs='This array/preparation CPU must be charged by run_job; zeroPDE does not mean zeroCPU.')
    write_csv('PARENT_METHOD_MEANS',parentrows);write_csv('METHOD_STATS',methodstats);write_csv('PAIRED_BOOTSTRAP',bootstrap)
    (out/'RESULTS.json').write_text(json.dumps(dict(summary=summary,rows=rows,candidate_only_pairs=candidates,
        parent_method_means=parentrows,paired_bootstrap=bootstrap),indent=2,allow_nan=False,default=json_default)+'\n')
    print(json.dumps(summary,allow_nan=False,default=json_default),flush=True)


if __name__=='__main__':main()
