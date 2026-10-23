"""Frozen P1B OFFLINE_ONLY truth scorer; no inference, physics or Book imports.

Only main may open the pre-existing OFFLINE_TRUTH.npz, after inference freeze
and point-validation completion markers have been checked. It never creates
truth/observations or changes predictions. Statistical output is not a GO vote.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
ARMS=('direct_point','progressive_point','direct_SMC_mean','progressive_SMC_mean')
BAYES='continuous_quadrature_Bayes_mean'
POINT_NAMES={'direct_point':'direct joint point','progressive_point':'progressive block point'}
SMC_NAMES={'direct_SMC_mean':'direct','progressive_SMC_mean':'progressive'}
CONTRASTS={
    'multi_hypothesis':(('direct_point','direct_SMC_mean'),('progressive_point','progressive_SMC_mean')),
    'progressive_schedule':(('direct_SMC_mean','progressive_SMC_mean'),('direct_point','progressive_point'))}


def _safe(value):
    """Preserve proper infinite log scores explicitly; never clip probabilities."""
    if isinstance(value,np.ndarray):return _safe(value.tolist())
    if isinstance(value,np.generic):return _safe(value.item())
    if isinstance(value,dict):return {str(k):_safe(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [_safe(v) for v in value]
    if isinstance(value,float):
        if math.isnan(value):raise ValueError('NaN cannot be serialized as a score')
        if math.isinf(value):return '+Infinity' if value>0 else '-Infinity'
    return value


def _write(path,value):
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(_safe(value),indent=2,allow_nan=False)+'\n')
    temp.replace(path)


def _csv(path,rows):
    fields=list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields)
        writer.writeheader()
        for row in rows:
            converted={key:(json.dumps(_safe(value),allow_nan=False) if isinstance(value,(dict,list,tuple))
                            else _safe(value)) for key,value in row.items()}
            writer.writerow(converted)


def _number(value,name):
    if isinstance(value,(bool,np.bool_)) or not isinstance(value,(int,float,np.number)) or not np.isfinite(value):
        raise ValueError(name+' must be finite numeric')
    return float(value)


def _probabilities(values):
    if np.iscomplexobj(values):raise ValueError('branch masses must be real')
    p=np.asarray(values,dtype=float)
    if p.shape!=(2,) or not np.all(np.isfinite(p)) or np.any(p<0):
        raise ValueError('branch mass requires two finite nonnegative values')
    total=float(p.sum())
    if total<=0 or abs(total-1)>1e-8:
        raise ValueError('branch masses not normalized within floating-point tolerance')
    # A recorded sum normalization handles arithmetic roundoff only. No floor,
    # clipping, missing-branch repair or truth-dependent change is introduced.
    return p/total,total


def proper_branch_scores(masses,true_branch):
    p,total=_probabilities(masses)
    if true_branch not in (0,1):raise ValueError('binary true branch required')
    probability=float(p[true_branch])
    log_score=math.log(probability) if probability>0 else -math.inf
    return dict(branch_Brier=float((p[1]-true_branch)**2),
                proper_branch_log_score=log_score,branch_log_loss=-log_score,
                true_branch_probability=probability,true_branch_probability_zero=probability==0,
                scored_branch_mass=p.tolist(),raw_branch_mass_sum=total,
                roundoff_sum_normalized=total!=1.)


def _cdf_values(mapping,cuts,continuous=False):
    values={}
    for cut in cuts:
        key=str(cut)
        value=mapping[key]['order9'] if continuous else mapping[key]
        values[key]=_number(value,'CDF '+key)
    return values


def score_observation(inference,points,truth_parameters,cuts):
    """Pure scoring function. Caller supplies truth only AFTER freeze gate."""
    if np.iscomplexobj(truth_parameters):raise ValueError('truth parameters must be real')
    truth=np.asarray(truth_parameters,dtype=float)
    if truth.shape!=(5,) or not np.all(np.isfinite(truth)) or truth[0] not in (0,1):
        raise ValueError('truth parameters must be [branch,u,contrast,amplitude,phase]')
    branch=int(truth[0]);target=(-.9 if branch==0 else -.8)+float(truth[1])
    continuous=inference['continuous_ref'];finite=inference['finite_ref']
    cq=continuous.get('qualified') is True;fq=finite.get('qualified') is True
    ccdf=_cdf_values(continuous['CDF'],cuts,True)
    fcdf=_cdf_values(finite['CDF'],cuts)
    rows=[]
    common=dict(observation=int(inference['observation']),parent=int(inference['parent']),
                GT_role='OFFLINE_ONLY',true_branch=branch,true_target_x=target,
                continuous_reference_qualified=cq,finite_reference_qualified=fq,
                reference_qualification_population_retained=True,new_rhs=0)
    for arm in (*ARMS,BAYES):
        row=dict(**common,arm=arm,estimate_x=None,estimate_sd=None,
                 position_squared_error=None,position_absolute_error=None,
                 scoring_status='SCORED',probability_status='NOT_A_PROBABILITY_ARM')
        prediction=None
        if arm in POINT_NAMES:
            prediction=points['methods'][POINT_NAMES[arm]]
            row.update(full_validation_status=prediction.get('full_validation_status'),
                       fallback_used=prediction.get('fallback_used'),
                       response_error_white=prediction.get('response_error_white'),
                       potential_error=prediction.get('potential_error'),
                       Phi=prediction.get('Phi'),
                       point_parameters={k:prediction.get(k) for k in
                                         ('branch','u','contrast','amplitude','phase','x')})
            estimate=prediction.get('x')
        elif arm in SMC_NAMES:
            prediction=inference['ensemble'][SMC_NAMES[arm]]
            estimate=prediction.get('average_x_mean')
            row['estimate_sd']=prediction.get('ensemble_mixture_x_sd')
        else:
            estimate=continuous.get('x_mean')
            row['estimate_sd']=continuous.get('x_sd')
            row['reference_estimate_x_retained']=estimate
            if not cq:
                row['scoring_status']='UNQUALIFIED_CONTINUOUS_REFERENCE_NOT_SCORED'
                rows.append(row)
                continue
        try:
            x=_number(estimate,'position estimate')
            row.update(estimate_x=x,position_squared_error=(x-target)**2,
                       position_absolute_error=abs(x-target))
        except ValueError as exc:
            row.update(scoring_status='INVALID_POSITION_PREDICTION_RETAINED',score_error=str(exc))
        if arm in SMC_NAMES:
            row['probability_status']=('QUALIFIED_REFERENCE_CONTROLLED_FAMILY' if fq
                                       else 'REFERENCE_UNQUALIFIED_EMPIRICAL_SCORES_ONLY')
            try:
                row.update(proper_branch_scores(prediction['average_branch_mass'],branch))
                cdf=_cdf_values(prediction['average_CDF'],cuts)
                ordered=[cdf[str(cut)] for cut in cuts]
                coherent=(all(0<=v<=1 for v in ordered)
                          and all(a<=b+1e-12 for a,b in zip(ordered[:-1],ordered[1:])))
                ec={key:abs(cdf[key]-ccdf[key]) for key in cdf}
                ef={key:abs(cdf[key]-fcdf[key]) for key in cdf}
                row.update(CDF=cdf,CDF_coherent=coherent,
                           CDF_absolute_errors_vs_continuous=ec,CDF_max_error_vs_continuous=max(ec.values()),
                           CDF_absolute_errors_vs_finite=ef,CDF_max_error_vs_finite=max(ef.values()),
                           CDF_checkpoint_truth_Brier=({str(cut):(cdf[str(cut)]-float(target<=cut))**2
                                                        for cut in cuts} if coherent else None),
                           supplied_branch_MISS_run_rate=prediction.get('branch_MISS_run_rate'),
                           supplied_branch_MISS_eligible_rate=prediction.get('branch_MISS_eligible_rate'))
                ref_mass,_=_probabilities(finite['branch_mass'])
                row['missed_significant_branches_ensemble']=[b for b in (0,1)
                    if ref_mass[b]>.05 and row['scored_branch_mass'][b]==0]
                if not coherent:row['probability_status']='INVALID_CDF_RETAINED_NO_PROBABILITY_SUCCESS_CLAIM'
            except (ValueError,KeyError,TypeError) as exc:
                row.update(probability_status='INVALID_PROBABILITY_PREDICTION_RETAINED',
                           probability_score_error=str(exc))
        rows.append(row)
    return rows


def parent_scores(scores,replicas=2):
    """Keep incomplete/unqualified parents visible; full mean needs both noises."""
    parents=sorted({r['parent'] for r in scores});result=[]
    for parent in parents:
        for arm in (*ARMS,BAYES):
            rows=[r for r in scores if r['parent']==parent and r['arm']==arm]
            row=dict(parent=parent,arm=arm,expected_replicas=replicas,observations=len(rows),
                     continuous_qualified_replicas=sum(r['continuous_reference_qualified'] for r in rows),
                     finite_qualified_replicas=sum(r['finite_reference_qualified'] for r in rows))
            for metric in ('position_squared_error','position_absolute_error','branch_Brier',
                           'proper_branch_log_score','branch_log_loss',
                           'CDF_max_error_vs_continuous','CDF_max_error_vs_finite'):
                values=[r.get(metric) for r in rows]
                valid=[v for v in values if v is not None and not math.isnan(float(v))]
                complete=len(rows)==replicas and len(valid)==replicas
                row[metric]=float(np.mean(valid)) if complete else None
                row[metric+'_replicas_scored']=len(valid)
                if valid and not complete:
                    row[metric+'_partial_qualified_mean']=float(np.mean(valid))
            result.append(row)
    return result


def paired_bootstrap_gain(a,b,replicates=2000,seed=2026102317):
    """Pure paired independent-parent percentile bootstrap; positive favors b."""
    a,b=np.asarray(a,dtype=float),np.asarray(b,dtype=float)
    if a.ndim!=1 or a.shape!=b.shape or not len(a) or not np.all(np.isfinite(a-b)):
        raise ValueError('finite complete paired parent loss vectors required')
    if replicates<1:raise ValueError('positive bootstrap count required')
    differences=a-b
    draw=np.random.default_rng(seed).integers(0,len(a),size=(replicates,len(a)))
    means=differences[draw].mean(axis=1)
    lower,upper=np.quantile(means,[.025,.975])
    return dict(parent_count=len(a),mean_loss_a=float(a.mean()),mean_loss_b=float(b.mean()),
                gain_a_minus_b=float(differences.mean()),paired_percentile_CI95=[float(lower),float(upper)],
                replicates=int(replicates),seed=int(seed),unit='independent parent; noise replicas averaged first',
                positive_gain_favors='b',bootstrap_method='paired parent resampling, percentile 2.5/97.5',
                hypothesis_decision='NOT_ASSIGNED')


def summarize(scores,parents,replicates=2000,seed=2026102317):
    ids=sorted({r['parent'] for r in parents})
    lookup={(r['parent'],r['arm']):r for r in parents}
    summary=dict(parent_count=len(ids),score_rows=len(scores),new_rhs=0,GT_role='OFFLINE_ONLY',
                 decision='NO_AUTOMATIC_GO_PARENT_REVIEW_REQUIRED',arm_summaries={},contrasts={})
    for arm in (*ARMS,BAYES):
        rows=[r for r in parents if r['arm']==arm]
        coverage={}
        metrics={}
        for metric in ('position_squared_error','position_absolute_error','branch_Brier',
                       'proper_branch_log_score','branch_log_loss',
                       'CDF_max_error_vs_continuous','CDF_max_error_vs_finite'):
            values=[r[metric] for r in rows if r.get(metric) is not None]
            metrics[metric]=float(np.mean(values)) if values else None
            coverage[metric]=len(values)
        summary['arm_summaries'][arm]=dict(parent_means=metrics,complete_parent_coverage=coverage,
            all_parent_count=len(ids),qualified_subset_only=arm==BAYES,
            partial_parents_retained_in_PARENT_SCORES=True,
            infinite_branch_log_score_parents=sum(r.get('proper_branch_log_score')==-math.inf for r in rows))
    for family,comparisons in CONTRASTS.items():
        summary['contrasts'][family]=[]
        for a,b in comparisons:
            for metric in ('position_squared_error','position_absolute_error'):
                missing=[parent for parent in ids
                         if lookup[(parent,a)].get(metric) is None or lookup[(parent,b)].get(metric) is None]
                row=dict(arm_a=a,arm_b=b,loss=metric,primary=metric=='position_squared_error',
                         all_parent_ids=ids,missing_parent_ids=missing)
                if missing:
                    row.update(status='INCOMPLETE_FULL_POPULATION_NO_CI',
                               eligible_parent_count=len(ids)-len(missing),no_silent_subset_bootstrap=True)
                else:
                    row.update(status='PAIRED_PARENT_DESCRIPTIVE_INTERVAL',**paired_bootstrap_gain(
                        [lookup[(parent,a)][metric] for parent in ids],
                        [lookup[(parent,b)][metric] for parent in ids],replicates,seed))
                summary['contrasts'][family].append(row)
    # One observation contributes five rows; count qualification once only.
    objects={r['observation']:r for r in scores}
    summary['reference_qualification']=dict(
        object_count=len(objects),continuous_qualified=sum(r['continuous_reference_qualified'] for r in objects.values()),
        finite_qualified=sum(r['finite_reference_qualified'] for r in objects.values()),
        continuous_unqualified_observations=[k for k,r in objects.items() if not r['continuous_reference_qualified']],
        finite_unqualified_observations=[k for k,r in objects.items() if not r['finite_reference_qualified']],
        population_objects_not_dropped=True)
    summary['definitions']=dict(branch_Brier='binary event branch1: (p1-true_branch)^2',
        proper_branch_log_score='log p(true_branch), higher better; exact zero yields -Infinity, never floored',
        branch_log_loss='negative proper log score, lower better',
        CDF_errors='against both declared continuous order9 and finite references; qualification flags retained',
        position_mean='Bayes action under squared position loss within controlled family, not a new physical world',
        contrast_families='multi-hypothesis gain and progressive-schedule gain reported separately',
        CI_scope='descriptive paired parent bootstrap; no automatic GO or multiplicity-adjusted declaration')
    return summary


def _index(rows,name,count):
    if not isinstance(rows,list) or len(rows)!=count:raise ValueError(name+' must contain every observation')
    result={}
    for row in rows:
        index=row['observation']
        if isinstance(index,bool) or not isinstance(index,int) or index in result:
            raise ValueError('duplicate/invalid observation in '+name)
        if isinstance(row.get('parent'),bool) or not isinstance(row.get('parent'),int):
            raise ValueError('invalid parent ID in '+name)
        result[index]=row
    if set(result)!=set(range(count)):raise ValueError(name+' observation IDs must cover frozen population')
    return result


def _completion_record(record,name,count):
    if not isinstance(record,dict):raise ValueError(name+' marker must be an object')
    if (record.get('complete') is False or record.get('completed') is False
            or record.get('inference_frozen') is False
            or str(record.get('status','')).upper() in ('RUNNING','INCOMPLETE','FAILED','NOT_RUN','PREPARED')):
        raise RuntimeError(name+' marker does not establish a frozen completed process')
    if 'observations' in record and record['observations']!=count:
        raise RuntimeError(name+' completion count differs from frozen population')


def _cost_fields(row):
    return {key:value for key,value in row.items() if any(word in key.lower()
            for word in ('cost','cpu','wall_seconds','rhs','maxwell'))}


def main():
    """Sole truth-reading scoring entry. Never call before both processes freeze."""
    base=ROOT/'results/P1B'
    inference_dir=base/'inference';validation_dir=base/'point_validation'
    freeze=next((inference_dir/name for name in ('COMPLETE.json','GATE.json')
                 if (inference_dir/name).exists()),None)
    validation_marker=validation_dir/'COMPLETE.json'
    if freeze is None or not validation_marker.exists():
        raise RuntimeError('OFFLINE truth access blocked: inference freeze and point-validation completion required')
    freeze_record=json.loads(freeze.read_text())
    validation_record=json.loads(validation_marker.read_text())
    protocol=json.loads((ROOT/'configs/P1_PROBABILITY_PROTOCOL.json').read_text())
    scoring=protocol['P1B']['scoring']
    count=int(protocol['independent_parents'])*int(protocol['noise_replicas'])
    if (protocol['independent_parents']!=16 or protocol['noise_replicas']!=2
            or scoring['paired_bootstrap_replicates']!=2000 or scoring['bootstrap_seed']!=2026102317):
        raise RuntimeError('frozen 16-parent/two-noise/2000-bootstrap scoring protocol differs')
    _completion_record(freeze_record,'inference freeze',count)
    _completion_record(validation_record,'point validation',count)
    inference=_index(json.loads((inference_dir/'RESULTS.json').read_text()),'inference',count)
    points=_index(json.loads((validation_dir/'RESULTS.json').read_text()),'point validation',count)
    for observation in range(count):
        if inference[observation]['parent']!=points[observation]['parent']:
            raise ValueError('parent mismatch between frozen inference/validation')
    # No truth path is opened above this line. Only two existing archive keys
    # are loaded; the noiseless F array is deliberately not read or regenerated.
    truth_path=base/'data/OFFLINE_TRUTH.npz'
    with np.load(truth_path,allow_pickle=False) as archive:
        parameters=archive['parameters'].copy()
        parent_ids=archive['parent_ids'].copy()
    if (parameters.shape!=(count,5) or parent_ids.shape!=(count,)
            or np.iscomplexobj(parameters) or parent_ids.dtype.kind not in 'iu'
            or not np.all(np.isfinite(parameters))):
        raise ValueError('offline truth shape/finite-value mismatch')
    parent_counts=Counter(parent_ids.tolist())
    if len(parent_counts)!=16 or set(parent_counts.values())!={2}:
        raise ValueError('offline truth must contain 16 independent parents with two noise replicas')
    for parent in parent_counts:
        values=parameters[parent_ids==parent]
        if not np.array_equal(values[0],values[1]):
            raise ValueError('two noise replicas have different parent truth parameters')
    start=time.process_time()
    scores=[];cost_provenance=[]
    for observation in range(count):
        if inference[observation]['parent']!=parent_ids[observation]:
            raise ValueError('offline parent ordering does not match frozen observations')
        scores.extend(score_observation(inference[observation],points[observation],
                                        parameters[observation],protocol['cdf_checkpoints']))
        cost_provenance.append(dict(observation=observation,parent=int(parent_ids[observation]),
            inference_cost_fields=_cost_fields(inference[observation]),
            point_validation_cost_fields=_cost_fields(points[observation]),
            validated_point_methods=points[observation]['methods']))
    parents=parent_scores(scores,2)
    summary=summarize(scores,parents,2000,2026102317)
    summary.update(identity='P1B_FROZEN_OFFLINE_ONLY_SCORE',new_rhs=0,new_truth_generated=False,
                   cpu_seconds=time.process_time()-start,
                   point_full_validation_status_counts=dict(Counter(
                       str(r.get('full_validation_status')) for r in scores if r['arm'] in POINT_NAMES)),
                   actual_population=dict(observations=count,parents=16,noise_replicas=2),
                   all_setup_correction_validation_costs='retained as source provenance; shared costs not double-summed or credited free',
                   outside_controlled_family_mass=None,original_minimum_feature='UNRESOLVED')
    out=base/'offline_scores';out.mkdir(parents=True,exist_ok=False)
    _write(out/'TRUTH_ACCESS.json',dict(role='OFFLINE_ONLY',path=str(truth_path),
        keys_read=['parameters','parent_ids'],inference_freeze_marker=str(freeze),
        point_validation_complete_marker=str(validation_marker),new_rhs=0,
        new_truth_generated=False,online_inference_called=False))
    _write(out/'SCORES.json',scores);_csv(out/'SCORES.csv',scores)
    _write(out/'PARENT_SCORES.json',parents);_csv(out/'PARENT_SCORES.csv',parents)
    _write(out/'COST_PROVENANCE.json',dict(inference_marker=freeze_record,
        point_validation_marker=validation_record,objects=cost_provenance,
        source_bank_setup_costs_remain_paid=True,scorer_new_rhs=0))
    _write(out/'PAIRED_CI.json',summary['contrasts'])
    _write(out/'SUMMARY.json',summary)
    _write(out/'COMPLETE.json',dict(status='OFFLINE_SCORES_READY_FOR_PARENT_REVIEW',
        observations=count,parent_count=16,new_rhs=0,GT_role='OFFLINE_ONLY',automatic_GO=False))
    return summary


if __name__=='__main__':
    argparse.ArgumentParser(description=__doc__).parse_args()
    main()
