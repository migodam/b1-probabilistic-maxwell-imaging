"""Offline uniform 4-micro-wall x 12-target physical finite reference.

Requires the public-macro online seal.  Reads already paid response arrays only;
never calls a forward solver and never changes the six-wall online prior.
"""
from pathlib import Path
import argparse
import csv
import json
import os
import sys
import time

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
sys.path.insert(0,str(ROOT.parent/"research_w_gq/src"))
from gq.contrast import logsumexp
from r2.accounting import BudgetBook

FREQUENCIES=(.85,1.,1.15)
CONTRACTS=("A0","A1","A2")
SIGMA=.01
WALL_IDS=((0,0),(0,1),(1,0),(1,1))
SCOPE="OFFLINE_PHYSICAL_FINITE_REFERENCE_NOT_GENERAL_ONLINE_OR_UNIFORM_CONTINUOUS_PDE"


def plain(value):
    if isinstance(value,Path):return str(value)
    if isinstance(value,np.ndarray):return plain(value.tolist())
    if isinstance(value,np.generic):return plain(value.item())
    if isinstance(value,dict):return {str(k):plain(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [plain(v) for v in value]
    if isinstance(value,float) and not np.isfinite(value):return None
    return value


def write_json(path,value):
    with Path(path).open("x") as stream:
        json.dump(plain(value),stream,indent=2,allow_nan=False);stream.write("\n")


def write_csv(path,rows):
    with Path(path).open("x",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def load(path,keys):
    with np.load(path,allow_pickle=False) as archive:
        return {key:np.array(archive[key],copy=True) for key in keys}


def entropy(probability):
    p=np.asarray(probability);positive=p>0
    return float(-np.sum(p[positive]*np.log(p[positive])))


def nested(value,contract):
    return value[...,1,:2,:] if contract=="A0" else value[...,1,:,:] if contract=="A1" else value


def normalized(logmass,axis=None):
    p=np.exp(logmass-logsumexp(logmass,axis=axis) if axis is None else
             logmass-np.expand_dims(logsumexp(logmass,axis=axis),axis))
    return p/p.sum() if axis is None else p/p.sum(axis=axis,keepdims=True)


def weighted_variance(values,mass):
    mean=mass@values
    return mean,mass@((values-mean)**2)


def correlation(a,b,joint):
    """A is a scalar wall descriptor, B is a scalar target descriptor."""
    wa=joint.sum(axis=1);tb=joint.sum(axis=0)
    am=float(wa@a);bm=float(tb@b)
    av=float(wa@((a-am)**2));bv=float(tb@((b-bm)**2))
    if av<=1e-30 or bv<=1e-30:
        return None
    return float(np.sum(joint*(a-am)[:,None]*(b-bm)[None,:])/np.sqrt(av*bv))


def hpd95(q):
    order=np.lexsort((np.arange(len(q)),-q))
    count=min(len(q),int(np.searchsorted(np.cumsum(q[order]),.95,side="left"))+1)
    selected=np.zeros(len(q),dtype=bool);selected[order[:count]]=True
    return selected


def bootstrap(parent_rows,output,seed=28102931,count=2000):
    """Paired macro-stratified target-parent cluster bootstrap after noise mean.

    The two macros are fixed strata.  Six target-parent scenarios are resampled
    within each; the same 2,000 draws compare every nested contract.  This is
    descriptive, and provides no population inference from two wall parents.
    """
    metrics=("position_MAP_squared_error","position_mean_squared_error","posterior_position_risk",
        "radius_MAP_squared_error","radius_mean_squared_error","contrast_MAP_squared_error",
        "contrast_mean_squared_error","loss_MAP_squared_error","loss_mean_squared_error",
        "shape_marginal_MAP_correct","HPD95_contains_true","HPD95_size",
        "texture_probability_true","texture_entropy_bits","macro_probability_true",
        "conditional_texture_entropy_true_macro_bits",
        "conditional_texture_probability_true_given_true_macro",
        "global_wall_pixel_variance_ratio_to_prior","truth_macro_pixel_variance_ratio_to_prior",
        "joint_wall_target_MI_nats")
    arrays=np.array([[[next(r for r in parent_rows if r["parent"]==p and r["contract"]==c)[m]
                      for m in metrics] for c in CONTRACTS] for p in range(12)])
    macro=np.array([next(r["wall_parent"] for r in parent_rows if r["parent"]==p) for p in range(12)])
    rng=np.random.default_rng(seed)
    draws=np.concatenate([rng.choice(np.flatnonzero(macro==p),size=(count,int(np.sum(macro==p))),replace=True)
                          for p in range(2)],axis=1)
    samples=arrays[draws].mean(axis=1)
    point=arrays.mean(axis=0)
    np.savez_compressed(output/"BOOTSTRAP_2000.npz",parent_indices=draws,
        metrics=np.array(metrics),contracts=np.array(CONTRACTS),bootstrap_parent_means=samples,
        parent_metric_values=arrays,seed=seed)
    rows=[]
    statistic={}
    for k,metric in enumerate(metrics):
        rmse=metric.endswith("_squared_error")
        label=metric.replace("_squared_error","_RMSE") if rmse else metric
        statistic[label]=(np.sqrt(point[:,k]) if rmse else point[:,k],
                          np.sqrt(samples[:,:,k]) if rmse else samples[:,:,k])
    for metric,(p,s) in statistic.items():
        for i,contract in enumerate(CONTRACTS):
            lo,hi=np.quantile(s[:,i],[.025,.975])
            rows.append(dict(contract=contract,statistic=metric,point=p[i],
                descriptive_p025=lo,descriptive_p975=hi,parent_units=12,macro_strata=2,
                noise_repeats_averaged_before_resampling=2,bootstrap_draws=count))
    write_csv(output/"BOOTSTRAP_DESCRIPTIVE.csv",rows)
    paired=[]
    for start,end in ((0,1),(1,2),(0,2)):
        for metric in ("position_MAP_RMSE","position_mean_RMSE","posterior_position_risk",
                       "texture_entropy_bits","conditional_texture_entropy_true_macro_bits",
                       "truth_macro_pixel_variance_ratio_to_prior","global_wall_pixel_variance_ratio_to_prior"):
            p,s=statistic[metric];delta=s[:,end]-s[:,start]
            lo,hi=np.quantile(delta,[.025,.975])
            paired.append(dict(comparison=CONTRACTS[end]+"-"+CONTRACTS[start],statistic=metric,
                point_difference=p[end]-p[start],descriptive_p025=lo,descriptive_p975=hi,
                paired_draws=count,scope="DESCRIPTIVE_FIXED_TWO_MACRO_STRATA_NOT_POPULATION_EVIDENCE"))
    write_csv(output/"PAIRED_ADDED_DATA_DESCRIPTIVE.csv",paired)
    return rows,paired


def micro_consensus_offline(out):
    """Posthoc complete-table distances in the declared physical micro class.

    All 4x12x3 responses and all fixed finite likelihoods are already available.
    No statement here licenses an online stop, a bound from forward proximity,
    a missing-row certificate or a cold-PDE step.
    """
    imaging=ROOT/"results/imaging"
    if not (imaging/"ONLINE_SEALED.json").is_file() or not (out/"SUMMARY.json").is_file():
        raise RuntimeError("sealed online output and completed offline finite reference required")
    b=BudgetBook()
    prior=load(out/"DECLARED_FINITE_PRIOR.npz",("target_positions","target_shape_index",
               "target_radius","target_real_contrast","target_loss"))
    saved=load(out/"RAW_JOINT_POSTERIORS.npz",("joint","log_likelihood","truth_parent",
        "truth_noise","truth_wall_parent","truth_texture","truth_target_id"))
    y=np.load(imaging/"OBSERVATIONS.npy",allow_pickle=False)
    if saved["log_likelihood"].shape!=(24,3,4,12) or y.shape!=(24,3,4,16):
        raise ValueError("completed finite reference schema differs")
    data=np.empty((4,12,3,4,16),dtype=complex)
    for wall,(p,t) in enumerate(WALL_IDS):
        for j in range(12):
            for fi,f in enumerate(FREQUENCIES):
                path=(ROOT/f"results/micro/p{p}_n16_t{t}_f1_q{j+2}.npz" if f==1
                    else imaging/f"micro_p{p}_t{t}_f{f:g}_q{j}.npz")
                data[wall,j,fi]=load(path,("data",))["data"]
    # Fixed legal target-label anchor, same for every background and observation.
    anchor=0
    left,right=np.triu_indices(4,1)
    same_macro=(left//2)==(right//2)
    pair_index={tuple(pair):k for k,pair in enumerate(zip(left,right))}
    positions=prior["target_positions"]
    def realify(value,leading):
        z=np.asarray(value);shape=z.shape[:leading]+(-1,)
        return np.concatenate((z.real.reshape(shape),z.imag.reshape(shape)),axis=-1)/SIGMA
    static=[];f_all=[];delta_all=[];basis_all=[]
    for ci,contract in enumerate(CONTRACTS):
        f=realify(nested(data,contract),2);delta=f-f[:,[anchor]]
        _,singular,vh=np.linalg.svd(delta.reshape(48,-1),full_matrices=False)
        rank=int(np.sum(singular>singular[0]*1e-12)) if singular[0]>0 else 0
        basis=vh[:rank].T
        pair_static=[]
        for k,(a,c) in enumerate(zip(left,right)):
            df=f[a]-f[c];dd=delta[a]-delta[c];df0=f[a,anchor]-f[c,anchor]
            pair_delta=np.vstack((delta[a],delta[c]))
            _,sv,pvh=np.linalg.svd(pair_delta,full_matrices=False)
            pair_rank=int(np.sum(sv>sv[0]*1e-12)) if sv[0]>0 else 0
            pair_basis=pvh[:pair_rank].T
            pair_static.append(dict(contract=contract,pair_index=k,wall_a=int(a),wall_b=int(c),
                macro_a=int(a//2),macro_b=int(c//2),texture_a=int(a%2),texture_b=int(c%2),
                within_same_macro=bool(same_macro[k]),real_observation_dimension=f.shape[-1],
                all_background_increment_span_rank=rank,pair_increment_span_rank=pair_rank,
                full_response_whitened_real_RMS=float(np.sqrt(np.mean(df**2))),
                full_response_targetwise_white_norm_min=float(np.min(np.linalg.norm(df,axis=1))),
                full_response_targetwise_white_norm_max=float(np.max(np.linalg.norm(df,axis=1))),
                full_response_targetwise_white_norm_mean=float(np.mean(np.linalg.norm(df,axis=1))),
                increment_whitened_real_RMS=float(np.sqrt(np.mean(dd**2))),
                increment_targetwise_white_norm_max=float(np.max(np.linalg.norm(dd,axis=1))),
                anchor_background_difference_white_norm=float(np.linalg.norm(df0)),
                anchor_background_difference_projection_global_increment_span_norm=float(np.linalg.norm(df0@basis)),
                anchor_background_difference_projection_pair_increment_span_norm=float(np.linalg.norm(df0@pair_basis))))
        static.append(pair_static);f_all.append(f);delta_all.append(delta);basis_all.append(basis)
    pair_rows=[];case_rows=[];conditional_all=[];contrast_all=[];pair_tv_all=[];pair_d_all=[]
    pair_osc_all=[];group_equal_error_all=[];residual_projection_norm_all=[]
    max_single_identity=0.;max_pair_identity=0.;max_cached_likelihood_error=0.
    for case in range(24):
        parent=int(saved["truth_parent"][case]);noise=int(saved["truth_noise"][case])
        true_macro=int(saved["truth_wall_parent"][case]);true_j=int(saved["truth_target_id"][case])
        for ci,contract in enumerate(CONTRACTS):
            ll=saved["log_likelihood"][case,ci];conditional=normalized(ll,axis=1)
            d=ll-ll[:,[anchor]];f=f_all[ci];delta=delta_all[ci]
            yw=realify(nested(y[case],contract),0)
            ll_check=-.5*np.sum((f-yw)**2,axis=-1)
            max_cached_likelihood_error=max(max_cached_likelihood_error,float(np.max(np.abs(ll-ll_check))))
            r0=yw-f[:,anchor]
            calculated=np.einsum("wm,wxm->wx",r0,delta)-.5*np.sum(delta**2,axis=-1)
            max_single_identity=max(max_single_identity,float(np.max(np.abs(calculated-d))))
            tv=[];dinfs=[];oscs=[];agreements=[];position_agreements=[]
            maps=np.argmax(conditional,axis=1);means=conditional@positions
            for k,(a,c) in enumerate(zip(left,right)):
                diff=d[a]-d[c]
                bar_delta=.5*(delta[a]+delta[c]);delta_delta=delta[a]-delta[c]
                rbar_x=yw-.5*(f[a]+f[c]);df0=f[a,anchor]-f[c,anchor]
                increment_term=np.sum(rbar_x*delta_delta,axis=-1)
                background_projection_term=-bar_delta@df0
                error=float(np.max(np.abs(increment_term+background_projection_term-diff)))
                max_pair_identity=max(max_pair_identity,error)
                distance=float(.5*np.sum(np.abs(conditional[a]-conditional[c])))
                map_agreement=bool(maps[a]==maps[c])
                map_position_distance=float(np.linalg.norm(positions[maps[a]]-positions[maps[c]]))
                row=dict(case=case,parent=parent,noise=noise,true_macro=true_macro,
                    true_texture=int(saved["truth_texture"][case]),true_target_id=true_j,
                    **static[ci][k],conditional_target_TV=distance,
                    contrast_difference_Linf=float(np.max(np.abs(diff))),
                    contrast_difference_oscillation=float(np.ptp(diff)),
                    increment_difference_residual_term_Linf=float(np.max(np.abs(increment_term))),
                    projected_background_residual_term_Linf=float(np.max(np.abs(background_projection_term))),
                    pair_contrast_identity_max_abs_error=error,
                    conditional_target_MAP_a=int(maps[a]),conditional_target_MAP_b=int(maps[c]),
                    conditional_target_MAP_label_agreement=map_agreement,
                    conditional_target_MAP_position_distance=map_position_distance,
                    conditional_target_MAP_position_agreement=map_position_distance<=1e-12,
                    conditional_target_MAP_shape_agreement=bool(prior["target_shape_index"][maps[a]]==prior["target_shape_index"][maps[c]]),
                    conditional_target_MAP_real_contrast_difference=float(abs(prior["target_real_contrast"][maps[a]]-prior["target_real_contrast"][maps[c]])),
                    conditional_target_MAP_loss_difference=float(abs(prior["target_loss"][maps[a]]-prior["target_loss"][maps[c]])),
                    conditional_target_MAP_radius_difference=float(abs(prior["target_radius"][maps[a]]-prior["target_radius"][maps[c]])),
                    conditional_target_mean_position_distance=float(np.linalg.norm(means[a]-means[c])),
                    conditional_target_MAP_a_true_position_squared_error=float(np.sum((positions[maps[a]]-positions[true_j])**2)),
                    conditional_target_MAP_b_true_position_squared_error=float(np.sum((positions[maps[c]]-positions[true_j])**2)),
                    conditional_target_TV_above_002=distance>.02,
                    scope="POSTHOC_ALL_RESPONSES_QUERIED_NO_EARLY_STOP_NO_QUERY_SAVINGS")
                pair_rows.append(row);tv.append(distance);dinfs.append(row["contrast_difference_Linf"])
                oscs.append(row["contrast_difference_oscillation"]);agreements.append(map_agreement)
                position_agreements.append(map_position_distance<=1e-12)
            equal_errors=[]
            for p in range(2):
                equal=.5*(conditional[2*p]+conditional[2*p+1])
                exact=normalized(logsumexp(ll[2*p:2*p+2],axis=0))
                equal_errors.append(float(.5*np.abs(equal-exact).sum()))
            within_indices=np.flatnonzero(same_macro)
            true_pair=pair_index[(2*true_macro,2*true_macro+1)]
            case_rows.append(dict(case=case,parent=parent,noise=noise,contract=contract,
                true_macro=true_macro,all_four_background_conditional_TV_diameter=max(tv),
                both_macro_within_texture_TV_max=max(tv[k] for k in within_indices),
                true_macro_within_texture_TV=tv[true_pair],
                within_macro0_TV=tv[pair_index[(0,1)]],within_macro1_TV=tv[pair_index[(2,3)]],
                all_four_MAP_label_agreement=bool(np.all(maps==maps[0])),
                all_four_MAP_position_agreement=bool(np.all(np.linalg.norm(positions[maps]-positions[maps[0]],axis=1)<=1e-12)),
                both_macro_within_MAP_label_agreement=all(agreements[k] for k in within_indices),
                true_macro_MAP_label_agreement=agreements[true_pair],
                true_macro_equal_texture_mixture_TV_to_exact_macro_conditional=equal_errors[true_macro],
                macro0_equal_texture_mixture_TV_to_exact_macro_conditional=equal_errors[0],
                macro1_equal_texture_mixture_TV_to_exact_macro_conditional=equal_errors[1]))
            conditional_all.append(conditional);contrast_all.append(d);pair_tv_all.append(tv)
            pair_d_all.append(dinfs);pair_osc_all.append(oscs);group_equal_error_all.append(equal_errors)
            residual_projection_norm_all.append(np.linalg.norm(r0@basis_all[ci],axis=1))
    if max_cached_likelihood_error>1e-9 or max_single_identity>1e-9 or max_pair_identity>1e-9:
        raise AssertionError("saved physical micro finite contrast identity failed")
    write_csv(out/"MICRO_CONSENSUS_OFFLINE.csv",pair_rows)
    write_csv(out/"MICRO_CONSENSUS_CASES.csv",case_rows)
    parent_rows=[]
    for parent in range(12):
        for contract in CONTRACTS:
            rr=[r for r in case_rows if r["parent"]==parent and r["contract"]==contract]
            parent_rows.append(dict(parent=parent,contract=contract,noise_repeats=2,true_macro=rr[0]["true_macro"],
                all_four_conditional_diameter_noise_mean=float(np.mean([r["all_four_background_conditional_TV_diameter"] for r in rr])),
                all_four_conditional_diameter_noise_max=max(r["all_four_background_conditional_TV_diameter"] for r in rr),
                within_both_macro_TV_noise_max=max(r["both_macro_within_texture_TV_max"] for r in rr),
                true_macro_within_texture_TV_noise_mean=float(np.mean([r["true_macro_within_texture_TV"] for r in rr])),
                true_macro_within_texture_TV_noise_max=max(r["true_macro_within_texture_TV"] for r in rr),
                true_macro_MAP_label_agreement_noise_fraction=float(np.mean([r["true_macro_MAP_label_agreement"] for r in rr])),
                true_macro_equal_texture_mixture_target_TV_noise_mean=float(np.mean([
                    r["true_macro_equal_texture_mixture_TV_to_exact_macro_conditional"] for r in rr]))))
    write_csv(out/"MICRO_CONSENSUS_PARENT.csv",parent_rows)
    np.savez_compressed(out/"MICRO_CONSENSUS_RAW.npz",conditional_target=np.array(conditional_all).reshape(24,3,4,12),
        contrast=np.array(contrast_all).reshape(24,3,4,12),
        pair_target_TV=np.array(pair_tv_all).reshape(24,3,6),
        pair_contrast_Linf=np.array(pair_d_all).reshape(24,3,6),
        pair_contrast_oscillation=np.array(pair_osc_all).reshape(24,3,6),
        macro_equal_texture_mixture_target_TV=np.array(group_equal_error_all).reshape(24,3,2),
        background_residual_projection_norm=np.array(residual_projection_norm_all).reshape(24,3,4),
        pair_left=left,pair_right=right,within_same_macro=same_macro,contracts=np.array(CONTRACTS),anchor_target_label=anchor)
    def stats(values):
        a=np.asarray(values)
        return dict(minimum=float(np.min(a)),median=float(np.median(a)),p90=float(np.quantile(a,.9)),
                    maximum=float(np.max(a)),mean=float(np.mean(a)))
    summaries=[]
    for contract in CONTRACTS:
        rr=[r for r in pair_rows if r["contract"]==contract and r["within_same_macro"]]
        cc=[r for r in case_rows if r["contract"]==contract]
        summary=dict(contract=contract,within_macro_pair_observations=len(rr),
            within_macro_conditional_TV=stats([r["conditional_target_TV"] for r in rr]),
            true_macro_conditional_TV=stats([r["true_macro_within_texture_TV"] for r in cc]),
            all_four_background_conditional_diameter=stats([r["all_four_background_conditional_TV_diameter"] for r in cc]),
            within_macro_TV_above_002_count=sum(r["conditional_target_TV_above_002"] for r in rr),
            true_macro_TV_above_002_count=sum(r["true_macro_within_texture_TV"]>.02 for r in cc),
            within_macro_MAP_label_agreement_fraction=float(np.mean([r["conditional_target_MAP_label_agreement"] for r in rr])),
            within_macro_MAP_position_agreement_fraction=float(np.mean([r["conditional_target_MAP_position_agreement"] for r in rr])),
            all_four_MAP_label_agreement_fraction=float(np.mean([r["all_four_MAP_label_agreement"] for r in cc])),
            true_macro_equal_texture_mixture_TV_to_exact_macro_conditional=stats([
                r["true_macro_equal_texture_mixture_TV_to_exact_macro_conditional"] for r in cc]),
            within_macro_full_response_targetwise_white_norm_max=max(r["full_response_targetwise_white_norm_max"] for r in rr),
            within_macro_increment_targetwise_white_norm_max=max(r["increment_targetwise_white_norm_max"] for r in rr),
            within_macro_contrast_difference_Linf=stats([r["contrast_difference_Linf"] for r in rr]),
            within_macro_projected_background_residual_term_Linf=stats([r["projected_background_residual_term_Linf"] for r in rr]))
        summaries.append(summary)
    failed=[{k:r[k] for k in ("case","parent","noise","contract","macro_a","conditional_target_TV",
        "conditional_target_MAP_label_agreement")} for r in pair_rows if r["within_same_macro"] and r["conditional_target_TV_above_002"]]
    result=dict(scope="POSTHOC_FULL_QUERY_PHYSICAL_MICRO_FINITE_CLASS_ONLY",
        status="DESCRIPTIVE_NUMBERS_PARENT_SCIENTIFIC_JUDGMENT_REQUIRED",summary=summaries,
        all_failures_above_002=failed,case_count=24,target_parent_scenarios=12,macro_parents=2,
        texture_alternatives_per_macro=2,proper_prior="same uniform 4 micro-wall x12 target generation family, offline only",
        conditional_target_prior="common uniform twelve joint target labels given each fixed wall",
        target_anchor_label=anchor,full_response_arrays_read=144,
        cached_likelihood_comparison_max_abs_error=max_cached_likelihood_error,
        residual_projection_contrast_identity_max_abs_error=max_single_identity,
        two_background_residual_projection_contrast_identity_max_abs_error=max_pair_identity,
        identity="D=<y-f0,Delta>-.5||Delta||^2; D1-D2=<rbar_X,deltaDelta>-<deltaf0,barDelta>",
        distance_units="realified /sigma(.01); response and increment RMS over target labels and real observation coordinates; norms over observation coordinates",
        projections="empirical complete finite increment spans, relative SVD cutoff1e-12, all teacher responses used",
        all_static_pair_metrics=static,actual_new_full_FE=0,actual_new_RHS=0,
        online_early_stop=False,online_query_savings=0,TV_certificate_from_forward_distance=False,
        unqueried_background_bound_available=False,cold_PDE_proposed=False,
        no_uniform_continuous_PDE_claim=True,online_prior_changed=False,
        true_fixed_wall_parameters_used_online=False,
        limitation="all fixed walls/targets/frequencies queried before diagnostic; same-generation controlled finite reference only; empirical floating arithmetic, not uniform physical certificate")
    write_json(out/"MICRO_CONSENSUS_OFFLINE.json",result)
    b.event("E3","physical_micro_consensus_offline",actual_new_full_FE=0,actual_new_RHS=0,
        offline_paid_response_loads=144,query_savings_claim=False,TV_certificate=False,
        finite_class_walls=4,targets=12,observations=24)
    print(json.dumps(plain(dict(scope=result["scope"],summary=summaries,
        failure_rows=len(failed),actual_new_full_FE=0,actual_new_RHS=0))),flush=True)
    return result


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/"results/e2_micro_reference")
    parser.add_argument("--consensus-only",action="store_true")
    args=parser.parse_args(argv)
    out=args.output.resolve();imaging=ROOT/"results/imaging"
    seal_path=imaging/"ONLINE_SEALED.json"
    if not seal_path.is_file():
        raise RuntimeError("WAIT_FOR_PARENT: public-macro ONLINE_SEALED.json is required")
    seal=json.loads(seal_path.read_text())
    if seal.get("no_truth_input") is not True or not (imaging/"ONLINE_OUTPUTS.npz").is_file():
        raise ValueError("public-macro online output seal/schema unavailable")
    if args.consensus_only:
        if not out.is_relative_to(ROOT):raise ValueError("offline output must stay under R2")
        return micro_consensus_offline(out)
    if not out.is_relative_to(ROOT) or out.exists():
        raise ValueError("fresh offline reference directory under R2 required")
    b=BudgetBook();started=time.process_time();out.mkdir(parents=True)
    write_json(out/"PLAN_AND_BARRIER.json",dict(scope=SCOPE,online_seal=seal,
        seal_required_before_truth_or_micro_reference_access=True,
        proper_joint_prior="uniform 4 micro walls x 12 target labels, 48 atoms each mass 1/48",
        micro_wall_ids=WALL_IDS,wall_resolution="eta=1/16 lambda0, paid discretization",
        online_prior_unchanged="six public macro walls x12 targets uniform counting",
        online_claim=False,PDE_speed_claim=False,uniform_continuous_PDE_claim=False,
        actual_new_full_FE=0,actual_new_RHS=0,bootstrap_draws=2000,
        bootstrap_units="12 target-parent scenarios, two fixed macro strata, noise-average 2 first"))
    # Public target definitions are not inferred or fitted to these observations.
    from r2.worlds import target_defs
    targets=target_defs()[1]
    if len(targets)!=12:
        raise ValueError("predeclared twelve target labels are required")
    positions=np.array([t["pose"] for t in targets],dtype=float)
    radius=np.array([t["radius"] for t in targets]);contrast=np.array([t["contrast"] for t in targets])
    loss=np.array([t["loss"] for t in targets]);shapes=tuple(sorted({t["shape"] for t in targets}))
    shape_index=np.array([shapes.index(t["shape"]) for t in targets])
    if len(shapes)!=2:
        raise ValueError("two predeclared target shapes required")
    observations=np.load(imaging/"OBSERVATIONS.npy",allow_pickle=False)
    truth=json.loads((imaging/"TRUTH_OFFLINE.json").read_text())
    if observations.shape!=(24,3,4,16) or len(truth)!=24 or np.any(~np.isfinite(observations)):
        raise ValueError("schema unavailable: 24 observations x3f x4source x16RX")
    for parent in range(12):
        records=[t for t in truth if t["parent"]==parent]
        if len(records)!=2 or {t["noise"] for t in records}!={0,1} or any(t["target_id"]!=parent for t in records):
            raise ValueError("twelve target-parent scenarios with two noise repeats required")
    if {t["wall_parent"] for t in truth}!={0,1}:
        raise ValueError("two declared macro parents required")
    if any(sum(t["wall_parent"]==p and t["noise"]==0 for t in truth)!=6 for p in range(2)):
        raise ValueError("six target-parent scenarios per fixed macro stratum required")
    data=np.empty((4,12,3,4,16),dtype=complex);source=[]
    for wall,(parent,texture) in enumerate(WALL_IDS):
        for target in range(12):
            for fi,f in enumerate(FREQUENCIES):
                path=(ROOT/f"results/micro/p{parent}_n16_t{texture}_f1_q{target+2}.npz" if f==1
                      else imaging/f"micro_p{parent}_t{texture}_f{f:g}_q{target}.npz")
                value=load(path,("data",))["data"]
                if value.shape!=(4,16) or np.any(~np.isfinite(value)):
                    raise ValueError("paid response schema differs: "+str(path))
                data[wall,target,fi]=value
                source.append(dict(path=path.relative_to(ROOT),wall_id=wall,wall_parent=parent,
                                   texture=texture,target_id=target,frequency=f))
    wall_fields=[load(ROOT/f"results/micro/WALL_p{p}_n16_t{t}.npz",("x","y","eps","ell")) for p,t in WALL_IDS]
    if any(not np.array_equal(wall_fields[0][key],w[key]) for key in ("x","y") for w in wall_fields):
        raise ValueError("four paid walls must share a physical pixel grid")
    ny,nx=wall_fields[0]["eps"].shape
    eps=np.array([w["eps"].ravel() for w in wall_fields]);ell=np.array([w["ell"].ravel() for w in wall_fields])
    occupied=np.any(eps>1,axis=0)
    wall_prior=np.ones(4)/4;target_prior=np.ones(12)/12
    prior_eps_mean,prior_eps_var=weighted_variance(eps,wall_prior)
    prior_ell_mean,prior_ell_var=weighted_variance(ell,wall_prior)
    prior_cond_eps=[];prior_cond_ell=[]
    for p in range(2):
        prior_cond_eps.append(weighted_variance(eps[2*p:2*p+2],np.ones(2)/2)[1])
        prior_cond_ell.append(weighted_variance(ell[2*p:2*p+2],np.ones(2)/2)[1])
    prior_cond_eps=np.array(prior_cond_eps);prior_cond_ell=np.array(prior_cond_ell)
    prior_position_mean=target_prior@positions
    prior_position_risk=float(np.sum(target_prior*np.sum((positions-prior_position_mean)**2,axis=1)))
    np.savez_compressed(out/"DECLARED_FINITE_PRIOR.npz",joint_mass=np.ones((4,12))/48,
        wall_parent=np.array([p for p,t in WALL_IDS]),texture=np.array([t for p,t in WALL_IDS]),
        target_positions=positions,target_radius=radius,target_real_contrast=contrast,target_loss=loss,
        target_shape_index=shape_index,target_shapes=np.array(shapes),
        wall_x=wall_fields[0]["x"],wall_y=wall_fields[0]["y"],wall_epsilon=eps.reshape(4,ny,nx),
        wall_loss=ell.reshape(4,ny,nx),prior_pixel_epsilon_variance=prior_eps_var.reshape(ny,nx),
        prior_pixel_loss_variance=prior_ell_var.reshape(ny,nx),
        prior_conditional_macro_pixel_epsilon_variance=prior_cond_eps.reshape(2,ny,nx),
        prior_conditional_macro_pixel_loss_variance=prior_cond_ell.reshape(2,ny,nx),
        prior_position_risk=prior_position_risk)
    write_json(out/"SOURCE_ACCESS.json",dict(offline_after_seal=True,read_paid_response_atoms=144,
        paid_responses=source,truth_keys=list(truth[0]),frequency_order=FREQUENCIES,
        real_imag_noise_sigma=SIGMA,full_reference_used_online=False,
        source_files_modified=False,actual_new_full_FE=0,actual_new_RHS=0,
        physical_scope="fixed stored finite full-wave response predictions; discretization/model/roundoff error not certified"))
    joint_all=[];ll_all=[];target_all=[];wall_all=[];texture_cond_all=[];hpd_all=[]
    pixel_mean_all=[];pixel_var_all=[];pixel_loss_var_all=[];pixel_cond_var_all=[];pixel_cond_loss_var_all=[]
    pixel_within_all=[];pixel_between_all=[];pixel_corr_all=[];pixel_corr_valid_all=[]
    row_records=[];identity=dict(max_probability_sum_error=0.,max_marginal_identity_error=0.,
        max_total_pixel_variance_decomposition_error=0.,max_true_response_observation_noise_error=0.)
    for case,(y,tt) in enumerate(zip(observations,truth)):
        true_wall=2*tt["wall_parent"]+tt["texture"];true_j=tt["target_id"]
        if targets[true_j]!=tt["target"]:
            raise ValueError("truth target metadata differs from fixed target prior")
        for contract in CONTRACTS:
            pred=nested(data,contract);yy=nested(y,contract)
            axes=tuple(range(2,pred.ndim))
            ll=-.5*np.sum(np.abs(pred-yy)**2,axis=axes)/SIGMA**2
            joint=normalized(ll-np.log(48.))
            q=joint.sum(axis=0);wall=joint.sum(axis=1)
            log_wall=logsumexp(ll,axis=1)-np.log(48.)
            conditional_tex=np.vstack([normalized(log_wall[2*p:2*p+2]) for p in range(2)])
            macro=wall.reshape(2,2).sum(axis=1);texture=wall.reshape(2,2).sum(axis=0)
            hpd=hpd95(q)
            map_id=int(np.argmax(q));joint_wall,joint_target=np.unravel_index(np.argmax(joint),joint.shape)
            pos_mean=q@positions;posvar=np.sum(q[:,None]*(positions-pos_mean)**2,axis=0)
            risk=float(posvar.sum());map_pos=positions[map_id]
            posterior_radius=q@radius;posterior_contrast=q@contrast;posterior_loss=q@loss
            shape_q=np.array([q[shape_index==i].sum() for i in range(2)])
            position_distances=np.linalg.norm(positions-pos_mean,axis=1)
            distance_order=np.argsort(position_distances,kind="stable")
            cred_index=min(11,int(np.searchsorted(np.cumsum(q[distance_order]),.95,side="left")))
            credible_radius=float(position_distances[distance_order[cred_index]])
            eps_mean,eps_var=weighted_variance(eps,wall);ell_mean,ell_var=weighted_variance(ell,wall)
            cond_mean=[];cond_var=[];cond_ell_var=[];cond_entropy=[]
            for p in range(2):
                cm,cv=weighted_variance(eps[2*p:2*p+2],conditional_tex[p])
                cond_mean.append(cm);cond_var.append(cv)
                cond_ell_var.append(weighted_variance(ell[2*p:2*p+2],conditional_tex[p])[1])
                cond_entropy.append(entropy(conditional_tex[p])/np.log(2))
            cond_mean=np.array(cond_mean);cond_var=np.array(cond_var)
            within=macro@cond_var;between=macro@((cond_mean-eps_mean)**2)
            identity["max_total_pixel_variance_decomposition_error"]=max(identity["max_total_pixel_variance_decomposition_error"],
                float(np.max(np.abs(eps_var-within-between))))
            cov=(eps-eps_mean).T@(joint@(positions-pos_mean))
            denominator=np.sqrt(eps_var[:,None]*posvar[None,:]);valid=denominator>1e-15
            pixel_corr=np.full_like(cov,np.nan);np.divide(cov,denominator,out=pixel_corr,where=valid)
            joint_log=ll-logsumexp(ll)
            wall_log=logsumexp(joint_log,axis=1);target_log=logsumexp(joint_log,axis=0)
            mutual=float(np.sum(joint*(joint_log-wall_log[:,None]-target_log[None,:])))
            scalar=dict(case=case,parent=tt["parent"],noise=tt["noise"],wall_parent=tt["wall_parent"],
                true_texture=tt["texture"],true_target_id=true_j,contract=contract,
                target_MAP_id=map_id,joint_MAP_wall_id=int(joint_wall),joint_MAP_target_id=int(joint_target),
                target_MAP_x=float(map_pos[0]),target_MAP_y=float(map_pos[1]),
                posterior_mean_x=float(pos_mean[0]),posterior_mean_y=float(pos_mean[1]),
                position_MAP_squared_error=float(np.sum((map_pos-positions[true_j])**2)),
                position_mean_squared_error=float(np.sum((pos_mean-positions[true_j])**2)),
                posterior_position_risk=risk,position_credible_radius95=credible_radius,
                target_MAP_radius=float(radius[map_id]),posterior_mean_radius=float(posterior_radius),
                radius_MAP_squared_error=float((radius[map_id]-radius[true_j])**2),
                radius_mean_squared_error=float((posterior_radius-radius[true_j])**2),
                posterior_radius_risk=float(q@((radius-posterior_radius)**2)),
                target_MAP_real_contrast=float(contrast[map_id]),posterior_mean_real_contrast=float(posterior_contrast),
                contrast_MAP_squared_error=float((contrast[map_id]-contrast[true_j])**2),
                contrast_mean_squared_error=float((posterior_contrast-contrast[true_j])**2),
                posterior_contrast_risk=float(q@((contrast-posterior_contrast)**2)),
                target_MAP_loss=float(loss[map_id]),posterior_mean_loss=float(posterior_loss),
                loss_MAP_squared_error=float((loss[map_id]-loss[true_j])**2),
                loss_mean_squared_error=float((posterior_loss-loss[true_j])**2),
                posterior_loss_risk=float(q@((loss-posterior_loss)**2)),
                target_MAP_shape=targets[map_id]["shape"],
                target_label_MAP_shape_correct=int(shape_index[map_id]==shape_index[true_j]),
                shape_marginal_MAP_correct=int(int(np.argmax(shape_q))==shape_index[true_j]),
                shape_probability_true=float(shape_q[shape_index[true_j]]),
                shape_entropy_bits=entropy(shape_q)/np.log(2),
                target_probability_true=float(q[true_j]),target_entropy_nats=entropy(q),
                HPD95_contains_true=int(hpd[true_j]),HPD95_size=int(hpd.sum()),HPD95_mass=float(q[hpd].sum()),
                wall_probability_true=float(wall[true_wall]),macro_probability_true=float(macro[tt["wall_parent"]]),
                texture_probability_true=float(texture[tt["texture"]]),texture_entropy_bits=entropy(texture)/np.log(2),
                macro_entropy_bits=entropy(macro)/np.log(2),wall_entropy_bits=entropy(wall)/np.log(2),
                conditional_texture_entropy_macro0_bits=float(cond_entropy[0]),
                conditional_texture_entropy_macro1_bits=float(cond_entropy[1]),
                conditional_texture_entropy_true_macro_bits=float(cond_entropy[tt["wall_parent"]]),
                conditional_texture_entropy_posterior_macro_average_bits=float(macro@cond_entropy),
                conditional_texture_probability_true_given_true_macro=float(conditional_tex[tt["wall_parent"],tt["texture"]]),
                global_wall_pixel_epsilon_variance_public=float(eps_var.mean()),
                global_wall_pixel_epsilon_variance_occupied=float(eps_var[occupied].mean()),
                global_wall_pixel_variance_ratio_to_prior=float(eps_var.mean()/prior_eps_var.mean()),
                global_wall_pixel_loss_variance_public=float(ell_var.mean()),
                truth_macro_pixel_epsilon_variance_public=float(cond_var[tt["wall_parent"]].mean()),
                truth_macro_pixel_loss_variance_public=float(np.asarray(cond_ell_var)[tt["wall_parent"]].mean()),
                truth_macro_pixel_variance_ratio_to_prior=float(cond_var[tt["wall_parent"]].mean()/prior_cond_eps[tt["wall_parent"]].mean()),
                pixel_within_macro_variance_public=float(within.mean()),pixel_between_macro_variance_public=float(between.mean()),
                joint_wall_target_MI_nats=max(0.,mutual),
                wall_macro_position_corr_x=correlation(np.array([0.,0.,1.,1.]),positions[:,0],joint),
                wall_macro_position_corr_y=correlation(np.array([0.,0.,1.,1.]),positions[:,1],joint),
                texture_position_corr_x=correlation(np.array([0.,1.,0.,1.]),positions[:,0],joint),
                texture_position_corr_y=correlation(np.array([0.,1.,0.,1.]),positions[:,1],joint),
                wall_epsilon_mean_position_corr_x=correlation(eps.mean(axis=1),positions[:,0],joint),
                wall_epsilon_mean_position_corr_y=correlation(eps.mean(axis=1),positions[:,1],joint),
                actual_new_full_FE=0,actual_new_RHS=0)
            for w,(p,t) in enumerate(WALL_IDS):scalar[f"wall_p{p}_t{t}_probability"]=float(wall[w])
            for p in range(2):
                for t in range(2):scalar[f"texture{t}_probability_given_macro{p}"]=float(conditional_tex[p,t])
            row_records.append(scalar)
            joint_all.append(joint);ll_all.append(ll);target_all.append(q);wall_all.append(wall)
            texture_cond_all.append(conditional_tex);hpd_all.append(hpd)
            pixel_mean_all.append(eps_mean.reshape(ny,nx));pixel_var_all.append(eps_var.reshape(ny,nx))
            pixel_loss_var_all.append(ell_var.reshape(ny,nx));pixel_cond_var_all.append(cond_var.reshape(2,ny,nx))
            pixel_cond_loss_var_all.append(np.asarray(cond_ell_var).reshape(2,ny,nx))
            pixel_within_all.append(within.reshape(ny,nx));pixel_between_all.append(between.reshape(ny,nx))
            pixel_corr_all.append(pixel_corr.reshape(ny,nx,2));pixel_corr_valid_all.append(valid.reshape(ny,nx,2))
            identity["max_probability_sum_error"]=max(identity["max_probability_sum_error"],abs(float(joint.sum())-1))
            alternate=normalized(logsumexp(ll,axis=0)-np.log(48.))
            identity["max_marginal_identity_error"]=max(identity["max_marginal_identity_error"],float(np.max(np.abs(q-alternate))))
    reshape=lambda values,*tail:np.asarray(values).reshape(24,3,*tail)
    np.savez_compressed(out/"RAW_JOINT_POSTERIORS.npz",joint=reshape(joint_all,4,12),
        log_likelihood=reshape(ll_all,4,12),target=reshape(target_all,12),wall=reshape(wall_all,4),
        texture_given_macro=reshape(texture_cond_all,2,2),HPD95=reshape(hpd_all,12),
        contracts=np.array(CONTRACTS),wall_ids=np.array(WALL_IDS),
        truth_parent=np.array([t["parent"] for t in truth]),truth_noise=np.array([t["noise"] for t in truth]),
        truth_wall_parent=np.array([t["wall_parent"] for t in truth]),truth_texture=np.array([t["texture"] for t in truth]),
        truth_target_id=np.array([t["target_id"] for t in truth]))
    np.savez_compressed(out/"POSTERIOR_WALL_PIXEL_DIAGNOSTICS.npz",
        posterior_pixel_epsilon_mean=reshape(pixel_mean_all,ny,nx),
        posterior_pixel_epsilon_variance=reshape(pixel_var_all,ny,nx),
        posterior_pixel_loss_variance=reshape(pixel_loss_var_all,ny,nx),
        posterior_pixel_epsilon_variance_given_macro=reshape(pixel_cond_var_all,2,ny,nx),
        posterior_pixel_loss_variance_given_macro=reshape(pixel_cond_loss_var_all,2,ny,nx),
        posterior_pixel_within_macro_variance=reshape(pixel_within_all,ny,nx),
        posterior_pixel_between_macro_variance=reshape(pixel_between_all,ny,nx),
        posterior_pixel_epsilon_position_correlation=reshape(pixel_corr_all,ny,nx,2),
        posterior_pixel_correlation_defined=reshape(pixel_corr_valid_all,ny,nx,2),
        prior_pixel_epsilon_variance=prior_eps_var.reshape(ny,nx),
        prior_pixel_epsilon_variance_given_macro=prior_cond_eps.reshape(2,ny,nx),
        x=wall_fields[0]["x"],y=wall_fields[0]["y"])
    write_csv(out/"CASE_RESULTS.csv",row_records)
    identity.update(status="PASS_FIXED_SAVED_FLOATING_REFERENCE_IDENTITIES_NOT_PDE_CERTIFICATION",
        no_new_PDE=True,online_prior_changed=False,reference_wall_count=4,reference_target_count=12,
        reference_prior_sum=1.,reference_atom_count=48,paid_response_arrays_read=144,
        probability_nonnegative=bool(np.all(np.asarray(joint_all)>=0)),
        HPD_mass_min=min(r["HPD95_mass"] for r in row_records))
    if (identity["max_probability_sum_error"]>1e-12 or identity["max_marginal_identity_error"]>1e-12 or
            identity["max_total_pixel_variance_decomposition_error"]>1e-10 or identity["HPD_mass_min"]<.95-1e-12):
        raise AssertionError("finite reference identity check failed")
    identity.pop("max_true_response_observation_noise_error")
    write_json(out/"REFERENCE_IDENTITY.json",identity)
    average_metrics=("position_MAP_squared_error","position_mean_squared_error","posterior_position_risk",
        "radius_MAP_squared_error","radius_mean_squared_error","contrast_MAP_squared_error",
        "contrast_mean_squared_error","loss_MAP_squared_error","loss_mean_squared_error",
        "target_label_MAP_shape_correct","shape_marginal_MAP_correct","shape_probability_true",
        "HPD95_contains_true","HPD95_size","position_credible_radius95","texture_probability_true",
        "texture_entropy_bits","macro_probability_true","conditional_texture_entropy_true_macro_bits",
        "conditional_texture_probability_true_given_true_macro","global_wall_pixel_variance_ratio_to_prior",
        "truth_macro_pixel_variance_ratio_to_prior","joint_wall_target_MI_nats")
    parent_rows=[]
    for parent in range(12):
        for contract in CONTRACTS:
            selected=[r for r in row_records if r["parent"]==parent and r["contract"]==contract]
            parent_rows.append(dict(parent=parent,wall_parent=selected[0]["wall_parent"],
                true_texture=selected[0]["true_texture"],contract=contract,noise_repeats=2,
                **{k:float(np.mean([r[k] for r in selected])) for k in average_metrics}))
    write_csv(out/"PARENT_RESULTS_NOISE_AVERAGED.csv",parent_rows)
    bootstrap_rows,paired=bootstrap(parent_rows,out)
    summary=[]
    for contract in CONTRACTS:
        selected=[r for r in parent_rows if r["contract"]==contract]
        record=dict(contract=contract,independent_target_scenarios=12,wall_macro_parents=2,
            stored_noise_observations=24,prior_position_risk=prior_position_risk,
            **{k:float(np.mean([r[k] for r in selected])) for k in average_metrics})
        for k in average_metrics:
            if k.endswith("_squared_error"):
                record[k.replace("_squared_error","_RMSE")]=float(np.sqrt(record[k]))
        record["posterior_position_risk_ratio_to_prior"]=record["posterior_position_risk"]/prior_position_risk
        summary.append(record)
    comparison=[]
    macro_results=imaging/"RESULTS.json"
    if macro_results.is_file():
        for row in json.loads(macro_results.read_text())["summary"]:
            if row["method"] in ("FINITE_MACRO_TARGET_MARGINAL","KNOWN_WALL_ORACLE"):
                comparison.append(dict(**row,reference_scope="MODEL_MISMATCH_DIAGNOSTIC" if row["method"]=="FINITE_MACRO_TARGET_MARGINAL" else "TRUTH_WALL_KNOWN_ORACLE_DIAGNOSTIC"))
    result=dict(status="OFFLINE_FINITE_REFERENCE_SCORED_PARENT_JUDGMENT_REQUIRED",scope=SCOPE,
        summary=summary,paired_added_data=paired,public_macro_comparison=comparison,
        independent_target_scenarios=12,wall_macro_parents=2,micro_texture_alternatives_per_macro=2,
        reference_wall_atoms=4,reference_target_atoms=12,proper_uniform_joint_atoms=48,
        target_label_prior="twelve fixed joint position/radius/real-contrast/loss/shape labels, not an independent-factorial continuous target prior",
        conditional_texture_loss_prior_variance_zero=bool(np.max(prior_cond_ell)<1e-24),
        noise_repeats_per_scenario=2,observation_count=24,
        target_and_wall_truth_assignment="12 target scenarios, paired fixed macro/texture assignment; not a crossed population sample",
        bootstrap=dict(draws=2000,seed=28102931,unit="target-parent scenario after averaging two noise repetitions",
            strata="two fixed macro parents, six target scenarios each",paired_acquisitions=True,
            scope="descriptive conditional on this finite design; no population coverage guarantee"),
        online_prior_unchanged=True,online_results_modified=False,actual_new_full_FE=0,actual_new_RHS=0,
        paid_response_arrays_read=144,finite_atom_log_likelihood_evaluations=24*3*48,
        numerical_identity=identity,process_CPU_seconds=time.process_time()-started,
        physical_scope="fixed paid eta=1/16 micro-wall finite response family; no uniform continuous-PDE, model-error or roundoff certificate",
        scientific_claim_judgment="RESERVED_FOR_PARENT")
    write_json(out/"SUMMARY.json",result)
    b.event("E2","offline_micro_finite_reference",actual_new_full_FE=0,actual_new_RHS=0,
        offline_paid_response_loads=144,finite_atom_log_likelihood_evaluations=24*3*48,
        online_prior_unchanged=True,scope=SCOPE)
    report=["# E2 四微观墙 × 十二目标的离线有限参考", "",
        "此分析在六宏观墙在线结果封存后进行。它使用已付全波响应、均匀 4×12 joint counting prior，给出该固定响应族的有限后验数值。在线 prior 仍为六宏观墙；此参考不追溯加入在线模型，不作在线加速或连续 PDE 一致性主张。科学证据裁决留给主线程。", "",
        "设计包含 12 个目标场景、2 个宏观墙 parent、每个 macro 的 2 种纹理候选；24 条观测是每个目标场景 2 次噪声重复，不是 24 个独立 parent。目标与墙/纹理的真值分配固定，没有全部交叉抽样。A0=f1 的前两源；A1=f1 的全部四源；A2=.85/1/1.15 的全部四源。", "",
        "| 合同 | MAP位置RMSE | 后验均值位置RMSE | 后验位置风险 | 目标半径均值RMSE | 实contrast均值RMSE | loss均值RMSE | shape正确率 | HPD95覆盖 |", "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in summary:
        report.append("| %s | %.6f | %.6f | %.6g | %.6f | %.6f | %.6f | %.3f | %.3f |" % (r["contract"],r["position_MAP_RMSE"],r["position_mean_RMSE"],r["posterior_position_risk"],r["radius_mean_RMSE"],r["contrast_mean_RMSE"],r["loss_mean_RMSE"],r["shape_marginal_MAP_correct"],r["HPD95_contains_true"]))
    report.extend(["", "| 合同 | 纹理真值概率 | 纹理熵bit | 真macro内纹理熵bit | 真macro内纹理真值概率 | 全局wall像素方差/prior | 真macro内纹理像素方差/prior |", "|---|---:|---:|---:|---:|---:|---:|"])
    for r in summary:
        report.append("| %s | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f |" % (r["contract"],r["texture_probability_true"],r["texture_entropy_bits"],r["conditional_texture_entropy_true_macro_bits"],r["conditional_texture_probability_true_given_true_macro"],r["global_wall_pixel_variance_ratio_to_prior"],r["truth_macro_pixel_variance_ratio_to_prior"]))
    report.extend(["",f"prior 位置风险为 {prior_position_risk:.6g}。全部指标先在每个目标 parent 内平均两次噪声重复，再等权汇总 12 个目标 parent。2,000 次配对 bootstrap 在两个固定 macro 分层内各重采样六个目标 parent；区间仅为此有限设计的描述性波动，不表示从两个墙 parent 得到的总体覆盖。", "",
        "PAIRED_ADDED_DATA_DESCRIPTIVE.csv 保存 A1-A0、A2-A1、A2-A0 的目标风险、实际 RMSE、纹理熵和像素方差差值及描述性区间。数值差值不直接替代科学主张判决。", "",
        "RAW_JOINT_POSTERIORS.npz 保存全部 24×3×4×12 joint 后验、target/wall 边缘和给定 macro 的纹理后验；CASE_RESULTS.csv 保存各候选墙质量、MAP/mean 目标参数、HPD、风险、纹理熵和 wall-position 相关。POSTERIOR_WALL_PIXEL_DIAGNOSTICS.npz 保存同一物理网格的 prior/posterior 像素方差、macro 条件纹理方差、全方差分解和像素-位置相关；零方差下的相关标为未定义。", "",
        "六宏观墙结果与此正确指定的有限微观响应族是两个不同模型类。public_macro_comparison 保留 MODEL_MISMATCH_DIAGNOSTIC 标签；更换数据合同后的宏观模型误差不能直接解释为新增数据本身损失信息。数值归一化与有限恒等式检查见 REFERENCE_IDENTITY.json；它们不认证 FEM 离散误差、模型误差、舍入误差或整个连续参数族。", "",
        "本分析仅加载 144 个已付响应数组，新增 FE=0、RHS=0；未调参、未改旧在线结果。", ""])
    (out/"REPORT_ZH.md").write_text("\n".join(report),encoding="utf-8")
    print(json.dumps(plain(dict(status=result["status"],scope=SCOPE,
        summaries=[{k:r[k] for k in ("contract","position_MAP_RMSE","position_mean_RMSE","posterior_position_risk","texture_entropy_bits","conditional_texture_entropy_true_macro_bits","truth_macro_pixel_variance_ratio_to_prior")} for r in summary],
        actual_new_full_FE=0,actual_new_RHS=0))),flush=True)
    return result


if __name__=="__main__":main()
