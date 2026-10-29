"""Preregistered possibility-cloud statistics; imports no Maxwell/online engine.

NOT EXECUTED by author. Prediction-only loader reads sealed public outputs.
The only disk read of private truth is inside score_sealed_predictions(), after
all32/allarms sealing, frozen-statistics and separate root approval checks.
No posterior probability interpretation is attached to finite-cloud weights.
"""
from __future__ import annotations
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_k]='1'
import argparse,csv,json,math,time
from pathlib import Path
from collections import Counter
import numpy as np
P=Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT=P/'results/possibility_volume'
DEFAULT_STATISTICS=P/'config/possibility_statistics_frozen.json'
DEFAULT_PROTOCOL=P/'config/possibility_volume_protocol.json'


def read_json(path):return json.loads(Path(path).read_text())
def clean(value):
    if isinstance(value,np.ndarray):return [clean(v) for v in value.tolist()]
    if isinstance(value,np.generic):return clean(value.item())
    if isinstance(value,float) and not math.isfinite(value):return None
    if isinstance(value,dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    return value

def write_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():raise RuntimeError('Existing statistical artifact: no implicit overwrite or post-GT refit')
    path.write_text(json.dumps(clean(value),ensure_ascii=False,allow_nan=False,indent=2)+'\n')

def _material(path):
    # Public prediction state only. Do not import the engine's load_npz: that
    # would also materialize large fields/J and could expose inappropriate data.
    path=Path(path)
    with np.load(path,allow_pickle=False) as z:
        if 'material' not in z.files:raise ValueError('prediction material missing: '+str(path))
        a=np.asarray(z['material'],dtype=np.complex128).copy()
    if a.shape!=(576,) or not np.isfinite(a).all():raise ValueError('finite original576 complex material required')
    return a

def _counter(mapping,key):
    v=mapping.get(key)
    return float(v) if v is not None else None

def _qualification(row,material,anchor,protocol,anchor_only=False):
    phi=float(row['Phi']);legal=bool(row.get('legal',True if anchor_only else False))
    legal=legal and bool(np.isfinite(material).all()) and bool(np.min(material.imag)>=-1e-12)
    noise=phi<=float(protocol['primary']['noise_ball'])
    if anchor_only:return legal and noise
    radius=float(row.get('target_RMS',row['radius']))
    actual=float(np.linalg.norm(material-anchor)/24)
    geometric=abs(actual-radius)<=float(protocol['candidate']['radius_relative_tolerance'])*radius
    progress=float(row['signed_progress'])>=float(protocol['candidate']['finite_progress_fraction'])*radius
    result=legal and noise and geometric and progress
    if abs(actual-float(row['material_RMS']))>1e-10:raise ValueError('sealed material/claimed RMS mismatch')
    if bool(row.get('primary_qualified',False))!=result:raise ValueError('sealed primary qualification disagrees with material/noise/progress contract')
    return result


def load_sealed_predictions(output,statistics_path=DEFAULT_STATISTICS,protocol_path=DEFAULT_PROTOCOL):
    """Prediction-only adapter to possibility_volume.py's native directory.

    Never reads private_output, generation truth, or any ground-truth array.
    Root can call this for readiness before authorizing offline scoring.
    """
    output=Path(output).resolve();stats=read_json(statistics_path);protocol=read_json(protocol_path)
    if stats.get('approved_by_root') is not True:raise RuntimeError('statistics not root-frozen')
    if protocol.get('approved_by_root') is not True:raise RuntimeError('online protocol not root-frozen')
    seal_path=output/'PREDICTIONS_SEALED.json';seal=read_json(seal_path)
    if seal.get('status')!='ALL32_ONLINE_PREDICTIONS_SEALED' or seal.get('truth_read') is not False:raise RuntimeError('all32 before-truth prediction seal required')
    if seal.get('statistics_config_id')!=stats['id']:raise RuntimeError('statistics/seal identity mismatch')
    if seal['protocol_id']!=protocol['id']:raise RuntimeError('protocol/seal identity mismatch')
    if stats.get('expected_protocol_id')!=protocol['id']:raise RuntimeError('statistical protocol identity mismatch')
    if float(stats['frozen_utc'])>float(seal['utc']):raise RuntimeError('statistics must be frozen before final prediction seal')
    parents=[r['parent'] for r in protocol['parents']];arms=list(protocol['arms'])
    if len(parents)!=32 or len(set(parents))!=32 or arms!=stats['methods']:raise RuntimeError('exact32 paired parent/arm contract')
    if stats['parent_ids']!=parents or [r['parent'] for r in seal['parents']]!=parents:raise RuntimeError('parent selection/order mismatch')
    if protocol['K_total_per_arm']!=16 or stats['nonanchor_direction_budget']!=16:raise RuntimeError('16 direction-budget contract')
    if list(protocol['radii'])!=stats['radii']:raise RuntimeError('radius contract mismatch')
    if float(protocol['candidate']['finite_progress_fraction'])!=.5:raise RuntimeError('signed progress must be at least half target radius')
    if float(protocol['primary']['noise_ball'])!=stats['primary_online_qualification']['absolute_original_noise_Phi_max']:raise RuntimeError('absolute original95% noise-ball mismatch')
    with np.load(protocol['shared_setup']['path'],allow_pickle=False) as z:roi=np.asarray(z['ROI'],int).copy()
    if roi.shape!=(256,) or len(set(roi.tolist()))!=256:raise RuntimeError('original256 ROI required')
    sealed_files={item['path']:item['bytes'] for p in seal['parents'] for item in p['files']}
    def public_path(path):
        q=Path(path).resolve()
        if not q.is_relative_to(output):raise RuntimeError('prediction path escapes public bundle')
        if str(q) not in sealed_files or q.stat().st_size!=sealed_files[str(q)]:raise RuntimeError('prediction artifact not in seal or size changed')
        return q
    bundle=dict(output=output,statistics=stats,protocol=protocol,seal=seal,seal_path=seal_path,parents=parents,arms=arms,roi=roi,groups={},anchors={})
    for parent in parents:
        folder=output/'parents'/parent;complete=read_json(public_path(folder/'COMPLETE.json'))
        if complete.get('truth_read') is not False or complete.get('data_only_selection') is not True:raise RuntimeError('parent predictions not data-only')
        anchor=_material(public_path(folder/'anchor/ANCHOR.npz'))
        astat=read_json(public_path(folder/'anchor/ANCHOR_STATUS.json'))
        if 'legal' not in astat:raise RuntimeError('anchor legality must be recorded by frozen original decoder; no assumed-true fallback')
        anchor_qualified=_qualification(astat,anchor,anchor,protocol,True)
        bundle['anchors'][parent]=dict(material=anchor,Phi=float(astat['Phi']),qualified=anchor_qualified,status=astat)
        rows=complete['rows']
        if len(rows)!=48:raise RuntimeError('all48 direction records including failures must be sealed')
        order={r['id']:i for i,r in enumerate(complete['proposal_manifest'])}
        if len(order)!=16:raise RuntimeError('immutable16 direction order required')
        for arm in arms:
            selected=[r for r in rows if r['arm']==arm]
            if len(selected)!=16:raise RuntimeError('method missing direction/failure rows')
            selected=sorted(selected,key=lambda r:r.get('slot',order[r['proposal']]))
            if [r.get('slot',order[r['proposal']]) for r in selected]!=list(range(16)):raise RuntimeError('slots0..15 required')
            candidates=[];raw=[]
            for slot,row in enumerate(selected):
                expected_pair=slot//8;expected_radius=stats['radii'][(slot%8)//2];expected_sign=1 if slot%2==0 else -1
                if row['pair']!=expected_pair or row['sign']!=expected_sign or row['radius']!=expected_radius:raise RuntimeError('radius/pair/antithetic online order mismatch')
                if row.get('GT_used',False):raise RuntimeError('truth-tagged candidate in online predictions')
                r=dict(row,slot=slot,material=None,eligible=False)
                path=row.get('selected_path') or str(folder/arm/row['proposal']/'SELECTED.npz')
                if row['status']=='EXECUTED':
                    r['material']=_material(public_path(path));r['eligible']=_qualification(row,r['material'],anchor,protocol)
                elif row.get('primary_qualified',False):raise RuntimeError('failed direction claims qualification')
                candidates.append(r)
                visited_path=folder/arm/row['proposal']/'VISITED.json'
                if str(visited_path.resolve()) in sealed_files:
                    visited=read_json(public_path(visited_path))['rows']
                    for visit_index,v in enumerate(visited):
                        vr=dict(v,slot=slot,visit_index=visit_index,material=_material(public_path(v['x_path'])))
                        vr['eligible']=bool(v.get('legal',False)) and float(v['Phi'])<=float(protocol['primary']['noise_ball'])
                        raw.append(vr)
                    # Independent audit of lowest-Phi online endpoint selection.
                    qualified_visits=[v for v in visited if v.get('primary_qualified',False)]
                    if qualified_visits and (not r['eligible'] or float(r['Phi'])>min(float(v['Phi']) for v in qualified_visits)+1e-10):raise RuntimeError('selected endpoint is not minimum dataPhi among qualified visited candidates')
            bundle['groups'][(parent,arm)]=dict(selected=candidates,raw=raw,complete=complete,arm=arm)
        for slot in range(16):
            rr=[bundle['groups'][(parent,a)]['selected'][slot] for a in arms]
            tt=[r.get('timing') for r in rr]
            if any(t is None for t in tt) or len({t['window_id'] for t in tt})!=1:raise RuntimeError('same-window paired arm timing metadata required')
            pair=slot//8;ri=(slot%8)//2;parent_index=parents.index(parent);rot=(parent_index+ri+pair)%3;expected_order=arms[rot:]+arms[:rot]
            if any(t['arm_order']!=expected_order or t['arm_order_index']!=expected_order.index(a) for a,t in zip(arms,tt)):raise RuntimeError('frozen rotated timing order mismatch')
    return bundle


def _distances(X):
    if len(X)==0:return np.zeros((0,0))
    g=X@X.conj().T;q=np.sum(abs(X)**2,axis=1)
    return np.sqrt(np.maximum(q[:,None]+q[None,:]-2*g.real,0)/576)

def _unique(records,threshold):
    # Approximate equality is deliberately order-dependent, never called an
    # exact quotient or posterior branch mass. No GT is involved in this order.
    if not records:return []
    X=np.stack([r['material'] for r in records]);d=_distances(X);keep=[]
    for i in range(len(records)):
        if not keep or np.min(d[i,keep])>threshold:keep.append(i)
    return [records[i] for i in keep]

def _cloud(bundle,parent,arm,K,radius=None,raw=False):
    a=bundle['anchors'][parent];members=[]
    if a['qualified']:members.append(dict(material=a['material'],Phi=a['Phi'],kind='ANCHOR'))
    group=bundle['groups'][(parent,arm)];records=group['raw'] if raw else group['selected'][:K]
    members += [r for r in records if r['eligible'] and (radius is None or r.get('radius',r.get('target_RMS'))==radius)]
    return _unique(members,float(bundle['statistics']['uniqueness']['material_RMS_threshold']))


def _binomial(hit,confidence=.95):
    # Small-n coverage intervals supplement, not replace, paired bootstrap.
    from scipy.stats import norm,beta
    n=len(hit);s=int(np.sum(hit));p=s/n;alpha=1-confidence;z=float(norm.ppf(1-alpha/2));den=1+z*z/n
    mid=(p+z*z/(2*n))/den;half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    lo=0. if s==0 else float(beta.ppf(alpha/2,s,n-s+1));hi=1. if s==n else float(beta.ppf(1-alpha/2,s+1,n-s))
    one=0. if s==0 else float(beta.ppf(alpha,s,n-s+1))
    return dict(successes=s,parents=n,rate=p,Wilson95=[max(0.,mid-half),min(1.,mid+half)],Clopper_Pearson95=[lo,hi],exact_one_sided95_lower=one,scope='fresh-parent binomial model; not posterior confidence')


def _components(d,threshold):
    seen=set();count=0
    for i in range(len(d)):
        if i in seen:continue
        count+=1;todo=[i];seen.add(i)
        while todo:
            j=todo.pop()
            for k in np.where(d[j]<=threshold)[0]:
                if int(k) not in seen:seen.add(int(k));todo.append(int(k))
    return count


def _cloud_score(members,truth,anchor,roi,cfg):
    eps=[cfg['epsilon_primary']]+cfg['epsilon_secondary'];result=dict(qualified_unique_total=len(members),qualified_unique_nonanchor=sum(r.get('kind')!='ANCHOR' for r in members),available=bool(members))
    if not members:
        result.update(best_material_RMS=None,best_material_MSE=None,nearest_error_right_censored=True,hits={str(e):False for e in eps},pair_diameter=None,mean_pair_distance=None,RMS_dispersion=None,max_anchor_distance=None,tube={str(e):dict(outer_diameter=None,components=0,packing_count=0) for e in eps},task_coverage=None)
        X=anchor[None];phis=np.array([0.]);fallback=True
    else:
        X=np.stack([m['material'] for m in members]);phis=np.array([m['Phi'] for m in members]);err=np.sum(abs(X-truth)**2,axis=1)/576;j=int(np.argmin(err));d=_distances(X);upper=d[np.triu_indices(len(X),1)];center=np.mean(X,axis=0)
        result.update(best_material_RMS=float(np.sqrt(err[j])),best_material_MSE=float(err[j]),offline_nearest_member_index=j,nearest_error_right_censored=False,hits={str(e):bool(err[j]<=e*e) for e in eps},pair_diameter=float(np.max(d)),mean_pair_distance=float(np.mean(upper)) if len(upper) else 0.,RMS_dispersion=float(np.sqrt(np.mean(np.sum(abs(X-center)**2,axis=1)/576))),max_anchor_distance=float(np.max(np.linalg.norm(X-anchor,axis=1)/24)),tube={})
        for e in eps:
            packing=[]
            for i in range(len(X)):
                if not packing or np.min(d[i,packing])>e:packing.append(i)
            result['tube'][str(e)]=dict(outer_diameter=float(np.max(d)+2*e),components=_components(d,2*e),packing_count=len(packing))
        fallback=False
    # Actions are deployable approximate ensemble means, not posterior means.
    w=np.exp(-(phis-np.min(phis)));w/=np.sum(w);uniform=np.mean(X,axis=0);weighted=w@X
    result['actions']=dict(anchor_fallback_used=fallback,uniform_MSE=float(np.mean(abs(uniform-truth)**2)),likelihood_weighted_MSE=float(np.mean(abs(weighted-truth)**2)),anchor_baseline_MSE=float(np.mean(abs(anchor-truth)**2)),weight_semantics='uniform or exp(-Phi) on deduplicated finite cloud; no prior/proposal/Jacobian correction, NOT posterior')
    if members:
        real=X[:,roi].real;gt=truth[roi].real;means=real.mean(axis=1);counts=(real>=cfg['tasks']['pixel_threshold']).sum(axis=1);gtcount=int(np.sum(gt>=cfg['tasks']['pixel_threshold']));events=counts>=cfg['tasks']['event_count_threshold'];gtevent=gtcount>=cfg['tasks']['event_count_threshold'];temperature=cfg['tasks']['soft_occupancy_temperature']
        from scipy.special import expit
        soft=expit((real-cfg['tasks']['pixel_threshold'])/temperature);softgt=expit((gt-cfg['tasks']['pixel_threshold'])/temperature);softmean=soft.mean(axis=1)
        result['task_coverage']=dict(ROI_mean_range=[float(means.min()),float(means.max())],ROI_mean_range_width=float(means.max()-means.min()),ROI_mean_covered=bool(means.min()<=gt.mean()<=means.max()),nearest_ROI_mean_error=float(np.min(abs(means-gt.mean()))),threshold_count_range=[int(counts.min()),int(counts.max())],nearest_threshold_count_error=int(np.min(abs(counts-gtcount))),hard_event_represented=bool(np.any(events==gtevent)),uniform_event_score=float(np.mean(events)),likelihood_weighted_event_score=float(w@events),soft_ROI_mean_range=[float(softmean.min()),float(softmean.max())],soft_ROI_mean_width=float(softmean.max()-softmean.min()),soft_ROI_mean_covered=bool(softmean.min()<=softgt.mean()<=softmean.max()),soft_pixel_range_mean_width=float(np.mean(soft.max(axis=0)-soft.min(axis=0))),soft_pixel_range_max_width=float(np.max(soft.max(axis=0)-soft.min(axis=0))),soft_occupancy_note='scalar task and marginal pixel diagnostics only; pixel ranges never define the primary material tube')
    return result


def _cost(bundle,parent,arm):
    group=bundle['groups'][(parent,arm)];done=group['complete'];rows=group['selected'];N=len(bundle['parents']);M=len(bundle['arms']);setup=bundle['protocol']['shared_setup'];anchor=done.get('anchor_cost');frame=done.get('fixed_frame_setup',done.get('frame_setup_cost',done.get('setupmeta')));runtime=done.get('runtime_setup',{});overhead=done.get('common_unallocated_overhead',{})
    counters=('reduced_RHS256','reduced_LU256','full_RHS576','full_LU576')
    totals={key:sum(float(r.get('cost',{}).get(key,0)) for r in rows) for key in counters}
    wall=sum(float(r.get('wall_seconds',0)) for r in rows);CPU=None if any(r.get('cost',{}).get('CPU_seconds',r.get('CPU_seconds')) is None for r in rows) else sum(float(r.get('cost',{}).get('CPU_seconds',r.get('CPU_seconds'))) for r in rows)
    anchor_wall=None if anchor is None else _counter(anchor,'wall_seconds');frame_wall=None if frame is None else _counter(frame,'wall_seconds')
    if anchor_wall is not None and frame_wall is not None:standalone_wall=wall+anchor_wall+frame_wall+float(runtime.get('wall_seconds',0))+float(overhead.get('wall_seconds',0))/M
    else:standalone_wall=None
    for key in counters:
        if anchor is not None:totals[key]+=float(anchor.get(key,0))
        if frame is not None:totals[key]+=float(frame.get(key,0))
    # Cached setup264 is explicitly allocated in a separate full576 cost view.
    oldsetup=float(setup.get('already_paid_setup_full_RHS576',0));totals['allocated_historical_setup_RHS576']=oldsetup/N
    totals['historical_setup_study_RHS576']=setup.get('historical_study_full_RHS576')
    eligible=[r for r in rows if r['eligible']];censored=not eligible;first=done.get('time_to_first_qualified_nonanchor',{}).get(arm)
    if first is None and eligible and standalone_wall is not None:
        first=anchor_wall+frame_wall+float(runtime.get('wall_seconds',0))+sum(float(r.get('wall_seconds',0)) for r in rows[:eligible[0]['slot']+1])
    return dict(standalone_counts_anchor_and_frame_included=totals,all_attempt_online_wall_seconds=wall,all_attempt_CPU_seconds=CPU,anchor_wall_seconds=anchor_wall,frame_wall_seconds=frame_wall,runtime_setup_wall_seconds=runtime.get('wall_seconds'),allocated_common_overhead_wall_seconds=float(overhead.get('wall_seconds',0))/M,standalone_wall_seconds=standalone_wall,study_allocated_anchor_frame_wall_seconds=None if standalone_wall is None else wall+(anchor_wall+frame_wall+float(runtime.get('wall_seconds',0))+float(overhead.get('wall_seconds',0)))/M,historical_setup_wall_seconds=setup.get('already_paid_setup_wall_seconds'),cold_wall_seconds=(standalone_wall+float(setup['already_paid_setup_wall_seconds'])/N) if standalone_wall is not None and setup.get('already_paid_setup_wall_seconds') is not None else None,cold_wall_note='historical setup wall unpriced unless separately frozen; common unknown setup cancels paired arm difference but absolute cold total is NOT measured',first_qualified_time=first,time_to_first_right_censored=censored,time_to_first_censor_duration=standalone_wall if censored else None,time_to_first_definition='returned selected qualified alternative in cumulative method-only clock; includes preceding failed direction attempts and common acquisition/setup. No within-trial first-hit time is invented.',failed_or_unavailable_directions=sum(r['status']!='EXECUTED' for r in rows),unqualified_directions=sum(not r['eligible'] for r in rows),timing_windows=[r.get('timing',r.get('timing_window')) for r in rows],cost_note='all directions including failures retained; standalone allocation duplicates common anchor/setup across arm comparison views, never sum these as campaign expenditure')


def _paired_interval(a,b,indices,confidence):
    a=np.asarray(a,float);b=np.asarray(b,float)
    if not np.isfinite(a).all() or not np.isfinite(b).all():return dict(status='NOT_ESTIMABLE_WITH_EMPTY_OR_MISSING_PARENTS',parents=len(a),note='No complete-case exclusion or invented finite nearest-error/cost penalty.')
    delta=a-b;boots=delta[indices].mean(axis=1);alpha=1-confidence
    return dict(status='ESTIMATED',parents=len(a),difference=float(delta.mean()),confidence=confidence,percentile_interval=[float(np.quantile(boots,alpha/2)),float(np.quantile(boots,1-alpha/2))],unit='paired parent',bootstrap_resamples=len(indices))


def score_sealed_predictions(output,statistics_path,protocol_path,root_approval_path):
    """ONLY entry allowed to open private truth; no call has been made yet."""
    bundle=load_sealed_predictions(output,statistics_path,protocol_path);cfg=bundle['statistics'];protocol=bundle['protocol'];approval=read_json(root_approval_path)
    if approval.get('approved_by_root') is not True or approval.get('action')!='OFFLINE_TRUTH_SCORING_AFTER_ALL32_SEALED':raise RuntimeError('separate root approval for offline truth scoring missing')
    if approval.get('statistics_id')!=cfg['id'] or approval.get('protocol_id')!=protocol['id']:raise RuntimeError('score approval identity mismatch')
    if float(approval['utc'])<float(bundle['seal']['utc']):raise RuntimeError('score approval must follow all32 prediction seal')
    score_dir=Path(output)/'offline_statistics'
    if score_dir.exists():raise RuntimeError('offline scoring already started; no post-GT correction/refit implicit replay')
    score_dir.mkdir(parents=True)
    write_json(score_dir/'TRUTH_ACCESS_BEGIN.json',dict(utc=time.time(),prediction_seal=str(bundle['seal_path']),statistics_id=cfg['id'],protocol_id=protocol['id'],root_approval=str(root_approval_path),role='OFFLINE_COVERAGE_CAPACITY_AND_ACTION_SCORE_ONLY',model_changes_after_truth='FORBIDDEN'))
    # Do not move this read into the loader, generator, online engine, or imports.
    truth_by_parent={}
    for parent in bundle['parents']:
        truth_path=Path(protocol['private_output'])/'truth'/f'{parent}.npz'
        with np.load(truth_path,allow_pickle=False) as z:
            if str(z['parent'])!=parent:raise RuntimeError('truth parent identity mismatch')
            truth=np.asarray(z['material'],np.complex128).copy()
        if truth.shape!=(576,) or not np.isfinite(truth).all():raise RuntimeError('truth original576 material required')
        truth_by_parent[parent]=truth
    allrows=[];full={};cost={};radius_rows=[];raw_rows=[];eps=[cfg['epsilon_primary']]+cfg['epsilon_secondary'];parents=bundle['parents'];arms=bundle['arms']
    for parent in parents:
        anchor=bundle['anchors'][parent]
        for arm in arms:
            group=bundle['groups'][(parent,arm)];cost[(parent,arm)]=_cost(bundle,parent,arm)
            for K in cfg['K_prefixes']:
                members=_cloud(bundle,parent,arm,K);row=dict(parent=parent,method=arm,K_direction_budget=K,shared_anchor_noise_qualified=anchor['qualified'],qualified_direction_count=sum(r['eligible'] for r in group['selected'][:K]),**_cloud_score(members,truth_by_parent[parent],anchor['material'],bundle['roi'],cfg))
                row['qualified_unique_nonanchor_yield']=row['qualified_unique_nonanchor']/K if K else None
                allrows.append(row)
                if K==16:full[(parent,arm)]=row
            for radius in cfg['radii']:
                members=_cloud(bundle,parent,arm,16,radius=radius);radius_rows.append(dict(parent=parent,method=arm,radius=radius,direction_budget=4,**_cloud_score(members,truth_by_parent[parent],anchor['material'],bundle['roi'],cfg)))
            members=_cloud(bundle,parent,arm,16,raw=True);raw_rows.append(dict(parent=parent,method=arm,raw_visited_states=len(group['raw']),variableK=True,diagnostic_only=True,**_cloud_score(members,truth_by_parent[parent],anchor['material'],bundle['roi'],cfg)))
    rng=np.random.default_rng(cfg['uncertainty']['statistical_rng_seed']);indices=rng.integers(0,len(parents),size=(cfg['uncertainty']['paired_bootstrap_resamples'],len(parents)))
    coverage=[]
    for arm in arms:
        for K in cfg['K_prefixes']:
            rr=[r for r in allrows if r['method']==arm and r['K_direction_budget']==K]
            for e in eps:coverage.append(dict(method=arm,K=K,epsilon=e,**_binomial([r['hits'][str(e)] for r in rr])))
    fixed=cfg['primary_comparison']['candidate'];dynamic=cfg['primary_comparison']['reference'];random=cfg['primary_comparison']['random_reference'];e=str(cfg['epsilon_primary']);conf=cfg['uncertainty']['familywise_two_sided_confidence']
    vec=lambda arm,key: [full[(p,arm)][key] if full[(p,arm)][key] is not None else np.nan for p in parents]
    hits=lambda arm:[float(full[(p,arm)]['hits'][e]) for p in parents]
    primary=dict(coverage_fixed_minus_dynamic=_paired_interval(hits(fixed),hits(dynamic),indices,conf),coverage_fixed_minus_random=_paired_interval(hits(fixed),hits(random),indices,conf),offline_nearest_MSE_fixed_minus_dynamic=_paired_interval(vec(fixed,'best_material_MSE'),vec(dynamic,'best_material_MSE'),indices,conf),all_attempt_standalone_wall_fixed_minus_dynamic=_paired_interval([cost[(p,fixed)]['standalone_wall_seconds'] if cost[(p,fixed)]['standalone_wall_seconds'] is not None else np.nan for p in parents],[cost[(p,dynamic)]['standalone_wall_seconds'] if cost[(p,dynamic)]['standalone_wall_seconds'] is not None else np.nan for p in parents],indices,conf))
    y=_paired_interval(vec(fixed,'qualified_unique_nonanchor_yield'),vec(dynamic,'qualified_unique_nonanchor_yield'),indices,.95)
    if y['status']=='ESTIMATED':y['secondary_noninferiority_margin']=-.10;y['lower_interval_above_margin']=y['percentile_interval'][0]>=-.10;y['interpretation']='secondary predeclared small-n bootstrap diagnostic; not a nominal noninferiority proof'
    actions={}
    for arm in arms:
        actions[arm]={k:dict(mean=float(np.mean([full[(p,arm)]['actions'][k] for p in parents])),parent_losses=[full[(p,arm)]['actions'][k] for p in parents],bootstrap95_percentile=[float(np.quantile(np.array([full[(p,arm)]['actions'][k] for p in parents])[indices].mean(axis=1),q)) for q in (.025,.975)]) for k in ('uniform_MSE','likelihood_weighted_MSE','anchor_baseline_MSE')}
    secondary_actions={}
    for k in ('uniform_MSE','likelihood_weighted_MSE'):
        secondary_actions[k]={}
        for ref in (dynamic,random):
            secondary_actions[k]['fixed_minus_'+ref]=_paired_interval([full[(p,fixed)]['actions'][k] for p in parents],[full[(p,ref)]['actions'][k] for p in parents],indices,.95)
    task_summary={};width_summary={}
    for arm in arms:
        rr=[full[(p,arm)] for p in parents];task_summary[arm]={}
        for name in ('ROI_mean_covered','hard_event_represented','soft_ROI_mean_covered'):
            task_summary[arm][name]=_binomial([False if r['task_coverage'] is None else bool(r['task_coverage'][name]) for r in rr])
        available=[r for r in rr if r['available']]
        width_summary[arm]=dict(parents=len(parents),available_parents=len(available),empty_parents=len(rr)-len(available),pair_diameter_mean=None if not available else float(np.mean([r['pair_diameter'] for r in available])),pair_diameter_median=None if not available else float(np.median([r['pair_diameter'] for r in available])),qualified_unique_nonanchor_mean=float(np.mean([r['qualified_unique_nonanchor'] for r in rr])),availability_note='width summaries conditional on nonempty cloud; coverage and risk denominators always all32')
    width_bins=[];edges=cfg['width_bins_for_descriptive_coverage']
    for arm in arms:
        for low,high in zip(edges[:-1],edges[1:]):
            rr=[full[(p,arm)] for p in parents if full[(p,arm)]['pair_diameter'] is not None and low<=full[(p,arm)]['pair_diameter']<high]
            width_bins.append(dict(method=arm,diameter_bin=[low,high],parents=len(rr),coverage=None if not rr else _binomial([r['hits'][e] for r in rr]),scope='predeclared descriptive width stratum, pointwise only; no joint box coverage or adaptive claim selection'))
    claimrow=next(r for r in coverage if r['method']==fixed and r['K']==16 and r['epsilon']==cfg['epsilon_primary']);probable=dict(**cfg['probable_GT_claim'],observed_parent_rate=claimrow['rate'],exact_one_sided95_lower=claimrow['exact_one_sided95_lower'],empirical_population_claim_supported=claimrow['exact_one_sided95_lower']>=1-cfg['probable_GT_claim']['delta'],scope='fresh-parent sampling population only; not posterior probability, not guaranteed coverage for a particular measurement')
    summary=dict(status='OFFLINE_SCORED_NO_POSTERIOR_CLAIM',parent_count=len(parents),bootstrap_unit='parent',bootstrap_resamples=len(indices),primary_four_contrasts=primary,secondary_yield_noninferiority=y,coverage=coverage,deployable_approximate_actions=actions,secondary_paired_action_differences=secondary_actions,material_task_coverage=task_summary,cloud_width_summary=width_summary,descriptive_width_coverage_bins=width_bins,probable_GT_empirical_claim=probable,anchor_unqualified_parent_count=sum(not bundle['anchors'][p]['qualified'] for p in parents),primary_empty_cloud_counts={a:sum(not full[(p,a)]['available'] for p in parents) for a in arms},claim_boundaries=['best-of-K is offline coverage capacity; nearest GT never selects an online endpoint, weight, stopping step, direction or model','C_epsilon is union of balls around actual material members, not a pixelwise box and not a 576D volume estimate','all32 parents retained, including zero-qualified clouds and failed directions','four primary contrasts use98.75% paired-parent percentile intervals; small-n bootstrap remains approximate','all curves/tasks/width strata are descriptive secondary results','original R9 failures remain unchanged; posterior accuracy is not a main gate'])
    write_json(score_dir/'SUMMARY.json',summary);write_json(score_dir/'PARENT_PREFIX_SCORES.json',allrows);write_json(score_dir/'PARENT_RADIUS_SCORES.json',radius_rows);write_json(score_dir/'RAW_NOISE_QUALIFIED_DIAGNOSTIC.json',raw_rows);write_json(score_dir/'FULL_SYSTEM_COSTS.json',[dict(parent=p,method=a,**cost[(p,a)]) for p in parents for a in arms])
    # Width/coverage/cost points retain parent identity for non-vacuity tradeoffs.
    points=[dict(parent=p,method=a,pair_diameter=full[(p,a)]['pair_diameter'],qualified_unique=full[(p,a)]['qualified_unique_total'],hit_epsilon075=full[(p,a)]['hits'][e],nearest_MSE=full[(p,a)]['best_material_MSE'],standalone_wall_seconds=cost[(p,a)]['standalone_wall_seconds'],actions=full[(p,a)]['actions']) for p in parents for a in arms]
    write_json(score_dir/'WIDTH_COVERAGE_COST_POINTS.json',points)
    accounting=dict(authoritative_campaign_resource_snapshot=read_json(Path(output)/'RESOURCE_SNAPSHOT.json'),generation_receipt=read_json(Path(output)/'GENERATION_SEALED.json'),parent_end_to_end_wall_seconds={p:bundle['groups'][(p,arms[0])]['complete'].get('parent_end_to_end_wall_seconds') for p in parents},cost_scope='campaign includes paid data generation, all anchor/correction attempts and failures. Deployment comparisons allocate common anchor/runtime/frame once per arm; shared per-parent reporting overhead equally across arms. Historical setup264RHS576 allocated1/32 per standalone parent, not recharged to fresh campaign. Neither reducedRHS256 nor LU256 are equated to full576 work.')
    write_json(score_dir/'COST_ACCOUNTING.json',accounting)
    write_json(score_dir/'SCORING_COMPLETE.json',dict(utc=time.time(),status='COMPLETE',statistics_id=cfg['id'],new_Maxwell_RHS=0,new_forward_calls=0,no_post_truth_model_changes=True))
    return summary


def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['validate_predictions','score']);parser.add_argument('--output',type=Path,default=DEFAULT_OUTPUT);parser.add_argument('--statistics',type=Path,default=DEFAULT_STATISTICS);parser.add_argument('--protocol',type=Path,default=DEFAULT_PROTOCOL);parser.add_argument('--root-approval',type=Path)
    a=parser.parse_args()
    if a.command=='validate_predictions':
        b=load_sealed_predictions(a.output,a.statistics,a.protocol);print(json.dumps(dict(status='SEALED_PREDICTIONS_VALIDATED_NO_TRUTH_READ',parents=len(b['parents']),groups=len(b['groups']))));return
    if a.root_approval is None:raise RuntimeError('--root-approval is required before any private truth read')
    print(json.dumps(clean(score_sealed_predictions(a.output,a.statistics,a.protocol,a.root_approval)),ensure_ascii=False,allow_nan=False))
if __name__=='__main__':main()
