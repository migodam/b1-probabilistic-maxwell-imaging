"""E0 read-only trajectory attribution. Never imports/constructs a Maxwell solver.

Only saved online arrays/history are read. NPZ __OFFLINE_TRUTH members are never
accessed. The independent field metrics below reproduce the historical formulas
without importing research_a.target (which imports physics).
"""
from __future__ import annotations
import csv
import json
import os
from pathlib import Path
import sys
import time
from collections import Counter, defaultdict
import numpy as np
from scipy.ndimage import gaussian_filter, label

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'research/research_g_20261015'
A = ROOT / 'research/research_a_20261015'
DELEG = OUT / 'research/delegated/g_trajectory'
PLAN = DELEG / 'ANALYSIS_PLAN.json'
START = time.perf_counter()
CPU = time.process_time()
READ = []
TAU_RADIAL = 6.796635897286134  # sqrt(chi2.ppf(.95, 32)), inherited exactly
# E0 sensitivity descriptor, not a new acceptance rule or tuning target.
SPEC = {
    'version': 'G-multiscale-support-v1 (2D) / g-e0-radial-support-v1',
    'primary_threshold': {'A30_2D': .08, 'radial': .025},
    'multipliers': [.5, 1., 2.],
    'radial': {'grid': '512 midpoint uniform radius samples, R=1', 'radii_smoothing_R': [.01, .025, .05], 'weights': 'exact spherical shell volume', 'support': 'decoded epsilon != vacuum; exact r/h decoding', 'boundary': 'nearest extension at r=0 and R'},
    'A30': {'grid': [24, 24], 'smoothing_pixels': [0., 1., 2.], 'contrast_amplitude': 'abs(chi-background)/1.2', 'occupancy': 'clip((gaussian_smoothed_amplitude-.20)/.10,0,1)', 'boundary': 'constant zero', 'additional_distinctness_requirement': 'unsmoothed abs(chi-background)/1.2>.25 masks XOR >=6 cells', 'foreground_connectivity': 4, 'background_connectivity': 8, 'minimum_feature_cells': 6, 'new_material_distance': 'sqrt(mean((delta_real/4.5)^2+(delta_imag/.5)^2)), epsilon .03'},
    'distance': '2D RMS of three occupancy views plus mandatory >=6 mask XOR; radial shell-volume-weighted RMS of smoothed support',
    'interpretation': 'E0 diagnostic separation only; finite sampled descriptor, not topology certification or G E1 discovery threshold',
    'connectivity': 'support descriptor does not count components/holes; no digital topology claim'
}

def load(path):
    READ.append(str(path.relative_to(ROOT)))
    return json.loads(path.read_text())


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def csvwrite(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, keys)
        w.writeheader()
        for row in rows:
            w.writerow({k: json.dumps(v) if isinstance(v, (list, dict)) else v for k, v in row.items()})


def radial_field(p, radius):
    r, h, er, ei = np.asarray(p, float)
    mask = (radius >= r) & ~((h > 0) & (radius >= .65-h) & (radius < .65+h))
    return np.where(mask, er + 1j*ei, 1.+0j)


def radial_breaks(*points):
    edges = {0., 1.}
    for p in points:
        r, h = np.asarray(p)[:2]
        if r > 0: edges.add(float(r))
        if h > 0: edges.update([float(.65-h), float(.65+h)])
    edges = np.array(sorted(edges))
    return (edges[:-1]+edges[1:])/2, np.diff(edges**3)


def radial_inner(a, b, c):
    mid, vol = radial_breaks(a, b, c)
    e = radial_field(b, mid) - radial_field(a, mid)
    u = radial_field(c, mid) - radial_field(a, mid)
    denom = float(np.sum(vol*abs(e)**2))
    return float(np.real(np.sum(vol*np.conj(e)*u))/denom) if denom > 1e-24 else None


def radial_distance(a, b):
    mid, vol = radial_breaks(a, b)
    return float(np.sqrt(np.sum(vol*abs(radial_field(a, mid)-radial_field(b, mid))**2))/abs(1.5+.1j))


def radial_geometry(a, b):
    return bool(abs(a[0]-b[0]) >= .05 or abs(a[1]-b[1]) >= .025)


def grid_distance(a, b):
    return float(np.sqrt(np.mean(abs(a-b)**2))/1.2)


def grid_primary_distance(a, b):
    d = a-b
    return float(np.sqrt(np.mean((d.real/4.5)**2+(d.imag/.5)**2)))


def grid_inner(a, b, c):
    e, u = b-a, c-a
    denom = float(np.sum(abs(e)**2))
    return float(np.real(np.vdot(e, u))/denom) if denom > 1e-24 else None


def geom_feature(value, platform, bg=None):
    if platform == 'radial':
        radius = (np.arange(512)+.5)/512
        mask = abs(radial_field(value, radius)-1.) > 1e-12
        return np.array([gaussian_filter(mask.astype(float), s*512, mode='nearest') for s in [.01,.025,.05]])
    amplitude = (abs(value-bg)/1.2).reshape(24,24)
    views = [np.clip(((gaussian_filter(amplitude,s,mode='constant',cval=0.) if s else amplitude)-.20)/.10,0.,1.) for s in [0.,1.,2.]]
    return np.array(views+[(amplitude>.25).astype(float)])


def geom_distance(a, b, platform):
    delta = (a-b)**2
    if platform == 'radial':
        edges = np.linspace(0.,1.,513)
        return float(np.sqrt(np.mean(np.sum(delta*np.diff(edges**3), axis=1))))
    if np.count_nonzero(a[3] != b[3]) < 6:
        return 0.
    return float(np.sqrt(np.mean(delta[:3])))


def pack_count(items, distance, threshold=.03):
    reps = []
    for item in items:
        if all(distance(item, r) >= threshold for r in reps): reps.append(item)
    return len(reps)


def geom_original_count(items, different):
    reps = []
    for item in items:
        if all(different(item, r) for r in reps): reps.append(item)
    return len(reps)


def stage_counts(items, distance, different, platform, bg=None):
    if not items:
        return {'original_material_distinct': 0, 'original_geometry_distinct': 0, 'G_primary_material_distinct': 0, 'multiscale_structure_distinct': {str(m): 0 for m in [.5,1.,2.]}, 'states': 0}
    feats = [geom_feature(x, platform, bg) for x in items]
    return dict(states=len(items), original_material_distinct=pack_count(items,distance),
                original_geometry_distinct=geom_original_count(items,different),
                G_primary_material_distinct=pack_count(items,grid_primary_distance if platform=='A30' else distance),
                multiscale_structure_distinct={str(m):pack_count(feats,lambda a,b:geom_distance(a,b,platform), (.08 if platform=='A30' else .025)*m) for m in [.5,1.,2.]})


def radial_audit():
    campaign_path = A/'results/CONFIRM_CAMPAIGN.json'
    campaign = load(campaign_path)
    rows, panels, improvements = [], [], []
    for parent in campaign:
        for method in parent['methods']:
            arm = method['arm']
            base = parent['base']
            archive = [{'id':'base','p':base['p'],'residual':base['residual']}] if base['valid'] else []
            raw_points, corrected_points, fullvalid_points = [], [], []
            for step, end in enumerate(method['attempts']):
                h = end.get('history',{})
                raw, point = np.asarray(h['proposal_p']), np.asarray(end['p'])
                par = np.asarray(h['parent_p'])
                prefix = archive.copy()
                d_raw = min([radial_distance(raw, c['p']) for c in prefix] or [np.inf])
                d_end = min([radial_distance(point, c['p']) for c in prefix] or [np.inf])
                nearest = min(prefix,key=lambda c:radial_distance(point,c['p'])) if prefix else None
                matches = [c for c in prefix if radial_distance(point,c['p']) < .03]
                duplicate = bool(end['valid'] and matches)
                raw_valid = h.get('raw_proxy_residual',np.inf) <= TAU_RADIAL
                stored = end.get('discovery_status')
                replay = 'DATA_INCOMPATIBLE' if not end['valid'] else 'DUPLICATE_MATERIAL' if duplicate else 'NEW_MATERIAL'
                if replay != stored: raise AssertionError(f'Radial replay mismatch: {parent["id"]}/{arm}/{step}: {replay}!={stored}')
                if duplicate and any(end['residual'] < c['residual']-1e-12 for c in matches):
                    improvements.append(dict(platform='radial',parent_id=parent['id'],arm=arm,endpoint_id=f'step_{step}',residual=end['residual'],
                        closest_retained_id=nearest['id'],closest_retained_residual=nearest['residual'],improvement_over_nearest=nearest['residual']-end['residual'],
                        improves_all_matching_representatives=all(end['residual'] < c['residual']-1e-12 for c in matches)))
                if replay == 'NEW_MATERIAL': archive.append({'id':f'step_{step}','p':point.tolist(),'residual':end['residual']})
                box=h.get('corrector',{}).get('protected_geometry')
                pbound=np.array([[0,.2925518253046813],[0,.11702073012187252],[2.2,2.8],[.04,.16]])
                bound_count=int(np.sum(np.minimum(abs(point-pbound[:,0]),abs(point-pbound[:,1])) < 1e-8))
                row = dict(platform='radial',parent_id=parent['id'],arm=arm,trajectory=step,
                    raw_source=h.get('source'),spatial_scale=None if box is None else [box[1]-box[0],box[3]-box[2]],
                    raw_full_valid=None,raw_full_valid_reason='MISSING: only raw proxy residual recorded',
                    raw_proxy_compatible=bool(raw_valid),endpoint_full_valid=end['valid'],archive_status=replay,
                    full_invalid_to_valid=None,proxy_invalid_to_full_valid=bool(not raw_valid and end['valid']),
                    proxy_valid_to_full_invalid=bool(raw_valid and not end['valid']),
                    raw_to_endpoint_material_distance=radial_distance(raw,point),
                    raw_to_parent_material_distance=radial_distance(raw,par),endpoint_to_parent_material_distance=radial_distance(point,par),
                    edit_direction_retention=radial_inner(par,raw,point),
                    raw_geometry_distinct_parent=radial_geometry(raw,par),endpoint_geometry_distinct_parent=radial_geometry(point,par),
                    structural_edit_retained=bool(not radial_geometry(raw,par) or radial_geometry(point,par)),
                    raw_novel_to_prefix=d_raw>=.03,endpoint_novel_to_prefix=d_end>=.03,
                    raw_novel_collapsed_to_prefix=bool(d_raw>=.03 and d_end<.03),
                    raw_distinct_parent_collapsed=bool(radial_distance(raw,par)>=.03 and radial_distance(point,par)<.03),
                    correction_converged=h.get('corrector',{}).get('converged'),bound_saturation_coordinates=bound_count,
                    protected_cells_unchanged=None,trust_region_exit=None,trust_region_reason='MISSING: no trajectory trust region recorded',
                    raw_residual=h.get('raw_proxy_residual'),raw_residual_kind='proxy',endpoint_residual=end['residual'],
                    source=f'{campaign_path.relative_to(ROOT)}#parent={parent["id"]}/method={arm}/attempt={step}')
                rows.append(row);raw_points.append(raw);corrected_points.append(point)
                if end['valid']: fullvalid_points.append(point)
            if len(archive)!=len(method['archive']): raise AssertionError('Archive length differs')
            for r,s in zip(archive,method['archive']):
                if radial_distance(r['p'],s['p'])>1e-10: raise AssertionError('Archive contents differ')
            panels.append(dict(platform='radial',parent_id=parent['id'],arm=arm,base_valid=base['valid'],
                raw_proposals=stage_counts(raw_points,radial_distance,radial_geometry,'radial'),
                corrected_endpoints=stage_counts(corrected_points,radial_distance,radial_geometry,'radial'),
                full_valid_endpoints=stage_counts(fullvalid_points,radial_distance,radial_geometry,'radial'),
                historical_archive=stage_counts([x['p'] for x in archive],radial_distance,radial_geometry,'radial'),
                raw_full_validation_count=None,raw_full_validation_reason='MISSING: ordinary proposals used proxy; see separately paid no-comp control',
                source=f'{campaign_path.relative_to(ROOT)}#parent={parent["id"]}/method={arm}'))
    attribution=load(A/'results/STRATIFICATION_ATTRIBUTION.json')
    pairwise=load(A/'results/PAIRWISE_AMBIGUITY.json')
    nocomp=load(A/'results/NO_COMPENSATION_ABLATION.json')
    summary=defaultdict(dict)
    for arm in dict.fromkeys(r['arm'] for r in rows):
        rr=[r for r in rows if r['arm']==arm];pp=[p for p in panels if p['arm']==arm]
        summary[arm]=summarize(rr)
        summary[arm]['mean_counts_by_stage']={stage:{key:float(np.mean([p[stage][key] for p in pp])) for key in ['states','original_material_distinct','original_geometry_distinct','G_primary_material_distinct']} for stage in ['raw_proposals','corrected_endpoints','full_valid_endpoints','historical_archive']}
        summary[arm]['historical_archive_total']=sum(p['historical_archive']['states'] for p in pp)
        summary[arm]['multiscale_archive_total_by_threshold']={str(m):sum(p['historical_archive']['multiscale_structure_distinct'][str(m)] for p in pp) for m in [.5,1.,2.]}
    return dict(parents=len(campaign),episodes=len(panels),summary=dict(summary),panels=panels,
        historical_attribution={'records':len(attribution['records']),'all_fullvalidated':attribution['all_fullvalidated'],'maximum_material_endpoint_difference':attribution['maximum_material_endpoint_difference'],'source':'research/research_a_20261015/results/STRATIFICATION_ATTRIBUTION.json'},
        pairwise_ambiguity={'structural_pairs':pairwise['total_pairs'],'pairs_at_whitened_distance_le1':pairwise['pairs_at_noise_distance_le1'],'source':'research/research_a_20261015/results/PAIRWISE_AMBIGUITY.json'},
        no_compensation_control=nocomp),rows,improvements


def summarize(rows):
    counts=Counter(r['archive_status'] for r in rows)
    return dict(trajectories=len(rows),archive_status_counts=dict(counts),valid_endpoints=sum(r['endpoint_full_valid'] for r in rows),
        invalid_to_valid=sum(bool(r.get('full_invalid_to_valid')) for r in rows) if any(r.get('raw_full_valid') is not None for r in rows) else None,
        proxy_invalid_to_full_valid=sum(bool(r.get('proxy_invalid_to_full_valid')) for r in rows),
        valid_to_invalid=sum(bool(r.get('full_valid_to_invalid')) for r in rows) if any(r.get('raw_full_valid') is not None for r in rows) else None,
        proxy_valid_to_full_invalid=sum(bool(r.get('proxy_valid_to_full_invalid')) for r in rows),
        raw_novel_collapsed_to_prefix=sum(bool(r.get('raw_novel_collapsed_to_prefix')) for r in rows),
        raw_distinct_parent_collapsed=sum(bool(r.get('raw_distinct_parent_collapsed')) for r in rows),
        structurally_edited=sum(bool(r['raw_geometry_distinct_parent']) for r in rows),
        structurally_edited_and_retained=sum(bool(r['raw_geometry_distinct_parent'] and r['endpoint_geometry_distinct_parent']) for r in rows),
        median_correction_material_distance=float(np.median([r['raw_to_endpoint_material_distance'] for r in rows])),
        bound_saturated_endpoints=sum(r['bound_saturation_coordinates']>0 for r in rows),
        max_protected_cells_defect=max([r['protected_defect_max'] for r in rows if r.get('protected_defect_max') is not None],default=None))


def digital_holes(image, background):
    # 2D diagnostic only: foreground 4 / complementary background 8.
    foreground=(abs(image-background)/1.2>.25).reshape(24,24)
    components,n=label(~foreground,structure=np.ones((3,3),int))
    boundary=set(np.r_[components[0],components[-1],components[:,0],components[:,-1]].tolist())
    return [components==i for i in range(1,n+1) if i not in boundary and np.count_nonzero(components==i)>=6]


def hole_transition(base, raw, endpoint, bg):
    bh,rh,eh=[digital_holes(x,bg) for x in [base,raw,endpoint]]
    born=[r for r in rh if not any(np.count_nonzero(r&b)/np.count_nonzero(r)>=.5 for b in bh)]
    retained=[r for r in born if any(np.count_nonzero(r&e)/np.count_nonzero(r)>=.5 for e in eh)]
    return dict(base_digital_holes=len(bh),raw_digital_holes=len(rh),endpoint_digital_holes=len(eh),raw_new_digital_holes=len(born),retained_new_digital_holes=len(retained),erased_new_digital_holes=len(born)-len(retained))


def grid_audit():
    history_path=A/'results/E2/SEARCH_HISTORY.json';history=load(history_path)
    protocol=load(A/'results/E2/PROTOCOL.json');metrics=load(A/'results/E2/METRICS.json')
    saved_summary=load(A/'results/E2/RESULT_SUMMARY.json')
    mm={(x['case_id'],x['arm']):x for x in metrics}
    arrays_path=A/'results/E2/MATERIAL_AND_PREDICTION_ARRAYS.npz'
    READ.append(str(arrays_path.relative_to(ROOT))+' (selected online keys only; no OFFLINE_TRUTH access)')
    arrays=np.load(arrays_path,allow_pickle=False)
    rows,panels,improvements,states=[],[],[],[]
    for case in protocol['cases']:
        bg=arrays[case+'__background'];base=arrays[case+'__base']
        base_row=next(r for r in history if r['case_id']==case and r['arm']=='base')
        for arm in protocol['arms']:
            rr=[r for r in history if r['case_id']==case and r['arm']==arm]
            archive=[dict(id='base',chi=base,residual=base_row['residual'])] if base_row['compatible'] else []
            trajectories=sorted(set(r['trajectory'] for r in rr))
            raw_images,endpoint_images,valid_images=[],[],[]
            prefix_for_traj={}
            for row in rr:
                image=arrays[row['array_prefix']+'__chi']
                valid=bool(row['compatible'])
                matches=[c for c in archive if grid_distance(image,c['chi'])<.03]
                status='DATA_INCOMPATIBLE' if not valid else 'DUPLICATE_MATERIAL' if matches else 'NEW_MATERIAL'
                if row['step']==0:prefix_for_traj[row['trajectory']]=[c['chi'].copy() for c in archive]
                if valid and matches:
                    nearest=min(matches,key=lambda c:grid_distance(image,c['chi']))
                    if any(row['residual']<c['residual']-1e-12 for c in matches):
                        improvements.append(dict(platform='A30_2D',parent_id=mm[case,arm]['parent_id'],case_id=case,arm=arm,
                            endpoint_id=row['array_prefix'],residual=row['residual'],closest_retained_id=nearest['id'],
                            closest_retained_residual=nearest['residual'],improvement_over_nearest=nearest['residual']-row['residual'],
                            improves_all_matching_representatives=all(row['residual']<c['residual']-1e-12 for c in matches)))
                if status=='NEW_MATERIAL':archive.append(dict(id=row['array_prefix'],chi=image,residual=row['residual']))
                states.append(dict(platform='A30_2D',parent_id=mm[case,arm]['parent_id'],case_id=case,arm=arm,trajectory=row['trajectory'],step=row['step'],
                    full_valid=valid,archive_status=status,residual=row['residual'],array_key=row['array_prefix']+'__chi',
                    source=f'{history_path.relative_to(ROOT)}#array_prefix={row["array_prefix"]}'))
            if len(archive)!=mm[case,arm]['valid_distinct_material_count']:raise AssertionError('E2 material archive replay mismatch')
            different=lambda a,b: np.count_nonzero((abs(a-bg)>.18)!=(abs(b-bg)>.18))>=6
            for t in trajectories:
                tr=sorted([r for r in rr if r['trajectory']==t],key=lambda r:r['step']);first,last=tr[0],tr[-1]
                raw=arrays[first['array_prefix']+'__chi'];end=arrays[last['array_prefix']+'__chi']
                prefix=prefix_for_traj[t]
                d_raw=min([grid_distance(raw,c) for c in prefix] or [np.inf]);d_end=min([grid_distance(end,c) for c in prefix] or [np.inf])
                states_last=next(s for s in reversed(states) if s['array_key']==last['array_prefix']+'__chi')
                pd=max(r.get('protected_defect',0.) for r in tr)
                actual_protected_defect=None
                edit=first.get('edit',{})
                if first.get('protected_cells',0)>0 and 'center_pixels' in edit:
                    yy,xx=np.mgrid[:24,:24];cy,cx=edit['center_pixels'];rad=edit['radius_pixels']
                    protected=((yy-cy)**2+(xx-cx)**2<=rad**2).ravel()
                    if int(protected.sum())!=first['protected_cells']:raise AssertionError('Protected mask replay count mismatch')
                    actual_protected_defect=max(float(np.max(abs(arrays[r['array_prefix']+'__chi'][protected]-bg[protected]))) for r in tr)
                    if actual_protected_defect>1e-12:raise AssertionError('Saved protected cells were changed')
                box=np.minimum.reduce([abs(end.real+.5),abs(end.real-4),abs(end.imag),abs(end.imag-.5)])
                row=dict(platform='A30_2D',parent_id=mm[case,arm]['parent_id'],case_id=case,arm=arm,trajectory=t,round=first.get('round'),real_parameters=first.get('real_parameters'),
                    raw_source=first.get('edit',{}).get('edit'),spatial_scale=first.get('edit',{}).get('radius_pixels'),
                    raw_full_valid=bool(first['compatible']),endpoint_full_valid=bool(last['compatible']),archive_status=states_last['archive_status'],
                    full_invalid_to_valid=bool(not first['compatible'] and last['compatible']),full_valid_to_invalid=bool(first['compatible'] and not last['compatible']),
                    raw_to_endpoint_material_distance=grid_distance(raw,end),raw_to_parent_material_distance=grid_distance(raw,base),endpoint_to_parent_material_distance=grid_distance(end,base),
                    G_primary_raw_to_endpoint_material_distance=grid_primary_distance(raw,end),G_primary_raw_to_parent_material_distance=grid_primary_distance(raw,base),G_primary_endpoint_to_parent_material_distance=grid_primary_distance(end,base),
                    edit_direction_retention=grid_inner(base,raw,end),raw_geometry_distinct_parent=bool(different(raw,base)),endpoint_geometry_distinct_parent=bool(different(end,base)),
                    raw_novel_to_prefix=bool(d_raw>=.03),endpoint_novel_to_prefix=bool(d_end>=.03),raw_novel_collapsed_to_prefix=bool(d_raw>=.03 and d_end<.03),
                    raw_distinct_parent_collapsed=bool(grid_distance(raw,base)>=.03 and grid_distance(end,base)<.03),
                    structural_edit_retained=bool(not different(raw,base) or different(end,base)),
                    protected_cells=first.get('protected_cells'),protected_defect_max=pd,protected_defect_recomputed_from_saved_arrays=actual_protected_defect,protected_cells_unchanged=bool(pd==0) if first.get('protected_cells',0)>0 else None,
                    raw_residual=first['residual'],raw_residual_kind='full-discrete',endpoint_residual=last['residual'],bound_saturation_coordinates=int(np.count_nonzero(box<1e-10)),
                    trust_region_exit=None,trust_region_reason='MISSING: no per-trajectory trust-exit field; QP trust=.35 is implementation-level',
                    source=f'{history_path.relative_to(ROOT)}#case={case}/arm={arm}/trajectory={t}',raw_array=first['array_prefix']+'__chi',endpoint_array=last['array_prefix']+'__chi')
                row.update(hole_transition(base,raw,end,bg))
                row['hole_attribution'] = 'NO_RAW_HOLE_BIRTH' if row['raw_new_digital_holes']==0 else 'ERASED_HOLE' if row['retained_new_digital_holes']==0 else 'RETAINED_BUT_INVALID' if not last['compatible'] else 'VALID_BUT_NOT_ARCHIVED_ENDPOINT' if states_last['archive_status']!='NEW_MATERIAL' else 'VALID_ARCHIVED_ENDPOINT'
                rows.append(row);raw_images.append(raw);endpoint_images.append(end)
                if last['compatible']: valid_images.append(end)
            panels.append(dict(platform='A30_2D',parent_id=mm[case,arm]['parent_id'],case_id=case,arm=arm,base_valid=bool(base_row['compatible']),
                raw_proposals=stage_counts(raw_images,grid_distance,different,'A30',bg),corrected_endpoints=stage_counts(endpoint_images,grid_distance,different,'A30',bg),
                full_valid_endpoints=stage_counts(valid_images,grid_distance,different,'A30',bg),historical_archive=stage_counts([c['chi'] for c in archive],grid_distance,different,'A30',bg),
                all_full_trial_states=len(rr),all_full_valid_trial_states=sum(r['compatible'] for r in rr),original_saved_best4_complex_MAE=mm[case,arm]['Best_of_4_complex_MAE'],
                original_saved_coverage4=mm[case,arm]['Coverage_at_4'],source=f'{history_path.relative_to(ROOT)}#case={case}/arm={arm}'))
    arrays.close()
    summary={}
    for arm in protocol['arms']:
        tr=[r for r in rows if r['arm']==arm];pp=[p for p in panels if p['arm']==arm]
        summary[arm]=summarize(tr)
        summary[arm]['digital_hole_attribution_counts']=dict(Counter(r['hole_attribution'] for r in tr))
        summary[arm].update(historical_archive_total=sum(p['historical_archive']['states'] for p in pp),
            no_valid_cases=sum(p['historical_archive']['states']==0 for p in pp),
            all_case_coverage4=float(np.mean([p['original_saved_coverage4'] for p in pp])),
            all_case_best4_MAE=None if any(p['original_saved_best4_complex_MAE'] is None for p in pp) else float(np.mean([p['original_saved_best4_complex_MAE'] for p in pp])),
            original_valid_structural_total=sum(p['historical_archive']['original_geometry_distinct'] for p in pp),
            multiscale_archive_total_by_threshold={str(m):sum(p['historical_archive']['multiscale_structure_distinct'][str(m)] for p in pp) for m in [.5,1.,2.]})
    return dict(parents=4,conditions=8,summary=summary,panels=panels,
                actual_historical_QP_diagnostics=saved_summary['QP_diagnostics'],
                reported_pre_dedup_offline_diagnostic='Original RESULT_SUMMARY/FAILURE_ANALYSIS records sub-1% within-neighborhood truth-MAE improvements. G did not read truth arrays or use them in replay.',
                max_protected_defect_recomputed_from_saved_arrays=max([r['protected_defect_recomputed_from_saved_arrays'] for r in rows if r.get('protected_defect_recomputed_from_saved_arrays') is not None],default=None),
                no_valid_base_conditions=sum(not x['compatible'] for x in history if x['arm']=='base'),
                protected_edit_conclusion='Recorded protected defects remain zero; cannot attribute failure to corrector filling protected cells'),rows,improvements,states


def main():
    dump(DELEG/'STRUCTURE_DIAGNOSTIC_SPEC.json',SPEC)
    status='FAILED'
    try:
        radial,rr,ri=radial_audit()
        grid,gr,gi,states=grid_audit()
        result=dict(status='COMPLETE_CACHED_ONLY',new_full_state_attempts=0,new_complex_full_solver_RHS=0,
            solver_instantiations=0,optimizer_reruns=0,truth_arrays_accessed=[],source_revision_declared='f1b674eb4893e39126701ed48f3d3dcef3098c3c',
            radial=radial,A30_2D=grid,
            first_member_dedup={'truth_free_residual_improvements_discarded':len(ri)+len(gi),
                'radial':len(ri),'A30_2D':len(gi),'interpretation':'Duplicate endpoints with lower residual than >=1 matching retained representative; not new world discoveries. Exact comparisons in CSV.',
                'all_matching_representatives_improved':sum(r['improves_all_matching_representatives'] for r in ri+gi)},
            structure_sensitivity_spec=SPEC,
            limitations=['Radial raw validity is proxy-only; no full raw->corrected validation claim outside recorded no-comp controls',
                          'E2 historical archive considers every full trajectory state, including raw step0, not only final endpoints',
                          'Structure diversity uses greedy pairwise separation, not connected components or posterior modes',
                          'Trust-region exits and full raw radial solver diagnostics unavailable',
                          'No grouped history predictor in this script; handled independently by integrator',
                          'No original source data overwritten or new physical calls'],
            evidence_files=READ)
        dump(OUT/'results/E0_TRAJECTORY_AUDIT.json',result)
        csvwrite(OUT/'results/E0_TRAJECTORY_TRANSITIONS.csv',rr+gr)
        csvwrite(OUT/'results/E0_FULLSTATE_ARCHIVE_TRANSITIONS.csv',states)
        csvwrite(OUT/'results/E0_FIRST_MEMBER_RESIDUAL_IMPROVEMENTS.csv',ri+gi)
        dump(DELEG/'ASSET_MAP.json',{'read_assets':READ,'output_root':str(OUT),'online_arrays_read_only':True,'offline_truth_arrays_read':False})
        status='COMPLETE'
        print(json.dumps({'status':status,'radial_attempts':len(rr),'A30_trajectories':len(gr),'A30_full_trial_states':len(states),'residual_improvements_discarded':len(ri)+len(gi)}))
    finally:
        dump(DELEG/'EXECUTION_RECEIPT.json',dict(status=status,controlled_execution_wall_seconds=time.perf_counter()-START,process_cpu_seconds=time.process_time()-CPU,
            new_full_state_attempts=0,new_complex_full_solver_RHS=0,GPU_seconds=0,solver_instantiations=0,optimizer_reruns=0,
            no_truth_members_accessed=True,process_id=os.getpid(),python=sys.executable,BLAS_threads={k:os.environ.get(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']},
            expected_wall_seconds_max=120,source_analysis_plan=str(PLAN.relative_to(ROOT))))

if __name__=='__main__':main()
