"""Join saved scalar reconstruction endpoints to the public certificate evidence."""

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np

HERE = Path(__file__).resolve().parent
VALIDATION = HERE.parent
WORKSPACE = HERE.parents[3]
DEFAULT_SRAP = WORKSPACE/'research/srap_cie_array_study_20261121'
DEFAULT_R14 = (WORKSPACE/'research/delegated/max_struct_r1/r14_worktree/research'
               /'research_max_struct_r14')
METHODS_SRAP = {'dual':'Dual-CIE','srap':'SRAP-CIE'}
METHODS_R14 = {
    'pair00':'00','pair01':'01 / original global A+-Arc','pair10':'10','pair11':'11',
    'pair01_to11':'01→11','cie':'CIE','dual_cie':'Dual-CIE','global_dual':'global-Dual',
    'g19':'G19','dr':'DR','direct':'Direct L-BFGS-B',
}


def json_read(path):
    return json.loads(path.read_text())


def by_csv(path):
    with path.open(newline='',encoding='utf-8') as f:
        yield from csv.DictReader(f)


def floats(x):
    v=float(x)
    if not np.isfinite(v):
        raise ValueError('A saved endpoint scalar is not finite.')
    return v


def max_caps(valid):
    f = valid/'paired_results.csv'
    if not f.exists():
        return {}
    best=defaultdict(float)
    for r in by_csv(f):
        if (r['arm'] in ('pure_joint','optical_joint','definite_joint','indefinite_joint')
                and r['passes_with_confirmed_coverage']=='True'):
            best[r['case']]=max(best[r['case']],float(r['cap']))
    return dict(best)


def common_row(cohort,name,code,method,view,budget,unit,source_grid,generation_grid,
               certificate_grid,peak,residual,norm,radius,cap,status,archived_path,archive_selection,paid_at):
    same_grid=source_grid==generation_grid==certificate_grid
    residual_abs=residual*norm
    peak_relation=('at_or_below_cap' if cap is not None and peak<=cap else
                   'above_cap' if cap is not None else 'certificate_pending')
    return dict(cohort=cohort,observation_id=name,method_code=code,method=method,
                budget=budget,budget_unit=unit,view=view,history_status=status,
                historical_selection=archive_selection,endpoint_available_under_historical_clock=paid_at,
                endpoint_grid_side=source_grid,generation_grid_side=generation_grid,
                certificate_grid_side=certificate_grid,same_native_contract=same_grid,
                peak_chi=peak,relative_residual_Frobenius=residual,
                numerical_absolute_residual_Frobenius=residual_abs,
                certificate_radius_Frobenius=radius,
                same_grid_residual_over_certificate_radius=residual_abs/radius if same_grid else None,
                same_grid_numerically_within_certificate_radius=bool(residual_abs<=radius) if same_grid else None,
                certified_joint_cap=cap,numeric_peak_relation=peak_relation,
                archived_prediction_present=archived_path.is_file(),
                archived_prediction=(archived_path.name if cohort=='SRAP' else
                                     f'{archived_path.parent.name}/{archived_path.name}'))


def srap(srap):
    valid=VALIDATION/'srap'
    proof=json_read(valid/'online.json')
    refs={r['id']:r for r in proof['objects']}
    with np.load(valid/'public_observations.npz',allow_pickle=False) as a:
        y={str(cid):(value,float(d)) for cid,value,d in zip(a['ids'],a['Y'],a['original_delta'],strict=True)}
    for cid,value in y.items():
        path=srap/'data/public'/f'{cid}.npz'
        # The release's data/ folder can be absent while the validation pack and
        # historical prediction receipts remain available.
        if not path.is_file():
            continue
        with np.load(path,allow_pickle=False) as f:
            if not (np.array_equal(f['Y'],value[0]) and float(f['delta'])==value[1]):
                raise ValueError(f'SRAP observation mismatch: {cid}.')
            meta=json.loads(f['case_json'].item())
            if (meta['fine'],meta['side'])!=(96,32):
                raise ValueError(f'Unexpected SRAP grid dimensions: {cid}.')
    selected={
        (r['case'],r['method'],int(r['budget']),r['view']):r
        for r in by_csv(srap/'results/confirmation_metrics.csv')
        if r['case'] in refs and r['method'] in METHODS_SRAP
        and r['view'] in ('raw','light') and r['budget'] in ('15','30')
    }
    results,missing=[],[]
    for cid,meta in refs.items():
        for code,method in METHODS_SRAP.items():
            file=srap/'results/predictions'/f'{cid}__{code}.npz'
            for budget in (15,30):
                for view in ('raw','light'):
                    key=cid,code,budget,view
                    row=selected.get(key)
                    if not row or not (row['peak_contrast'] and row['residual']):
                        missing.append(key)
                        continue
                    results.append(common_row('SRAP',cid,code,method,view,budget,'wall_seconds',32,96,96,
                                              floats(row['peak_contrast']),floats(row['residual']),
                                              np.linalg.norm(y[cid][0]),meta['padded_radius'],.5,
                                              row['status'],file,row['selected_stage'],floats(row['available_wall'])))
    return results,missing


def uniform_srap_audit(srap, endpoints):
    # Both archives use the registered disk-pulse cell-centre formula on the
    # same width. Array order is row-major (y before x).
    proof=json_read(VALIDATION/'srap/online.json')
    offline=json_read(VALIDATION/'srap/offline.json')
    if not all(offline['geometry_replay_binary64_equal'].values()):
        raise ValueError('The 96² proof grid has not been checked against the generator.')
    admitted={r['id'] for r in offline['rows'] if r['inside_padded_ball']}
    certified={r['id']:r for arm in proof['profiles'] if arm['label']=='global_real_050' for r in arm['rows']}
    rows=[]
    for e in endpoints:
        if e['numeric_peak_relation']!='at_or_below_cap':
            continue
        cid=e['observation_id']
        t,c=e['budget'],e['view']
        source=srap/'results/predictions'/e['archived_prediction']
        with np.load(source,allow_pickle=False) as f:
            material=f[f'm_{t}_{c}']
        with np.load(srap/'data/public'/f'{cid}.npz',allow_pickle=False) as f:
            coarse=f['centers'].reshape(32,32,2)
            width=float(json.loads(f['case_json'].item())['width'])
        if not (material.shape==(1024,) and np.isrealobj(material)
                and np.all(np.isfinite(material)) and np.min(material)>=0.):
            raise ValueError(f'The stored SRAP endpoint is not on the real nonnegative grid: {cid}/{t}/{c}.')
        if np.max(material)!=e['peak_chi'] or width!=2.:
            raise ValueError(f'SRAP scalar/physical grid disagrees with the retained output: {cid}/{t}/{c}.')
        fine_axis=(np.arange(96)+.5)*(width/96)-width/2
        yc,xc=np.meshgrid(fine_axis,fine_axis,indexing='ij')
        grouped=np.stack((xc,yc),axis=-1).reshape(32,3,32,3,2).mean(axis=(1,3))
        deviation=float(np.max(np.abs(coarse-grouped)))
        # Mathematical nesting is exact from the width and integer cell ratio.
        # This check also verifies the archived physical coordinate ordering.
        if deviation>1.e-14:
            raise ValueError(f'SRAP cell-centre nesting failed: {cid}.')
        lifted=np.repeat(np.repeat(material.reshape(32,32),3,axis=0),3,axis=1)
        if np.max(lifted)!=e['peak_chi'] or certified[cid]['joint_status']!='excluded' or cid not in admitted:
            raise ValueError(f'Insufficient evidence for the uniform lift of {cid}/{t}/{c}.')
        rows.append(dict(observation_id=cid,method=e['method'],budget_seconds=t,view=c,
                         coarse_min_chi=float(np.min(material)),coarse_max_chi=float(np.max(material)),
                         native_uniform_lift_min_chi=float(np.min(lifted)),
                         native_uniform_lift_max_chi=float(np.max(lifted)),
                         coarse_side=32,lifted_side=96,cells_per_coarse_pixel_axis=3,
                         max_cell_centre_nesting_error_metres=deviation,
                         historical_coarse_relative_residual=e['relative_residual_Frobenius'],
                         native_residual_lower_from_cap_proof=certified[cid]['joint_distance_lower'],
                         padded_radius=certified[cid]['padded_radius'],
                         uniform_lift_strictly_excluded=True))
    return rows


def r14(root):
    val=VALIDATION/'max_struct'
    caps=max_caps(val)
    with np.load(val/'public_inputs.npz',allow_pickle=False) as f:
        y={str(cid):(value,float(d)) for cid,value,d in zip(f['cases'],f['Y'],f['delta'],strict=True)}
        geo={k:f[k] for k in ('G','E0','S')}
        inverse_side=tuple(f['grid_shape'])
    if set(caps)!=set(y):
        raise ValueError('Complete same-grid strict certificates for all seven R14 materials first.')
    if inverse_side!=(24,24):
        raise ValueError('The R14 validation must use the matched 24² grid.')
    for cid,(data,d) in y.items():
        path=root/'data_public/independent24'/f'{cid}.npz'
        with np.load(path,allow_pickle=False) as f:
            if not (np.array_equal(f['Y'],data) and float(f['delta'])==d
                    and (int(f['forward_side']),int(f['inverse_side']))==(24,24)
                    and all(np.array_equal(f[k],a) for k,a in geo.items())):
                raise ValueError(f'R14 generator/public/validation contract differs: {cid}.')
    original={
        (r['case_id'],r['solver']):r for r in by_csv(root/'analysis/independent24_v1/ALL_VIEWS.csv')
        if r['view']=='primary' and r['case_id'] in y and r['solver'] in METHODS_R14
    }
    results,missing=[],[]
    for cid,(data,d) in y.items():
        for code,method in METHODS_R14.items():
            row=original.get((cid,code))
            if not row or row['available']!='True' or not row['material_max'] or not row['relative_residual']:
                missing.append((cid,code,'primary'))
                continue
            f=root/'RESULTS/independent24_v1'/row['task_id']/'primary.npz'
            r=common_row('MAX-STRUCT R14',cid,code,method,'primary',120,'whole_process_CPU_seconds',
                         24,24,24,floats(row['material_max']),floats(row['relative_residual']),
                         np.linalg.norm(data),d,caps.get(cid),row['run_status'],f,
                         row['source_view'],floats(row['CPU_seconds']))
            if row['residual_over_delta']:
                historical_ratio=floats(row['residual_over_delta'])
                if not np.isclose(historical_ratio,r['same_grid_residual_over_certificate_radius'],
                                  rtol=2e-12,atol=1e-14):
                    raise ValueError(f'R14 residual normalization mismatch at {cid}/{code}.')
            results.append(r)
    return results,missing


def describe(rows,missing,expected):
    counts=defaultdict(list)
    by_object=defaultdict(list)
    for r in rows:
        counts[(r['method_code'],r['budget'],r['view'])].append(r)
        by_object[r['observation_id']].append(r)
    native=all(row['same_native_contract'] for row in rows)
    methods=[]
    for (code,budget,view),records in counts.items():
        peaks=[r['peak_chi'] for r in records]
        res=[r['relative_residual_Frobenius'] for r in records]
        methods.append(dict(method_code=code,method=records[0]['method'],budget=budget,view=view,
                            materials=len({r['observation_id'] for r in records}),endpoints=len(records),
                            median_peak=float(np.median(peaks)),peak_span=[min(peaks),max(peaks)],
                            median_relative_residual=float(np.median(res)),
                            min_relative_residual=min(res),max_relative_residual=max(res),
                            numeric_peak_at_or_below_cap=sum(r['numeric_peak_relation']=='at_or_below_cap'
                                                             for r in records),
                            same_grid_numerically_in_ball=sum(r['same_grid_numerically_within_certificate_radius'] is True
                                                              for r in records) if native else None))
    expected_labels=Counter(r['numeric_peak_relation'] for r in rows)
    below=[r for r in rows if r['numeric_peak_relation']=='at_or_below_cap']
    categories={
        'materials':len(by_object),'expected_endpoints':expected,'endpoints':len(rows),
        'missing_scalar_endpoints':missing,
        'all_endpoints_use_certificate_grid':native,
        'archived_prediction_files':len({r['archived_prediction'] for r in rows
                                         if r['archived_prediction_present']}),
        'archived_prediction_missing_endpoints':sum(not r['archived_prediction_present'] for r in rows),
        'historical_views_status':dict(Counter(r['history_status'] for r in rows)),
        'peak_span':[min(r['peak_chi'] for r in rows),max(r['peak_chi'] for r in rows)],
        'numeric_peak_relation_rows':dict(expected_labels),
        'materials_with_any_endpoint_at_or_below_cap':len({r['observation_id'] for r in below}),
        'same_grid_inside_radius_numerical_endpoints':sum(r['same_grid_numerically_within_certificate_radius'] is True
                                                         for r in rows) if native else None,
        'same_grid_below_cap_and_inside_radius_numerical_endpoints':sum(
            r['same_grid_numerically_within_certificate_radius'] is True for r in below) if native else None,
        'methods':methods,
        'by_material':[
            dict(observation_id=cid,endpoints=len(group),certified_joint_cap=group[0]['certified_joint_cap'],
                 historical_peak_span=[min(r['peak_chi'] for r in group),max(r['peak_chi'] for r in group)],
                 numerical_endpoints_within_same_grid_radius=sum(
                     r['same_grid_numerically_within_certificate_radius'] is True for r in group) if native else None)
            for cid,group in by_object.items()],
    }
    return categories


def main(args):
    start,cpu=perf_counter(),process_time()
    a,miss_a=srap(args.srap)
    b,miss_b=r14(args.r14)
    out=a+b
    lift=uniform_srap_audit(args.srap,a)
    with (HERE/'srap_uniform_lift_audit.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=lift[0].keys())
        w.writeheader()
        w.writerows(lift)
    with (HERE/'endpoints.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=out[0].keys())
        w.writeheader()
        w.writerows(out)
    summary={
        'scope':'Historical scalar peak and residual columns only; no historical solver execution or GT',
        'srap_scalar_source':'SRAP results/confirmation_metrics.csv, peak_contrast and residual',
        'max_struct_scalar_source':'MAX-STRUCT R14 analysis/independent24_v1/ALL_VIEWS.csv, primary rows and no truth fields',
        'srap':describe(a,miss_a,176),
        'max_struct':describe(b,miss_b,77),
        'srap_uniform_lift':{'excluded_endpoints':len(lift),'materials':len({r['observation_id'] for r in lift}),
                             'maximum_cell_centre_nesting_error_metres':max(
                                 r['max_cell_centre_nesting_error_metres'] for r in lift)},
        'wall_seconds':perf_counter()-start,
        'cpu_seconds':process_time()-cpu,
    }
    (HERE/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--srap',type=Path,default=DEFAULT_SRAP)
    parser.add_argument('--r14',type=Path,default=DEFAULT_R14)
    main(parser.parse_args())
