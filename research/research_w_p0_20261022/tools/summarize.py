"""Cache-only decision/metrics; no imported physical solver."""
from pathlib import Path
import sys,json,csv,time
import numpy as np
from scipy.ndimage import label,binary_fill_holes
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'src'))
from research_w.accounting import Book

def write_json(path,value):path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False))
def csv_write(path,rows):
    columns=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(rows)

def main():
    out=root/'results';cfg=json.loads((root/'w_config.json').read_text());book=Book(out)
    old_events=[json.loads(l) for l in book.path.read_text().splitlines()]
    if not any(r.get('label')=='preinstrumentation_allowance' for r in old_events):
        book.event(event='cost',label='preinstrumentation_allowance',charged_cpu_seconds=600.,status='CONSERVATIVE_ALLOCATION_NOT_MEASUREMENT',
                   phase='validation',measurement='Early unzip/source audit, worker text/metadata extraction and unwrapped environment failures charged ten CPU minutes conservatively; actual physical reservations unaffected')
    p0a=json.loads((out/'p0a/results.json').read_text())
    corrected=json.loads((out/'qualification_source_average/WALL_GRIDS.json').read_text())
    original=json.loads((out/'qualification/WALL_GRIDS.json').read_text())
    derivative=json.loads((out/'qualification_source_average/DERIVATIVES.json').read_text())
    reciprocity=json.loads((out/'qualification/RECIPROCITY_REPAIR.json').read_text())
    four=json.loads((out/'four_blocks/RESULTS.json').read_text())
    rows=[]
    for r in p0a['bandwidth_sweep']:
        rows.append(dict(record='P0A_1D',status='ANALYTIC_REGRESSION',**r))
    for version,rr in [('original',original),('source_average',corrected)]:
        for r in rr:
            for key,value in r['mesh'].items():
                rows.append(dict(record='P0B_GRID_QUALIFICATION',backend_revision=version,seed=r['seed'],nx=r['nx'],previous_nx=r['previous_nx'],metric=key,value=value,status='PASS' if value<=.1 else 'FAIL_NUMERICAL_QUALIFICATION'))
    pair=[]
    with np.load(out/'qualification_source_average/WALL_seed2026102201_nx192.npz') as a, np.load(out/'qualification_source_average/WALL_seed2026102202_nx192.npz') as b:
        material_distance=float(np.linalg.norm(a['epsilon_cal']-b['epsilon_cal'])*(4/192))
        for band,frequencies in cfg['frequencies'].items():
            cal1=np.stack([a[f'data_{f:.8f}_cal'] for f in frequencies]);cal2=np.stack([b[f'data_{f:.8f}_cal'] for f in frequencies])
            task1=np.stack([a[f'data_{f:.8f}_task'] for f in frequencies]);task2=np.stack([b[f'data_{f:.8f}_task'] for f in frequencies])
            r=dict(record='REFERENCE_PAIR_NOT_A_SEARCH',seed1=2026102201,seed2=2026102202,band=band,nx=192,
                   d_cal=float(np.linalg.norm(cal1-cal2)/.01),d_fixed_target_total=float(np.linalg.norm(task1-task2)/.01),
                   d_fixed_target_induced=float(np.linalg.norm((task1-cal1)-(task2-cal2))/.01),
                   target_centroid_difference=0.,material_L2_distance=material_distance,
                   status='UNQUALIFIED_REFERENCE_DIAGNOSTIC_NOT_A_WITNESS')
            pair.append(r);rows.append(r)
        topology=[]
        for seed,z in [(2026102201,a),(2026102202,b)]:
            eps=z['epsilon_cal'];mask=eps.real>1.1
            holes=binary_fill_holes(mask)&~mask
            topology.append(dict(seed=seed,support_components=int(label(mask,np.ones((3,3)))[1]),air_holes=int(label(holes)[1]),
                                 support_area=float(mask.sum()*(4/192)**2),minimum_real_epsilon=float(eps.real.min()),maximum_real_epsilon=float(eps.real.max()),maximum_imaginary_epsilon=float(eps.imag.max()),
                                 threshold='Re epsilon>1.1; foreground8/background4; descriptive only'))
    write_json(out/'REFERENCE_PAIR.json',pair);write_json(out/'MATERIAL_LEGALITY.json',topology)
    for r in four:rows.append(dict(record='FOURBLOCK_DISCRETE_IDENTITY',**r))
    for seed in cfg['truth_seeds']:
        for band in cfg['frequencies']:
            for method in cfg['search']['methods']:
                for contract in ['CALIBRATION','JOINT']:
                    rows.append(dict(record='PLANNED_SEARCH',seed=seed,band=band,method=method,contract=contract,status='NOT_RUN_UNQUALIFIED_SOLVER',d_cal=None,d_task=None,d_joint=None,target_centroid_difference=None))
    csv_write(root/'RESULTS.csv',rows)
    latest=[r for r in corrected if r['nx']==192]
    numerical_max=max(v for r in latest for v in r['mesh'].values())
    numerical_min=min(v for r in latest for v in r['mesh'].values())
    check=dict(gate0='PASS',gate1='UNQUALIFIED_SOLVER',gate2='NOT_RUN_UNQUALIFIED_SOLVER',gate3='NOT_RUN_UNQUALIFIED_SOLVER',gate4='UNRESOLVED_2D',gate5='HOLD_P1_NOT_AUTHORIZED',
               overall='P0_B_UNRESOLVED_NUMERICAL_QUALIFICATION',search_executed=False,candidate_count=0,
               grid_whitened_min=numerical_min,grid_whitened_max=numerical_max,required=.1,
               observed_finest_grid=[192,144],matched_source_average=True,
               budget_exhausted=False,reason='Pre-registered qualification ladder does not meet full-data 0.1 tolerance. No search or physical identifiability conclusion is valid.',
               P1_A='HOLD; no qualified joint ambiguity reference',P1_B='HOLD; no qualified ambiguity or same-cost progressive evidence',
               NN='NOT_RUN_NOT_AUTHORIZED',three_D='NOT_RUN_NOT_AUTHORIZED')
    write_json(out/'GATE_DECISION.json',check)
    events=[json.loads(l) for l in book.path.read_text().splitlines()]
    reserves=[r for r in events if r.get('event')=='reserve']
    counts={key:sum(r.get(key,0) for r in reserves) for key in ['rhs','factorizations']}
    counts['by_phase']={phase:{key:sum(r.get(key,0) for r in reserves if r.get('phase')==phase) for key in ['rhs','factorizations']} for phase in Book.PHASE}
    counts['by_rhs_role']={name:sum(r.get('rhs',0) for r in reserves if r.get('label')==name) for name in sorted(set(r.get('label') for r in reserves if r.get('rhs')))}
    counts['native_state_evaluations']=sum(r.get('event')=='backend_state' for r in events)
    counts['job_wall_seconds']=sum(r.get('wall_seconds',0) for r in events if r.get('event')=='cost')
    counts['job_process_cpu_seconds']=sum(r.get('process_cpu_seconds',0) for r in events if r.get('event')=='cost')
    counts['charged_cpu_seconds']=sum(r.get('charged_cpu_seconds',0) for r in events)
    counts['unprofiled_setup_conservative_charge_seconds']=600.
    counts['maximum_job_peak_rss_bytes']=max((r.get('peak_rss_bytes',0) for r in events),default=0)
    counts['full_domain_green_matrix_materialized']=False
    counts['gpu_seconds']=0;counts['search_rhs']=counts['by_phase']['search']['rhs']
    counts['endpoint_search_rhs']=counts['by_phase']['endpoint']['rhs']
    write_json(out/'COST_SUMMARY.json',counts)
    failure=[dict(type='ENVIRONMENT_DEPENDENCY',stage='initial unwrapped freeze',error='System Python lacks NumPy; no config produced or physical call',retained='No change to original assets',cost_included='preinstrumentation allowance'),
             dict(type='MEMORY_ENGINEERING',stage='wall-grid-001',peak_rss_bytes=2054520832,target_bytes=1500000000,decision='Retain original results; chunked in-place assembly parity verified; corrected job peak 572342272 bytes',cost_included=True)]
    for r in events:
        if r.get('event')=='cost' and r.get('status')=='FAILED':failure.append(dict(type='FAILED_JOB',**r))
    for r in latest:
        failure.append(dict(type='SOLVER_QUALIFICATION_FAILURE',seed=r['seed'],nx=r['nx'],mesh=r['mesh'],tolerance=.1,decision='No candidate search, no GO, no uniqueness conclusion'))
    failure_path=out/'FAILURE_LEDGER.jsonl'
    existing=set(failure_path.read_text().splitlines()) if failure_path.exists() else set()
    with failure_path.open('a') as stream:
        for r in failure:
            line=json.dumps(r,ensure_ascii=False)
            if line not in existing:stream.write(line+'\n');existing.add(line)
    print(json.dumps(dict(gates=check,cost=counts),ensure_ascii=False))

if __name__=='__main__':main()
