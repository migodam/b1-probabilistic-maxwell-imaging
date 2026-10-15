"""Frozen-pool archive-only controlled replay; zero physics and no truth access."""
from pathlib import Path
import csv
import json
import time
import os
import numpy as np
import e0_trajectory_audit as audit

START,CPU=time.perf_counter(),time.process_time()
OUT,ROOT,A=audit.OUT,audit.ROOT,audit.A
DELEG=audit.DELEG


def packing_ok(items,distance):
    return all(distance(a['field'],b['field'])>=.03 for i,a in enumerate(items) for b in items[i+1:])


def insert(archive,item,distance,*,quality=False):
    matching=[(j,distance(item['field'],a['field'])) for j,a in enumerate(archive) if distance(item['field'],a['field'])<.03]
    if not matching:
        archive.append(item)
        return 'NEW_TO_ACTIVE_POOL',None
    if not quality:return 'DUPLICATE_KEPT_FIRST',None
    nearest=min(matching,key=lambda pair:pair[1])[0]
    old=archive[nearest]
    if item['residual']>=old['residual']-1e-12:return 'DUPLICATE_NO_RESIDUAL_GAIN',None
    if not all(distance(item['field'],a['field'])>=.03 for j,a in enumerate(archive) if j!=nearest):return 'REPLACEMENT_REJECTED_PACKING',None
    archive[nearest]=item
    return 'RESIDUAL_IMPROVEMENT_REPLACEMENT',dict(old_id=old['id'],new_id=item['id'],old_residual=old['residual'],new_residual=item['residual'],residual_reduction=old['residual']-item['residual'])


def replay(platform,parent_id,case_id,arm,items,distance,different,bg=None):
    first,quality,changes=[],[],[]
    immutable=[x for x in items if x['valid']]
    for item in immutable:
        insert(first,item,distance)
        status,record=insert(quality,item,distance,quality=True)
        if record:changes.append(record)
    if not packing_ok(first,distance) or not packing_ok(quality,distance):raise AssertionError('Archive pairwise packing failed')
    first_fields=[x['field'] for x in first];q_fields=[x['field'] for x in quality]
    original_lost=[x['id'] for x in first if all(different(x['field'],q['field']) for q in quality)]
    ffeatures=[audit.geom_feature(x['field'],'radial' if platform=='radial' else 'A30',bg) for x in first]
    qfeatures=[audit.geom_feature(x['field'],'radial' if platform=='radial' else 'A30',bg) for x in quality]
    scale=.025 if platform=='radial' else .08
    multiscale_lost=[first[i]['id'] for i,f in enumerate(ffeatures) if all(audit.geom_distance(f,q,'radial' if platform=='radial' else 'A30')>=scale for q in qfeatures)]
    r=dict(platform=platform,parent_id=parent_id,case_id=case_id,arm=arm,immutable_valid_endpoint_count=len(immutable),
        first_active_count=len(first),quality_active_count=len(quality),replacements=len(changes),
        first_min_full_residual=min([x['residual'] for x in first],default=None),quality_min_full_residual=min([x['residual'] for x in quality],default=None),
        first_mean_full_residual=float(np.mean([x['residual'] for x in first])) if first else None,quality_mean_full_residual=float(np.mean([x['residual'] for x in quality])) if quality else None,
        paired_replacement_residual_reduction_sum=sum(c['residual_reduction'] for c in changes),
        original_structure_count_first=audit.geom_original_count(first_fields,different),original_structure_count_quality=audit.geom_original_count(q_fields,different),
        multiscale_structure_count_first=audit.pack_count(ffeatures,lambda a,b:audit.geom_distance(a,b,'radial' if platform=='radial' else 'A30'),scale),
        multiscale_structure_count_quality=audit.pack_count(qfeatures,lambda a,b:audit.geom_distance(a,b,'radial' if platform=='radial' else 'A30'),scale),
        original_structure_coverage_lost_count=len(original_lost),multiscale_structure_coverage_lost_count=len(multiscale_lost),
        original_structure_lost_ids=original_lost,multiscale_structure_lost_ids=multiscale_lost,
        first_ids=[x['id'] for x in first],quality_ids=[x['id'] for x in quality],replacement_records=changes,
        first_base_retained=any(x['id']=='base' for x in first),quality_base_retained=any(x['id']=='base' for x in quality),
        new_physics=0,new_RHS=0,new_world_discovery_credit=0)
    return r


def main():
    status='FAILED'
    try:
        rows=[]
        campaign=json.loads((A/'results/CONFIRM_CAMPAIGN.json').read_text())
        for parent in campaign:
            for method in parent['methods']:
                source='research/research_a_20261015/results/CONFIRM_CAMPAIGN.json'
                items=[dict(id='base',field=np.array(parent['base']['p']),valid=parent['base']['valid'],residual=parent['base']['residual'])]
                items += [dict(id=f'step_{i}',field=np.array(x['p']),valid=x['valid'],residual=x['residual']) for i,x in enumerate(method['attempts'])]
                r=replay('radial',parent['id'],parent['id'],method['arm'],items,audit.radial_distance,audit.radial_geometry)
                r['source']=f'{source}#parent={parent["id"]}/arm={method["arm"]}'
                if r['first_active_count']!=len(method['archive']):raise AssertionError('Historical radial archive not reproduced')
                rows.append(r)
        history=json.loads((A/'results/E2/SEARCH_HISTORY.json').read_text())
        protocol=json.loads((A/'results/E2/PROTOCOL.json').read_text())
        with np.load(A/'results/E2/MATERIAL_AND_PREDICTION_ARRAYS.npz',allow_pickle=False) as arrays:
            for case in protocol['cases']:
                base_row=next(x for x in history if x['case_id']==case and x['arm']=='base')
                bg=arrays[case+'__background']
                for arm in protocol['arms']:
                    items=[dict(id='base',field=arrays[case+'__base'],valid=base_row['compatible'],residual=base_row['residual'])]
                    items += [dict(id=x['array_prefix'],field=arrays[x['array_prefix']+'__chi'],valid=x['compatible'],residual=x['residual']) for x in history if x['case_id']==case and x['arm']==arm]
                    different=lambda a,b: np.count_nonzero((abs(a-bg)>.18)!=(abs(b-bg)>.18))>=6
                    r=replay('A30_2D',case.split('_')[0],case,arm,items,audit.grid_distance,different,bg)
                    r['source']=f'research/research_a_20261015/results/E2/SEARCH_HISTORY.json#case={case}/arm={arm}'
                    rows.append(r)
        summary={}
        for platform in ['radial','A30_2D']:
            pp=[r for r in rows if r['platform']==platform]
            summary[platform]={arm:{'episodes':len([r for r in pp if r['arm']==arm]),
                **{k:sum(r[k] for r in pp if r['arm']==arm) for k in ['first_active_count','quality_active_count','replacements','original_structure_count_first','original_structure_count_quality','multiscale_structure_count_first','multiscale_structure_count_quality','original_structure_coverage_lost_count','multiscale_structure_coverage_lost_count','paired_replacement_residual_reduction_sum']}}
                for arm in dict.fromkeys(r['arm'] for r in pp)}
        audit.dump(OUT/'results/E0_ARCHIVE_REPAIR_REPLAY.json',dict(status='COMPLETE_CACHED_POOL_REPLAY',summary=summary,episodes=rows,
            plan=json.loads((DELEG/'ARCHIVE_REPAIR_PLAN.json').read_text()),no_new_world_discovery=True,truth_read=False,
            interpretation='Improves retained representative fit within fixed saved pools. Any lost structural coverage is reported explicitly. This does not reproduce counterfactual online adaptive proposals.'))
        audit.csvwrite(OUT/'results/E0_ARCHIVE_REPAIR_REPLAY.csv',rows)
        status='COMPLETE'
        print(json.dumps({'status':status,'episodes':len(rows),'replacements':sum(r['replacements'] for r in rows),'old_structure_neighborhoods_lost':sum(r['original_structure_coverage_lost_count'] for r in rows)}))
    finally:
        audit.dump(DELEG/'ARCHIVE_REPAIR_EXECUTION_RECEIPT.json',dict(status=status,controlled_execution_wall_seconds=time.perf_counter()-START,process_cpu_seconds=time.process_time()-CPU,
            new_full_state_attempts=0,new_complex_full_solver_RHS=0,GPU_seconds=0,solver_instantiations=0,optimizer_reruns=0,truth_members_accessed=[],expected_wall_seconds_max=30,
            BLAS_threads={k:os.environ.get(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']}))

if __name__=='__main__':main()
