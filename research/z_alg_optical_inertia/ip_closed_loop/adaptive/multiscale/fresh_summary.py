"""Cost and parity ledger for the fresh matched adaptive/exhaustive replays."""
import json
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
MAPS=HERE.parents[1]/'maps'
FRESH=HERE/'fresh'


def read(path):
    return json.loads(path.read_text())


def object_record(path):
    r=read(path)
    exact=r['fresh_singleton_replay']
    costs=exact['costs']
    adaptive=costs['adaptive_fresh_leaf_wall_s']
    exhaustive=costs['exhaustive_fresh_leaf_wall_s']
    parents=r['parent_stage_total']['wall_s']
    bank=costs['one_common_state_bank_wall_s']
    clip=costs['clip_partition_wall_s']
    adaptive_state=exact['adaptive_alone_chargeable_state_wall_s']
    model=MAPS/'results'/f"{r['oid']}__b0.03125__score_map.npz"
    with np.load(model,allow_pickle=False) as d:
        archive_directional_additions=int(np.count_nonzero((d['neumann']>0)&(d['sector']<=0)))
    return dict(oid=r['oid'],side=r['side'],tested_centres=r['total_centres'],
                parent_exclusions=r['pruned_by_positive_parent'],
                parent_witnesses=r['strict_parent_certificates'],
                strict_indefinite_parents=sum(x.get('strict_negative_inertia_parents',0) for x in r['counts_by_tile'].values()),
                fresh_singleton_inputs=r['fallback_singletons'],
                fresh_adaptive_total_exclusions=r['final_exclusions'],
                fresh_exhaustive_exclusions=exact['fresh_exhaustive_exclusions'],
                adaptive_and_exhaustive_identical=exact['same_as_fresh_exhaustive'],
                exhaustive_and_archive_identical=exact['fresh_exhaustive_parity_with_archive'],
                directional_additions_in_archive=archive_directional_additions,
                retained_physical_witness=r['old_physical_witness'],
                parent_stage_wall_s=parents,
                parent_factor_wall_s=r['new_parent_costs'].get('parent_local_wall_s',0.),
                shared_clip_partition_wall_s=clip,
                shared_singleton_state_proof_wall_s=bank,
                adaptive_leaf_online_wall_s=adaptive,exhaustive_leaf_online_wall_s=exhaustive,
                adaptive_charged_singleton_state_wall_s=adaptive_state,
                adaptive_branch_and_leaf_online_wall_s=parents+adaptive,
                exhaustive_online_wall_s=exhaustive,
                adaptive_charged_stage_wall_s=parents+adaptive_state+adaptive+clip,
                exhaustive_charged_stage_wall_s=bank+exhaustive+clip,
                both_routes_measured_wall_s=r['wall_s'],
                clip_types=exact['clip_types'],
                singleton_state_attempts=exact['common_state_bank_attempts'],
                centre_weight_evaluations=exact['route_centre_weight_evaluations'])


def main():
    files=sorted(FRESH.glob('ipcl-*.json'))
    if len(files)!=6:
        raise ValueError('The matched replay must contain all six independent observations.')
    rows={d['oid'].split('__')[1]:d for d in map(object_record,files)}
    for d in rows.values():
        if not d['adaptive_and_exhaustive_identical'] or not d['exhaustive_and_archive_identical']:
            print(d['oid'],'has a decision discrepancy. See raw margins.')
    summaries={}
    for side in (32,96):
        records=[r for r in rows.values() if r['side']==side]
        def add(key):
            return sum(r[key] for r in records)
        summaries[str(side)]=dict(total_centres=add('tested_centres'),
                                  centres_pruned_before_leaf=add('parent_exclusions'),
                                  parent_witnesses=add('parent_witnesses'),
                                  strict_indefinite_parent_witnesses=add('strict_indefinite_parents'),
                                  remaining_fresh_singleton_inputs=add('fresh_singleton_inputs'),
                                  final_adaptive_exclusions=add('fresh_adaptive_total_exclusions'),
                                  adaptive_singleton_checks=sum(r['centre_weight_evaluations']['adaptive'] for r in records),
                                  exhaustive_singleton_checks=sum(r['centre_weight_evaluations']['exhaustive'] for r in records),
                                  adaptive_online_wall_s=add('adaptive_branch_and_leaf_online_wall_s'),
                                  exhaustive_online_wall_s=add('exhaustive_online_wall_s'),
                                  common_singleton_state_wall_s=add('shared_singleton_state_proof_wall_s'),
                                  adaptive_charged_stage_wall_s=add('adaptive_charged_stage_wall_s'),
                                  exhaustive_charged_stage_wall_s=add('exhaustive_charged_stage_wall_s'))
    (FRESH/'summary.json').write_text(json.dumps(dict(summary=summaries,rows=rows,
            timing_scope='Recorded interleaved route-specific online times, plus independently metered common geometry and cached factor construction. Charged times are stage sums, accounting a state factor to each route that needed it. Neither route reads archived singleton decisions before solving. Separate archived directional labels are a comparison only.'),indent=2,allow_nan=False)+'\n')
    print(summaries)


if __name__=='__main__':main()
