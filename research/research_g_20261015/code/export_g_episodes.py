"""Export completed G endpoint assets and descriptive, truth-free statistics.

No physical solver/corrector or scorer import. Run only after the physical job's
JOB_END and final RESULTS.json; existing fields/receipts remain authoritative.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from g_archive import Archive, Candidate
from g_budget import GBudget, dump
from g_contract import validate_episode, VERSION, TARGET, PACKING, GAIN

ROOT = Path(__file__).resolve().parents[1]
MASS = dict(posterior_weight=None, represented_mass=None, unmodeled_mass=None)


def ref(path, prefix, component):
    return dict(path=path, array_key=prefix+'__'+component)


def diagnostic(value, key):
    number = value.get(key)
    return dict(value=number, reason=None if number is not None else 'Not exposed by original validation receipt: '+key)


def clean_online(value):
    # The original action records an unused explicit null label. Remove only
    # that null key; unexpected non-null labels must fail the contract.
    if isinstance(value, dict):
        return {key: clean_online(item) for key, item in value.items()
                if not (key == 'topology_label' and item is None)}
    if isinstance(value, list):
        return [clean_online(item) for item in value]
    return value


def episode_sidecar(episode):
    path = episode['array_path']
    validation = episode.get('validation', {})
    return dict(version=VERSION, episode_id=episode['endpoint_id'],
        native_package_refs=['results/E1_LEGACY/RESULTS.json', 'E1_PROTOCOL_FROZEN.json'],
        target=dict(id=TARGET, packing_version=PACKING, gain_version=GAIN, gain=1),
        observations=[dict(observation_id=episode['case_id'], data_path=dict(path=path, array_key='__observation'))],
        candidate=dict(candidate_id=episode['endpoint_id'], raw_field_path=ref(path,'raw','chi'),
                       corrected_field_path=ref(path,'correction2','chi'),
                       verified_field_path=ref(path,'validation','chi'),
                       prediction_path=ref(path,'validation','prediction')),
        validation=dict(solver_diagnostics=diagnostic(validation,'solver_backward_residual'),
                        constitutive_diagnostics=diagnostic(validation,'constitutive_relative_defect'),
                        full_residual=diagnostic(validation,'full_residual'),
                        material_legal=validation.get('material_legal',False),
                        full_validation_receipt='results/E1_LEGACY/trajectories/'+episode['endpoint_id']+'.json'),
        structure={key:episode.get(key) for key in ('new_structural_world','new_material_world',
                        'combined_new_world_reward','raw_structural_edit','retained_structural_edit',
                        'correction_collapse_to_prefix','retention')},
        history=clean_online(dict(round_id=episode['round_id'], case_id=episode['case_id'],
                    method=episode['method'], trajectory=episode['trajectory'], accepted=bool(episode.get('valid')),
                    action=episode.get('action'), prefix=episode.get('steps',[]),
                    validation_receipt_id=episode.get('validation_receipt'),
                    offer_status=episode.get('offer_status'), replaced_id=episode.get('replaced_id'))),
        probability=MASS.copy(), cost=dict(episode.get('cost_delta',{}),wall_seconds=episode.get('wall_seconds'),GPU=0))


def material_entry(item, case, validation=None, sidecar=None):
    prefix=item['array_prefix'];path=item['array_path']
    return dict(candidate_id=item['candidate_id'], case_id=case, is_base=bool(item.get('is_base',False)),
        full_residual=item.get('full_residual', (validation or {}).get('full_residual')),
        validation_receipt_id=item['validation_receipt'], sidecar_path=sidecar,
        arrays={key:ref(path,prefix,key) for key in ('chi','current','field','prediction')},
        probability=MASS.copy())


def verify_material_assets(item, root):
    for asset in item['arrays'].values():
        path=(root/asset['path']).resolve()
        if not path.is_relative_to(root.resolve()):raise ValueError('MATERIAL_PATH_OUTSIDE_ROOT')
        with np.load(path,allow_pickle=False) as arrays:
            if asset['array_key'] not in arrays.files:raise ValueError('MATERIAL_ASSET_KEY_MISSING:'+asset['array_key'])


def summarize_arm(method):
    episodes=sorted(method['episodes'],key=lambda row:row['trajectory'])
    rows=[]
    for count in range(1,len(episodes)+1):
        prefix=episodes[:count]; latest=prefix[-1].get('archive_after',{})
        rows.append(dict(round_id=method['round_id'],case_id=method['case_id'],method=method['method'],prefix_t=count,
            attempted_states=sum(item.get('cost_delta',{}).get('states',0) for item in prefix),
            complex_RHS=sum(item.get('cost_delta',{}).get('rhs',0) for item in prefix),
            actual_wall_seconds=sum(item.get('wall_seconds',0.) for item in prefix),
            raw_valid=sum(bool(item.get('steps') and item['steps'][0].get('compatible')) for item in prefix),
            endpoint_valid=sum(bool(item.get('valid')) for item in prefix),
            new_material_events=sum(bool(item.get('new_material_world')) for item in prefix),
            new_structure_events=sum(bool(item.get('new_structural_world')) for item in prefix),
            combined_material_OR_structure_reward=sum(item.get('combined_new_world_reward',0) for item in prefix),
            active_representatives=latest.get('active_combined_representatives'),
            material_representatives=latest.get('active_material_representatives'),
            structural_representatives=latest.get('active_structural_representatives'),
            raw_structural_proposals=sum(bool(item.get('raw_structural_edit')) for item in prefix),
            retained_structural_edits=sum(bool(item.get('retained_structural_edit')) for item in prefix),
            retained_from_raw_structural_proposals=sum(bool(item.get('raw_structural_edit') and item.get('retained_structural_edit')) for item in prefix),
            correction_collapses=sum(bool(item.get('correction_collapse_to_prefix')) for item in prefix),
            valid_duplicate_events=sum(item.get('offer_status')=='VALID_DUPLICATE_RETAINED_IN_LEDGER' for item in prefix),
            quality_replacements=sum(item.get('offer_status')=='REPLACED_LOWER_FULL_RESIDUAL' for item in prefix),
            mean_signed_projection_retention=float(np.mean([item['retention'] for item in prefix if item.get('retention') is not None])) if any(item.get('retention') is not None for item in prefix) else None))
    return rows


def sensitivity(method, root, threshold):
    # Repack only the immutable, already validated success pool. This is a
    # descriptive replay, never a counterfactual adaptive search or retuning.
    successes=method['success_ledger']
    if not successes:return None
    with np.load(root/successes[0]['array_path'],allow_pickle=False) as arrays:
        background=arrays['__background'].copy()
    archive=Archive(background,data_threshold=23.785449004057032,
                    structure_threshold=threshold,replacement_mode=method['round_id'])
    for item in successes:
        with np.load(root/item['array_path'],allow_pickle=False) as arrays:
            archive.offer(Candidate(item['candidate_id'],arrays[item['array_prefix']+'__chi'],
                arrays[item['array_prefix']+'__prediction'],item['full_residual'],True,True,
                item['validation_receipt'],is_base=item['is_base']))
    return dict(round_id=method['round_id'],case_id=method['case_id'],method=method['method'],
                deltaS=threshold,minimum_XOR_cells=6,scope='DESCRIPTIVE_SAME_VERIFIED_POOL_REPACK_ONLY',
                **archive.summary())


def pooled_packings(candidates, root):
    """Each observation has one ordered pool across all methods and rounds."""
    output=[]
    cases=sorted({item['case_id'] for item in candidates.values()})
    for case in cases:
        items=sorted((item for item in candidates.values() if item['case_id']==case),
                     key=lambda item:(not item['is_base'],item['candidate_id']))
        base=items[0]
        with np.load(root/base['arrays']['chi']['path'],allow_pickle=False) as arrays:
            background=arrays['__background'].copy()
        archive=Archive(background,data_threshold=23.785449004057032,replacement_mode='first')
        full=[]
        for item in items:
            with np.load(root/item['arrays']['chi']['path'],allow_pickle=False) as arrays:
                candidate=Candidate(item['candidate_id'],arrays[item['arrays']['chi']['array_key']],
                    arrays[item['arrays']['prediction']['array_key']],item['full_residual'],True,True,
                    item['validation_receipt_id'],is_base=item['is_base'])
            archive.offer(candidate);full.append(candidate)
        representatives={}
        for kind,criterion in [('material','material_distinct'),('structure','structural_distinct')]:
            selected=[]
            for candidate in full:
                if all(archive.pair(candidate,old)[criterion] for old in selected):selected.append(candidate)
            representatives[kind]=[item.candidate_id for item in selected]
        combined=[item.candidate_id for item in archive.active_archive]
        output.append(dict(case_id=case,observation_id=case,full_valid_unique_endpoint_ids=len(items),
            base_ids=[item['candidate_id'] for item in items if item['is_base']],
            ordering='Base first, then endpoint ID across both rounds/all methods; first-member fixed greedy pairwise packing',
            material_representative_ids=representatives['material'],
            structural_representative_ids=representatives['structure'],
            combined_representative_ids=combined,
            material_base_plus_unique=len(representatives['material']),
            material_additional_unique=max(0,len(representatives['material'])-1),
            structure_base_plus_unique=len(representatives['structure']),
            structure_additional_unique=max(0,len(representatives['structure'])-1),
            combined_base_plus_unique=len(combined),combined_additional_unique=max(0,len(combined)-1)))
    return dict(scope='Same-observation pooled descriptive packing; cross-arm/round event counts never summed as global worlds',
        material_threshold=.03,structure_threshold=.08,minimum_XOR_cells=6,
        separate_material_structure_packings='Both use all immutable full-valid pool, independently of combined active packing',
        per_case=output,
        per_parent_mean={key:float(np.mean([row[key] for row in output])) for key in
            ('material_base_plus_unique','structure_base_plus_unique','combined_base_plus_unique')},
        total_with_scope={key:sum(row[key] for row in output) for key in
            ('material_base_plus_unique','structure_base_plus_unique','combined_base_plus_unique')},
        total_scope='Sum of separate per-observation representative counts; not worlds compatible with one shared observation')


def export(root=ROOT, job='g-e1-episode-export-v1'):
    root=Path(root);book=GBudget(root)
    # job() refuses any unfinished earlier physical job, even if RESULTS exists.
    with book.job(job,stage='ANALYSIS'):
        result=json.loads((root/'results/E1_LEGACY/RESULTS.json').read_text())
        if result.get('status')!='COMPLETE':raise ValueError('COMPLETE_PHYSICAL_RESULTS_REQUIRED')
        ledger=book.rows();validation_rows={item['candidate_id']:item for item in ledger
                     if item.get('event') in ('G_ENDPOINT_VALIDATION','G_BASE_VALIDATION')}
        candidates={};exported=[];stats=[];sensitivity_rows=[]
        directory=root/'GSearchEpisode';directory.mkdir(parents=True,exist_ok=True)
        for base in result['bases']:
            identifier=base['candidate_id'];receipt=validation_rows.get(identifier,{})
            if receipt.get('receipt_id')!=base['validation_receipt'] or not base['validation'].get('compatible'):
                raise ValueError('BASE_VALIDATION_RECEIPT_MISMATCH:'+identifier)
            candidates[identifier]=material_entry(dict(base,is_base=True),base['case_id'],base['validation'])
        for method in result['methods']:
            for episode in method['episodes']:
                identifier=episode['endpoint_id'];receipt=validation_rows.get(identifier,{})
                if receipt.get('receipt_id')!=episode.get('validation_receipt'):
                    raise ValueError('ENDPOINT_VALIDATION_RECEIPT_MISMATCH:'+identifier)
                if receipt.get('validation_passed') is not episode.get('valid'):
                    raise ValueError('ENDPOINT_VALIDATION_STATUS_MISMATCH:'+identifier)
                sidecar=episode_sidecar(episode)
                validate_episode(sidecar,asset_root=root,check_assets=True)
                relative='GSearchEpisode/'+identifier+'.json';dump(root/relative,sidecar)
                exported.append(dict(episode_id=identifier,path=relative,accepted=bool(episode.get('valid'))))
            for item in method['success_ledger']:
                identifier=item['candidate_id']
                if identifier not in candidates:
                    receipt=validation_rows.get(identifier,{})
                    if not receipt.get('validation_passed'):raise ValueError('UNVERIFIED_SUCCESS_LEDGER_ENTRY:'+identifier)
                    candidates[identifier]=material_entry(item,method['case_id'],sidecar='GSearchEpisode/'+identifier+'.json')
            stats.extend(summarize_arm(method))
            sensitivity_rows.extend(sensitivity(method,root,value) for value in (.04,.08,.16))
        for item in candidates.values():verify_material_assets(item,root)
        dump(root/'CANDIDATE_MATERIALS/index.json',dict(status='VALIDATED_CANDIDATE_ONLY',
             definition='All immutable full-valid success entries; unique endpoint ID, bases counted once; includes physical duplicates',
             count=len(candidates),candidates=list(candidates.values()),**MASS))
        summary=dict(scope='DESCRIPTIVE_EXPOSED_DEVELOPMENT_NO_SCIENTIFIC_GATE',rows=stats,
             pooled_same_observation=pooled_packings(candidates,root),
             base_validation_count=len(result['bases']),
             base_costs_not_allocated_to_arm=True,
             definitions=dict(prefix_t='Cumulative first 1/2/3 trajectories within each arm; zero-based stored trajectory converted to count',
                 endpoint_valid='Fresh full validation including data threshold, material box, solver/constitutive and prediction replay',
                 raw_valid='Raw step compatible flag; not independent endpoint acceptance',
                 new_material_events='Valid endpoint material-distinct from every prior immutable success; base not counted',
                 new_structure_events='Valid endpoint structural-distinct from every prior immutable success; base not counted',
                 combined_material_OR_structure_reward='Original archive separated-from-all success reward; never labeled structure discovery',
                 representatives='After-prefix active packing and its ID-ordered material/structure subpackings; includes base where selected',
                 retention='Signed projection onto raw edit, relative to parent; not probability and may exceed [0,1]',
                 correction_collapses='Valid endpoint; raw differs structurally from every prefix success but endpoint does not',
                 actual_wall_seconds='Sum original trajectory wall; base validation and surrounding campaign setup excluded'),
             sensitivity=sensitivity_rows)
        dump(root/'results/E1_DISCOVERY_SUMMARY.json',summary)
        with (root/'results/E1_DISCOVERY_SUMMARY.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(stats[0]) if stats else [])
            writer.writeheader();writer.writerows(stats)
        payload=dict(status='COMPLETE',episodes=len(exported),accepted_episodes=sum(item['accepted'] for item in exported),
                     valid_unique_candidates_including_bases=len(candidates),sidecars=exported,
                     materials_index='CANDIDATE_MATERIALS/index.json',native_refs_preserved=True,
                     fields_copied_or_invented=False,truth_loaded=False,new_physics_states=0,new_physics_RHS=0,
                     statistics='results/E1_DISCOVERY_SUMMARY.json')
        dump(root/'results/G_EPISODE_EXPORT.json',payload)
    return payload


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--job',default='g-e1-episode-export-v1');args=parser.parse_args()
    print(json.dumps(export(args.root,args.job),indent=2))
