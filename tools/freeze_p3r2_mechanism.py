"""Freeze a separately preregistered21-model mechanism confirmation before data."""
from pathlib import Path
import json,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'tools'))
from p3r2.core import OUT,CAPS,TOTAL_RHS,write
from p3r2_mechanism.contract import (VERSION,KIND,PARENTS,FAMILIES,SEEDS,PRIMARY,PRIMARY_CONFIDENCE,
    REQUIRED_SOURCE_FILES,file_record,checkpoint_metadata,validate_freeze)
from infer_p3r2_locked import frozen_contract


def main():
    path=OUT/'MECHANISM_CONFIRMATION_FROZEN.json'
    if path.exists():raise FileExistsError('Confirmation already frozen; no replacement')
    old,_,_,_=frozen_contract(OUT,OUT/'METHODS_FROZEN_FOR_LOCKED.json')
    if json.loads((OUT/'FRESH_CHECKPOINT.json').read_text())['status']!='COMPLETE':raise RuntimeError('Preserve finished first cohort before confirmation')
    for parent in PARENTS:
        for p in (OUT/'inputs/online'/f'{parent}.npz',OUT/'inputs/offline_only'/f'{parent}.npz',OUT/'L/contexts'/f'{parent}.npz'):
            if p.exists():raise RuntimeError('A new-cohort artifact exists before freeze: '+str(p))
    for group in ('mechanism_controls','random_sensitivity_controls'):
        if json.loads((OUT/'A'/group/'COMPLETE.json').read_text()).get('status')!='COMPLETE':raise RuntimeError('All added controls must finish before new observations')
    source_files=list(REQUIRED_SOURCE_FILES)+['tools/freeze_p3r2_mechanism.py','tools/run_p3r2_mechanism_sequence.py',
        'src/p3r2/core.py','src/p3r2/anchors.py','src/p3r2/inversion.py','src/p3r2/fusion.py']
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    subprocess.run(['git','diff','--exit-code',revision,'--',*source_files],cwd=ROOT,check=True)
    if subprocess.check_output(['git','ls-files','--others','--',*source_files],cwd=ROOT,text=True).strip():raise RuntimeError('All executed confirmation source must be committed')
    original={'full_geometry':'local_geometry_cnn','local_anchor':'local_anchor_cnn',
        'capacity':'local_capacity_cnn','permuted':'local_permuted_geometry_cnn'}
    models=[]
    for family in FAMILIES:
        for seed in SEEDS:
            if family in original:
                name=f'{original[family]}_n512_s{seed}'
                matches=[e for e in old['ablation_checkpoints'] if e['name']==name]
                if len(matches)!=1:raise RuntimeError('Original immutable checkpoint missing: '+name)
                checkpoint=Path(matches[0]['checkpoint']);loader='original'
            else:
                group='random_sensitivity_controls' if family=='decoder_random_sensitivity' else 'mechanism_controls'
                checkpoint=OUT/'A'/group/'n512'/f'{family}_s{seed}'/'BEST.pt'
                loader='random_sensitivity' if group=='random_sensitivity_controls' else 'mechanism'
            if not checkpoint.is_absolute():checkpoint=OUT/checkpoint
            entry=dict(name=f'{family}_s{seed}',family=family,seed=seed,checkpoint=str(checkpoint.resolve()),loader=loader)
            checkpoint_metadata(checkpoint,entry);models.append(entry)
    target=OUT/'configs/TARGET.json';geometry=ROOT/'research/p1_material_target_20261018/inputs/geometry.json'
    artifacts=[ROOT/p for p in source_files]+[Path(e['checkpoint']) for e in models]+[target,geometry,
        OUT/'configs/MECHANISM_CONFIRMATION_001.json',OUT/'configs/MECHANISM_CONFIRMATION_002.json']
    records=[file_record(p) for p in artifacts]
    current=json.loads((OUT/'RESOURCE_LEDGER.json').read_text())
    if current['stages']['L']['RHS']+1640448>CAPS['L'] or sum(v['RHS'] for v in current['stages'].values())+1640448>TOTAL_RHS:
        raise RuntimeError('Worst-case confirmation no longer fits original caps')
    manifest=dict(version=VERSION,protocol_kind=KIND,status='FROZEN',source_revision=revision,frozen_unix=time.time(),
        parents=PARENTS,expected_count=256,models=models,frozen_source_files=source_files,frozen_files=records,
        primary_comparisons=[list(p) for p in PRIMARY],bootstrap_replicates=2000,bootstrap_seed=2026102393,
        primary_CI=PRIMARY_CONFIDENCE,primary_metric='complex_fullgrid_MSE',seed_averaging='Average three per-parent losses, not ensemble images',
        estimator_action='raw_unprojected_no_blend',target_config=str(target),geometry_config=str(geometry),
        generator_closed_json='L/jobs/data_locked_256_256/CLOSED.json',
        context_closed_jsons=['L/jobs/anchors_locked_256_256_96_64/CLOSED.json'],
        online_inference_closure='L/jobs/mechanism_confirmation_inference/CLOSED.json',
        old256_unchanged=True,old_scores_not_used_for_training=True,stop_after_this_confirmation=True,
        design_amendments=['configs/MECHANISM_CONFIRMATION_001.json','configs/MECHANISM_CONFIRMATION_002.json'],
        maximum_added_physics_RHS=1640448,authorization='Original P3LR2 L3000000/total20000000 caps unchanged',
        prefreeze_budget_sequence=current['sequence'],no_hash_checks=True)
    write(path,manifest);validate_freeze(OUT,path)
    print('21 models and 256 new IDs frozen before first new observation; source '+revision)


if __name__=='__main__':main()
