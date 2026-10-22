"""Collect scored P3L-R2 tables and explicit public JSON metadata, without physics.

Main runs this under its stage D analysis wrapper after both evaluations close:
  python tools/collect_p3r2_public_results.py
No Book import, truth/context/LU/checkpoint access, fitting, Git, or upload.
Existing generated targets are never replaced. Source JSON origins are recorded.
"""
from pathlib import Path
import argparse
import csv
from datetime import datetime, timezone
import json
import shutil
import time
import zipfile
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT/'research/p3l_r2_risk_posterior_exploration'
LIMIT = 64*1024*1024
TRAINING = ('training','training32768','residual_training','residual_controls','residual_local_controls',
            'mechanism_controls','random_sensitivity_controls')
JSON_NAMES = ('RESULT.json','RESULTS.json','PROTOCOL.json','CLOSED.json')
STATUSES = ('status','scientific_status','posterior_qualification','reference_qualification',
            'target_interpretation','scope','interpretation','calibrated_posterior','calibrated_density',
            'primary_baseline_qualification','original_FULL_posterior')


class Collector:
    def __init__(self, root):
        self.root=root
        self.files=[]
        self.missing=[]
        self.records=[]
        self.confirmation=False

    def source(self, relative):
        path=self.root/relative
        resolved=path.resolve()
        if self.root not in resolved.parents or path.is_symlink():
            raise ValueError('External source refused: '+str(relative))
        if any(p in {'inputs','offline_only','contexts','endpoint_geometry','source','checkpoints'} for p in path.relative_to(self.root).parts):
            raise ValueError('Forbidden data tree: '+str(relative))
        if not path.is_file():return None
        if path.stat().st_size>LIMIT:raise ValueError('Metadata input exceeds bounded read limit: '+str(relative))
        return path

    def read(self, relative):
        path=self.source(relative)
        if path is None:return None
        return json.loads(path.read_text())

    def absent(self, group, source):
        self.missing.append(dict(group=group,status='NOT_RUN',source=str(source)))

    def register(self, path, origins, statuses=None):
        item=dict(path=path.relative_to(self.root).as_posix(),size=path.stat().st_size,
                  source_origin_paths=[str((self.root/p).absolute()) for p in origins])
        if statuses:item['scientific_status_verbatim']=statuses
        self.files.append(item)
        return item

    def json(self, relative, value, origins=(), statuses=None):
        path=self.root/relative
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('x') as stream:
            json.dump(value,stream,indent=2,allow_nan=False,ensure_ascii=False)
            stream.write('\n')
        return self.register(path,origins,statuses)

    def table(self, name, columns, rows, origins):
        # Delay creation until the first actual observation: no empty placeholders.
        iterator=iter(rows)
        try:first=next(iterator)
        except StopIteration:
            self.absent(name,origins)
            return
        path=self.root/'tables'/name
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('x',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=columns,extrasaction='raise')
            writer.writeheader();writer.writerow(first)
            for row in iterator:writer.writerow(row)
        self.register(path,origins)

    def copy_json(self, relative):
        value=self.read(relative)
        if value is None:return False
        status={key:value[key] for key in STATUSES if key in value}
        item=self.json('public_metrics/'+relative,value,[relative],status)
        self.records.append(item)
        return True

    def raw_risks(self):
        relative='L/evaluation/RAW_RISKS.npz'
        path=self.source(relative)
        if path is None:
            self.absent('LOCKED_PARENT_RISKS',relative);return
        # NPZ file itself may be larger than summaries in other contexts; this
        # scored-only archive is tiny. Read only explicitly scored metric keys.
        wanted=('parent','method','fullgrid_MSE','fullgrid_MAE','ROI_MSE','ROI_MAE','method_used_fallback')
        with zipfile.ZipFile(path) as z:
            for key in wanted:
                member=key+'.npy'
                if member in z.namelist() and z.getinfo(member).file_size>LIMIT:
                    raise ValueError('Oversize scored key: '+key)
        with np.load(path,allow_pickle=False) as z:
            parents=z['parent'].astype(str);methods=z['method'].astype(str)
            if len(parents)!=256 or len(set(parents))!=256:
                raise ValueError('Scored parent risks must retain all256')
            if list(parents)!=self.parents:raise ValueError('Scored parent order differs from COMPLETE')
            metrics={k:z[k] for k in wanted[2:] if k in z.files}
            for key,value in metrics.items():
                if value.shape!=(256,len(methods)):raise ValueError('Scored risk shape mismatch: '+key)
            if not {'fullgrid_MSE','fullgrid_MAE'}<=metrics.keys():raise ValueError('Missing full-grid scored risks')
            def records():
                for i,parent in enumerate(parents):
                    for j,method in enumerate(methods):
                        row=dict(parent=parent,method=method)
                        for key in wanted[2:]:
                            row[key]=(bool(metrics[key][i,j]) if key=='method_used_fallback' else float(metrics[key][i,j])) if key in metrics else ''
                        yield row
            self.table('LOCKED_PARENT_RISKS.csv',wanted,records(),[relative])

    def validation(self):
        origins=[];rows=[];density=[]
        for folder in TRAINING:
            rel=f'A/{folder}/RESULTS.json';value=self.read(rel)
            if value is None:self.absent('validation/'+folder,rel);continue
            origins.append(rel)
            for run in value.get('results',[]):
                validation=run.get('validation',{})
                if 'complex_roi_mse' not in validation:
                    self.absent('validation/'+folder+'/'+str(run.get('kind')),rel);continue
                rows.append(dict(family=folder,kind=run['kind'],seed=run['seed'],N=run['train_size'],
                    validation_parents=run.get('validation_size',''),selected_update=run.get('selected_update',''),
                    parameter_count=run.get('parameter_count',''),val_fullgrid_MSE=validation['complex_roi_mse']*256/576,
                    cpu_seconds=run.get('cpu_seconds',''),wall_seconds=run.get('wall_seconds',''),source=rel))
        self.table('VALIDATION_RUNS.csv',('family','kind','seed','N','validation_parents','selected_update','parameter_count',
                   'val_fullgrid_MSE','cpu_seconds','wall_seconds','source'),rows,origins)
        density_origins=[]
        # Proper coordinate NLL is never mixed into material-MSE columns.
        for folder in ('density_training','weak_density_training','density','weak_density'):
            rel=f'A/{folder}/RESULTS.json';value=self.read(rel)
            if value is None:continue
            density_origins.append(rel)
            for run in value.get('results',[]):
                if 'validation_nll' in run:
                    density.append(dict(family=folder,kind=run.get('kind',''),seed=run.get('seed',''),N=run.get('train_size',''),
                        validation_NLL=run['validation_nll'],selected_update=run.get('selected_update',''),
                        parameter_count=run.get('parameter_count',''),cpu_seconds=run.get('cpu_seconds',''),
                        wall_seconds=run.get('wall_seconds',''),source=rel))
        self.table('DENSITY_VALIDATION_RUNS.csv',('family','kind','seed','N','validation_NLL','selected_update',
                   'parameter_count','cpu_seconds','wall_seconds','source'),density,density_origins)

    def trajectories(self):
        rel='L/trajectory_evaluation/RESULTS.json';result=self.read(rel)
        if result is None:self.absent('LOCKED_TRAJECTORY_RISKS',rel);return
        if result.get('parents')!=self.parents:raise ValueError('Trajectory result parent order mismatch')
        risks=result['per_parent_risks'];core=result['primary_arm'];summaries=result['summaries']
        def records():
            for arm,methods in summaries.items():
                cohort=self.parents if arm==core else self.parents[:32]
                for method,summary in methods.items():
                    if summary['parents']!=len(cohort):raise ValueError('Trajectory cohort count mismatch')
                    fallback=set(summary.get('fallback_parents',[]))
                    diagnostic=result.get('actual_candidate_diagnostics',{}).get(arm,{})
                    for parent in cohort:
                        value=risks[parent][arm+'__'+method]
                        yield dict(arm=arm,method=method,parent=parent,cohort_parents=len(cohort),
                            fullgrid_MSE=value['fullgrid_MSE'],fullgrid_MAE=value['fullgrid_MAE'],
                            fallback_to_anchor=(parent in fallback) if 'fallback_parents' in summary else '',
                            parent_incremental_RHS='',
                            arm_cohort_incremental_RHS=diagnostic.get('incremental_step_RHS',''),
                            cost_scope='cohort total selected-step RHS; per-parent cost unavailable in scored summary',source=rel)
        self.table('LOCKED_TRAJECTORY_RISKS.csv',('arm','method','parent','cohort_parents','fullgrid_MSE','fullgrid_MAE',
                   'fallback_to_anchor','parent_incremental_RHS','arm_cohort_incremental_RHS','cost_scope','source'),records(),[rel])

    def public_json(self):
        count_by_group={}
        for group in ('A','B','C','D','L'):
            count=0
            # Explicit stage/job metadata filenames only; no recursive bulk copy.
            for name in JSON_NAMES:
                count+=self.copy_json(f'{group}/{name}')
            jobs=self.root/group/'jobs'
            if jobs.exists() and not jobs.is_symlink():
                for job in jobs.iterdir():
                    if job.is_dir() and not job.is_symlink():
                        for name in JSON_NAMES:
                            count+=self.copy_json(f'{group}/jobs/{job.name}/{name}')
                        # Two meaningful MALA mechanism summaries; never chains.
                        for arm in ('frozen','dynamic'):
                            count+=self.copy_json(f'{group}/jobs/{job.name}/{arm}/RESULT.json')
            count_by_group[group]=count
        for folder in TRAINING+('density_training','weak_density_training','density','weak_density'):
            for name in ('RESULTS.json','RUN_CONTRACT.json','SELECTED.json','COMPLETE.json'):
                self.copy_json(f'A/{folder}/{name}')
            if folder in ('mechanism_controls','random_sensitivity_controls'):
                aggregate=self.read(f'A/{folder}/RESULTS.json')
                for run in (aggregate or {}).get('results',[]):
                    if 'path' not in run:continue
                    run_path=Path(run['path'])
                    if not run_path.is_absolute():run_path=self.root/run_path
                    base=(self.root/'A'/folder).resolve()
                    if base not in run_path.resolve().parents:raise ValueError('New run metadata outside its declared training group')
                    self.copy_json((run_path/'RESULT.json').relative_to(self.root).as_posix())
            # The selected run is a bounded pointer originating in that training
            # directory; checkpoints and training arrays are never opened.
            selected=self.read(f'A/{folder}/SELECTED.json')
            if selected:
                pointers=[selected[k] for k in ('path','selected_path') if k in selected]
                if isinstance(selected.get('result'),dict) and 'path' in selected['result']:
                    pointers.append(selected['result']['path'])
                for pointer in dict.fromkeys(pointers):
                    selected_path=Path(pointer)
                    if not selected_path.is_absolute():selected_path=self.root/selected_path
                    base=(self.root/'A'/folder).resolve()
                    if base not in selected_path.resolve().parents:continue
                    for name in ('RESULT.json','RESULTS.json'):
                        candidate=selected_path/name
                        self.copy_json(candidate.relative_to(self.root).as_posix())
        for folder in ('L/evaluation','L/trajectory_evaluation'):
            for name in ('RESULTS.json','COMPLETE.json'):
                self.copy_json(folder+'/'+name)
        for name in ('TARGET.json','PROTOCOL.json','AUTHORIZATION.json','A_DATA_SCALE_REVISION_002.json',
                     'LOCAL_RISK_CONTROLS_002.json','DENSITY_FRAME_ABLATION.json','EXECUTION_MEMORY_OVERLAY.json',
                     'MECHANISM_CONFIRMATION_001.json','MECHANISM_CONFIRMATION_002.json'):
            self.copy_json('configs/'+name)
        self.copy_json('METHODS_FROZEN_FOR_LOCKED.json')
        if self.confirmation:
            for relative in ('MECHANISM_CHECKPOINT.json','MECHANISM_CONFIRMATION_FROZEN.json',
                    'L/mechanism_confirmation/evaluation/RESULT.json',
                    'L/mechanism_confirmation/evaluation/COMPLETE.json'):
                if not self.copy_json(relative):raise ValueError('Required confirmation metadata missing: '+relative)
            relative='L/mechanism_confirmation/evaluation/PER_PARENT_RISKS.csv'
            source=self.source(relative)
            if source is None:raise ValueError('Confirmation scored CSV missing')
            # Preserve scored values/columns verbatim; never recompute scores or
            # merge the confirmation cohort into the original256 risk table.
            target=self.root/'public_metrics'/relative
            target.parent.mkdir(parents=True,exist_ok=True)
            with source.open('rb') as src,target.open('xb') as dst:shutil.copyfileobj(src,dst,1024*1024)
            self.records.append(self.register(target,[relative]))
        for group in count_by_group:
            count_by_group[group]=sum(item['path'].startswith('public_metrics/'+group+'/') for item in self.records)
            if not count_by_group[group]:self.absent('public_metrics/'+group,'No whitelisted existing metadata')
        return count_by_group

    def resource_snapshot(self):
        rel='RESOURCE_LEDGER.json';value=self.read(rel)
        if value is None:self.absent('resource_snapshot',rel);return
        sequence=value['sequence'];stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        self.json(f'public_metrics/resources/SNAPSHOT_seq{sequence}_{stamp}.json',
            dict(snapshot_kind='SEQUENCE_TIMESTAMPED_SNAPSHOT_NOT_FINAL_CLOSURE',captured_utc=stamp,
                 source_origin_path=str(self.root/rel),sequence=sequence,
                 may_include_active_packaging_or_analysis_job=True,ledger_snapshot=value),[rel])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=DEFAULT_OUT)
    args=parser.parse_args();out=args.root.resolve();c=Collector(out)
    markers=[c.read(f'L/{folder}/COMPLETE.json') for folder in ('evaluation','trajectory_evaluation')]
    if any(not m or m.get('status')!='COMPLETE' or m.get('count')!=256 or len(set(m.get('parents',[])))!=256 for m in markers):
        raise ValueError('Both exact256 scored evaluation COMPLETE gates required')
    for key in ('parents','source_revision','frozen_unix'):
        if key not in markers[0] or markers[0][key]!=markers[1].get(key):raise ValueError('Evaluation freeze mismatch: '+key)
    c.parents=markers[0]['parents']
    if len(c.parents)!=256:raise ValueError('COMPLETE parent list must contain exactly256 entries')
    checkpoint=c.read('MECHANISM_CHECKPOINT.json')
    if checkpoint is not None:
        if checkpoint.get('status')!='COMPLETE':raise ValueError('Mechanism confirmation is registered but not COMPLETE')
        confirm=c.read('L/mechanism_confirmation/evaluation/COMPLETE.json')
        expected=[f'R2_locked_{i:05d}' for i in range(256,512)]
        if not confirm or confirm.get('status')!='COMPLETE' or confirm.get('count')!=256 or confirm.get('parents')!=expected:
            raise ValueError('Separate second256 scored COMPLETE gate required')
        result=c.read('L/mechanism_confirmation/evaluation/RESULT.json')
        freeze=c.read('MECHANISM_CONFIRMATION_FROZEN.json')
        if not result or result.get('status')!='COMPLETE' or result.get('parents')!=expected or result.get('model_count')!=21:
            raise ValueError('Second cohort must preserve complete21-model scores')
        if not freeze or freeze.get('status')!='FROZEN' or freeze.get('parents')!=expected or len(freeze.get('models',[]))!=21:
            raise ValueError('Separate21-model mechanism freeze required')
        for amendment in ('MECHANISM_CONFIRMATION_001.json','MECHANISM_CONFIRMATION_002.json'):
            if c.read('configs/'+amendment) is None:raise ValueError('Confirmation amendment missing: '+amendment)
        for folder,count in (('mechanism_controls',6),('random_sensitivity_controls',3)):
            complete=c.read(f'A/{folder}/COMPLETE.json');runs=c.read(f'A/{folder}/RESULTS.json')
            if not complete or complete.get('status')!='COMPLETE' or complete.get('runs')!=count or not runs or len(runs.get('results',[]))!=count:
                raise ValueError('All new mechanism training run metadata required: '+folder)
        c.confirmation=True
    # Require an unused publication namespace to preserve earlier derived output.
    for path in (out/'public_metrics',out/'PUBLICATION_FILES.json'):
        if path.exists():raise FileExistsError('Preserve existing public collection: '+str(path))
    for name in ('VALIDATION_RUNS.csv','DENSITY_VALIDATION_RUNS.csv','LOCKED_PARENT_RISKS.csv','LOCKED_TRAJECTORY_RISKS.csv'):
        if (out/'tables'/name).exists():raise FileExistsError('Preserve existing table: '+name)
    start=time.monotonic();complete=False
    try:
        c.validation();c.raw_risks();c.trajectories();groups=c.public_json();c.resource_snapshot()
        for name in ('REPRODUCE.md','ENVIRONMENT.json'):
            path=c.source(name)
            if path is None:c.absent('root_reproduction_artifact',name)
            else:
                entry=c.register(path,[name]);entry['publication_role']='EXISTING_ROOT_REPRODUCTION_DOCUMENT'
        dataset=dict(independent_training_parents=32768,independent_validation_parents=2048,first_fresh_locked_parents=256,
            second_independent_confirmation_parents=256 if c.confirmation else 0,
            total_independent_parents=35328 if c.confirmation else 35072,anchor_training_subset=512,anchor_validation_subset=128,
            model_scope=dict(original_training_runs_declared=77,new_mechanism_runs_declared=9 if c.confirmation else 0,
                total_training_runs_declared=86 if c.confirmation else 77,
                second_cohort_evaluated_models=21 if c.confirmation else None,second_cohort_families=7 if c.confirmation else None,
                seeds_per_confirmation_family=3 if c.confirmation else None,
                confirmation_action='separate full-grid MSE; mean of three seed-specific per-parent risks, not ensemble images',
                multiplicity='three co-primary contrasts;98.33333333333333% CI; published scored values copied without recalculation',
                original256_winner_reselection=False,cohorts_not_pooled=True),
            anchor_subsets_are_additional_independent_parents=False,
            interpretation='512 training and128 validation anchors are subsets of the32768/2048 cohorts, never additional independent parents')
        c.json('public_metrics/DATASET_METADATA.json',dataset)
        c.json('public_metrics/INDEX.json',dict(dataset=dataset,groups=groups,files=list(c.files),
            missing_groups=c.missing,source_revision=markers[0]['source_revision'],frozen_unix=markers[0]['frozen_unix'],
            status='DERIVED_METADATA_COLLECTION',truth_context_LU_checkpoint_reads=False,
            scientific_status_policy='Existing source status/qualification fields copied verbatim; no new scientific verdict',
            resource_policy='Sequence-timestamped snapshot only; final campaign closure is written separately'))
        complete=True
    finally:
        # Even failed collections expose exactly which files were created;
        # main can inspect this manifest before any directed staging/upload.
        manifest=dict(status='COMPLETE' if complete else 'FAILED_PARTIAL_OUTPUT_RETAINED',
            wall_seconds=time.monotonic()-start,files=list(c.files),missing_groups=c.missing,
            paths_relative_to=str(out),git_paths=[str((out/item['path']).relative_to(ROOT)) for item in c.files],
            no_hashes=True,no_git_or_upload_performed=True,
            root_reproduction_documents=[item['path'] for item in c.files if item.get('publication_role')=='EXISTING_ROOT_REPRODUCTION_DOCUMENT'],
            manifest_path=str(out/'PUBLICATION_FILES.json'))
        path=out/'PUBLICATION_FILES.json'
        manifest['git_paths'].append(str(path.relative_to(ROOT)))
        with path.open('x') as stream:json.dump(manifest,stream,indent=2,allow_nan=False)
    print(json.dumps(dict(status='COMPLETE',manifest=str(out/'PUBLICATION_FILES.json'),generated_files=len(c.files))))


if __name__=='__main__':main()
