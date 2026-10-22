"""Build one explicitly selected GPT Pro handoff ZIP, never run research.

Main accounts execution externally. Example:
  python tools/package_p3r2_pro_handoff.py --output-directory /NEW/pro-handoff
Only stdlib; no Book, Git, network, GT, NPZ/model interpretation, or scoring.
The default150MiB cap applies to both selected uncompressed input+metadata and
actual ZIP bytes. Oversize/missing inputs fail; partial outputs remain intact.
"""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import argparse
import json
import os
import resource
import stat
import sys
import time
import zipfile

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'research/p3l_r2_risk_posterior_exploration'
MIB=1024*1024
MAX_FILES=5000
DEFAULT_RELEASE='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr2-risk-posterior-20261023'


class Selection:
    def __init__(self,small_limit):
        self.small_limit=small_limit;self.files={};self.exclusions=[];self.missing=[]

    def add(self,path,group,required=False,small=False):
        path=Path(path)
        relative=path.relative_to(ROOT).as_posix()
        if path.is_symlink():
            self.exclusions.append(dict(path=relative,reason='symlink provenance only',link_target=os.readlink(path)))
            if required:self.missing.append(relative+' (required symlink is not bundled)')
            return
        if not path.is_file():
            if required:self.missing.append(relative)
            return
        if not path.resolve().is_relative_to(ROOT):raise ValueError('External path refused: '+relative)
        info=path.lstat()
        if not stat.S_ISREG(info.st_mode):raise ValueError('Special file refused: '+relative)
        if any(p in {'offline_only','endpoints','attempts','cache','endpoint_geometry','__pycache__','.git'} for p in path.relative_to(ROOT).parts):
            raise ValueError('Forbidden raw/GT/cache tree: '+relative)
        if small and info.st_size>self.small_limit:
            self.exclusions.append(dict(path=relative,size=info.st_size,reason='Public metadata above explicit small-file threshold; full evidence Release retains source'))
            return
        if relative not in self.files:
            if len(self.files)>=MAX_FILES:raise ValueError('Selected file-count bound exceeded; no truncation')
            self.files[relative]=dict(path=relative,size=info.st_size,mtime_ns=info.st_mtime_ns,
                                      inode=info.st_ino,mode=stat.S_IMODE(info.st_mode),groups=[group])
        elif group not in self.files[relative]['groups']:self.files[relative]['groups'].append(group)

    def tree(self,folder,group,suffixes,required=False,small=False):
        if not folder.is_dir() or folder.is_symlink():
            if required:self.missing.append(folder.relative_to(ROOT).as_posix()+' (directory missing or symlink)')
            return
        found=False
        for directory,folders,files in os.walk(folder,followlinks=False):
            folders[:]=sorted(x for x in folders if x not in {'__pycache__','.git','source','endpoints','attempts','cache'} and not (Path(directory)/x).is_symlink())
            for name in sorted(files):
                path=Path(directory)/name
                if path.suffix.lower() in suffixes:
                    found=True;self.add(path,group,small=small)
        if required and not found:self.missing.append(folder.relative_to(ROOT).as_posix()+' (no selected files)')

    def metadata(self,path):
        # Only fixed JSON configs/manifests are interpreted; never binary assets.
        if path.suffix!='.json' or path.is_symlink() or not path.is_file() or path.stat().st_size>8*MIB:
            raise ValueError('Required small real JSON manifest missing/oversize: '+str(path))
        return json.loads(path.read_text())

    def relocate(self,value):
        text=str(value);marker='/research/p3l_r2_risk_posterior_exploration/'
        path=OUT/text.split(marker,1)[1] if marker in text else OUT/Path(text)
        if not path.resolve().is_relative_to(OUT):raise ValueError('Replay artifact outside OUT')
        return path


def select(small_limit):
    plan=Selection(small_limit)
    for path in sorted(OUT.glob('*.md')):plan.add(path,'root_documents')
    for name in ('START_HERE.md','THEORY_AND_METHODS.md','REPRODUCE.md'):
        plan.add(OUT/name,'root_documents',required=True)
    plan.tree(OUT/'THEORY','theory', {'.md','.json','.tex','.pdf'})
    plan.tree(OUT/'GPT_PRO_HANDOFF','pro_questions',{'.md','.json','.csv'},required=True)
    for path in sorted((OUT/'configs').glob('*.json')):plan.add(path,'configs_and_freezes')
    for path in sorted(OUT.glob('*FROZEN*.json')):plan.add(path,'configs_and_freezes')
    for name in ('ENVIRONMENT.json','SOURCE_INVENTORY.json','FRESH_CHECKPOINT.json',
                 'MECHANISM_CHECKPOINT.json','CAMPAIGN_CLOSED.json'):
        plan.add(OUT/name,'provenance')
    for name in ('execution_review','risk_and_lock_review','locked_statistics','mechanism_confirmation','resources'):
        plan.tree(OUT/'AUDIT'/name,'audit/'+name,{'.md','.json'},required=True)
    plan.tree(OUT/'public_metrics','public_small_metadata',{'.json','.csv'},required=True,small=True)
    plan.tree(OUT/'tables','scored_tables',{'.csv'},required=True)
    for name in ('main_v2','confirmation'):
        plan.tree(OUT/'figures'/name,'figures/'+name,{'.png','.pdf','.json'},required=True)
    # Explicit mechanism arrays only. No arbitrary walk of the paid job tree.
    bracket=OUT/'C/jobs/C_projector_brackets'
    for name in ('RESULTS.json','PROTOCOL.json','CLOSED.json'):plan.add(bracket/name,'projector_brackets',required=True)
    for parent in ('DEV_000','DEV_030'):
        plan.add(bracket/parent/'CENTER.npz','projector_bracket_raw',required=True)
        for h in ('h0p001','h0p002'):
            labels=['BRACKET','Xminus','Xplus','Yminus','Yplus']
            labels += [f'{order}_leg{i}' for order in ('XY','YX') for i in range(4)]
            for label in labels:plan.add(bracket/parent/f'{h}_{label}.npz','projector_bracket_raw',required=True)
    # Selected final nonGaussian pilot arrays, not intermediate chains/LU.
    for parent in ('DEV_000','DEV_030'):
        folder=OUT/'B/jobs'/('B_source_'+parent)
        for name in ('PARTICLES.npz','EXPLICIT_TARGET_MATERIAL_MOMENTS.npz','RESULT.json','LANDSCAPE_PROGRESS.json'):
            plan.add(folder/name,'nonGaussian_small_raw',required=True)
    for arm in ('frozen','dynamic'):
        folder=OUT/'B/jobs/B_metric_pair_DEV030'/arm
        for name in ('MATERIAL_DIAGNOSTICS.npz','RESULT.json'):
            plan.add(folder/name,'MALA_diagnostics',required=True)
    # Chosen from metadata-only path/size inventory: direct GL8/12/20 moments
    # about22kB each; frame about492kB; GH97 moments about21kB. Exclude the
    # 16.6MB GH97_RAW node material array and all GL patch RAW/paid endpoints.
    direct=OUT/'B/jobs/B_field2_direct_CDF_DEV030_r0'
    field=OUT/'B/jobs/B_field2_slice_DEV030_r0'
    for name in ('FROZEN_PARTITION.npz','GL_8_MOMENTS.npz','GL_12_MOMENTS.npz','GL_20_MOMENTS.npz',
                 'RESULT.json','PROTOCOL.json','CLOSED.json','REFINEMENT_PROGRESS.json',
                 'INDEPENDENT_SAME_TARGET_COMPARISON.json','UNIT_SQUARE_DIRECT_CDF_CONTROLS.json'):
        plan.add(direct/name,'field2_direct_CDF',required=True)
    for name in ('FRAME.npz','GH_9_MOMENTS.npz','GH_17_MOMENTS.npz','GH_33_MOMENTS.npz',
                 'GH_65_MOMENTS.npz','GH_97_MOMENTS.npz','CHAIN_MOMENTS_AND_JENSEN.npz',
                 'COORDINATE_DIAGNOSTICS.npz','MATERIAL_DIAGNOSTICS.npz','RESULT.json','PROTOCOL.json'):
        plan.add(field/name,'field2_small_arrays',required=True)
    for name in ('p3r2','p3r2_mechanism','p1_material'):
        plan.tree(ROOT/'src'/name,'source/'+name,{'.py'},required=True)
    for path in sorted((ROOT/'tools').glob('*p3r2*.py')):plan.add(path,'research_tools')
    for name in ('requirements.txt','README.md','AGENTS.md'):
        plan.add(ROOT/name,'source_metadata')
    plan.add(ROOT/'research/p1_material_target_20261018/inputs/geometry.json','P1_geometry',required=True)
    # Preserve first-cohort direct/residual/fused replay, without loading weights.
    manifest_path=OUT/'METHODS_FROZEN_FOR_LOCKED.json'
    plan.add(manifest_path,'first_cohort_replay',required=True)
    manifest=plan.metadata(manifest_path)
    if manifest.get('status')!='FROZEN':raise ValueError('Original first-cohort method manifest must be FROZEN')
    for key in ('direct_checkpoint','residual_checkpoint'):
        model=plan.relocate(manifest[key])
        if model.name!='BEST.pt':raise ValueError('Unexpected replay checkpoint name')
        plan.add(model,'selected_replay_models',required=True)
    fusion_path=plan.relocate(manifest['fusion_json'])
    plan.add(fusion_path,'first_cohort_replay',required=True)
    fusion=plan.metadata(fusion_path)
    chosen=plan.relocate(fusion['selected']['checkpoint'])
    if chosen.name!='BEST.pt':raise ValueError('Unexpected selected fusion checkpoint')
    plan.add(chosen,'selected_replay_models',required=True)
    for name in ('L/contexts/R2_locked_00000.npz','L/contexts/R2_locked_00000.json',
                 'L/estimates/R2_locked_00000.npz','reproduction/REPLAY_00000.json'):
        plan.add(OUT/name,'first_cohort_replay',required=True)
    plan.provenance=dict(first_cohort_frozen_source_revision=manifest.get('source_revision'),
                         first_cohort_frozen_unix=manifest.get('frozen_unix'))
    return plan


def write_json(path,value):
    with path.open('x') as stream:json.dump(value,stream,indent=2,ensure_ascii=False,allow_nan=False);stream.write('\n')


class CappedFile:
    def __init__(self,file,limit):self.file=file;self.limit=limit
    def write(self,data):
        if self.file.tell()+len(data)>self.limit:raise ValueError('ZIP byte cap exceeded; partial archive retained')
        return self.file.write(data)
    def tell(self):return self.file.tell()
    def seek(self,*args):return self.file.seek(*args)
    def flush(self):return self.file.flush()


def readme(release):
    return f'''# GPT Pro research handoff — P3L-R2

Start with research/p3l_r2_risk_posterior_exploration/GPT_PRO_HANDOFF/START_HERE.md,
then THEORY_AND_METHODS.md and GPT_PRO_HANDOFF/THREE_THEORY_QUESTIONS.md.
HANDOFF_CONTENTS.json enumerates every selected repository-relative file, size,
source origin and explicit omissions. Preserve this tree when extracting.

This is a deliberately selected upload bundle. The approximately93GB full raw
evidence is delivered separately as multipart Release archives. Full endpoints, histories,
LU/cache artifacts, all training checkpoints and final full ledger are outside
this ZIP. Full evidence Release entry: {release}
This tool does not check remote availability and does not upload anything.
FULL posterior/reference qualification has NOT been established. Field2 SLICE
CDF evidence and projector-bracket arrays have their stated conditional/local
scope; they do not upgrade FULL. The first256 and the independent second256
confirmation cohorts must remain separate: do not pool them, refit a winner,
or call the confirmation21 models a new winner over the original256 study.
The confirmation uses7 families×3 seeds, mean of per-seed per-parent full-grid
MSE losses (not an ensemble image), with3 co-primary98.33333333333333% intervals.
Original training77+new9 runs and confirmation21 evaluated models are different
counts.512/128 anchor contexts are training/validation subsets, not new parents.

All source/array/model assets are copied as bytes. No truth archive was selected,
no NPZ/checkpoint was interpreted and no research computation was executed.
Binary arrays are explicitly selected:54 projector-bracket files, final small
source-SMC particle/moment arrays, MALA diagnostics, chosen small field2
CDF/Jensen summaries and one first-cohort saved online context/estimate.
The two selected BEST models and frozen fusion rule support the documented
zero-Maxwell direct/residual/fused replay. See REPRODUCE.md and ENVIRONMENT.json.
External historical DEV symlink targets and external A30 source are not bundled.
The source snapshot may include later delivery tools; frozen revision provenance
is preserved in the original manifests, not inferred by running Git.

Unzip into a NEW empty directory. A successful package has RECEIPT.json status
COMPLETE alongside pro_handoff.zip. Failed or oversize .partial archives are
retained for diagnosis, are not upload-ready, and require a new output directory
for another attempt. This packaging task does not certify scientific conclusions.
'''


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-directory',required=True,type=Path)
    parser.add_argument('--max-bytes',type=int,default=150*MIB,help='cap on selected uncompressed bytes+metadata AND actual ZIP bytes')
    parser.add_argument('--small-metadata-bytes',type=int,default=MIB,help='explicit public_metrics JSON/CSV selection threshold')
    parser.add_argument('--release-url',default=DEFAULT_RELEASE,help='existing release entry; recorded, not fetched')
    args=parser.parse_args();destination=args.output_directory.expanduser().resolve()
    if destination==OUT.resolve() or OUT.resolve() in destination.parents:raise ValueError('Output must be outside OUT')
    if destination.exists():raise FileExistsError('Choose a NEW output directory; preserve prior handoff evidence')
    if args.max_bytes<=0 or args.small_metadata_bytes<=0:raise ValueError('Byte limits must be positive')
    destination.mkdir(parents=True,exist_ok=False)
    start=time.perf_counter();receipt=dict(status='FAILED',started_utc=datetime.now(timezone.utc).isoformat(),
        max_bytes=args.max_bytes,source_root=str(ROOT),output_directory=str(destination),hashes_used=False,
        physics_calls=0,truth_read=False,NPZ_or_checkpoint_interpretation=False)
    try:
        plan=select(args.small_metadata_bytes)
        rows=sorted(plan.files.values(),key=lambda r:r['path'])
        groups=Counter()
        for row in rows:
            for group in row['groups']:groups[group]+=row['size']
        contents=dict(status='EXPLICIT_SELECTED_CONTENTS_NOT_COMPLETION_RECEIPT',files=rows,
            file_count=len(rows),selected_input_bytes=sum(r['size'] for r in rows),
            source_root=str(ROOT),version_provenance=plan.provenance,group_bytes=dict(groups),
            public_metadata_small_file_limit=args.small_metadata_bytes,explicit_exclusions=plan.exclusions,
            missing_required=plan.missing,
            selection_policy='All root Markdown; explicit review/config/source/figure trees; fixed mechanism arrays; selected first-cohort replay only',
            never_selected=['inputs/offline_only','full RESOURCE_LEDGER.jsonl','LU/cache directories','paid endpoints/attempts',
                'GL patch RAW arrays','GH97_RAW node arrays','all historical chains','all training checkpoints except selected replay BEST models',
                'external symlink targets','older studies'],release_entry_not_verified=args.release_url)
        inventory_path=destination/'HANDOFF_CONTENTS.json';write_json(inventory_path,contents)
        readme_path=destination/'HANDOFF_README.md'
        with readme_path.open('x') as stream:stream.write(readme(args.release_url))
        total=contents['selected_input_bytes']+inventory_path.stat().st_size+readme_path.stat().st_size
        receipt.update(selected_file_count=len(rows),selected_input_bytes=contents['selected_input_bytes'],
                       total_uncompressed_bytes_including_metadata=total,missing_required=plan.missing,
                       excluded_public_metadata_count=len(plan.exclusions),group_bytes=dict(groups))
        if plan.missing:raise ValueError('Required handoff assets missing; see HANDOFF_CONTENTS.json')
        if total>args.max_bytes:raise ValueError(f'Selected bytes {total} exceed cap {args.max_bytes}; no truncation or automatic downselection')
        partial=destination/'pro_handoff.zip.partial'
        with partial.open('xb') as raw:
            capped=CappedFile(raw,args.max_bytes)
            with zipfile.ZipFile(capped,mode='w',compression=zipfile.ZIP_DEFLATED,compresslevel=1,allowZip64=True) as archive:
                for row in rows:
                    source=ROOT/row['path'];before=source.lstat()
                    if source.is_symlink() or (before.st_size,before.st_mtime_ns,before.st_ino)!=(row['size'],row['mtime_ns'],row['inode']):
                        raise ValueError('Source changed before copy: '+row['path'])
                    info=zipfile.ZipInfo(row['path'],date_time=time.localtime(max(before.st_mtime,315532800))[:6])
                    info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=row['mode']<<16;info.file_size=row['size']
                    fd=os.open(source,os.O_RDONLY|os.O_NOFOLLOW)
                    with os.fdopen(fd,'rb') as stream,archive.open(info,'w',force_zip64=True) as target:
                        opened=os.fstat(stream.fileno())
                        if (opened.st_size,opened.st_mtime_ns,opened.st_ino)!=(row['size'],row['mtime_ns'],row['inode']):raise ValueError('Source changed at open')
                        while True:
                            chunk=stream.read(MIB)
                            if not chunk:break
                            target.write(chunk)
                    after=source.lstat()
                    if (after.st_size,after.st_mtime_ns,after.st_ino)!=(row['size'],row['mtime_ns'],row['inode']):raise ValueError('Source changed during copy: '+row['path'])
                # Small generated metadata, streamed like every other asset.
                for path in (readme_path,inventory_path):
                    with path.open('rb') as source,archive.open(path.name,'w',force_zip64=True) as target:
                        while True:
                            chunk=source.read(MIB)
                            if not chunk:break
                            target.write(chunk)
            raw.flush();os.fsync(raw.fileno())
        final=destination/'pro_handoff.zip';partial.rename(final)
        receipt.update(status='COMPLETE',archive=final.name,archive_bytes=final.stat().st_size)
    except BaseException as error:
        receipt['error']=repr(error)
        raise
    finally:
        partial=destination/'pro_handoff.zip.partial'
        if partial.exists():receipt['retained_partial_bytes']=partial.stat().st_size
        receipt.update(wall_seconds=time.perf_counter()-start,
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024))
        write_json(destination/'RECEIPT.json',receipt)
    print(json.dumps(dict(status=receipt['status'],archive=str(destination/'pro_handoff.zip'),bytes=receipt['archive_bytes'])))


if __name__=='__main__':main()
