"""Transport-only relocation into a NEW tree; original archives stay unchanged.
Raw shards are repository-relative; dependencies may be Gaussian-relative.
Synthetic blinded GT is copied only with explicit offline-evidence opt-in, after
completed seal+score. No array-value parsing, scientific imports or physics.
"""
from pathlib import Path
import argparse,hashlib,json,os,shutil,time
TEXT={'.py','.json','.jsonl','.md','.csv','.html','.txt','.toml','.yml','.yaml'}
ALWAYS_EXCLUDED={'.git','__pycache__','.venv','.venv_nn'}
BLINDED={'private','offline_truth','offline_only'}
REPO_PREFIX=Path('Gaussian/B1/b1-p3l-r9-stabilized-subspace')
STUDY_REL=Path('research/p3l_r9_stabilized_subspace')
STUDY=REPO_PREFIX/STUDY_REL

def digest(path):
    h=hashlib.blake2b()
    with path.open('rb') as f:
        while block:=f.read(4*1024*1024):h.update(block)
    return h.hexdigest()
def dump(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,ensure_ascii=False,allow_nan=False,indent=2)+'\n')
def read(path):return json.loads(path.read_text())
def layout(archive):
    if (archive/STUDY).is_dir():return archive/STUDY,Path()
    if (archive/STUDY_REL).is_dir():return archive/STUDY_REL,REPO_PREFIX
    raise RuntimeError('Expected repository-relative research/study or Gaussian-root layout')
def relocate(archive,replay,original_root,dependency_root=None,include_offline=False,prepare_score_replay=False):
    archive=Path(archive).absolute();replay=Path(replay).absolute()
    old=str(Path(original_root).absolute()).rstrip('/');new=str(replay).rstrip('/')
    source_study,prefix=layout(archive);repo_layout=prefix!=Path()
    if replay.exists() or replay.is_relative_to(archive):raise RuntimeError('Replay must be NEW and outside archive')
    original_protocol=read(source_study/'config/possibility_volume_protocol.json')
    original_seal=read(source_study/'results/possibility_volume/PREDICTIONS_SEALED.json')
    if original_seal.get('status')!='ALL32_ONLINE_PREDICTIONS_SEALED' or len(original_seal['parents'])!=32:
        raise RuntimeError('Completed all32 public prediction seal required')
    if include_offline:
        complete=read(source_study/'results/possibility_volume/offline_statistics/SCORING_COMPLETE.json')
        if complete.get('status')!='COMPLETE':raise RuntimeError('Synthetic OFFLINE_ONLY copying requires completed scoring')
    if prepare_score_replay and not include_offline:
        raise RuntimeError('--prepare-score-replay requires explicit --include-offline-synthetic-truth')
    old_private=Path(original_protocol['private_output'])
    if not old_private.is_relative_to(Path(old)):raise RuntimeError('Unexpected blinded synthetic evidence root')
    allowed_private=old_private.relative_to(Path(old))
    if not allowed_private.is_relative_to(STUDY):raise RuntimeError('Only this E5 synthetic evidence namespace is permitted')
    if include_offline:
        source_private=archive/(allowed_private.relative_to(prefix) if repo_layout else allowed_private)
        missing=[p['parent'] for p in original_protocol['parents'] if not (source_private/'truth'/(p['parent']+'.npz')).exists()]
        if missing:raise RuntimeError('Completed raw archive lacks opted-in OFFLINE_ONLY synthetic truth: '+str(missing))
    replay.mkdir(parents=True,exist_ok=False)
    inventories={};disabled=[];sources=[(archive,prefix,'REPOSITORY_RELATIVE_SHARDS' if repo_layout else 'GAUSSIAN_LAYOUT_ARCHIVE')]
    if dependency_root is not None:
        dep=Path(dependency_root).absolute()
        if not (dep/'Gaussian').is_dir():raise RuntimeError('Dependency overlay must contain Gaussian/')
        sources.append((dep,Path(),'GAUSSIAN_RELATIVE_DEPENDENCIES'))
    def may_enter(rel):
        if any(p in ALWAYS_EXCLUDED for p in rel.parts):return False
        if not any(p in BLINDED for p in rel.parts):return True
        return include_offline and (rel.is_relative_to(allowed_private) or allowed_private.is_relative_to(rel))
    for srcroot,destprefix,source_role in sources:
        for base,dirs,files in os.walk(srcroot,followlinks=False):
            relbase=Path(base).relative_to(srcroot);mappedbase=destprefix/relbase
            dirs[:]=[d for d in dirs if not (Path(base)/d).is_symlink() and may_enter(mappedbase/d)]
            for name in files:
                src=Path(base)/name;rel=mappedbase/name
                if not may_enter(rel) or name.endswith(('.lock','.tmp')):continue
                if src.is_symlink():raise RuntimeError('Symlink not copied: '+str(src))
                offline=rel.is_relative_to(allowed_private)
                # Offline synthetic files are byte-copied, never parsed or text-rebased.
                dst=replay/rel;dst.parent.mkdir(parents=True,exist_ok=True);oldhash=digest(src)
                if dst.exists():
                    if str(rel) not in inventories or inventories[str(rel)]['original_blake2b']!=oldhash:raise RuntimeError('Conflicting archive/dependency overlay: '+str(rel))
                    inventories[str(rel)]['also_present_in'].append(source_role);continue
                replacements=0
                if not offline and src.suffix.lower() in TEXT:
                    try:
                        with src.open('r',encoding='utf-8') as fi,dst.open('w',encoding='utf-8') as fo:
                            for line in fi:
                                replacements+=line.count(old);fo.write(line.replace(old,new))
                    except UnicodeDecodeError:shutil.copy2(src,dst);replacements=0
                else:shutil.copy2(src,dst)
                if not offline and 'config' in rel.parts and 'enable' in name.lower() and src.suffix=='.json':
                    data=read(dst)
                    if isinstance(data,dict) and data.get('enabled') is True:
                        data.update(enabled=False,transport_only_disabled=True,new_root_enable_required_for_physics=True)
                        dump(dst,data);disabled.append(str(rel))
                h=digest(dst)
                inventories[str(rel)]=dict(archive_path=str(rel),runtime_path=str(rel),
                    original_bytes=src.stat().st_size,original_blake2b=oldhash,
                    replay_bytes=dst.stat().st_size,replay_blake2b=h,byte_identical=oldhash==h,
                    root_prefix_replacements=replacements,evidence_role='OFFLINE_ONLY_SYNTHETIC_TRUTH_OR_GENERATOR_STATE' if offline else 'PUBLIC_SOURCE_OR_PREDICTION',
                    array_values_parsed=False,also_present_in=[source_role])
    study=replay/STUDY;protocol_path=study/'config/possibility_volume_protocol.json';protocol=read(protocol_path)
    original_source_identities=original_protocol['source_files'];runtime_sources=[]
    for row in protocol['source_files']:
        path=Path(row['path'])
        if not path.is_relative_to(replay) or not path.exists():raise RuntimeError('Missing transported dependency: '+str(path))
        st=path.stat();runtime_sources.append(dict(path=str(path),bytes=st.st_size,mtime_ns=st.st_mtime_ns,blake2b=digest(path)))
    protocol['source_files']=runtime_sources;stats_path=study/'config/possibility_statistics_frozen.json';st=stats_path.stat()
    protocol['statistics_config_identity']=dict(path=str(stats_path),bytes=st.st_size,mtime_ns=st.st_mtime_ns,blake2b=digest(stats_path))
    protocol['transport_only_relocation']=dict(original_protocol_blake2b=digest(source_study/'config/possibility_volume_protocol.json'),
        original_source_identities=original_source_identities,original_root=old,replay_root=new,
        scientific_values_and_seeds_unchanged=True,archived_approval_retained_for_deterministic_evaluation=True,
        physics_enablement_disabled=True,offline_synthetic_evidence_opt_in=include_offline)
    dump(protocol_path,protocol)
    seal_path=study/'results/possibility_volume/PREDICTIONS_SEALED.json';seal=read(seal_path)
    for parent in seal['parents']:
        for row in parent['files']:
            f=Path(row['path'])
            if not f.is_relative_to(replay) or not f.exists():raise RuntimeError('Missing replay prediction: '+str(f))
            row['bytes']=f.stat().st_size
    seal['source_files']=runtime_sources
    seal['transport_only_relocation']=dict(original_prediction_seal_blake2b=digest(source_study/'results/possibility_volume/PREDICTIONS_SEALED.json'),
        original_source_identities=original_source_identities,original_scientific_freeze_utc_preserved=True,
        new_prediction_or_fit=False,note='Transport derivative of existing before-truth seal')
    dump(seal_path,seal)
    replay_approval=None;archived_scoring=None
    if prepare_score_replay:
        score=study/'results/possibility_volume/offline_statistics'
        archived_scoring=score.with_name('offline_statistics_ARCHIVED_ORIGINAL_SCORING')
        if not score.is_dir() or archived_scoring.exists():raise RuntimeError('Cannot preserve copied original scoring safely')
        score.rename(archived_scoring)  # Only a NEW derivative is renamed; archive remains untouched.
        oldscore_rel=score.relative_to(replay);newscore_rel=archived_scoring.relative_to(replay)
        for item in inventories.values():
            rp=Path(item['runtime_path'])
            if rp.is_relative_to(oldscore_rel):item['runtime_path']=str(newscore_rel/rp.relative_to(oldscore_rel))
        approval=read(study/'config/possibility_volume_offline_score_approval.json')
        if approval.get('approved_by_root') is not True or approval.get('action')!='OFFLINE_TRUTH_SCORING_AFTER_ALL32_SEALED':
            raise RuntimeError('Original approved score authority missing')
        approval.update(utc=time.time(),role='OFFLINE_ONLY_DETERMINISTIC_SCORE_REPLAY',
            transport_replay_opt_in=True,original_approval_preserved_in_archive=True,new_physics=False,
            archived_scoring_preserved_at=str(archived_scoring),
            authorization_scope='Explicit relocation flags under standing authorization for full synthetic research evidence replay; no new model or fitting')
        replay_approval=study/'config/possibility_volume_offline_score_REPLAY_approval.json';dump(replay_approval,approval)
    refreshed={str(protocol_path.relative_to(replay)),str(seal_path.relative_to(replay))}
    for item in inventories.values():
        if item['runtime_path'] in refreshed:
            f=replay/item['runtime_path'];h=digest(f);item.update(replay_bytes=f.stat().st_size,replay_blake2b=h,byte_identical=item['original_blake2b']==h)
    inventory=list(inventories.values());offline=[r for r in inventory if r['evidence_role'].startswith('OFFLINE_ONLY')]
    if include_offline:
        truth_dir=replay/allowed_private/'truth'
        missing=[p['parent'] for p in original_protocol['parents'] if not (truth_dir/(p['parent']+'.npz')).exists()]
        if missing:raise RuntimeError('Opted-in archive lacks synthetic truth: '+str(missing))
    meta=dict(status='TRANSPORT_DERIVATIVE_COMPLETE_NOT_NEW_EXPERIMENT',utc=time.time(),
        archive_root=str(archive),archive_layout='REPOSITORY_RELATIVE' if repo_layout else 'GAUSSIAN_RELATIVE',
        dependency_root=None if dependency_root is None else str(Path(dependency_root).absolute()),
        original_host_root=old,replay_root=new,original_archive_mutated=False,
        blinded_synthetic_truth_opt_in=include_offline,offline_only_files_copied=len(offline),
        array_values_parsed=False,physics_calls=0,mapping=[dict(original_prefix=old,runtime_prefix=new)],
        original_protocol_blake2b=digest(source_study/'config/possibility_volume_protocol.json'),
        runtime_protocol_blake2b=digest(protocol_path),disabled_enable_files=disabled,
        archived_scoring_preserved_at=None if archived_scoring is None else str(archived_scoring),
        score_replay_approval=None if replay_approval is None else str(replay_approval),
        changed_file_count=sum(not r['byte_identical'] for r in inventory),files=inventory,
        scope='private/ here means blinded synthetic experiment evidence, not personal/user-private data. Offline role isolation remains; online code never opens it.',
        limitations=['No original file or scoring output is deleted. Only a new derivative scoring directory may be renamed.',
                     'Physics remains disabled; offline deterministic score replay needs no paid regeneration.',
                     'Original and runtime protocol identities intentionally differ and are both retained.'])
    dump(replay/'TRANSPORT_RELOCATION_MANIFEST.json',meta)
    dump(replay/'TRANSPORT_CHANGED_FILES.json',[r for r in inventory if not r['byte_identical']])
    dump(replay/'OFFLINE_ONLY_SYNTHETIC_EVIDENCE.json',offline)
    print(json.dumps({k:meta[k] for k in ['status','archive_layout','offline_only_files_copied','changed_file_count','score_replay_approval','physics_calls']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--archive-root',required=True,type=Path)
    p.add_argument('--dependency-root',type=Path);p.add_argument('--replay-root',required=True,type=Path)
    p.add_argument('--original-root',required=True);p.add_argument('--include-offline-synthetic-truth',action='store_true')
    p.add_argument('--prepare-score-replay',action='store_true');a=p.parse_args()
    relocate(a.archive_root,a.replay_root,a.original_root,a.dependency_root,a.include_offline_synthetic_truth,a.prepare_score_replay)
