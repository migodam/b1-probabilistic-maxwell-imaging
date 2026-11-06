"""Plan or build R8 release evidence; no upload, gitignore edits, solver, GT or SHA.

Default plans using stat only. After all runs finish and parent reviews seals:
  python delegated/cost_package/package_release.py --execute --finalized \
    --archive-dir /absolute/outside/R8/release
Archive members use repository-relative research/research_max_struct_r8 paths.
Existing source files are NEVER changed; compact endpoints are independent files.
"""
import argparse
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timezone
import json
import os
import fnmatch
import zipfile
import time

ROOT=Path(__file__).resolve().parents[2]
WORK=ROOT/'delegated/cost_package'
PREFIX='research/research_max_struct_r8/'
LIMIT=480_000_000  # 20 MB room under 500 MB for ZIP member metadata/overhead.
MAX_ZIP=500_000_000
GENERATED={'ARCHIVE_INDEX.json','ARCHIVE_INDEX.md','IGNORE_RULES.txt','PACKAGING_EXECUTION_RECEIPT.json','COMPACT_INDEX.json'}
KEYS=('m','material','chi','prediction','Y_prediction','Y','gradient','grad','physical_gradient','data_gradient','prediction6','prediction8')


def private_path(path):
    return any(any(w in part.lower() for w in ['truth','private','sealed_gt']) for part in Path(path).parts)


def exclude_reason(rel,include_scored_truth=False):
    p=Path(rel);parts=[v.lower() for v in p.parts];lo=rel.lower();name=p.name.lower()
    if private_path(rel) and not include_scored_truth:return 'HOLD_GT: no GT payload is opened or copied by this script'
    if any(x in parts for x in ['data_local','inputs']):return 'Original chat or third-party original measurements/input'
    if '.local.' in lo:return 'Local third-party fulltext/download'
    if any(x in parts for x in ['.git','__pycache__','remote_vendor','site-packages','.venv','node_modules','.ssh']):return 'Local dependency/cache'
    if name.endswith('.pdf') and parts[0]!='figures':return 'Third-party/local PDF outside owned figure outputs'
    if name.endswith(('.pyc','.pyo','.pem','.key','.pub','.ppk','.p12','.download','.epub','.eml','.exp')):return 'Third-party download or identity'
    if p.suffix.lower()=='.pdf' and parts[0]!='figures':return 'Third-party fulltext PDF; source-owned FIGURES PDFs are allowed'
    if name.startswith(('.env','id_rsa','id_ed25519','known_hosts','authorized_keys')) or any(x in lo for x in ['ssh_identity','ssh_config','credentials','credential.json','token.json']):return 'Credential/SSH artifact'
    if any(x in lo for x in ['rawchat','raw_chat','source_chat','chat_export','chat_history']):return 'Original chat'
    if parts[0]=='audit' and name.endswith(('.zip','.tgz','.tar','.gz')):return 'Temporary vendor/transfer/timing bundle'
    if p.suffix.lower()=='.html' and parts[0]=='delegated':return 'Downloaded third-party fulltext HTML'
    return None


def classify(rel,size,include_scored_truth=False):
    p=Path(rel);parts=[v.lower() for v in p.parts];lo=rel.lower();name=p.name.lower()
    denied=exclude_reason(rel,include_scored_truth)
    if denied:return ('hold' if denied.startswith('HOLD_GT') else 'exclude'),denied
    if rel=='delegated/cost_package/COST_LEDGER_FINAL.json':return 'git','Authoritative final cost metadata; may refresh after packaging; no sealed scientific source changed'
    if parts[:2]==['audit','remote_imports']:
        if p.suffix.lower() in ['.json','.jsonl','.log','.py','.md','.txt'] and size<=256000:return 'git','Small copied remote audit/source provenance mirror only; never new compute or duplicated scientific array'
        return 'exclude','Remote audit mirror payload outside small-text whitelist; canonical scientific RESULTS originals preserved separately'
    if private_path(rel):return 'release','SCORING_GATE_PASSED: own scored truth, separate reproducibility archive; public development after release'
    if parts[0]=='figures' and p.suffix.lower()=='.pdf' and size<=20_000_000:return 'git','Own scientific figure/gallery PDF; not a third-party paper'
    if parts[:3]==['delegated','cost_package','compact_endpoints']:return 'git','Independent compact result (original unchanged)'
    if p.suffix.lower() in ['.npz','.npy','.h5','.hdf5','.mat']:
        if parts[0]=='data_public' and 'operators' not in parts and size<=256_000:return 'git','Small public synthetic observation/initializer; original measured Y remains excluded'
        return 'release','Complete arrays/operators/training data and original state retained with original path'
    if p.suffix.lower() in ['.pt','.pth']:
        return ('git','Frozen learned weights') if size<=5_000_000 else ('release','Large checkpoint')
    if parts[0]=='data_training':return 'release','Complete supervised training dataset; manifests/receipts retained in Git by override below'
    if any(t in name for t in ['history','budget','corrector']) and p.suffix.lower() in ['.json','.jsonl','.csv','.tsv']:return 'release','Full iteration/budget/corrector trajectory'
    if size>2_000_000 and p.suffix.lower() not in ['.png','.jpg','.jpeg','.svg']:return 'release','Large original evidence; preserve complete bytes in release'
    return 'git','Code/config/receipts/scalar scores/figures/weights/compact public evidence'


def branch(rel):
    if private_path(rel):return 'scored_truth'
    parts=Path(rel).parts
    if parts[0]=='RESULTS' and len(parts)>1:return parts[1]
    if parts[0]=='data_training':return 'training_data'
    if parts[0]=='data_public':return 'physics_operators' if 'operators' in parts else 'public_data'
    if parts[:2]==('delegated','paper_baselines'):return 'paper_baselines'
    if parts[0]=='delegated':return 'delegated_'+parts[1]
    return parts[0]


def endpoint(rel):
    p=Path(rel);name=p.name.lower()
    return p.suffix.lower()=='.npz' and p.parts[0]=='RESULTS' and (
        name in ['prediction.npz','derived.npz','paid_incumbent.npz'] or name.endswith('_prediction.npz')) and 'auxiliary' not in name and not private_path(rel)


def inventory(root,include_scored_truth=False):
    rows=[]
    for base,dirs,names in os.walk(root,followlinks=False):
        dirs[:]=sorted(d for d in dirs if d not in ['.git','release_archives'])
        for name in sorted(names):
            p=Path(base)/name;rel=str(p.relative_to(root))
            if p.is_symlink():
                rows.append(dict(path=rel,bytes=p.lstat().st_size,category='exclude',reason='Symlink not followed'));continue
            # Current generated index/receipt are included in Git handoff later,
            # never recursively archived into themselves or counted as science.
            if p.parent==WORK and name in GENERATED:continue
            size=p.stat().st_size;cat,reason=classify(rel,size,include_scored_truth)
            if p.suffix.lower()=='.json' and p.parts[-2] in ['data_training','paper_domain'] and 'manifest' in name:
                cat,reason='git','Supervised data manifest and provenance'
            rows.append(dict(path=rel,bytes=size,category=cat,reason=reason,branch=branch(rel),compact_candidate=endpoint(rel)))
    return rows


def manifests(root):
    # Explicit public task sources. Do not inspect private/scoring metadata.
    result=[]
    paths=sorted((root/'configs').glob('warm_*manifest.json'))+sorted((root/'experiments').glob('history*manifest.json'))
    for p in paths:
        d=json.loads(p.read_text());tasks=d if isinstance(d,list) else d.get('tasks',[])
        if isinstance(tasks,dict):tasks=list(tasks.values())
        tasks=[dict(task_id=t.get('task_id'),output_dir=t.get('output_dir')) for t in tasks]
        result.append(dict(path=str(p.relative_to(root)),tasks=len(tasks),entries=tasks))
    # Additional branch evidence discovered by explicit result root; counts are
    # snapshot metadata, never assumed completion of a fixed historic number.
    branches={}
    for name in ['warm_primary','warm_development','warm_anchor_development','warm_anchor_development_extended','warm_anchor_confirmation','history_warm','paper_benchmarks','vector3d','fresnel']:
        base=root/'RESULTS'/name
        branches[name]=dict(receipts=len(list(base.rglob('*RECEIPT.json'))) if base.exists() else 0,
                            array_files=len(list(base.rglob('*.npz'))) if base.exists() else 0)
    return result,branches


def packs(rows):
    output=[]
    for group in sorted(set(r.get('branch') for r in rows if r['category']=='release')):
        items=sorted((r for r in rows if r['category']=='release' and r['branch']==group),key=lambda r:r['path'])
        batch=[];size=0;idx=1
        for r in items:
            if r['bytes']>LIMIT:raise ValueError('Single source exceeds 480 MB; retain original and implement explicit chunk/restoration contract before packaging: '+r['path'])
            if batch and size+r['bytes']>LIMIT:
                output.append(dict(name=f'max-struct-r8-{group}-part{idx:03d}.zip',uncompressed_bytes=size,members=batch));idx+=1;batch=[];size=0
            batch.append(dict(path=r['path'],archive_member=PREFIX+r['path'],bytes=r['bytes']));size+=r['bytes']
        if batch:output.append(dict(name=f'max-struct-r8-{group}-part{idx:03d}.zip',uncompressed_bytes=size,members=batch))
    return output


def write_plan(root,rows,include_scored_truth=False):
    ms,bs=manifests(root);planned=packs(rows);groups=defaultdict(lambda:dict(files=0,bytes=0))
    for r in rows:groups[r['category']]['files']+=1;groups[r['category']]['bytes']+=r['bytes']
    plan=dict(status='PLAN_ONLY_NOT_PACKAGED_NOT_UPLOADED',created_utc=datetime.now(timezone.utc).isoformat(),
        root_prefix=PREFIX,all_sources_unchanged=True,GT_payload_reads=0,hash_checks='NONE',include_scored_truth=include_scored_truth,
        truth_release_gate='Original750/42 seal+scoring+750scores AND separate126/42 anchor seal+confirmation scoring+126scores; all complete. DefaultHOLD.',
        archive_limit_bytes=MAX_ZIP,planned_source_limit_bytes=LIMIT,groups=dict(groups),
        source_task_manifests=ms,explicit_result_roots=bs,files=rows,archives=planned,
        restore='Download all archives, extract at repository root in a fresh copy. Every original release file retains its exact repository-relative member path. Do not overwrite existing sealed files.',
        compact_policy=dict(location='delegated/cost_package/compact_endpoints/<original-relative-path>',keys=list(KEYS),
            timing='Created only by --execute --finalized; plan mode does not open arrays',
            caveat='Predictions/material/gradients only; gradients copied only if already present. No new physics evaluation or invented zero gradient. Fresnel includes predicted Y, never measured raw Y.'),
        excluded_not_claimed_recoverable='Original chats, third-party downloads/raw measurements, credentials/local dependencies and HOLD truth are not release members. Official download/retrieval code is retained.')
    WORK.mkdir(parents=True,exist_ok=True)
    (WORK/'ARCHIVE_INDEX.json').write_text(json.dumps(plan,indent=2,ensure_ascii=False)+'\n')
    # Anchored literal file paths can be pasted into Git ignore after root review;
    # do not modify an actual .gitignore here. Explicit paths match this snapshot.
    ignore=['# PROPOSAL ONLY; apply after final refresh/root review. Paths relative to R8 directory.',
            '# R1-R7 unaffected. Every release-only original is listed in ARCHIVE_INDEX.json.',
            '# Compact endpoints, code, receipts, scalar scores, figures and frozen weights stay in Git.']
    ignore += ['/'+r['path'] for r in rows if r['category'] in ['release','exclude','hold']]
    (WORK/'IGNORE_RULES.txt').write_text('\n'.join(ignore)+'\n')
    lines=['# R8发布档案计划','',f"Git候选 {groups['git']['files']} 文件 / {groups['git']['bytes']} bytes；release {groups['release']['files']} 文件 / {groups['release']['bytes']} bytes。",'',
        '默认计划只stat原文件与读取公开manifest；未读数组/GT。执行结果以PACKAGING_EXECUTION_RECEIPT为准；所有数值都是当前快照，最终完成后重跑刷新。','',
        '|包|成员数|非压缩字节|','|---|---:|---:|']
    lines += [f"|{p['name']}|{len(p['members'])}|{p['uncompressed_bytes']}|" for p in planned]
    lines += ['','完整原prediction(j/Q等)、轨迹/预算/checkpoints、训练数据及physics算子保留于release原路径。Git的compact为独立新文件，不改写源文件，不参与封存。默认所有private/truth路径HOLD且不读取。显式--include-scored-truth仅在750/42与126/42两套完整seal、各自评分回执和750/126行SCORES都通过后把自有GT单独分包；原Chat和第三方raw永远排除。Fresnel原Y只从官方重取代码取得，不分发原测量文件。','',
        '最终执行：`python delegated/cost_package/package_release.py --execute --finalized --archive-dir /absolute/outside/R8/release`。',
        '执行器拒绝覆盖原档案/compact；zip生成后只核对成员路径、文件数、尺寸和CRC元数据，不执行SHA。恢复时在新目录按repo相对路径解压。',
        'IGNORE_RULES.txt为本轮明确文件路径建议；不是已修改的gitignore。新实验出现后必须刷新。']
    (WORK/'ARCHIVE_INDEX.md').write_text('\n'.join(lines)+'\n')
    return plan


def execute(root,plan,archive_dir):
    import numpy as np
    archive_dir=archive_dir.resolve()
    if archive_dir==root or root in archive_dir.parents:raise ValueError('Archive directory must be outside R8 to avoid recursive packaging')
    archive_dir.mkdir(parents=True,exist_ok=True)
    for p in plan['archives']:
        if (archive_dir/p['name']).exists():raise FileExistsError('Archive exists; choose fresh destination: '+p['name'])
    targets=[r for r in plan['files'] if r.get('compact_candidate') and r['category'] not in ['hold','exclude']]
    for r in targets:
        target=WORK/'compact_endpoints'/r['path']
        if target.exists():raise FileExistsError('Compact exists, original untouched; choose fresh copy or preserve existing compact: '+str(target))
    begin,cpu=time.perf_counter(),time.process_time();compacts=[]
    for r in targets:
        source=root/r['path'];target=WORK/'compact_endpoints'/r['path'];before=source.stat()
        if private_path(r['path']):raise ValueError('GT prohibited')
        with np.load(source,allow_pickle=False) as z:
            if any(any(w in key.lower() for w in ['truth','ground_truth','m_true','m_truth']) for key in z.files):
                raise ValueError('Unexpected truth-key in endpoint; no payload read: '+r['path'])
            arrays={key:z[key] for key in KEYS if key in z.files}
        if not arrays:continue
        target.parent.mkdir(parents=True,exist_ok=True)
        with target.open('xb') as stream:np.savez_compressed(stream,**arrays)
        after=source.stat()
        if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):raise RuntimeError('Source mutated during compact creation: '+r['path'])
        compacts.append(dict(source=r['path'],compact=str(target.relative_to(root)),keys=list(arrays),bytes=target.stat().st_size,shape={k:list(v.shape) for k,v in arrays.items()},
            scope='Copied existing endpoint values only; no GT/physical/scoring calls; source file never opened for write'))
    (WORK/'COMPACT_INDEX.json').write_text(json.dumps(dict(status='INDEPENDENT_COMPACTS_CREATED',records=compacts,solver_invocations=0,GT_payload_reads=0),indent=2)+'\n')
    built=[]
    for p in plan['archives']:
        dest=archive_dir/p['name'];member_expected={m['archive_member']:m['bytes'] for m in p['members']}
        with zipfile.ZipFile(dest,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=3,allowZip64=True) as z:
            for member in p['members']:
                if (private_path(member['path']) and not plan['include_scored_truth']) or exclude_reason(member['path'],plan['include_scored_truth']):raise ValueError('Excluded payload attempted: '+member['path'])
                source=root/member['path'];before=source.stat()
                if before.st_size!=member['bytes']:raise RuntimeError('Source changed since plan; rerun after completion: '+member['path'])
                z.write(source,member['archive_member'])
                after=source.stat()
                if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):raise RuntimeError('Source mutated during packaging: '+member['path'])
        if dest.stat().st_size>MAX_ZIP:raise RuntimeError('Archive exceeds 500MB: '+str(dest))
        with zipfile.ZipFile(dest) as z:
            actual={i.filename:i.file_size for i in z.infolist()}
            if actual!=member_expected:raise RuntimeError('Archive member/size mismatch '+str(dest))
        built.append(dict(name=p['name'],path=str(dest),bytes=dest.stat().st_size,members=len(member_expected),uncompressed_bytes=p['uncompressed_bytes']))
    receipt=dict(status='PACKAGED_NOT_UPLOADED',finished_utc=datetime.now(timezone.utc).isoformat(),archives=built,compact_endpoints=len(compacts),cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-begin,
        original_files_unchanged=True,GT_payload_reads='ARCHIVE_BYTE_COPY_ONLY_AFTER_SEAL_AND_SCORE' if plan['include_scored_truth'] else 0,include_scored_truth=plan['include_scored_truth'],hash_checks='NONE',physical_calls=0,scoring_calls=0,gitignore_changed=False,
        validation='ZIP member paths and uncompressed file sizes checked; no source SHA or content equality claim; source stat unchanged before/after read',
        excluded=[dict(path=r['path'],reason=r['reason']) for r in plan['files'] if r['category'] in ['exclude','hold']])
    (WORK/'PACKAGING_EXECUTION_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    # Add new compact files to refreshed Git plan; original archive partition is
    # stable because generated receipt/index files are not recursive evidence.
    refreshed=write_plan(root,inventory(root,plan['include_scored_truth']),plan['include_scored_truth']);refreshed['status']='PACKAGED_NOT_UPLOADED';refreshed['built_archives']=built;refreshed['GT_payload_reads']=receipt['GT_payload_reads']
    (WORK/'ARCHIVE_INDEX.json').write_text(json.dumps(refreshed,indent=2,ensure_ascii=False)+'\n')
    return receipt


def scored_truth_gate(root):
    seal=root/'audit/WARM_ALL_PUBLIC_SEAL.json';receipt=root/'RESULTS/warm_scoring/RECEIPT.json';scores=root/'RESULTS/warm_scoring/SCORES.json'
    for p in [seal,receipt,scores]:
        if not p.is_file():raise ValueError('--include-scored-truth denied: missing public seal/scoring file '+str(p))
    se=json.loads(seal.read_text());re=json.loads(receipt.read_text());sc=json.loads(scores.read_text())
    if se.get('status')!='ALL_FROZEN_WARM_TASKS_SEALED' or se.get('tasks')!=750 or se.get('objects')!=42 or len(se.get('entries',[]))!=750:
        raise ValueError('--include-scored-truth denied: all750/42 public prediction seal incomplete')
    if re.get('status')!='SCORED_ALL_FROZEN_PREDICTIONS' or re.get('tasks')!=750 or re.get('objects')!=42:
        raise ValueError('--include-scored-truth denied: all750/42 scoring receipt incomplete')
    if not isinstance(sc,dict) or not isinstance(sc.get('rows'),list) or len(sc['rows'])!=750:
        raise ValueError('--include-scored-truth denied: SCORES.json must contain750 rows')
    anchor_seal=root/'audit/WARM_ANCHOR_PUBLIC_SEAL.json'
    anchor_receipt=root/'delegated/warm_anchor_analysis/CONFIRMATION_RECEIPT.json'
    anchor_scores=root/'delegated/warm_anchor_analysis/CONFIRMATION_SCORES.json'
    for p in [anchor_seal,anchor_receipt,anchor_scores]:
        if not p.is_file():raise ValueError('--include-scored-truth denied: missing separate126 anchor public seal/scoring '+str(p))
    ase=json.loads(anchor_seal.read_text());are=json.loads(anchor_receipt.read_text());asc=json.loads(anchor_scores.read_text())
    if ase.get('status')!='ALL_ANCHOR_TASKS_SEALED' or ase.get('tasks')!=126 or ase.get('objects')!=42 or len(ase.get('entries',[]))!=126:
        raise ValueError('--include-scored-truth denied: separate126/42 anchor prediction seal incomplete')
    if are.get('status')!='SEPARATE126_ANCHOR_OFFLINE_SCORED' or are.get('tasks')!=126 or are.get('objects')!=42:
        raise ValueError('--include-scored-truth denied: separate126/42 anchor confirmation scoring incomplete')
    if not isinstance(asc,dict) or not isinstance(asc.get('rows'),list) or len(asc['rows'])!=126:
        raise ValueError('--include-scored-truth denied: anchor CONFIRMATION_SCORES must contain126 rows')


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--include-scored-truth',action='store_true',help='After both complete750/42 and126/42 seals and scorings only, publish own truth as separate reproducibility package; defaultHOLD');ap.add_argument('--execute',action='store_true');ap.add_argument('--finalized',action='store_true');ap.add_argument('--archive-dir',type=Path);a=ap.parse_args()
    if a.execute and (not a.finalized or a.archive_dir is None):ap.error('--execute requires --finalized and --archive-dir after root review of completed runs')
    if a.include_scored_truth:scored_truth_gate(ROOT)
    plan=write_plan(ROOT,inventory(ROOT,a.include_scored_truth),a.include_scored_truth)
    if a.execute:
        receipt=execute(ROOT,plan,a.archive_dir);print(json.dumps(dict(status=receipt['status'],archives=len(receipt['archives']),compact_endpoints=receipt['compact_endpoints'])))
    else:print(json.dumps(dict(status=plan['status'],groups=plan['groups'],archives=len(plan['archives']),compact_candidates=sum(r.get('compact_candidate',False) for r in plan['files']))))

if __name__=='__main__':main()
