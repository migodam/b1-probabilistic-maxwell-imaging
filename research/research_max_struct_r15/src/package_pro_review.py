#!/usr/bin/env python3
"""Plan/build finite R15 Pro review ZIP, complementary to full raw releases.

Run --plan first. --build is for the parent only after final docs/key states are
ready. Stdlib only; stat/JSON/source imports are inspected, scientific modules
and numeric NPZ/GT payloads are never loaded. ZIP build copies original bytes.
"""
from __future__ import annotations
import argparse,ast,csv,json,os,re,time,zipfile
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TEXT={'.json','.csv','.md','.py','.txt','.yml','.yaml','.html'}
PUBLIC_EXT=TEXT|{'.png','.pdf','.svg','.npz','.npy'}
BLOCKED={'private_operational','secrets','credentials','chat_readback','__pycache__','.git','fulltexts','fulltext'}
EXCLUDED_NAMES={'TASK_RESUME.md','CANDIDATES.npz','CANDIDATE_POOL.npz','HISTORY.json','NATIVE_HISTORY.json','CYCLES.jsonl','PAID_HISTORY.csv','PAID_HISTORY.json','PAID_CANDIDATE_HISTORY.csv','PAID_CANDIDATE_HISTORY.json','NATIVE_STAGE_HISTORY.csv','NATIVE_STAGE_HISTORY.json','NATIVE_STEPS.csv','NATIVE_STEPS.json'}
REQUIRED_DOCS=('RESEARCH_REPORT_ZH.md','THEORY_DERIVATIONS.md','THEORY_STATUS.md','PRO_HANDOFF.md','REPRODUCE.md','COST_LEDGER.md','CONFIG_PRECEDENCE.md')
SCORING_MEMBERS=('research/research_max_struct_r14/runtime/scoring_history/src/score_results.py',
 'research/research_max_struct_r14/runtime/scoring_history/src/seal_results.py',
 'research/research_max_struct_r14/src/score_dualmetric.py','research/research_max_struct_r14/src/score_extensions.py',
 'research/research_max_struct_r14/supplement_observable19/src/analyze_results.py')
PORTABLE_READER='''"""Single-threaded OFFLINE reader; normalizes paths in memory, never NPZ bytes.
No science runs at import. Inspecting/restoring states is an explicit user action.
"""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'runtime'),str(ROOT.parents[1]),str(ROOT.parents[1]/'vendor')]

def load_complete_state(path):
    import dualmetric_core as codec
    state=codec._load_state_npz(Path(path))
    for key in ('fixed_cache_relative','background_cache_relative'):
        if isinstance(state.get('r15',{}).get(key),str):state['r15'][key]=state['r15'][key].replace('\\\\','/')
    return state

def restore_portable(solver_class,model,Y,config,path,**kwargs):
    # Temporary reader-side adapter for the original classmethod. Single-threaded
    # offline use only; restores the codec immediately and changes no source/state.
    import dualmetric_core as codec
    original=codec._load_state_npz
    def reader(p):
        state=original(p)
        for key in ('fixed_cache_relative','background_cache_relative'):
            if isinstance(state.get('r15',{}).get(key),str):state['r15'][key]=state['r15'][key].replace('\\\\','/')
        return state
    codec._load_state_npz=reader
    try:return solver_class.from_checkpoint(model,Y,config,path,**kwargs)
    finally:codec._load_state_npz=original
'''


def dump(path,obj):path.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n',encoding='utf8')


def walk(folder):
    if not folder.is_dir() or folder.is_symlink():return
    for parent,dirs,names in os.walk(folder,followlinks=False):
        here=Path(parent);dirs[:]=[n for n in dirs if n.lower() not in BLOCKED and not (here/n).is_symlink()]
        for name in sorted(names):
            p=here/name
            if not p.is_symlink() and p.is_file():yield p


def collect(root):
    W=root.parents[1];prefix=str(root.relative_to(W))+'/'
    selected={};missing=[];excluded=[];closure_issues=[]
    truth_allow={W/f'research/research_max_struct_r14/sealed_truth/independent24/r14_new_{i:04d}.npz' for i in range(24)}
    def add(path,role,member=None,required=False):
        original=Path(path)
        if original.is_symlink():return
        p=original.resolve()
        if not p.is_file():
            if required:missing.append(dict(path=str(p),role=role))
            return
        if not p.is_relative_to(W):raise ValueError('Outside external worktree: '+str(p))
        relative=p.relative_to(W)
        if any(part.lower() in BLOCKED for part in relative.parts):return
        if 'sealed_truth' in [x.lower() for x in relative.parts] and p not in truth_allow:return
        if p.is_symlink():return
        if p.name in EXCLUDED_NAMES or p.suffix.lower() not in PUBLIC_EXT:return
        if p.suffix.lower() in TEXT and p.stat().st_size>32_000_000:
            excluded.append(dict(path=str(p),bytes=p.stat().st_size,reason='Noncompact text >32MB; full release preserves it'));return
        target=member or str(relative)
        if target in selected:
            if selected[target][0]!=p:raise ValueError('Conflicting archive source mapping: '+target)
            return
        selected[target]=(p,role)
    for p in root.glob('*.md'):
        if p.name!='TASK_RESUME.md':add(p,'final author document snapshot')
    for name in ('ORIGINAL_REQUEST_ZH.txt','DOWNLOADS.json'):add(root/name,'request/download metadata snapshot; publication verification not asserted')
    for name in REQUIRED_DOCS:
        if not (root/name).is_file():missing.append(dict(path=str(root/name),role='required final review document'))
    for folder,role in [('src','R15 source'),('configs','frozen configuration'),('runtime','frozen R15 runtime'),('data_public/independent24','public observations24'),('audit/rank_all24','saved public real/complex geometry; no new SVD'),('predictions','original and approved prediction seals')]:
        for p in walk(root/folder):add(p,role)
    for p in walk(root/'analysis'):
        if p.suffix.lower() in TEXT|{'.png','.pdf','.svg'}:add(p,'compact scored analysis/diagnostic evidence')
    for p in walk(root/'FIGURES'):
        if p.suffix.lower() in TEXT|{'.png','.pdf','.svg'}:add(p,'final or explicitly archived figure/receipt')
    audit_folders=('old_endpoint_geometry_v1','paid_endpoint_diagnostic_v1','readout_counterfactual_v1','denominator_wiring','fullstate_protocol','packaging_failures','publication','analysis_attempts')
    for folder in audit_folders:
        for p in walk(root/'audit'/folder):
            # Same public inputs already included at canonical data_public path.
            if folder=='old_endpoint_geometry_v1' and 'public_inputs' in p.relative_to(root/'audit'/folder).parts:continue
            role='OFFLINE_GEOMETRY_WITH_AUTHORIZED_GT_DERIVED_ARRAYS' if folder=='old_endpoint_geometry_v1' and p.suffix=='.npz' else 'scoped audit/diagnostic evidence'
            add(p,role)
    for name in ('OLD_ENDPOINT_GEOMETRY_AUTHORIZATION.json','MANIFEST_FREEZE_twofold_exposed24_v1.json','PUBLIC_ACQUISITION_SUMMARY.json','ACTUAL_LEGACY_SOURCE_COMPARISON.json','INHERITED_SOURCES_AND_PUBLIC_INPUTS.json','PLUMBING_SOURCE_COPIES.json','LOCAL_ENVIRONMENT.json','RUNNER_IMPLEMENTATION.json'):
        add(root/'audit'/name,'public provenance/contract receipt')
    # Full selected endpoint views, not full paid candidate pools or native cycles.
    for task in sorted((root/'RESULTS/twofold_exposed24_v1').glob('*')):
        if not task.is_dir() or task.is_symlink():continue
        for view in ('initial','raw','corrected','primary','secondary'):
            for ext in ('.json','.npz'):add(task/(view+ext),'selected saved material/data view; not necessarily full J/Q/PR')
        for name in ('RECEIPT.json','ACTUAL_CONFIG.json','TWOFOLD_SETUP.json','WRAPPER_PROVENANCE.json'):
            add(task/name,'actual run/setup/paid-cost receipt')
    keydir=root/'archives/pro_key_states_v1'
    manifest_path=keydir/'R15_PRO_KEY_FULLSTATES_V1_MANIFEST.json'
    manifest=json.loads(manifest_path.read_text()) if manifest_path.is_file() else None
    if manifest is None:missing.append(dict(path=str(manifest_path),role='selected complete key-state extraction manifest'))
    else:
        for item in manifest['members']:
            p=root/Path(item['member'].replace('\\','/'))
            add(p,'selected complete key state/cache: '+', '.join(item.get('roles',[])),required=True)
            if p.is_file() and p.stat().st_size!=item['bytes']:closure_issues.append(dict(path=str(p),issue='KEY_MEMBER_BYTES_DIFFER_FROM_RESTORED_MANIFEST',expected=item['bytes'],actual=p.stat().st_size))
    for name in ('SUMMARY.json','RESTORE.md','MEMBERS.csv','R15_PRO_KEY_FULLSTATES_V1_MANIFEST.json','R15_PRO_KEY_FULLSTATES_V1_RECEIPT.json','LOCAL_PORTABLE_DEPENDENCY_RECEIPT.json','LOCAL_RESTORE_RECEIPT.json'):
        add(keydir/name,'key-state subset/cache closure provenance; outer ZIP intentionally omitted')
    portable=keydir/'LOCAL_PORTABLE_DEPENDENCY_RECEIPT.json'
    references=json.loads(portable.read_text()).get('references',[]) if portable.is_file() else []
    for item in references:
        checkpoint=root/Path(item['checkpoint'].replace('\\','/'))
        normalized=item['original_relative'].replace('\\','/')
        target=(checkpoint.parent/normalized).resolve()
        archive_path=str(target.relative_to(W)) if target.is_relative_to(W) else None
        if not target.is_file() or archive_path not in selected:closure_issues.append(dict(checkpoint=str(checkpoint),original_relative=item['original_relative'],normalized_target=str(target),issue='CACHE_NOT_SELECTED_OR_MISSING'))
    # Freeze19-source runtime snapshot both in provenance location and canonical
    # import path. All18 project sources matched full checkout in prior receipt.
    overrides={}
    snapshot=root/'sources/legacy_runtime_snapshot/legacy_sources'
    for p in walk(snapshot):
        if p.suffix=='.py':
            canonical=str(p.relative_to(snapshot));overrides[canonical]=p
            add(p,'actual inherited runtime snapshot')
            add(p,'actual inherited runtime mounted at canonical import path',member=canonical)
    for p in (root/'sources/legacy_runtime_snapshot').glob('*.json'):add(p,'legacy runtime provenance')
    # R15 analyzers dynamically import these original sibling helpers.
    for name in SCORING_MEMBERS:add(W/name,'required sibling scoring/analyzer dependency',required=True)
    for p in (W/'research/research_max_struct_r14/runtime/scoring_history/src').glob('*.py'):add(p,'original scoring-history helper')
    add(W/'research/research_max_struct_r14/runtime/scoring_history/PROVENANCE.json','original scoring provenance')
    add(W/'research/research_max_struct_r14/configs/independent24_v1.json','frozen parent config reference')
    # Source-only AST closure. No scientific Python import/execution occurs.
    parsed=set();pending=[(member,p) for member,(p,role) in selected.items() if p.suffix=='.py'];unresolved=set()
    def module_source(module):
        canonical='/'.join(module.split('.'))+'.py'
        if canonical in overrides:return overrides[canonical],canonical
        q=W/canonical
        if q.is_file():return q,canonical
        init=W/Path(*module.split('.'))/'__init__.py'
        if init.is_file():return init,str(init.relative_to(W))
        return None,None
    while pending:
        member,p=pending.pop()
        if p in parsed:continue
        parsed.add(p)
        tree=ast.parse(p.read_text(encoding='utf8'))
        for node in ast.walk(tree):
            modules=[]
            if isinstance(node,ast.Import):modules=[a.name for a in node.names]
            elif isinstance(node,ast.ImportFrom) and node.module:
                modules=[node.module]+[node.module+'.'+a.name for a in node.names if a.name!='*']
                if node.level:
                    q=p.parent/(node.module.replace('.','/')+'.py')
                    if q.is_file():
                        add(q,'relative Python source dependency');pending.append((str(q.relative_to(W)),q))
                    continue
            for module in modules:
                q,target=module_source(module) if module.startswith('research.') else (None,None)
                if q is None and '.' not in module:
                    candidates=[p.parent/(module+'.py'),root/'src'/(module+'.py'),root/'runtime'/(module+'.py'),W/'research/research_max_struct_r14/src'/(module+'.py'),W/'research/research_max_struct_r14/runtime/scoring_history/src'/(module+'.py')]
                    q=next((x for x in candidates if x.is_file()),None);target=str(q.relative_to(W)) if q else None
                if q:
                    add(q,'required source import dependency',member=target);pending.append((target,q))
                elif module.startswith('research.'):
                    parent_module=module.rsplit('.',1)[0]
                    parent_source,_=module_source(parent_module)
                    if not (W/Path(*module.split('.'))).is_dir() and parent_source is None:unresolved.add(module)
    # Explicit release-only JSON list overrides Git publication layout. No git
    # ignore query is used, and old failed/partial analysis cost records remain.
    release_list=root/'audit/publication/RELEASE_ONLY_TABLES.json'
    release_items=json.loads(release_list.read_text()).get('items',[]) if release_list.is_file() else []
    for item in release_items:
        target=root/item['path'];add(target,'explicit release-only JSON; authoritative offline supplement',required=True)
        if target.is_file() and target.stat().st_size!=item['bytes']:closure_issues.append(dict(path=str(target),issue='RELEASE_ONLY_LIST_SIZE_CHANGED',expected=item['bytes'],actual=target.stat().st_size))
    release_missing=[item['path'] for item in release_items if str((root/item['path']).relative_to(W)) not in selected]
    closure_issues.extend(dict(path=x,issue='REQUIRED_RELEASE_ONLY_JSON_NOT_SELECTED') for x in release_missing)
    # Parent explicitly authorized opaque24 GT inclusion only for offline scoring.
    for p in sorted(truth_allow):add(p,'OFFLINE_SCORING_ONLY: already-sealed exposed24 material truth; never online input',required=True)
    # Actual canonical full-checkout path; do not invent nonexistent analyze_independent.py.
    unavailable_optional=[]
    candidate=W/'research/research_max_struct_r14/src/analyze_independent.py'
    if candidate.is_file():add(candidate,'optional actual parent analyzer')
    else:unavailable_optional.append(str(candidate))
    entries=[]
    for member,(p,role) in sorted(selected.items()):
        stat=p.stat();entries.append(dict(member=member,source=str(p),bytes=stat.st_size,mtime_ns=stat.st_mtime_ns,role=role,
            compression='ZIP_DEFLATED_level1' if p.suffix.lower() in TEXT else 'ZIP_STORED'))
    return entries,dict(missing=missing,cache_closure_issues=closure_issues,unresolved_project_imports=sorted(unresolved),unavailable_optional_paths=unavailable_optional,
        selected_key_members=len(manifest['members']) if manifest else 0,cache_reference_records=len(references),
        release_only_JSON_required=len(release_items),release_only_JSON_selected=len(release_items)-len(release_missing),git_ignore_policy='Ignored for selected package payloads; Git/Pages layout never filters review JSON',offlinescoring_truth_members=24,excluded=excluded,required_legacy_mounts=sorted(overrides),
        source_validation='AST and path/byte-count only; no module import, physics, scoring, numeric NPZ/GT reads or hashes')


def start_here(plan):
    return '''# R15 Pro finite review packet

This is a selected review snapshot, complementary to the16 full raw release archives. It is not the full candidate pool, every checkpoint, every cycle or an executed replay validation. Source and every original selected NPZ byte are preserved. Extract the ZIP into an empty directory: research/research_max_struct_r15/ is the study root; required legacy/sibling research/... and vendor paths preserve the original import tree.

Read PRO_HANDOFF.md, RESEARCH_REPORT_ZH.md, THEORY_DERIVATIONS.md, THEORY_STATUS.md and CONFIG_PRECEDENCE.md first. Current file availability and exact bytes are in MANIFEST.json/CSV. No missing document is represented as completed. Original paid algorithms, priors, selectors, predictions and costs are not changed by packaging.

Selected material/data views are not automatically full optimizer states. The80 key-state members selected by archives/pro_key_states_v1/R15_PRO_KEY_FULLSTATES_V1_MANIFEST.json supply complete initial/final fork and main-final/checkpoint512 state/cache subsets only. Paid-endpoint diagnostic arrays separately provide selected physical m/J/Q/gradients. Old endpoint geometry arrays contain explicitly authorized offline GT-derived materials/errors; they are NOT online inputs.

OFFLINE_SCORING_ONLY: the24 original sibling R14 sealed_truth/independent24 NPZ files are copied opaquely under their original paths after current predictions were sealed/scored. They are included solely to permit authorized offline scoring/inspection. Never expose them or cached GT/reference arrays to an online inverse method, candidate selector, new R15 prediction or tuning. Geometry reference/material errors are likewise offline-only. The packager only copies truth files opaquely; it never loads or interprets numeric truth arrays.

Keep original Windows prediction seals intact. Absolute Windows/source-volume receipt paths describe provenance; local replay needs an explicit relocation map/new local seal if the scorer requires local inventory. Copying this packet does not make those remote absolute paths portable or authorize new scoring. See REPRODUCE.md for actual command/environment requirements and the key-state RESTORE.md for subset choices.

Windows relative cache strings inside original checkpoint NPZs are unchanged. On POSIX normalize backslashes to slash only in the reader. portable_state_reader.py provides a single-threaded in-memory adapter and complete-state loader; it never rewrites NPZ bytes or source files. All nonnull key-state cache targets are checked in the selected-file plan. This package does not construct missing currents/PR history from m-only views, regenerate SVDs, or guarantee independent replay of unselected historical tests. Full analyze_campaign paid-pool/common-prefix recomputation needs the omitted original HISTORY/CANDIDATES and full release; this review packet provides selected-view metric imports and saved analysis, not all-pool replay. The historical wiring check's smaller rank_characterization cache is omitted in favor of the actual all24 caches; full raw release retains the former.

Restoration/source import closure was inspected statically, not executed in a fresh environment. Third-party installed libraries (NumPy/SciPy/Matplotlib/psutil etc.) remain environment requirements. Required frozen19-source legacy files are included in provenance form and at canonical original import paths. Large full paid pools and whole raw archives are omitted intentionally; get their actual URLs from the owner's final DOWNLOADS.json/release inventory. Publication links will be supplied later by the owner through Git/release metadata; this snapshot makes no verified-download claim. No release URLs are fabricated and no upload/publication is performed.
'''


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT)
    mode=p.add_mutually_exclusive_group(required=True);mode.add_argument('--plan',action='store_true');mode.add_argument('--build',action='store_true',help='Parent only, after final docs/key-state readiness')
    a=p.parse_args();root=a.root.resolve();W=root.parents[1];cpu,wall=time.process_time(),time.perf_counter()
    if not str(W).startswith('/Volumes/'):raise ValueError('External-volume workspace/output required; no internal-disk staging')
    out=root/'archives/pro_review_v1';out.mkdir(parents=True,exist_ok=True)
    entries,info=collect(root)
    stats=Counter();roles=Counter()
    for e in entries:stats[Path(e['member']).suffix.lower()]+=e['bytes'];roles[e['role']]+=e['bytes']
    plan=dict(schema='R15_PRO_REVIEW_EXPLICIT_FILES_V1',status='PLAN_ONLY_NO_ZIP_WRITTEN' if a.plan else 'BUILD_REQUESTED',snapshot_utc=datetime.now(timezone.utc).isoformat(),
        root=str(root),files=len(entries),source_bytes=sum(e['bytes'] for e in entries),binary_stored_bytes=sum(e['bytes'] for e in entries if e['compression']=='ZIP_STORED'),
        bytes_by_suffix=dict(stats),bytes_by_role=dict(roles),members=entries,**info,
        generated_members=[str((out/name).relative_to(W)) for name in ('MANIFEST.json','MANIFEST.csv','START_HERE.md','portable_state_reader.py')],
        completeness='Finite selected review subset plus offline scoring truth, not full raw archives or all candidates/states',
        numeric_GT_array_loads=0,science_runs=0,hash_operations='NONE',archive_target='~200-400MB intent; actual size reported only after build, never fabricated')
    dump(out/'MANIFEST.json',plan)
    with (out/'MANIFEST.csv').open('w',newline='',encoding='utf8') as f:
        writer=csv.DictWriter(f,fieldnames=list(entries[0]));writer.writeheader();writer.writerows(entries)
    (out/'START_HERE.md').write_text(start_here(plan),encoding='utf8')
    (out/'portable_state_reader.py').write_text(PORTABLE_READER,encoding='utf8')
    if a.plan:
        plan.update(CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
        dump(out/'PLAN_RECEIPT.json',{k:v for k,v in plan.items() if k!='members'})
        with (out/'PLAN_ATTEMPTS.jsonl').open('a') as f:f.write(json.dumps({k:v for k,v in plan.items() if k not in ('members','bytes_by_role')},ensure_ascii=False)+'\n')
    else:
        if info['missing'] or info['cache_closure_issues'] or info['unresolved_project_imports']:raise ValueError('Incomplete final packet inputs; inspect MANIFEST/plan. No ZIP written.')
        archive=out/'R15_PRO_REVIEW_v1.zip'
        if archive.exists():raise FileExistsError('Existing review archive retained; choose new version in source before rebuilding')
        # Contents frozen in this build inventory. Never re-encode numeric arrays.
        with zipfile.ZipFile(archive,'x',allowZip64=True) as z:
            for e in entries:
                source=Path(e['source']);method=zipfile.ZIP_DEFLATED if e['compression'].startswith('ZIP_DEFLATED') else zipfile.ZIP_STORED
                z.write(source,e['member'],compress_type=method,compresslevel=1)
                stat=source.stat()
                if (stat.st_size,stat.st_mtime_ns)!=(e['bytes'],e['mtime_ns']):raise RuntimeError('Source changed during build: '+str(source))
            for name in ('MANIFEST.json','MANIFEST.csv','START_HERE.md','portable_state_reader.py'):
                z.write(out/name,str((root/'archives/pro_review_v1'/name).relative_to(W)),compress_type=zipfile.ZIP_DEFLATED,compresslevel=1)
        with zipfile.ZipFile(archive) as z:
            bad=z.testzip()
            if bad:raise RuntimeError('Ordinary CRC failure: '+bad)
            actual={i.filename:i for i in z.infolist()}
            expected={e['member'] for e in entries}|{str((out/name).relative_to(W)) for name in ('MANIFEST.json','MANIFEST.csv','START_HERE.md','portable_state_reader.py')}
            if set(actual)!=expected:raise RuntimeError('Member set mismatch')
            for e in entries:
                i=actual[e['member']]
                if i.file_size!=e['bytes']:raise RuntimeError('Member byte count mismatch')
                if i.compress_type!=(zipfile.ZIP_DEFLATED if e['compression'].startswith('ZIP_DEFLATED') else zipfile.ZIP_STORED):raise RuntimeError('Compression method mismatch')
        directory=[dict(member=i.filename,bytes=i.file_size,compressed_bytes=i.compress_size,compression=i.compress_type,CRC32=f'{i.CRC:08x}') for i in actual.values()]
        receipt=dict(status='BUILT_CRC_SIZE_MEMBER_CHECKED_NOT_PUBLISHED',archive=str(archive),archive_bytes=archive.stat().st_size,members=len(actual),source_bytes=plan['source_bytes'],ordinary_CRC='PASS',
            archive_member_directory=directory,opaque_OFFLINE_SCORING_ONLY_truth_files=24,subset_scope=plan['completeness'],numeric_GT_array_loads=0,science_runs=0,hash_operations='NONE',CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
        receipt['archive_under_parent_500MB_limit']=archive.stat().st_size<500_000_000
        dump(out/'BUILD_RECEIPT.json',receipt)
    print(json.dumps(dict(mode='plan' if a.plan else 'build',files=len(entries),source_bytes=plan['source_bytes'],binary_stored_bytes=plan['binary_stored_bytes'],missing=info['missing'],cache_closure_issues=info['cache_closure_issues'],unresolved_project_imports=info['unresolved_project_imports'],out=str(out)),indent=2))


if __name__=='__main__':main()
