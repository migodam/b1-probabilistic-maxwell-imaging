#!/usr/bin/env python3
"""Build two separate lossless evidence ZIPs; stdlib only, no science execution."""
import ast
import csv
import json
from pathlib import Path
import time
import zipfile
from datetime import datetime, timezone

R = Path(__file__).resolve().parents[1]
W = R.parents[1]
PREFIX = 'research/research_max_struct_r14/'
TEXT = {'.py','.md','.csv','.json','.txt','.yml','.yaml'}


def dump(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2)+'\n', encoding='utf8')


def collect_pro():
    chosen = {}
    missing = []
    def add(p, role, required=True):
        p = Path(p)
        if not p.is_file():
            if required: missing.append(str(p.relative_to(R)))
            return
        member = str(p.relative_to(W))
        if any(x in member.lower() for x in ['sealed_truth','private_operational','chat_readback','fulltext']):
            raise ValueError('Forbidden Pro member: '+member)
        chosen[member] = (p, role)
    docs = 'RESEARCH_REPORT_ZH THEORY_DERIVATIONS THEORY_STATUS METHODS_AND_ABLATIONS BASELINE_AUDIT PRIOR_ART_AND_NOVELTY CONTROL_DERIVATION TRAJECTORY_ANALYSIS BENCHMARKS COST_LEDGER FAILURES REPRODUCE PRO_HANDOFF CHAT_HANDOFF PAPER_DRAFT TECH_STACK_MAP PRO_PACKET_CONTENTS EXPERIMENT_CONTRACT INDEPENDENT_CONTRACT CORE_IMPLEMENTATION_NOTES'.split()
    for n in docs: add(R/(n+'.md'), 'author report')
    for p in (R/'configs').glob('*.json'): add(p, 'frozen public config')
    for p in (R/'data_public').rglob('*'):
        if p.is_file() and p.suffix in {'.npz','.json'}: add(p, 'public observation/input')
    tables = {
      'exposed48_v2':'RECEIPT.json SUMMARY.csv ALL_VIEWS.csv PAIRED_METRIC_DIFFERENCES.csv RECOVERY_GAIN_LOSS.csv COMMON_CPU_PREFIXES.csv CORRECTOR_EFFECTS.csv RAW_INITIAL_AUDIT.csv STRICT48_PAID_POOL_BRACKETS.csv COST_LEDGER.csv',
      'independent24_v1':'RECEIPT.json SUMMARY.csv ALL_VIEWS.csv POSTHOC_GAIN_LOSS.csv COMMON_CPU_PREFIXES.csv DIRECT_PHYSICAL_SEARCH_VIEWS.csv PUBLIC_SELECTOR_REPLAY.csv PUBLIC_SELECTOR_COUNTS.json RAW_INITIAL_AUDIT.csv COST_LEDGER.csv ISSUES.json',
      'self_only16_v1':'RECEIPT.json SUMMARY.csv ALL_VIEWS.csv COST_LEDGER.csv',
      'stage_control16_v1':'RECEIPT.json SUMMARY.csv ALL_VIEWS.csv COST_LEDGER.csv',
      'shared_readout4_v1':'SCORING_RECEIPT.json METHOD_VIEW_SUMMARY.csv SCORES.csv CORRECTOR_DECOMPOSITION.csv MATERIAL_ARRAY_DIFFERENCES.csv COST_TOTALS.json RAW_INITIAL_AUDIT.csv',
      'vector4_v1':'SCORING_RECEIPT.json METHOD_VIEW_SUMMARY.csv SCORES.csv CORRECTOR_DECOMPOSITION.csv COST_LEDGER.csv'}
    for d,names in tables.items():
        for n in names.split(): add(R/'analysis'/d/n, 'compact scored evidence')
    old = [('0002',128),('0002',512),('0017',128),('0017',512),('0017',2500),('0037',128),('0037',512)]
    for c,k in old:
        add(R/f'RESULTS/dualmetric_pilot12_v1/r13_fresh_{c}__pair01/CHECKPOINT_{k:05}.npz','original full checkpoint')
        d=R/f'audit/interventions/r13_fresh_{c}_k{k:05}'
        for n in 'RECEIPT.json BRANCHES.json CONFIGURATION.json COST_LEDGER.json FAILURES.json'.split(): add(d/n,'intervention provenance')
    for c in ['0012','0021']:
        for pair in ['01','11']:
            for k in [128,512,1500]: add(R/f'RESULTS/independent24_v1/r14_new_{c}__pair{pair}/CHECKPOINT_{k:05}.npz','new independent full checkpoint')
    for n in ['INTERVENTION_BRANCHES.json','INTERVENTION_RECEIPTS.csv']: add(R/'audit/dynamics_all192_v1'/n,'intervention provenance')
    for n in 'FIRST_STEP_FORCE_ARRAYS.npz FIRST_STEP_MATERIAL_FORCES.csv MECHANICAL_GROUP_SUMMARY.csv MECHANICAL_BRIEF.json PROTOCOL_BEFORE_ACTIONS.json COUNTER_AND_SOURCE_SCOPE.json VALIDATION_COST_EVENTS.csv RECEIPT.json FAILURES.json SOURCE_AT_ACTIONS.py REPORT_ZH.md'.split(): add(R/'audit/material_force_decomposition_v1'/n,'force decomposition')
    for p in (R/'audit/local_closed_loop_v1').glob('*'):
        if p.suffix in {'.npz','.json','.csv','.md'}: add(p,'offline root diagnostic')
    for d in ['finite_mode_trajectory_v1','coupled_mode_trajectory_v1']:
        for n in 'FROZEN_CONFIG.json MODE_SELECTION.json TRAJECTORIES.csv FIRST_STEP_POLICY_EQUALITY.json RECEIPT.json COST_LEDGER.json FAILURES.json RESULTS_BRIEF.md INTERPRETATION_BOUNDARIES.md MODE_AND_IDEAL_PREDICTIONS.npz'.split(): add(R/'audit'/d/n,'finite/coupled diagnostic')
    for n in ['REVIEW.md','CLAIM_EVIDENCE.csv','SEARCH_LOG.md','VERIFIED_CITATIONS.json']: add(R/'audit/prior_art'/n,'bibliography evidence')
    for n in ['VECTOR_CONTRACT.md','PARENT_INDEPENDENT24_AUTHORIZATION.json']: add(R/'audit'/n,'public contract/authorization')
    figure_dirs=[R/'FIGURES/exposed48',R/'analysis/independent24_v1/FIGURES',R/'analysis/vector4_v1/FIGURES',R/'audit/material_force_decomposition_v1',R/'audit/finite_mode_trajectory_v1/FIGURES',R/'audit/coupled_mode_trajectory_v1/FIGURES']
    for d in figure_dirs:
        for p in d.rglob('*'):
            if p.suffix.lower() in {'.png','.pdf'} and 'FIRST_LAYOUT' not in p.name: add(p,'actual complete figure set')
    excluded = {'generate_independent24.py','launch_remote_upload.py','upload_draft_remote.py','package_evidence.py'}
    pending = [p for p in (R/'src').glob('*.py') if p.name not in excluded]
    pending += list((R/'runtime').rglob('*.py'))
    pending += list((W/'research/research_max_struct_r13/runtime').glob('*.py'))
    # Preserve compact provenance of the unchanged historical runtime.
    for p in (R/'runtime').rglob('PROVENANCE.json'): add(p,'historical runtime provenance')
    seen=set(); unresolved=set()
    # Parse imports without importing or executing any scientific module.
    while pending:
        p=pending.pop()
        if p in seen: continue
        seen.add(p);add(p,'project source dependency')
        tree=ast.parse(p.read_text(encoding='utf8'))
        for node in ast.walk(tree):
            modules=[]
            if isinstance(node,ast.Import): modules=[a.name for a in node.names]
            elif isinstance(node,ast.ImportFrom) and node.module:
                modules=[node.module]+[node.module+'.'+a.name for a in node.names if a.name!='*']
            for m in modules:
                if m.startswith('research.'):
                    base=W/Path(*m.split('.'))
                    q=base.with_suffix('.py')
                    if q.is_file(): pending.append(q)
                    elif (base/'__init__.py').is_file(): pending.append(base/'__init__.py')
                    elif not base.is_dir() and not any((W/Path(*m.rsplit('.',1)[0].split('.'))).with_suffix('.py').is_file() for _ in [0]): unresolved.add(m)
                elif '.' not in m:
                    # Frozen runners add the historical R13 src directory to sys.path.
                    q=W/'research/research_max_struct_r13/src'/(m+'.py')
                    if q.is_file(): pending.append(q)
    return chosen,missing,sorted(unresolved)


def collect_truth():
    chosen={}; mappings=[]
    for cohort in ['exposed48_v2','independent24_v1','vector4_v1']:
        receipt=R/'analysis'/cohort/'TRUTH_SOURCES.csv'
        with receipt.open(newline='',encoding='utf8') as f: rows=list(csv.DictReader(f))
        for row in rows:
            original=row['path']
            if cohort=='independent24_v1':
                name=original.replace('\\','/').split('/')[-1]
                p=R/'sealed_truth/independent24'/name
                member=PREFIX+'sealed_truth/independent24/'+name
                status='local original-generation counterpart; no remote byte-identity assertion'
            else:
                p=Path(original)
                suffix=str(p).split('/research/research_max_struct_r13/',1)[1]
                member='research/research_max_struct_r13/'+suffix
                status='exact scoring source path'
            if not p.is_file(): raise FileNotFoundError(p)
            chosen[member]=(p,'synthetic offline truth')
            mappings.append(dict(cohort=cohort,recorded_scoring_source=original,local_source=str(p),archive_path=member,provenance=status))
        chosen['provenance/'+cohort+'_TRUTH_SOURCES.csv']=(receipt,'source receipt')
    for n in ['GENERATOR_MANIFEST.json','GENERATION_DIAGNOSTICS.json']:
        p=R/'sealed_truth/independent24'/n
        if p.is_file(): chosen[PREFIX+'sealed_truth/independent24/'+n]=(p,'original local-generation scoring metadata')
    for n in ['GENERATION24_RECEIPT.json','PARENT_INDEPENDENT24_AUTHORIZATION.json']:
        chosen['provenance/'+n]=(R/'audit'/n,'generation/authorized scoring provenance')
    return chosen,mappings


def package(directory, filename, chosen, explanation, extra):
    cpu=time.process_time(); wall=time.perf_counter()
    directory.mkdir(parents=True,exist_ok=False)
    snapshot=datetime.now(timezone.utc).isoformat()
    entries=[]
    for member,(p,role) in sorted(chosen.items()):
        stat=p.stat()
        entries.append(dict(archive_path=member,source_path=str(p),bytes=stat.st_size,source_mtime_ns=stat.st_mtime_ns,role=role,compression='DEFLATED_level1' if p.suffix.lower() in TEXT else 'STORED'))
    manifest=dict(snapshot_utc=snapshot,source_files=len(entries),source_bytes=sum(e['bytes'] for e in entries),members=entries,**extra)
    dump(directory/'MANIFEST.json',manifest)
    with (directory/'MANIFEST.csv').open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=list(entries[0]));w.writeheader();w.writerows(entries)
    (directory/'START_HERE.md').write_text(explanation,encoding='utf8')
    archive=directory/filename
    with zipfile.ZipFile(archive,'w',allowZip64=True) as z:
        for e in entries:
            p=Path(e['source_path']); z.write(p,e['archive_path'],compress_type=zipfile.ZIP_DEFLATED if p.suffix.lower() in TEXT else zipfile.ZIP_STORED,compresslevel=1)
            stat=p.stat()
            if (stat.st_size,stat.st_mtime_ns)!=(e['bytes'],e['source_mtime_ns']): raise RuntimeError('Source changed during snapshot: '+str(p))
        for name in ['START_HERE.md','MANIFEST.json','MANIFEST.csv']: z.write(directory/name,name,compress_type=zipfile.ZIP_DEFLATED,compresslevel=1)
    with zipfile.ZipFile(archive) as z:
        bad=z.testzip()
        infos=z.infolist()
        if bad: raise RuntimeError('CRC failed: '+bad)
        actual={i.filename:i for i in infos}
        if set(actual)!={e['archive_path'] for e in entries}|{'START_HERE.md','MANIFEST.json','MANIFEST.csv'}: raise RuntimeError('Member set mismatch')
        for e in entries:
            i=actual[e['archive_path']]
            if i.file_size!=e['bytes']: raise RuntimeError('Size mismatch')
            expected=zipfile.ZIP_DEFLATED if e['compression'].startswith('DEFLATED') else zipfile.ZIP_STORED
            if i.compress_type!=expected: raise RuntimeError('Compression mismatch')
        verification=[dict(path=i.filename,bytes=i.file_size,compressed_bytes=i.compress_size,CRC32=f'{i.CRC:08x}',compression=i.compress_type) for i in infos]
    receipt=dict(archive=str(archive),archive_bytes=archive.stat().st_size,members=len(infos),uncompressed_bytes=sum(i.file_size for i in infos),snapshot_utc=snapshot,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,ordinary_CRC='PASS all members',member_paths_sizes_compression='PASS',hashes_used=False,science_executed=False,verification=verification,**extra)
    dump(directory/'RECEIPT.json',receipt)
    return receipt


def main():
    pro,missing,unresolved=collect_pro()
    if missing: raise FileNotFoundError('Required packet files missing: '+repr(missing))
    explanation='''# START HERE — R14 Pro finite evidence packet\n\nRead research/research_max_struct_r14/PRO_PACKET_CONTENTS.md, then RESEARCH_REPORT_ZH.md, THEORY_DERIVATIONS.md and THEORY_STATUS.md. Reports are the author snapshot; packaging is no scientific acceptance. MANIFEST lists each original source path, byte count and role. Extract into an empty directory to preserve research/... source layout. Frozen C2 uses unqualified historical runtime imports: the original field_current.py is supplied under research/research_max_struct_r13/runtime; expose that directory and historical R13/src on PYTHONPATH when restoring a minimal packet runtime. Python environment and historical Git requirements are in REPRODUCE.md; dependency imports were parsed without executing code. No environment or scientific reproducibility test was run by this packaging step.\n\nIncludes 7 original full checkpoints, 12 new independent full checkpoints (0012/0021 × 01/11 × 128/512/1500), the original force NPZ, offline diagnostic arrays, compact scored tables, all 48+24 public observations and 4 public vector inputs, original source dependencies and complete selected figure directories. All binary files retain original bytes. Root artificial known-model diagnostic arrays are declared offline mechanism inputs, not cohort truth. Scored CSVs and published scored material figures contain authorized outcome summaries; cohort truth arrays are excluded.\n\nThe large native histories, all paid candidates and full intervention trajectories are omitted from this curated packet and belong to full raw release. DOWNLOADS.json (when published by the owner) is the authoritative actual asset URL index; this build does not invent release URLs or upload anything. The synthetic offline truth archive is deliberately separate and must never be added to this packet. No credentials, chat caches or third-party fulltexts are included.\n\n19 full checkpoint NPZs are original inputs, not selected-m scoring exports. Optional downstream physics/scoring consumes budget and is not performed by opening this package. Ordinary ZIP CRC and member sizes are checked; no cryptographic hashes are used.\n'''
    pr=package(R/'archives/pro_packet_v1','R14_PRO_PACKET_v1.zip',pro,explanation,dict(missing_required=missing,unresolved_project_imports=unresolved,source_dependency_validation='AST import closure; no module import/execution'))
    truth,mappings=collect_truth()
    explanation='''# SYNTHETIC OFFLINE TRUTH — separate authorized scoring archive\n\nThis archive contains original bytes for 48 exposed scalar cases, 24 independent scalar cases and 4 vector cases. It is not part of the Pro evidence packet. Do not expose these arrays to inverse methods, public selectors or prospective prediction choices. Existing independent24 scoring was parent-authorized only after all 264 predictions were sealed; the authorization receipt is included. Exposed48 is an exposed development cohort, not a new blind test. Vector qualification remains the scope declared in the author reports.\n\nUnzip into an empty repository-layout directory; consult the per-member manifest and source mapping. Exposed/vector sources are the exact paths recorded by scoring receipts. Independent sources are the local original-generation NPZ counterparts of the Windows scoring paths; their remote byte identity was not checked or asserted. Original GENERATOR_MANIFEST and GENERATION_DIAGNOSTICS are included separately as private offline scoring metadata, never in the Pro packet. The provenance map records this distinction. No truth arrays were loaded by packaging, and no scorer/physics was run.\n\nPreserve original remote seals. If relocating, follow REPRODUCE.md post-restore mapping/seal instructions; never rewrite historical Windows seals as Mac originals. New scoring eligibility requires a fixed prediction seal and explicit authorization before truth access. Ordinary ZIP CRC verifies transport only, not scientific validity.\n'''
    tr=package(R/'archives/offline_truth_v1','R14_SYNTHETIC_OFFLINE_TRUTH_v1.zip',truth,explanation,dict(source_mappings=mappings,truth_arrays_loaded=False,missing_sources=[]))
    out=R/'audit/publication';out.mkdir(exist_ok=True)
    dump(out/'PRO_PACKET_BUILD_RECEIPT.json',pr)
    dump(out/'PRO_PACKET_OFFLINE_TRUTH_RECEIPT.json',tr)
    print(json.dumps({k:{n:v[n] for n in ['archive','archive_bytes','members','CPU_seconds','wall_seconds']} for k,v in [('pro',pr),('offline_truth',tr)]},indent=2))

if __name__=='__main__': main()
