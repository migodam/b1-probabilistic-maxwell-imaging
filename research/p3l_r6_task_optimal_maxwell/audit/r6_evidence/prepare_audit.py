from pathlib import Path
import shutil,json,hashlib,subprocess,re
root=Path.cwd(); out=root/'research/p3l_r6_task_optimal_maxwell/audit/r6_evidence'; inp=root/'research/p3l_r6_task_optimal_maxwell/inputs/Research_P3LR5T'
a=out/'analytical_replay';a.mkdir(exist_ok=True);shutil.copy2(inp/'analytical_units.py',a/'analytical_units.py')
d=out/'decoder_replay';d.mkdir(exist_ok=True);shutil.copytree(inp/'decoder_sources',d/'decoder_sources',dirs_exist_ok=True)
s=(inp/'decoder_atom_check.py').read_text().replace('hashlib.sha256((SRC / name).read_bytes()).hexdigest()','hashlib.blake2b((SRC / name).read_bytes()).hexdigest()').replace('"local_sha256"','"local_blake2b"')
s=s.replace('b = mean + z_rest @ B.T\n            host.append(float(ndtr((-tau - b.max(axis=1)) / a0).mean()))\n            pure.append(float(ndtr((b.min(axis=1) - tau) / a0).mean()))', '''# Bounded decoder-only audit adaptation; largest b is 4096x256 float64.
            host_sum = pure_sum = 0.0
            for start_batch in range(0, len(z_rest), 4096):
                b = mean + z_rest[start_batch:start_batch+4096] @ B.T
                host_sum += ndtr((-tau - b.max(axis=1)) / a0).sum()
                pure_sum += ndtr((b.min(axis=1) - tau) / a0).sum()
            host.append(float(host_sum / len(z_rest)))
            pure.append(float(pure_sum / len(z_rest)))''')
(d/'decoder_atom_check.py').write_text(s)
# compare byte identities to fixed git objects, never SHA256.
checks=[]
for name,path in {'prior.py':'src/p1_material/prior.py','ONLINE_TARGET.json':'research/p1_material_target_20261018/configs/ONLINE_TARGET.json','geometry.json':'research/p1_material_target_20261018/inputs/geometry.json'}.items():
    old=subprocess.check_output(['git','show','ad71a86ef3c28bf48049199f733080a0f2b95e33:'+path]); new=(inp/'decoder_sources'/name).read_bytes();checks.append(dict(path=path,bytes=len(new),equals_fixed_git=old==new,blake2b=hashlib.blake2b(new).hexdigest()))
# Every absolute checkpoint/source reference in original frozen selection, not just finalists.
frozen=root/'research/p3l_r5_parallel_ab/stage_b_material_risk/SELECTION_FROZEN.json';f=json.loads(frozen.read_text());paths=[x['path'] for x in f['checkpoints']]+list(f['source_assets'].values());assets=[]
for p in paths:
    q=Path(p); assets.append(dict(path=p,exists=q.exists(),bytes=q.stat().st_size if q.exists() else None,role='historical_read_only'))
# Inventory actual R5B raw/checkpoint directories without opening truth.
b=root.parent/'b1-p3l-r5-b-confirm/research/p3l_r5_b_independent_confirmation'
raw=[]
for q in b.rglob('*'):
    if q.is_file() and q.suffix in ('.npz','.pt','.jsonl'):
        raw.append(dict(path=str(q),bytes=q.stat().st_size))
headings=[]
for lineno,line in enumerate((inp/'output/Research_P3LR5T_Source.txt').read_text().splitlines(),1):
    if re.match(r'^#{1,3} ',line):headings.append(dict(line=lineno,heading=line,status='READ_COVERAGE_ONLY_NOT_PROOF_APPROVAL'))
(out/'SOURCE_AND_ASSET_AUDIT.json').write_text(json.dumps(dict(decoder_source_checks=checks,original_R5_selection_assets=assets,r5b_raw_checkpoint_inventory=raw,r5b_raw_inventory_count=len(raw),missing_original_selection_paths=[x['path'] for x in assets if not x['exists']],theory_section_coverage=headings,identity_method='byte comparison plus BLAKE2b for three small decoder inputs; no SHA256; no large raw files hashed',physics_RHS=0,training=0),indent=2))
