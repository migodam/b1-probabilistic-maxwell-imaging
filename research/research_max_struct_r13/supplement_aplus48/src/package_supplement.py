"""Package existing public evidence; no computation of experiment results."""
from pathlib import Path
import json,time,zipfile
C=time.process_time();T=time.perf_counter();S=Path(__file__).resolve().parents[1];R=S.parent;REPO=R.parents[1]
D=S/'release_preparation';D.mkdir(exist_ok=True)
def files(root):return [p for p in root.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc']
def archive(name,paths,base):
    c=time.process_time();t=time.perf_counter();target=D/name;rows=[]
    with zipfile.ZipFile(target,'x',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(set(paths)):
            member=p.relative_to(base).as_posix();z.write(p,member);rows.append(dict(member=member,bytes=p.stat().st_size))
    return dict(name=name,local_path=str(target),bytes=target.stat().st_size,uncompressed_bytes=sum(x['bytes'] for x in rows),files=len(rows),CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t,inventory=rows)
analysis=files(S/'analysis')
code=files(S/'src')+files(S/'configs')+files(S/'audit')
code += [p for p in S.iterdir() if p.is_file() and p.suffix in ('.md','.json','.html')]
public=files(R/'data_public/fresh48')+files(R/'data_public/fresh16_nn')
out=[archive('max-struct-r13-aplus48-analysis.zip',analysis,S),archive('max-struct-r13-aplus48-public-inputs-and-sources.zip',code+public,REPO)]
raw=['aplus48_extension32_full_results.zip','aplus48_dynamics8_full_results.zip','aplus_stack112_full_results.zip','aplus_stack112_execution_sources.zip']
for name in raw:
    p=S/'archives'/name
    with zipfile.ZipFile(p) as z:
        names=z.namelist()
        if any('private_operational' in n or n.startswith('D:') or '/D:' in n for n in names):raise ValueError('Private/alias member')
    out.append(dict(name=name,local_path=str(p),bytes=p.stat().st_size,files=len(names),scope='Actual completed campaign/execution archive; unchanged'))
manifest=dict(status='LOCAL_REVIEWABLE_ARTIFACTS_PREPARED_NOT_YET_PUBLISHED',assets=out,CPU_seconds=time.process_time()-C,wall_seconds=time.perf_counter()-T,new_physics=0,new_scoring=0,policy='No private SSH scripts/logs/credentials; no local D: hardlink alias trees; original data/weights remain separate references; no SHA checks')
(S/'audit/RELEASE_PACKAGE_RECEIPT.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps([{k:r[k] for k in ('name','bytes','files')} for r in out]))
