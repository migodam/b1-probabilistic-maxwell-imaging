"""Finite R4 public-only SSH jobs; existing transport, isolated own remote root."""
from pathlib import Path
import argparse,importlib.util,os,subprocess,sys,zipfile,json,time
STUDY=Path(__file__).resolve().parents[1];REPO=STUDY.parents[1];WORKSPACE=REPO.parents[3]
REMOTE=os.environ.get('GN_NN_REMOTE_ROOT','D:/AI/GN_NN_R4_20261107');REL='research/research_gn_nn_preservation_r4'
MASK=20560 # Pcores4,6,12,14; historical85 overlapped another session's0,2.
def transport():
    spec=importlib.util.spec_from_file_location('existing_ssh',WORKSPACE/'Gaussian/A9/code/remote.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
def ps(script,timeout=120,*,echo=True):
    r=transport().powershell(script,timeout=timeout)
    out=r.stdout.decode('utf8','replace');err=r.stderr.decode('utf8','replace')
    if echo:print(out,end='',flush=True)
    if err:print(err,file=sys.stderr,flush=True)
    if r.returncode:raise RuntimeError(f'Remote returncode{r.returncode}: {err[-1000:]}')
    return out
def frozen_test_contract():
    spec=importlib.util.spec_from_file_location('r4_public_test_partition',STUDY/'experiments/run_test_partition.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module.load_contract()

def transfer_receipt(action,started,**fields):
    receipt=dict(action=action,wall_seconds=time.perf_counter()-started,unix=time.time(),truth_access='NONE',**fields)
    path=STUDY/f'audit/REMOTE_TRANSFER_{action}_{time.time_ns()}.json';path.write_text(json.dumps(receipt,indent=2)+'\n');return receipt

def upload(include_data=True):
    started=time.perf_counter();pack_started=started;status='FAILED';error=None;pack_seconds=scp_seconds=extract_seconds=None
    archive=STUDY/'audit/REMOTE_public.zip'
    files=list((STUDY/'src').glob('*.py'))+list((STUDY/'experiments').glob('*.py'))
    files += [f for f in (STUDY/'vendor').rglob('*') if f.is_file() and '__pycache__' not in f.parts and f.suffix in ('.py','.json','.md','.txt')]
    files += [STUDY/'audit'/name for name in ('METHOD_FREEZE.json','FINAL_RUNS.json','TEST_RESOURCE_PARTITION.json') if (STUDY/'audit'/name).exists()]
    for folder in ('research_z_alg_r4_diversity','research_gn_nn_completion_r3','research_gn_search_space_r2'):
        files+=list((REPO/'research'/folder/'src').glob('*.py'))
    if include_data:files+=[STUDY/'data/operators.npz',STUDY/'data/manifests/online.json']+list((STUDY/'data/public').glob('*.npz'))
    try:
        with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=1) as z:
            for f in files:
                if f.exists():z.write(f,f.relative_to(REPO))
        pack_seconds=time.perf_counter()-pack_started;ps(f"New-Item -ItemType Directory -Force -Path '{REMOTE}' | Out-Null")
        m=transport();tick=time.perf_counter();subprocess.run(['scp','-i',m.KEY,'-o','BatchMode=yes',str(archive),m.HOST+':'+REMOTE+'/public.zip'],check=True);scp_seconds=time.perf_counter()-tick
        tick=time.perf_counter()
        check="""import pathlib,zipfile,json
root=pathlib.Path(r'REMOTE_ROOT');rel=pathlib.Path('REL_ROOT')
with zipfile.ZipFile(root/'public.zip') as z:
    for name in ('METHOD_FREEZE.json','FINAL_RUNS.json','TEST_RESOURCE_PARTITION.json'):
        entry=(rel/'audit'/name).as_posix();dst=root/'project'/rel/'audit'/name
        if entry in z.namelist() and dst.exists() and json.loads(dst.read_text(encoding='utf-8'))!=json.loads(z.read(entry)):raise RuntimeError('Preserve mismatched existing freeze: '+name)
""".replace('REMOTE_ROOT',REMOTE).replace('REL_ROOT',REL)
        ps("& D:/python/python.exe -c '"+check.replace("'","''")+"'; if ($LASTEXITCODE -ne 0) {exit $LASTEXITCODE}",300)
        ps(f"Expand-Archive -LiteralPath '{REMOTE}/public.zip' -DestinationPath '{REMOTE}/project' -Force",300);extract_seconds=time.perf_counter()-tick
        status='PUBLIC_ONLY_UPLOAD_COMPLETE';print(status+'; NO_TRUTH')
    except Exception as exc:error=repr(exc);raise
    finally:transfer_receipt('source_upload',started,status=status,error=error,archive=str(archive),archive_bytes=archive.stat().st_size if archive.exists() else None,packaging_seconds=pack_seconds,scp_seconds=scp_seconds,remote_extract_seconds=extract_seconds,included_public_observations=include_data,includes_weights_or_truth=False)

def upload_predictions():
    started=time.perf_counter();freeze,final,allocation,_=frozen_test_contract();ids=allocation['remote']['obs_ids'];files=[];missing=[]
    for run in final['nn_runs']:
        for obs in ids:
            for ext in ('.npz','.json'):
                f=STUDY/'predictions'/run/(obs+ext)
                if not f.exists():missing.append(str(f.relative_to(STUDY)))
                else:files.append(f)
    if missing:
        transfer_receipt('prediction_upload',started,status='HOLD_MISSING_PUBLIC_PREDICTIONS',missing=missing,new_physics_counts={});raise FileNotFoundError(f'{len(missing)} missing remote-assigned prediction files')
    archive=STUDY/f'audit/REMOTE_predictions_{time.time_ns()}.zip';status='FAILED';error=None;report=None;scp_seconds=None;packaging_seconds=None;remote_import_seconds=None
    try:
        tick=time.perf_counter()
        with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=1) as z:
            for f in files:z.write(f,f.relative_to(REPO))
            for name in ('METHOD_FREEZE.json','FINAL_RUNS.json','TEST_RESOURCE_PARTITION.json'):z.write(STUDY/'audit'/name,(STUDY/'audit'/name).relative_to(REPO))
        packaging_seconds=time.perf_counter()-tick;m=transport();tick=time.perf_counter();subprocess.run(['scp','-i',m.KEY,'-o','BatchMode=yes',str(archive),m.HOST+':'+REMOTE+'/predictions.zip'],check=True);scp_seconds=time.perf_counter()-tick
        code="""import pathlib,zipfile,json
root=pathlib.Path(r'REMOTE_ROOT')/'project';rel=pathlib.Path('REL_ROOT');imported=[];existing=[];held=[]
with zipfile.ZipFile(pathlib.Path(r'REMOTE_ROOT')/'predictions.zip') as z:
    for name in ('METHOD_FREEZE.json','FINAL_RUNS.json','TEST_RESOURCE_PARTITION.json'):
        entry=(rel/'audit'/name).as_posix();dest=root/rel/'audit'/name;content=z.read(entry)
        if dest.exists() and json.loads(dest.read_text(encoding='utf-8'))!=json.loads(content):raise RuntimeError('Frozen remote metadata mismatch: '+name)
        if not dest.exists():dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(content)
    for entry in z.namelist():
        if '/predictions/' not in entry or not entry.endswith('.npz'):continue
        path=pathlib.PurePosixPath(entry)
        if '..' in path.parts or path.is_absolute():raise ValueError('Archive path escape')
        dest=root/path;j=dest.with_suffix('.json');meta=str(path.with_suffix('.json'))
        if dest.exists() and j.exists():existing.append(entry);continue
        if dest.exists() or j.exists():held.append(entry);continue
        dest.parent.mkdir(parents=True,exist_ok=True)
        with dest.open('xb') as out:out.write(z.read(entry))
        with j.open('xb') as out:out.write(z.read(meta))
        imported.append(entry)
print(json.dumps(dict(imported=imported,existing=existing,held=held)))
""".replace('REMOTE_ROOT',REMOTE).replace('REL_ROOT',REL)
        tick=time.perf_counter();report=json.loads(ps("& D:/python/python.exe -c '"+code.replace("'","''")+"'",300,echo=False).strip().splitlines()[-1]);remote_import_seconds=time.perf_counter()-tick;status='PUBLIC_REMOTE_PREDICTIONS_UPLOADED' if not report['held'] else 'HOLD_PARTIAL_CANONICAL_PREDICTIONS'
        print(json.dumps(dict(status=status,imported_pairs=len(report['imported']),existing_pairs=len(report['existing']),held_pairs=len(report['held']))),flush=True)
    except Exception as exc:error=repr(exc);raise
    finally:transfer_receipt('prediction_upload',started,status=status,error=error,archive=str(archive),nn_runs=final['nn_runs'],observations=len(ids),prediction_files=len(files),archive_bytes=archive.stat().st_size if archive.exists() else None,packaging_seconds=packaging_seconds,scp_seconds=scp_seconds,remote_import_seconds=remote_import_seconds,report=report,includes_weights_or_truth=False)

def require_new_remote_allocation_free():
    # Read only our public completion files and our process tree. Never stop it.
    code="""import pathlib,json
s=pathlib.Path(r'REMOTE_ROOT/project/REL_ROOT');records=json.loads((s/'data/manifests/online.json').read_text())['records'];ids=[r['obs_id'] for r in records if r['split'] in ('val','test')];n=sum((s/'anchors'/(obs+suffix+'.npz')).exists() and (s/'anchors'/(obs+suffix+'.json')).exists() for obs in ids for suffix in ('_discrepancy','_weakreg1e-5'));print(n)
if n!=744:raise RuntimeError('Wait for all744 existing baseline pairs; no new queue started')
""".replace('REMOTE_ROOT',REMOTE).replace('REL_ROOT',REL)
    ps("& D:/python/python.exe -c '"+code.replace("'","''")+"'; if ($LASTEXITCODE -ne 0) {exit $LASTEXITCODE}")
    script="""$all=@(Get-CimInstance Win32_Process); $own=@{};
foreach($p in $all){if($p.CommandLine -and $p.CommandLine.Contains('REMOTE_ROOT')){$own[[int]$p.ProcessId]=$true}}
for($i=0;$i -lt 5;$i++){foreach($p in $all){if($own.ContainsKey([int]$p.ParentProcessId)){$own[[int]$p.ProcessId]=$true}}}
$busy=@($all | Where-Object {$own.ContainsKey([int]$_.ProcessId) -and $_.Name -eq 'python.exe' -and $_.CommandLine -match 'run_anchors.py|run_baselines.py|run_test_partition.py|run_cie_reference.py|run_evaluation.py'});
if($busy.Count -gt 0){$busy | Select-Object ProcessId,CommandLine | Format-List; throw 'Own R4 numerical allocation still busy; no new queue started'}
""".replace('REMOTE_ROOT',REMOTE)
    ps(script)

def start(kind='anchors'):
    log=REMOTE+'/'+kind+'.log';err=REMOTE+'/'+kind+'.err';script_path=REMOTE+'/'+kind+'.ps1'
    args='--split train --split val --split test --workers 4 --core-mask '+str(MASK) if kind=='anchors' else '--split val --split test --workers 4 --core-mask '+str(MASK)
    target='run_anchors.py' if kind=='anchors' else 'run_baselines.py'
    body=f"Set-Location '{REMOTE}/project/{REL}'\n$env:OMP_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OPENBLAS_NUM_THREADS='1'; $env:VECLIB_MAXIMUM_THREADS='1'\n$env:PYTHONPATH='D:/AI/ZALG_R4_DIVERSITY_20261030/project/vendor'\n[System.Diagnostics.Process]::GetCurrentProcess().ProcessorAffinity=[IntPtr]{MASK}\n& D:/python/python.exe -B -u experiments/{target} {args}\nexit $LASTEXITCODE\n"
    escaped=body.replace("'","''")
    result=ps(f"Set-Content -LiteralPath '{script_path}' -Value '{escaped}' -Encoding UTF8; $p=Start-Process powershell.exe -ArgumentList '-NoProfile','-File','{script_path}' -RedirectStandardOutput '{log}' -RedirectStandardError '{err}' -PassThru; Write-Output $p.Id")
    (STUDY/f'audit/REMOTE_{kind}_start.json').write_text(json.dumps(dict(kind=kind,pid=result.strip(),remote=REMOTE,affinity_mask=MASK,cores=[4,6,12,14],historical_mask85_not_used='avoid existing unrelated CPU tasks on0,2',public_truth_access='NONE',start_unix=time.time()),indent=2)+'\n')
def run(kind='anchors'):
    started=time.perf_counter();new_kind=kind in ('test','cie','discrepancy-all-scales')
    try:
        if new_kind:
            from datetime import datetime,timezone
            remaining=datetime(2026,11,7,10,44,27,tzinfo=timezone.utc).timestamp()-time.time()
            if remaining<=0:raise RuntimeError('Fixed absolute deadline reached; no new work')
            require_new_remote_allocation_free()
        if kind in ('test','cie'):frozen_test_contract()
    except Exception as exc:
        if new_kind:transfer_receipt('queue_preflight_'+kind,started,status='HOLD_BEFORE_LAUNCH',error=repr(exc),new_physics_counts={})
        raise
    if kind=='test':
        args='--execution-partition remote --workers 4 --execute';target='run_test_partition.py'
    elif kind=='cie':
        args='--split test --scope all108 --methods cie_dual --workers 4 --execution-partition remote --execute --test-methods-frozen';target='run_cie_reference.py'
    elif kind=='discrepancy-all-scales':
        args='--split val --split test --workers 4 --core-mask '+str(MASK)+' --kind discrepancy-all-scales';target='run_baselines.py'
    else:
        # Existing anchor/baseline command branches are retained verbatim.
        args='--split train --split val --split test --workers 4 --core-mask '+str(MASK) if kind=='anchors' else '--split val --split test --workers 4 --core-mask '+str(MASK)
        target='run_anchors.py' if kind=='anchors' else 'run_baselines.py'
    body=f"Set-Location '{REMOTE}/project/{REL}'\n$env:OMP_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OPENBLAS_NUM_THREADS='1'; $env:VECLIB_MAXIMUM_THREADS='1'\n$env:PYTHONPATH='D:/AI/ZALG_R4_DIVERSITY_20261030/project/vendor'\n[System.Diagnostics.Process]::GetCurrentProcess().ProcessorAffinity=[IntPtr]{MASK}\n& D:/python/python.exe -B -u experiments/{target} {args}\nexit $LASTEXITCODE\n"
    ps("Set-Content -LiteralPath '"+REMOTE+'/'+kind+".ps1' -Value '"+body.replace("'","''")+"' -Encoding UTF8")
    redirection='*>>' if new_kind else '*>'
    setup_wall_seconds=time.perf_counter()-started;command_started=time.perf_counter()
    status='FAILED_OR_INTERRUPTED';error=None
    try:
        if kind=='discrepancy-all-scales':
            # Wrapper deadline owns this exact newly launched child tree only.
            script=f"$p=Start-Process powershell.exe -ArgumentList '-NoProfile','-File','{REMOTE}/{kind}.ps1' -PassThru -NoNewWindow; $end=[DateTime]::Parse('2026-11-07T10:44:27Z').ToUniversalTime(); while(-not $p.HasExited -and [DateTime]::UtcNow -lt $end){{Start-Sleep -Milliseconds 250; $p.Refresh()}}; if(-not $p.HasExited){{taskkill /PID $p.Id /T /F; Write-Output 'OWN_ALLSCALES_CHILD_INTERRUPTED_COST_UNKNOWN_NOT_ZERO'; exit 124}} else {{exit $p.ExitCode}}"
            ps(script,max(30,int(remaining)+15))
        else:ps(f"& powershell.exe -NoProfile -File '{REMOTE}/{kind}.ps1' {redirection} '{REMOTE}/{kind}.log'; exit $LASTEXITCODE",max(30,int(remaining)+15) if new_kind else 21600)
        status='REMOTE_COMMAND_RETURNED'
    except Exception as exc:error=repr(exc);raise
    finally:
        if new_kind:transfer_receipt('queue_'+kind,started,status=status,error=error,setup_and_resource_preflight_wall_seconds=setup_wall_seconds,foreground_command_wall_seconds=time.perf_counter()-command_started,affinity_mask=MASK,scientific_cost_scope='read canonical physical case receipts, not transport wall',interrupted_cost_status='UNKNOWN_NOT_ZERO if connection or own child interrupted; never infer zero',automatic_retry=False)

def status():
    ps(f"$d='{REMOTE}/project/{REL}/anchors'; (Get-ChildItem $d -Filter 'obs_*.npz' -ErrorAction SilentlyContinue | Measure-Object).Count; Get-Content '{REMOTE}/anchors.log' -Tail 2 -ErrorAction SilentlyContinue; Get-Content '{REMOTE}/anchors.err' -Tail 4 -ErrorAction SilentlyContinue")

def fetch_evaluation():
    """Completed remote-assigned frozen evaluation pairs and immutable receipt snapshots."""
    started=time.perf_counter();_,final,allocation,_=frozen_test_contract();stamp=time.time_ns();remote_zip=REMOTE+'/test_evaluation.zip';archive=STUDY/f'audit/REMOTE_evaluation_{stamp}.zip';status='FAILED';error=None;scp_seconds=None;fetched=0;held=[]
    code="""import pathlib,zipfile,json
s=pathlib.Path(r'REMOTE_ROOT/project/REL_ROOT');allocation=json.loads((s/'audit/TEST_RESOURCE_PARTITION.json').read_text(encoding='utf-8'));final=json.loads((s/'audit/FINAL_RUNS.json').read_text(encoding='utf-8'));ids=allocation['remote']['obs_ids'];pairs=[]
with zipfile.ZipFile(r'REMOTE_ROOT/test_evaluation.zip','w',zipfile.ZIP_STORED) as z:
    for run in final['runs']:
        for obs in ids:
            f=s/'evaluation'/run/(obs+'.npz');j=f.with_suffix('.json')
            if f.exists() and j.exists():z.write(f,str(f.relative_to(s)));z.write(j,str(j.relative_to(s)));pairs.append((run,obs))
    for pattern in ('evaluation/BATCH_*.json','delegated/cie/*.json','delegated/data/test_queues/*.json'):
        for f in s.glob(pattern):z.write(f,str(f.relative_to(s)))
print(len(pairs))
""".replace('REMOTE_ROOT',REMOTE).replace('REL_ROOT',REL)
    try:
        ps("& D:/python/python.exe -c '"+code.replace("'","''")+"'",300)
        m=transport();tick=time.perf_counter();subprocess.run(['scp','-i',m.KEY,'-o','BatchMode=yes',m.HOST+':'+remote_zip,str(archive)],check=True);scp_seconds=time.perf_counter()-tick
        with zipfile.ZipFile(archive) as z:
            members={e.filename.replace('\\','/'):e for e in z.infolist() if not e.is_dir()}
            for name,entry in members.items():
                path=Path(name)
                if path.is_absolute() or '..' in path.parts:raise ValueError('Archive path escape')
                target=STUDY/path
                if path.parts[0]=='evaluation' and name.endswith('.npz'):
                    mate=str(path.with_suffix('.json'));j=target.with_suffix('.json')
                    if mate not in members:continue
                    metadata=json.loads(z.read(members[mate]))
                    if metadata.get('execution',{}).get('execution_partition')!='remote':held.append(name+' HOST_PROVENANCE');continue
                    if target.exists() and j.exists():continue
                    if target.exists() or j.exists():held.append(name);continue
                    target.parent.mkdir(parents=True,exist_ok=True)
                    with target.open('xb') as out:out.write(z.read(entry))
                    with j.open('xb') as out:out.write(z.read(members[mate]))
                    fetched+=2
                elif path.parts[0]=='delegated' or path.name.startswith('BATCH_'):
                    # Keep every transferred version even if a running queue's
                    # progress receipt changes; canonical case files stay immutable.
                    content=z.read(entry);snapshot=STUDY/'audit/REMOTE_EVALUATION_SNAPSHOTS'/str(stamp)/path;snapshot.parent.mkdir(parents=True,exist_ok=True);snapshot.write_bytes(content)
                    target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(content);fetched+=1
        status='FETCHED_COMPLETE_PUBLIC_EVALUATIONS'
    except Exception as exc:error=repr(exc);raise
    finally:transfer_receipt('evaluation_fetch',started,status=status,error=error,archive=str(archive),new_files=fetched,held_partial_canonical=held,scp_seconds=scp_seconds,weights_or_truth=False)
    print(json.dumps(dict(status=status,new_files=fetched,held=held,archive=str(archive))))
def fetch(split=None):
    name='anchors'+('_'+split if split else '');remote_zip=REMOTE+'/'+name+'.zip';archive=STUDY/f'audit/REMOTE_{name}.zip'
    # zip snapshots only completed NPZ+receipt pairs, no lock or interruption.
    code="import pathlib,json,zipfile; s=pathlib.Path(r'"+REMOTE+"/project/"+REL+"'); rec=json.loads((s/'data/manifests/online.json').read_text(encoding='utf-8'))['records']; ids={r['obs_id'] for r in rec if "+("r['split']=='"+split+"'" if split else 'True')+"}; z=zipfile.ZipFile(r'"+remote_zip+"','w',zipfile.ZIP_STORED); files=[f for f in (s/'anchors').glob('obs_*.npz') if f.stem.split('_weakreg')[0].split('_fulltrf')[0].split('_discrepancy')[0] in ids and f.with_suffix('.json').exists()]; [(z.write(f,'anchors/'+f.name),z.write(f.with_suffix('.json'),'anchors/'+f.with_suffix('.json').name)) for f in files]; [z.write(f,'anchors/'+f.name) for f in (s/'anchors').glob('BATCH_*.json')]; z.close(); print(len(files))"
    ps("& D:/python/python.exe -c '"+code.replace("'","''")+"'",300)
    m=transport();subprocess.run(['scp','-i',m.KEY,'-o','BatchMode=yes',m.HOST+':'+remote_zip,str(archive)],check=True)
    count=0
    with zipfile.ZipFile(archive) as z:
        for ent in z.infolist():
            if ent.is_dir():continue
            path=Path(ent.filename.replace('\\','/'))
            if path.is_absolute() or '..' in path.parts:raise ValueError('Archive path escape')
            target=STUDY/path
            if target.exists():continue
            target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(z.read(ent));count+=1
    write_canonical_accounting()
    print(json.dumps(dict(fetched_files=count,split=split,archive=str(archive))))

# Called only by the finite local controller, not a daemon or service.
def upload_local_pairs(obs_ids,already_uploaded):
    selected=[]
    for obs_id in obs_ids:
        if obs_id in already_uploaded:continue
        p=STUDY/'anchors'/f'{obs_id}.npz';j=p.with_suffix('.json')
        if not p.exists() or not j.exists():continue
        r=json.loads(j.read_text())
        if r.get('execution_origin',{}).get('machine')!='arm64':continue
        selected.append(obs_id)
    if not selected:return dict(imported=[],existing=[],held=[])
    archive=STUDY/'audit/REMOTE_local_completed.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_STORED) as z:
        for obs_id in selected:
            for ext in ('.npz','.json'):
                p=STUDY/'anchors'/(obs_id+ext);z.write(p,'anchors/'+p.name)
    m=transport();subprocess.run(['scp','-i',m.KEY,'-o','BatchMode=yes',str(archive),m.HOST+':'+REMOTE+'/local_completed.zip'],check=True)
    code="""import pathlib,zipfile,json,shutil
root=pathlib.Path(r'REMOTE_ROOT'); dst=root/'project/REL_ROOT/anchors'; imported=[];existing=[];held=[]
with zipfile.ZipFile(root/'local_completed.zip') as z:
    for name in z.namelist():
        if not name.endswith('.json'):continue
        obs=pathlib.Path(name).stem;p=dst/(obs+'.npz');j=dst/(obs+'.json')
        if p.exists() and j.exists():existing.append(obs);continue
        if p.exists() or j.exists():held.append(obs);continue
        with p.open('xb') as out:out.write(z.read('anchors/'+obs+'.npz'))
        with j.open('xb') as out:out.write(z.read(name))
        imported.append(obs)
print(json.dumps(dict(imported=imported,existing=existing,held=held)))
""".replace('REMOTE_ROOT',REMOTE).replace('REL_ROOT',REL)
    result=ps("& D:/python/python.exe -c '"+code.replace("'","''")+"'",300)
    return json.loads(result.strip().splitlines()[-1])


def write_canonical_accounting():
    """One actual canonical physics cost per observation; no reused-job inflation."""
    from collections import Counter
    records=json.loads((STUDY/'data/manifests/online.json').read_text())['records']
    costs={};cases=Counter();states=Counter();seconds=Counter();all_cases=[]
    for record in records:
        obs=record['obs_id'];p=STUDY/'anchors'/f'{obs}.json'
        if not p.exists() or not p.with_suffix('.npz').exists():continue
        r=json.loads(p.read_text())
        if r.get('config',{}).get('discrepancy_stop',False):continue
        origin=r.get('execution_origin',{})
        hardware='M4_CPU' if origin.get('machine')=='arm64' else 'i9_CPU' if origin else 'i9_CPU_LEGACY_INFERRED_FROM_OBSERVED_REMOTE_EXECUTION'
        costs.setdefault(hardware,Counter()).update(r['counts']);cases[hardware]+=1;states[r['status']]+=1;seconds[hardware]+=r['solve_wall_seconds']
        all_cases.append(dict(obs_id=obs,split=record['split'],origin=hardware,counts=r['counts'],solve_seconds=r['solve_wall_seconds']))
    actions=Counter();batches=[]
    for path in sorted((STUDY/'anchors').glob('BATCH_*.json')):
        batch=json.loads(path.read_text());batch_actions=Counter()
        # Parse actual child actions even from the live remote coordinator,
        # whose old Python objects cannot reload new aggregation logic.
        for job in batch.get('jobs',[]):
            action=job.get('worker_action')
            if not action:
                try:action=json.loads(job.get('stdout','').strip().splitlines()[-1]).get('worker_action','LEGACY_COMPUTED')
                except (ValueError,IndexError):action='LEGACY_UNKNOWN'
            batch_actions[action]+=1
        actions.update(batch_actions)
        actions['PREFLIGHT_REUSED_COMPLETE_PAIR']+=len(batch.get('skipped_existing',[]))
        batches.append(dict(path=str(path),actions=dict(batch_actions),preflight_reused=len(batch.get('skipped_existing',[])),legacy_aggregate_scope='If coordinator was loaded before dedup patch, its aggregate_model_counts is a canonical-case sum including imported M4 costs, not fresh remote physics; use this ledger instead'))
    misconfigured=STUDY/'audit/ANCHOR_PROTOCOL_CORRECTION.json'
    correction=json.loads(misconfigured.read_text()) if misconfigured.exists() else {}
    value=dict(schema='gn-nn-r4-canonical-physics-costs-v1',updated_unix=time.time(),canonical_completed=len(all_cases),canonical_expected=len(records),case_counts_by_origin=dict(cases),canonical_counts_by_origin={k:dict(v) for k,v in costs.items()},canonical_solve_seconds_by_origin=dict(seconds),canonical_statuses=dict(states),job_actions_from_completed_batches=dict(actions),reused_jobs_new_physics_counts={},batch_action_status='PARTIAL_UNTIL_REMOTE_BATCH_RECEIPTS_COMPLETE',batches=batches,misconfigured_attempts=dict(completed=correction.get('completed'),counts=correction.get('completed_counts'),inflight_cost_status=correction.get('inflight_at_stop')),truth_access='NONE',test_truth_status='SEALED_NOT_SCORED',accounting_rule='sum each canonical observation receipt once, separate failed preparation and overhead; never add local and remote batch aggregate totals')
    p=STUDY/'audit/data_canonical_costs.json';tmp=p.with_suffix('.accounting.tmp');tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(p)
    return value

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['upload','upload-code','upload-predictions','run','status','fetch','fetch-evaluation']);p.add_argument('--kind',default='anchors',choices=['anchors','baselines','cie','test','discrepancy-all-scales']);p.add_argument('--split',choices=['train','val','test']);a=p.parse_args()
    if a.action=='upload-predictions':upload_predictions()
    elif a.action.startswith('upload'):upload(a.action=='upload')
    elif a.action=='start':start(a.kind)
    elif a.action=='run':run(a.kind)
    elif a.action=='status':status()
    elif a.action=='fetch-evaluation':fetch_evaluation()
    else:fetch(a.split)

