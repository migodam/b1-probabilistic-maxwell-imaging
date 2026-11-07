"""Finite R4 public-only SSH jobs; existing transport, isolated own remote root."""
from pathlib import Path
import argparse,importlib.util,os,subprocess,sys,zipfile,json,time
STUDY=Path(__file__).resolve().parents[1];REPO=STUDY.parents[1];WORKSPACE=REPO.parents[3]
REMOTE=os.environ.get('GN_NN_REMOTE_ROOT','D:/AI/GN_NN_R4_20261107');REL='research/research_gn_nn_preservation_r4'
MASK=20560 # Pcores4,6,12,14; historical85 overlapped another session's0,2.
def transport():
    spec=importlib.util.spec_from_file_location('existing_ssh',WORKSPACE/'Gaussian/A9/code/remote.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
def ps(script,timeout=120):
    r=transport().powershell(script,timeout=timeout)
    out=r.stdout.decode('utf8','replace');err=r.stderr.decode('utf8','replace');print(out,end='',flush=True)
    if err:print(err,file=sys.stderr,flush=True)
    if r.returncode:raise RuntimeError(f'Remote returncode{r.returncode}: {err[-1000:]}')
    return out
def upload(include_data=True):
    archive=STUDY/'audit/REMOTE_public.zip';files=list((STUDY/'src').glob('*.py'))+list((STUDY/'experiments').glob('*.py'))
    for folder in ('research_z_alg_r4_diversity','research_gn_nn_completion_r3','research_gn_search_space_r2'):
        files+=list((REPO/'research'/folder/'src').glob('*.py'))
    if include_data:files+=[STUDY/'data/operators.npz',STUDY/'data/manifests/online.json']+list((STUDY/'data/public').glob('*.npz'))
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=1) as z:
        for f in files:
            if f.exists():z.write(f,f.relative_to(REPO))
    ps(f"New-Item -ItemType Directory -Force -Path '{REMOTE}' | Out-Null")
    m=transport();subprocess.run(['scp','-i',m.KEY,'-o','BatchMode=yes',str(archive),m.HOST+':'+REMOTE+'/public.zip'],check=True)
    ps(f"Expand-Archive -LiteralPath '{REMOTE}/public.zip' -DestinationPath '{REMOTE}/project' -Force",300)
    print('PUBLIC_ONLY_UPLOAD_COMPLETE; NO_TRUTH')
def start(kind='anchors'):
    log=REMOTE+'/'+kind+'.log';err=REMOTE+'/'+kind+'.err';script_path=REMOTE+'/'+kind+'.ps1'
    args='--split train --split val --split test --workers 4 --core-mask '+str(MASK) if kind=='anchors' else '--split val --split test --workers 4 --core-mask '+str(MASK)
    target='run_anchors.py' if kind=='anchors' else 'run_baselines.py'
    body=f"Set-Location '{REMOTE}/project/{REL}'\n$env:OMP_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OPENBLAS_NUM_THREADS='1'; $env:VECLIB_MAXIMUM_THREADS='1'\n$env:PYTHONPATH='D:/AI/ZALG_R4_DIVERSITY_20261030/project/vendor'\n[System.Diagnostics.Process]::GetCurrentProcess().ProcessorAffinity=[IntPtr]{MASK}\n& D:/python/python.exe -B -u experiments/{target} {args}\nexit $LASTEXITCODE\n"
    escaped=body.replace("'","''")
    result=ps(f"Set-Content -LiteralPath '{script_path}' -Value '{escaped}' -Encoding UTF8; $p=Start-Process powershell.exe -ArgumentList '-NoProfile','-File','{script_path}' -RedirectStandardOutput '{log}' -RedirectStandardError '{err}' -PassThru; Write-Output $p.Id")
    (STUDY/f'audit/REMOTE_{kind}_start.json').write_text(json.dumps(dict(kind=kind,pid=result.strip(),remote=REMOTE,affinity_mask=MASK,cores=[4,6,12,14],historical_mask85_not_used='avoid existing unrelated CPU tasks on0,2',public_truth_access='NONE',start_unix=time.time()),indent=2)+'\n')
def run(kind='anchors'):
    # Windows OpenSSH tears down detached children on session close. Keep one
    # foreground SSH session alive; local exec runs it in the background.
    args='--split train --split val --split test --workers 4 --core-mask '+str(MASK) if kind=='anchors' else '--split val --split test --workers 4 --core-mask '+str(MASK)
    target='run_anchors.py' if kind=='anchors' else 'run_baselines.py'
    body=f"Set-Location '{REMOTE}/project/{REL}'\n$env:OMP_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OPENBLAS_NUM_THREADS='1'; $env:VECLIB_MAXIMUM_THREADS='1'\n$env:PYTHONPATH='D:/AI/ZALG_R4_DIVERSITY_20261030/project/vendor'\n[System.Diagnostics.Process]::GetCurrentProcess().ProcessorAffinity=[IntPtr]{MASK}\n& D:/python/python.exe -B -u experiments/{target} {args}\nexit $LASTEXITCODE\n"
    ps("Set-Content -LiteralPath '"+REMOTE+'/'+kind+".ps1' -Value '"+body.replace("'","''")+"' -Encoding UTF8")
    ps(f"& powershell.exe -NoProfile -File '{REMOTE}/{kind}.ps1' *> '{REMOTE}/{kind}.log'; exit $LASTEXITCODE",21600)

def status():
    ps(f"$d='{REMOTE}/project/{REL}/anchors'; (Get-ChildItem $d -Filter 'obs_*.npz' -ErrorAction SilentlyContinue | Measure-Object).Count; Get-Content '{REMOTE}/anchors.log' -Tail 2 -ErrorAction SilentlyContinue; Get-Content '{REMOTE}/anchors.err' -Tail 4 -ErrorAction SilentlyContinue")
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
    p=argparse.ArgumentParser();p.add_argument('action',choices=['upload','upload-code','run','status','fetch']);p.add_argument('--kind',default='anchors',choices=['anchors','baselines']);p.add_argument('--split',choices=['train','val','test']);a=p.parse_args()
    if a.action.startswith('upload'):upload(a.action=='upload')
    elif a.action=='start':start(a.kind)
    elif a.action=='run':run(a.kind)
    elif a.action=='status':status()
    else:fetch(a.split)

