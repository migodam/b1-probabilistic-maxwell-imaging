"""Parent-owned transport for the already installed finite Windows CPU queue."""
from pathlib import Path
import argparse, importlib.util, json, shutil, subprocess, time, zipfile
ROOT=Path(__file__).resolve().parents[1]
workspace=ROOT.parents[5]
spec=importlib.util.spec_from_file_location('existing_remote',workspace/'Gaussian/A9/code/remote.py')
remote=importlib.util.module_from_spec(spec);spec.loader.exec_module(remote)
REMOTE='D:/AI/MAX_STRUCT_R4_20261103'
AUDIT=ROOT/'audit/confirmation'

def save(name,entry):
    (AUDIT/name).write_text(json.dumps(entry,indent=2))

def ps(script,timeout=30):
    t=time.perf_counter();p=remote.powershell(script,timeout=timeout)
    return dict(returncode=p.returncode,stdout=p.stdout.decode('utf8','replace'),
        stderr=p.stderr.decode('utf8','replace'),wall_seconds=time.perf_counter()-t)

def transfer(source,target):
    t=time.perf_counter()
    p=subprocess.run(['scp','-i',remote.KEY,'-o','BatchMode=yes','-o','ConnectTimeout=8',source,target],capture_output=True,timeout=120)
    return dict(returncode=p.returncode,stdout=p.stdout.decode('utf8','replace'),stderr=p.stderr.decode('utf8','replace'),wall_seconds=time.perf_counter()-t)

def prepare():
    probe=ps(r"""$jobs=Get-Process -Name python,pythonw -ErrorAction SilentlyContinue | Select-Object Id,ProcessName,CPU,WorkingSet64
$gpu=& nvidia-smi --query-gpu=memory.used,utilization.gpu --format=csv,noheader
[ordered]@{python_processes=$jobs;gpu=$gpu;destination_exists=(Test-Path 'D:/AI/MAX_STRUCT_R4_20261103')}|ConvertTo-Json -Depth 4 -Compress""")
    save('PRELAUNCH_ACTIVITY.json',probe)
    print(json.dumps(probe),flush=True)
    if probe['returncode']:return
    state=json.loads(probe['stdout'])
    if state['python_processes'] or state['destination_exists']:
        print('Activity or destination requires parent inspection; no launch or overwrite.');return
    package=ROOT/'audit/remote_package'
    path=package/'configs/confirmation_manifest.json'
    manifest=json.loads(path.read_text());manifest.update(status='FROZEN_PARENT_EXECUTION_AUTHORIZED',execution_authorized=True)
    path.write_text(json.dumps(manifest,indent=2))
    save('EXECUTION_MANIFEST.json',manifest)
    archive=Path(shutil.make_archive(str(ROOT/'audit/confirmation_public_package'),'zip',package))
    entry=transfer(str(archive),remote.HOST+':D:/AI/MAX_STRUCT_R4_20261103_package.zip')
    save('UPLOAD_RECEIPT.json',dict(**entry,archive_bytes=archive.stat().st_size,truth_access='NONE_PUBLIC_ONLY'))
    if entry['returncode']:raise RuntimeError('Upload failed; receipt saved')
    entry=ps("Expand-Archive -LiteralPath 'D:/AI/MAX_STRUCT_R4_20261103_package.zip' -DestinationPath '"+REMOTE+"'; if ($?) { Write-Output 'PUBLIC_PACKAGE_READY' }",timeout=60)
    save('REMOTE_PREPARATION_RECEIPT.json',entry);print(json.dumps(entry),flush=True)

def run():
    entry=ps("Set-Location '"+REMOTE+"'; & 'D:/python/python.exe' remote_runner.py run configs/confirmation_manifest.json; exit $LASTEXITCODE",timeout=2400)
    save('REMOTE_RUN_RECEIPT.json',dict(**entry,actual_cpu_workers=4,actual_cuda_workers=0,
        scope='foreground SSH elapsed includes remote process startup, all finite queue tasks and shutdown; not sum of worker wall'))
    print(json.dumps(entry),flush=True)

def collect():
    entry=ps("Compress-Archive -Path '"+REMOTE+"/results','"+REMOTE+"/predictions' -DestinationPath 'D:/AI/MAX_STRUCT_R4_20261103_results.zip'",timeout=60)
    save('REMOTE_RESULT_ARCHIVE_RECEIPT.json',entry)
    if entry['returncode']:raise RuntimeError('Archive failed')
    target=ROOT/'audit/confirmation_result_package.zip'
    entry=transfer(remote.HOST+':D:/AI/MAX_STRUCT_R4_20261103_results.zip',str(target))
    save('DOWNLOAD_RECEIPT.json',entry)
    if entry['returncode']:raise RuntimeError('Download failed')
    destination=ROOT/'results/confirmation_remote'
    destination.mkdir()
    with zipfile.ZipFile(target) as archive:archive.extractall(destination)
    print(json.dumps(dict(status='COLLECTED_UNSCORED',path=str(destination),truth_access='NONE')),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run','collect']);args=p.parse_args()
    globals()[args.action]()
