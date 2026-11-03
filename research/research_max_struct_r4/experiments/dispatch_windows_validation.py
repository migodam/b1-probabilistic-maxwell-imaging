"""Single public development validation after the18-task queue has finished."""
from pathlib import Path
import json,time,zipfile
from dispatch_confirmation import ROOT,REMOTE,remote,ps,transfer,save,extract_windows_archive

def main():
    tick=time.perf_counter()
    done=ps("Test-Path '"+REMOTE+"/results/confirmation/runner_summary.json'")
    if done['returncode'] or done['stdout'].strip()!='True':
        raise RuntimeError('Original queue not yet closed; do not share its allocation')
    files=[(ROOT/'experiments/validate_windows_development.py',REMOTE+'/research/research_max_struct_r4/experiments/validate_windows_development.py'),
        (ROOT.parent/'research_max_struct_r2/data_public/obs_b20d5a55677d4eca.npz',REMOTE+'/data_public/development_validation.npz')]
    uploaded=[]
    for source,dest in files:
        receipt=transfer(str(source),remote.HOST+':'+dest);uploaded.append(receipt)
        if receipt['returncode']:
            save('WINDOWS_VALIDATION_UPLOAD_FAILURE.json',dict(receipts=uploaded));raise RuntimeError('Upload failed')
    command="Set-Location '"+REMOTE+"'; & 'D:/python/python.exe' research/research_max_struct_r4/experiments/validate_windows_development.py --public data_public/development_validation.npz --config configs/confirmation_frozen_PUBLIC.json --output validation/development_windows; exit $LASTEXITCODE"
    run=ps(command,timeout=900)
    save('WINDOWS_VALIDATION_RUN.json',dict(upload=uploaded,execution=run,truth_access='NONE_PUBLIC_ONLY',reason='critical positive cross-platform validation; no algorithm retuning'))
    print(json.dumps(run),flush=True)
    if run['returncode']:return
    archived=ps("Compress-Archive -Path '"+REMOTE+"/validation/development_windows' -DestinationPath 'D:/AI/MAX_STRUCT_R4_20261103_validation.zip'",timeout=60)
    target=ROOT/'audit/windows_validation_package.zip'
    downloaded=transfer(remote.HOST+':D:/AI/MAX_STRUCT_R4_20261103_validation.zip',str(target))
    if archived['returncode'] or downloaded['returncode']:
        save('WINDOWS_VALIDATION_COLLECTION_FAILURE.json',dict(archive=archived,download=downloaded));return
    destination=ROOT/'results/cross_platform_remote';destination.mkdir()
    extract_windows_archive(target,destination)
    save('WINDOWS_VALIDATION_TRANSPORT.json',dict(archive=archived,download=downloaded,total_transport_and_run_wall_seconds=time.perf_counter()-tick,physical_cost_in_nested_RUN_RECEIPT=True))
    print(json.dumps(dict(status='COLLECTED_DEVELOPMENT_VALIDATION_UNSCORED',path=str(destination))),flush=True)
if __name__=='__main__':main()
