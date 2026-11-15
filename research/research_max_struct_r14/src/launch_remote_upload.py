"""Capture gh token in memory and pass only via SSH stdin to Windows helper.

No archive data travels back to this host. No implicit helper copy, release
creation or publication. Run only when the parent explicitly starts uploading.
"""
from __future__ import annotations
import argparse,base64,json,subprocess,time
from pathlib import Path


def ps_quote(value):return "'"+str(value).replace("'","''")+"'"


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--host',required=True);p.add_argument('--identity-file',type=Path)
    p.add_argument('--remote-python',required=True);p.add_argument('--remote-script',required=True)
    p.add_argument('--archive-dir',required=True);p.add_argument('--manifest',required=True);p.add_argument('--remote-receipt',required=True)
    p.add_argument('--release-id',required=True,type=int);p.add_argument('--local-receipt',type=Path,required=True)
    a=p.parse_args()
    if a.host.startswith('-') or any(c.isspace() for c in a.host):raise ValueError('Invalid SSH host')
    if a.local_receipt.exists():raise FileExistsError('Use a new local launcher receipt')
    arguments=[a.remote_script,'--release-id',str(a.release_id),'--archive-dir',a.archive_dir,
               '--manifest',a.manifest,'--receipt',a.remote_receipt]
    command='& '+ps_quote(a.remote_python)+' '+ ' '.join(ps_quote(x) for x in arguments)+'; exit $LASTEXITCODE'
    encoded=base64.b64encode(command.encode('utf-16-le')).decode('ascii')
    ssh=['ssh','-T','-o','BatchMode=yes']
    if a.identity_file is not None:ssh+=['-i',str(a.identity_file)]
    ssh += [a.host,'powershell.exe','-NoProfile','-NonInteractive','-EncodedCommand',encoded]
    # gh auth token stdout remains captured only; neither stdout nor stderr of
    # this auth command is written or printed, including on a failure.
    start=time.perf_counter();cpu=time.process_time()
    authentication=subprocess.run(['gh','auth','token','--hostname','github.com'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
    if authentication.returncode or not authentication.stdout.strip():raise RuntimeError('Local GitHub authentication command failed; credential output suppressed')
    token=authentication.stdout.strip();authentication=None
    try:
        remote=subprocess.run(ssh,input=token+b'\n',stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
        # Defense in depth: even unexpected SSH/remote output cannot persist
        # the exact in-memory token in a receipt or terminal output.
        safeout=remote.stdout.replace(token,b'[REDACTED]').decode('utf8',errors='replace')
        safeerr=remote.stderr.replace(token,b'[REDACTED]').decode('utf8',errors='replace')
        record=dict(status='REMOTE_HELPER_EXITED' if remote.returncode==0 else 'REMOTE_HELPER_FAILED',
            returncode=remote.returncode,host=a.host,remote_release_id=a.release_id,
            remote_receipt=a.remote_receipt,remote_stdout=safeout,remote_stderr=safeerr,
            controller_CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start,
            token_transport='ghcapturememory to SSHstdin only;notargv/environment/disk',
            archive_bytes_transferred_to_local=0,release_created=False,publication_performed=False)
        a.local_receipt.parent.mkdir(parents=True,exist_ok=True)
        with a.local_receipt.open('x',encoding='utf8') as stream:json.dump(record,stream,indent=2);stream.write('\n')
        print(json.dumps({k:record[k] for k in ('status','returncode','wall_seconds','archive_bytes_transferred_to_local')}))
        return remote.returncode
    finally:token=None


if __name__=='__main__':raise SystemExit(main())
