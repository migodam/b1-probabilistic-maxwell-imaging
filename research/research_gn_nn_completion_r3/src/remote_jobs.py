"""Finite remote jobs using the project's existing SSH transport; no service."""
from pathlib import Path
import argparse, importlib.util, os, subprocess, sys, zipfile

STUDY=Path(__file__).resolve().parents[1]
REPO=STUDY.parents[1]
WORKSPACE=REPO.parents[3]
REMOTE=os.environ.get('GN_NN_REMOTE_ROOT','D:/AI/GN_NN_R3_20261107')
REL='research/research_gn_nn_completion_r3'

def transport():
    path=WORKSPACE/'Gaussian/A9/code/remote.py'
    spec=importlib.util.spec_from_file_location('existing_ssh',path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m

def ps(script,timeout=120):
    r=transport().powershell(script,timeout=timeout)
    sys.stdout.write(r.stdout.decode('utf8','replace'));sys.stdout.flush()
    sys.stderr.write(r.stderr.decode('utf8','replace'));sys.stderr.flush()
    if r.returncode:raise SystemExit(r.returncode)

def upload(include_data=True):
    archive=STUDY/'audit/remote_public.zip'
    files=list((STUDY/'src').glob('*.py'))+list((STUDY/'experiments').glob('*.py'))
    files+=list((REPO/'research/research_z_alg_r4_diversity/src').glob('*.py'))
    files+=list((REPO/'research/research_gn_search_space_r2/src').glob('*.py'))
    if include_data:
        files+=[STUDY/'data/operators.npz',STUDY/'data/manifests/online.json']
        files+=list((STUDY/'data/manifests').glob('train_expansion_*.json'))
        files+=list((STUDY/'data/public').glob('*.npz'))
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=1) as z:
        for f in files:
            if f.exists():z.write(f,f.relative_to(REPO))
    ps("New-Item -ItemType Directory -Force -Path '"+REMOTE+"' | Out-Null")
    m=transport()
    subprocess.run(['scp','-i',m.KEY,'-o','BatchMode=yes',str(archive),m.HOST+':'+REMOTE+'/public.zip'],check=True)
    ps("Expand-Archive -LiteralPath '"+REMOTE+"/public.zip' -DestinationPath '"+REMOTE+"/project' -Force",120)
    print('Uploaded only public observations/operators and source; no truth files.')

def run(args):
    script=f"""
Set-Location '{REMOTE}/project/{REL}'
$env:OMP_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OPENBLAS_NUM_THREADS='1'; $env:VECLIB_MAXIMUM_THREADS='1'
$env:PYTHONPATH='D:/AI/ZALG_R4_DIVERSITY_20261030/project/vendor'
[System.Diagnostics.Process]::GetCurrentProcess().ProcessorAffinity=[IntPtr]85
& D:/python/python.exe -B -u experiments/run_anchors.py {args}
exit $LASTEXITCODE
"""
    ps(script,21600)

def fetch():
    ps(f"Compress-Archive -Path '{REMOTE}/project/{REL}/anchors' -DestinationPath '{REMOTE}/anchors.zip' -Force",300)
    archive=STUDY/'audit/remote_anchors.zip';m=transport()
    subprocess.run(['scp','-i',m.KEY,'-o','BatchMode=yes',m.HOST+':'+REMOTE+'/anchors.zip',str(archive)],check=True)
    with zipfile.ZipFile(archive) as z:
        for entry in z.infolist():
            if entry.is_dir():continue
            target=STUDY/entry.filename.replace('\\','/')
            target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(z.read(entry))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['upload','upload-code','run','fetch']);p.add_argument('--args',default='--workers 4');a=p.parse_args()
    if a.action.startswith('upload'):upload(a.action=='upload')
    elif a.action=='run':run(a.args)
    else:fetch()
