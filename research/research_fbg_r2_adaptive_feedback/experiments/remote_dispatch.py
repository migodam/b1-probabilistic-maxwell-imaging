"""Existing SSH transport; finite foreground CUDA process, no service."""
from pathlib import Path
import argparse,importlib.util,json,subprocess,sys,zipfile
HERE=Path(__file__).resolve().parents[1]
WORKSPACE=HERE.parents[5]
REMOTE='D:/AI/FBG_R2_20261030'
HOST='jayzh@100.121.97.39'
KEY=str(Path.home()/'.ssh/heatcryst_xinan_ed25519')

def transport(script,timeout):
    path=WORKSPACE/'Gaussian/A9/code/remote.py'
    spec=importlib.util.spec_from_file_location('existing_remote',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    result=module.powershell(script,timeout=timeout)
    sys.stdout.write(result.stdout.decode('utf-8','replace'))
    sys.stderr.write(result.stderr.decode('utf-8','replace'))
    if result.returncode:raise SystemExit(result.returncode)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['upload','run','fetch'])
    ap.add_argument('name',nargs='?',default='stage1');a=ap.parse_args()
    if a.action=='upload':
        archive=HERE/'audit/remote_source.zip'
        with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as out:
            for folder in ('src','experiments'):
                for file in (HERE/folder).glob('*.py'):out.write(file,file.relative_to(HERE))
            for file in (HERE/'data_contract').glob('*remote*.json'):out.write(file,file.name)
        subprocess.run(['scp','-i',KEY,'-o','BatchMode=yes',str(archive),HOST+':'+REMOTE+'/source.zip'],check=True)
        transport("Expand-Archive -LiteralPath 'D:\\AI\\FBG_R2_20261030\\source.zip' -DestinationPath 'D:\\AI\\FBG_R2_20261030' -Force",30)
    elif a.action=='run':
        script="""
Set-Location D:\AI\FBG_R2_20261030
$env:CUDA_PIPELINE_ROOT='D:\AI\ZALG_CUDA_PIPELINE_20261030'
$env:OMP_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OPENBLAS_NUM_THREADS='1'
[System.Diagnostics.Process]::GetCurrentProcess().ProcessorAffinity=[IntPtr]21845
& D:\python\python.exe -u -m src.adaptive_solver NAME_remote.json
""".replace('NAME',a.name)
        transport(script,7200)
    else:
        transport("Compress-Archive -Path 'D:\\AI\\FBG_R2_20261030\\results\\NAME' -DestinationPath 'D:\\AI\\FBG_R2_20261030\\NAME.zip' -Force".replace('NAME',a.name),30)
        target=HERE/'audit'/f'{a.name}_receipts.zip'
        subprocess.run(['scp','-i',KEY,'-o','BatchMode=yes',HOST+':'+REMOTE+f'/{a.name}.zip',str(target)],check=True)
        with zipfile.ZipFile(target) as archive:
            for item in archive.infolist():
                if not item.is_dir():
                    dest=HERE/'results'/item.filename.replace('\\','/')
                    dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(archive.read(item))

if __name__=='__main__':main()
