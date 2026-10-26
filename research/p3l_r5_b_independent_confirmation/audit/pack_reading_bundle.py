"""Small scientific reading bundle; complete raw evidence stays in four archives."""
from pathlib import Path
import zipfile,subprocess,json,argparse
P=Path(__file__).resolve().parents[1];REPO=P.parents[1]
D=REPO.parents[2]/'research/delegated/r5bconfirm_delivery_20261027'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=D/'Research_P3LR5_B_Confirm_Results_and_Pro_Handoff.zip');a=ap.parse_args()
    source=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip()
    files=[]
    for f in P.rglob('*'):
        if not f.is_file() or 'raw' in f.relative_to(P).parts or '__pycache__' in f.parts:continue
        if f.suffix in ['.md','.py','.png','.json','.csv'] and f.stat().st_size<2500000 and f.name not in ['PUBLICATION_RECEIPT.json','WEB_READBACK.json','RELEASE_UPLOAD_RECEIPT.json']:
            files.append(f)
    files+=list((P/'figures/confirmation').glob('*SOURCE.npz'))
    files+=[P/'raw/offline_scores/PARENT_LOSSES.npz',P/'raw/offline_scores/PARENT_METRICS.npz']
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(a.output,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        z.writestr('R5B_Confirm_Results/READ_ME_FIRST.txt',
          'Completed independent512-parent confirmation. Scientific/delivery source commit '+source+'\n'+
          '80% retention statistical gate HOLD; no new training. Fullraw/models in four lossless Release archives.\n'+
          'Start with START_HERE.md and GPT_PRO_HANDOFF/START_HERE.md. Existing validation was never relabeled as new TEST.\n')
        for f in sorted(set(files)):z.write(f,'R5B_Confirm_Results/'+str(f.relative_to(P)))
    with zipfile.ZipFile(a.output) as z:assert z.testzip() is None
    print(json.dumps(dict(bytes=a.output.stat().st_size,files=len(set(files))+1,source=source,new_RHS=0)))
if __name__=='__main__':main()
