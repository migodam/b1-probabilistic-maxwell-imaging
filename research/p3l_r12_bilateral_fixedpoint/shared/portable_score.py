"""Score published clouds with bundled reuse inputs, without a physical call."""
from pathlib import Path
import argparse,json,subprocess,sys
P=Path(__file__).resolve().parents[1]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--dataset',choices=['P1','N2'],required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
    out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True)
    source=P/'audit'/('P1_PREDICTION_MANIFEST.json' if a.dataset=='P1' else 'N2_PREDICTION_MANIFEST.json')
    m=json.loads(source.read_text());truth=[]
    for g in m['groups']:
        oid=g['parent'];data='P1' if a.dataset=='P1' else 'N'
        g['cloud_path']=str(P/'results'/a.dataset/oid/g['method']/'CLOUD.npz')
        g['anchor_path']=str(P/'inputs/reuse'/data/oid/'ANCHOR.npz')
        g['roi_path']=str(P/'inputs/FIXED_HOST_SETUP.npz')
        if g.get('initial_frame_path'):g['initial_frame_path']=str(P/'inputs/reuse'/data/oid/'ANCHOR_FRAME.npz')
        if not any(t['parent']==oid for t in truth):truth.append(dict(dataset=data,parent=oid,path=str(P/'inputs/reuse'/data/oid/'truth.npz')))
    (out/'manifest.json').write_text(json.dumps(m,indent=2));(out/'truth.json').write_text(json.dumps(dict(truth=truth),indent=2))
    cmd=[sys.executable,str(P/'agent_c/offline_score.py'),'score','--manifest',str(out/'manifest.json'),'--truth-manifest',str(out/'truth.json'),'--out',str(out/'scores')]
    if a.dataset=='N2':cmd+=['--blind-read',str(P/'audit/N2_MAIN_BLIND_READ.json')]
    subprocess.run(cmd,check=True)
if __name__=='__main__':main()
