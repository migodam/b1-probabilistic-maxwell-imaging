"""Predict an anonymous public cohort using a frozen checkpoint; no GT reads."""
import argparse,json,sys
from pathlib import Path
S=Path(__file__).resolve().parents[1];sys.path.insert(0,str(S/'src'))
from predict import predict_anchors
p=argparse.ArgumentParser();p.add_argument('--checkpoint',type=Path,required=True);p.add_argument('--run',required=True);p.add_argument('--split',choices=['dev','val','test'],required=True);p.add_argument('--device',default='auto');a=p.parse_args()
records=json.loads((S/'data/manifests/online.json').read_text())['records']
paths=[S/'anchors'/f"{r['obs_id']}.npz"for r in records if r['split']==a.split]
out=S/'predictions'/a.run
receipts=predict_anchors(a.checkpoint,paths,S/'data/operators.npz',out,a.device)
print(json.dumps(dict(status='PREDICTED_PUBLIC_INPUT_ONLY',observations=len(receipts),output=str(out))))
