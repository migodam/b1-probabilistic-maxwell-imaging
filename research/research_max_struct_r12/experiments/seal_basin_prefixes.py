"""Choose only paid original-Phi candidates, no ground truth access."""
from pathlib import Path
import argparse,json
import numpy as np
def main():
 p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args()
 records=[]
 for folder in sorted(a.root.glob('*__*')):
  if not (folder/'RECEIPT.json').exists():continue
  h=json.loads((folder/'history.json').read_text())
  with np.load(folder/'prediction.npz',allow_pickle=False) as z:
   selected={}
   for cap in [30,60,120,240]:
    eligible=[r for r in h if r['cpu']<=cap]
    if not eligible:continue
    r=min(eligible,key=lambda x:x['Phi'])
    selected[f'm_{cap}']=z['candidates'][r['candidate']].copy()
    records.append(dict(run=folder.name,cap=cap,candidate=r['candidate'],paid_cpu=r['cpu'],Phi=r['Phi']))
  np.savez_compressed(folder/'paid_prefix.npz',**selected)
 (a.root/'PREFIX_SEAL.json').write_text(json.dumps(dict(selection='minimum original Phi at paid body CPU prefix, no GT',records=records),indent=2))
if __name__=='__main__':main()
