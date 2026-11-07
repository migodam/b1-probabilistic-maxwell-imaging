"""Run finite VAL/TEST public-only full TRF and weakreg1e-5 control queues."""
import argparse,subprocess,sys
from pathlib import Path
STUDY=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=2);p.add_argument('--core-mask',type=int);p.add_argument('--split',action='append',choices=['val','test']);p.add_argument('--kind',action='append',choices=['fulltrf','weakreg1e-5','discrepancy','discrepancy-all-scales']);a=p.parse_args()
for kind in a.kind or ['discrepancy','weakreg1e-5']:
    command=[sys.executable,str(STUDY/'experiments/run_anchors.py'),'--workers',str(a.workers)]
    for split in a.split or ['val','test']:command+=['--split',split]
    if kind=='weakreg1e-5':command+=['--weakreg','1e-5']
    elif kind=='discrepancy':command+=['--discrepancy']
    elif kind=='discrepancy-all-scales':command+=['--discrepancy-all-scales']
    else:command+=['--full-trf']
    if a.core_mask:command+=['--core-mask',str(a.core_mask)]
    subprocess.run(command,check=True)
