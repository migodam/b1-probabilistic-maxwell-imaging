import argparse,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('action',choices=['audit','pilot','refine','infer','summarize']);p.add_argument('--root',default='.');p.add_argument('--observation');p.add_argument('--output');p.add_argument('--model',choices=['refined','full-grid'],default='refined');a=p.parse_args()
root=Path(a.root).resolve();config=json.loads((root/'configs/protocol.json').read_text())
if a.action=='audit':
 from .audit import run;run(root,config)
elif a.action=='pilot':
 from .pilot import run;run(root,config)
elif a.action=='refine':
 from .refine import run;run(root,config)
elif a.action=='infer':
 from .cli_infer import run
 if not a.observation or not a.output:raise SystemExit('--observation and --output required')
 run(root,config,a.observation,a.output,a.model)
else:
 from .budget import Ledger
 print(json.dumps(Ledger(root/'results',config['physics_new_cap']).summary(),indent=2))
