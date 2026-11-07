"""Attach authoritative external process CPU while retaining self timers as nested."""
import json
from pathlib import Path

def attach_external_cpu(directory,cost_path):
 directory=Path(directory);cost_path=Path(cost_path);receipt=directory/'RECEIPT.json'
 if not receipt.exists():return
 r=json.loads(receipt.read_text());c=json.loads(cost_path.read_text());external=c['whole_child_CPU_seconds']
 if 'CPU_seconds' in r and 'internal_self_CPU_seconds' not in r:r['internal_self_CPU_seconds']=r['CPU_seconds']
 if 'secondary_seal_CPU_seconds' in r:r['internal_secondary_seal_CPU_seconds']=r.pop('secondary_seal_CPU_seconds')
 r['CPU_seconds']=external;r['CPU_scope']='EXTERNAL full-process wait4/RUSAGE_CHILDREN CPU including imports and finalization; internal timers nested, never added';r['external_process_cost_receipt']=str(cost_path)
 tmp=receipt.with_suffix('.external.tmp');tmp.write_text(json.dumps(r,indent=2)+'\n');tmp.replace(receipt)
