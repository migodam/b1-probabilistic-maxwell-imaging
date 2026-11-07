"""One already exposed16 saved endpoint: exact/mismatch/missing-gradient logic only."""
from pathlib import Path
import json,time
import numpy as np
from endpoint_bookkeeping import endpoint,load_material_gradient,strict_match
HERE=Path(__file__).absolute().parent;ROOT=HERE.parents[2]
def main():
 wall=time.perf_counter();seal_path=ROOT/'delegated/diagnostics/bothblocks_confirmation/ALL80_TERMINAL_PREDICTION_SEAL.json';seal=json.loads(seal_path.read_text())
 if seal['status']!='ALL80_TERMINAL_BEFORE_ANY_CONFIRMATION_GT_SCORING':raise RuntimeError('Old exposed16 public seal required')
 e=seal['entries'][0];cap_path=Path(e['prediction_path']);result=endpoint(cap_path.parent)
 if not result['available']:raise RuntimeError('The single predefined saved record has no matched gradient; record this limitation without trying another object')
 cap,_=load_material_gradient(cap_path);m,g=load_material_gradient(Path(result['gradient_source']));positive=strict_match(cap,m,g);changed=m.copy();changed[0]=np.nextafter(changed[0],np.inf);negative=strict_match(cap,changed,g);missing=strict_match(cap,m,None)
 assert positive[0] and not negative[0] and not missing[0]
 receipt=dict(status='PASS_SINGLE_EXPOSED_RECORD_MATCH_LOGIC',task_id=e['task']['task_id'],source_seal=str(seal_path),source_cap=str(cap_path),endpoint_result=result,positive_exact_match=positive,one_ULP_material_difference_rejected=negative,exact_material_missing_gradient_rejected=missing,objects_examined=1,GT_arrays_read=0,physical_or_optimizer_modules_imported=False,new_forward=0,new_adjoint=0,new_inverse=0,CPU_seconds_near_write=time.process_time(),wall_seconds=time.perf_counter()-wall,scope='Bookkeeping matching logic on one existing public saved record only; no repeated identity or physical validation')
 (HERE/'ONE_RECORD_LOGIC_CHECK.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
