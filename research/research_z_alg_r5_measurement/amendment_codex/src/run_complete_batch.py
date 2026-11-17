from complete_quotient import *
import traceback
rows=[];start=time.process_time()
for case in ['T0','T1','T2','T3','T4','T5']:
 for arm in ('clean','noisy'):
  path=OUT/'results'/f'complete_{case}_{arm}_128.json'
  if path.exists():result=json.loads(path.read_text())
  else:
   try:result=run(case,arm)
   except Exception as ex:
    result=dict(case=case,arm=arm,status='FAILED',error=repr(ex),traceback=traceback.format_exc());save(OUT/'results'/f'complete_{case}_{arm}_failure.json',result);print(result,flush=True)
  rows.append({k:result.get(k) for k in ['case','arm','quotient_complete','material_box_complete','legal_count','CPU_seconds','error']})
save(OUT/'results/complete_batch.json',dict(rows=rows,CPU_seconds_new_batch=time.process_time()-start))
