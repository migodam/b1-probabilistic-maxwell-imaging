from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run import OUT,write

def main():
 job=json.loads((OUT/'TRAINING_JOB.json').read_text());rows=[]
 for r in job['records']:
  for stage in ['provider','head']:
   rows.append(dict(arm=r['arm'],seed=r['seed'],stage=stage,optimizer_updates=r[f'{stage}_updates'],selected_update=r[f'{stage}_selected_update'],fit_and_validation_wall_seconds=r[f'{stage}_training_seconds'],device=r['device'],CPU_threads=1,train_parents=1024,selection_validation_parents=256,physics_RHS=0,receipt=str(OUT/'jobs'/f"{r['arm']}_s{r['seed']}.json")))
 write(OUT/'COST_LEDGER.json',dict(status=job['status'],actual_training_jobs=len(rows),records=rows,sum_fit_and_validation_seconds=sum(r['fit_and_validation_wall_seconds'] for r in rows),conservative_MPS_reserved_job_wall_seconds=job['wall_seconds'],peak_rss_bytes=job['peak_rss_bytes'],peak_recorded_MPS_driver_bytes=max(r['MPS_driver_allocated_bytes'] for r in job['records']),upper_bound_peakRSS_plus_peak_recorded_MPS_driver_bytes=job['peak_rss_bytes']+max(r['MPS_driver_allocated_bytes'] for r in job['records']),physics_RHS=0,historical_costs='original Gaussian teacher, TRAIN32768, full1024/256 anchors and original head/student training remain sunk dependencies; reused and not repriced here',memory_scope='max sampled driver allocation plus process peak RSS; not continuous whole-system proof',new_joint_simulation_or_MCMC=False))
if __name__=='__main__':main()
