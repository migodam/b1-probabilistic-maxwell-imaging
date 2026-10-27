"""CPU1 deployment costs, fixed first64 VAL, excluding common historical anchor generation."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
from pathlib import Path
import sys,time,resource
sys.path.insert(0,str(Path(__file__).resolve().parent))
from inference import Action
from run import OUT,HIST,OLD,write
import torch,numpy as np

def main():
 torch.set_num_threads(1);va=dict(np.load(HIST/'prepared_validation.npz'));y=va['y'][:64]
 with np.load(OLD/'research/p3l_r4_conditional_transport/inputs/strong_context_1024_256/anchor_validation.npz') as z:online={k:z[k][:64] for k in ['observation','residual','q0','chi0','J']}
 records={}
 for arm in ['original_feature32768','feature_control1024','direct_gt','feature_action','gt_action']:
  t=time.perf_counter();a=Action(arm);cold=time.perf_counter()-t
  for _ in range(3):a.predict(y,online)
  rec={}
  for label,fn in [('provider',lambda:a.moments(y)),('provider_feature_head',lambda:a.predict(y,online))]:
   raw=[]
   for _ in range(20):
    t=time.perf_counter();fn();raw.append(time.perf_counter()-t)
   rec[label]=dict(raw_batch64_seconds=raw,median_seconds_per_parent=float(np.median(raw)/64),p10_seconds_per_parent=float(np.quantile(raw,.1)/64),p90_seconds_per_parent=float(np.quantile(raw,.9)/64))
  rec['cold_load_seconds']=cold;records[arm]=rec
 write(OUT/'DEPLOYMENT_COST.json',dict(status='COMPLETE',device='cpu',CPU_threads=1,batch=64,parents='historical VAL first64 fixed',warmups=3,repeats=20,records=records,shared_anchor_cost='excluded, historical physics still required and disclosed',physics_RHS=0,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
if __name__=='__main__':main()
