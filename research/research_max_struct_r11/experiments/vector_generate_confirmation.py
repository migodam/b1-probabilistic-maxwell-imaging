"""Root-frozen new3D observation generation only; never called by inverse."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
 os.environ[name]='1'
import time
ENTRY=time.perf_counter()
from pathlib import Path
import sys,json,traceback,gc
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
import numpy as np
import torch
from research.research_max_struct_r10.experiments.vector_run import initialize,dump,stamp
from research.research_max_struct_r10.src.vector_materials import material_from_design

def main():
 cfg=json.loads((ROOT/'configs/vector_confirmation_contract.json').read_text());design=json.loads((ROOT/'configs/vector_confirmation_material_design.json').read_text())
 if cfg['status']!='ROOT_FROZEN_BEFORE_NEW_GENERATION' or not cfg['private_generation_authorized'] or not design['private_generation_authorized']:raise RuntimeError('Missing root scientific freeze')
 out=ROOT/'delegated/vector/confirmation_generation';out.mkdir(parents=True,exist_ok=False)
 private=ROOT/'delegated/vector/private_confirmation_truth';public=ROOT/'delegated/vector/confirmation_observations'
 private.mkdir(parents=True,exist_ok=False);public.mkdir(parents=True,exist_ok=False)
 torch.set_num_threads(1);torch.set_num_interop_threads(1)
 cases={c['design_id']:c for c in design['cases']};ids=cfg['case_ids'];cross=cfg['crossgrid']['design_ids']
 observations=[(key,12) for key in ids]+[(key,18) for key in cross]
 streams=np.random.SeedSequence(design['seed']).spawn(len(observations));stream_ids={pair:i for i,pair in enumerate(observations)}
 config={k:cfg[k] for k in ('physics','acquisition')};config.update(execution_backend='cpu-dense',execution=dict(memory_limit_bytes=2_000_000_000))
 counts={};records=[];costs=[];model=None;failure=None;status='FAILED'
 try:
  for side,chosen in ((12,ids),(18,cross)):
   model=initialize(config,out/('grid'+str(side)),side)
   for key in chosen:
    wall,cpu=time.perf_counter(),time.process_time();before=dict(model.counts)
    m=material_from_design(cases[key],model.cpu.points,model.cpu.h)
    mt=torch.as_tensor(m,dtype=torch.float64,device=model.device)
    J=model.solve(mt,torch.repeat_interleave(mt,3)[:,None]*model.inc,label='new_confirmation_generation_'+key+'_grid'+str(side))
    Q=model.inc+model.green(J);clean=(model.S@J).cpu().numpy()
    folder=private/key;folder.mkdir(parents=True,exist_ok=True);truth=folder/('grid'+str(side)+'.npz')
    if truth.exists():raise RuntimeError('Write-once new truth already exists')
    np.savez_compressed(truth,material=m,J=J.cpu().numpy(),Q=Q.cpu().numpy(),Y_clean=clean)
    sigma=.03*np.linalg.norm(clean)/np.sqrt(clean.size);ordinal=stream_ids[(key,side)];rng=np.random.default_rng(streams[ordinal])
    noise=sigma/np.sqrt(2)*(rng.standard_normal(clean.shape)+1j*rng.standard_normal(clean.shape))
    public_id=key+('_matched12' if side==12 else '_cross18to12');dest=public/public_id;dest.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(dest/'observations.npz',Y=clean+noise)
    dump(dest/'PUBLIC.json',dict(public_id=public_id,model='iid_circular_complex',sigma=float(sigma),expected_relative_Frobenius=.03,
      convention='Eabs(n)^2=sigma^2',whitening=False,realization_renormalized=False,inverse_side=12,generation_side=side,
      source_positions_resolution_independent=True,truth_access_required_for_inverse=False))
    row=dict(public_id=public_id,private_file=stamp(truth),public_files=[stamp(dest/'PUBLIC.json'),stamp(dest/'observations.npz')],
      CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,counts_delta={k:int(v-before.get(k,0)) for k,v in model.counts.items() if v-before.get(k,0)},
      noise_private_child_ordinal=ordinal,truth_generation_only=True)
    costs.append(row);dump(out/'GENERATED_SO_FAR.json',dict(observations=costs,complete=len(costs)==len(observations)))
    print(json.dumps(dict(event='NEW_PUBLIC_OBSERVATION_GENERATED',public_id=public_id,CPU_seconds=row['CPU_seconds'],wall_seconds=row['wall_seconds'])),flush=True)
   for k,v in model.counts.items():counts[k]=counts.get(k,0)+int(v)
   records.extend(model.solve_records);del model,mt,J,Q,clean,noise;model=None;gc.collect()
  dump(out/'GENERATION_SEAL.json',dict(status='ALL10_PUBLIC_OBSERVATIONS_GENERATED_PRIVATE_FIELDS_SEALED',
    observations=costs,independent_material_designs=8,matched_observations=8,crossgrid_observations=2,
    source_contract=stamp(ROOT/'configs/vector_confirmation_contract.json'),design_contract=stamp(ROOT/'configs/vector_confirmation_material_design.json'),
    inverse_truth_access=False,private_arrays_not_inverse_inputs=True))
  status='COMPLETE'
 except BaseException as exc:failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc());print(json.dumps(failure),flush=True)
 finally:
  if model is not None:
   for k,v in model.counts.items():counts[k]=counts.get(k,0)+int(v)
   records.extend(model.solve_records)
  dump(out/'PROCESS_COST.json',dict(status=status,failure=failure,platform='Mac_CPU',CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-ENTRY,
   wholeprocess_includes_imports_setup_exports=True,counts=counts,actual_solve_records=records,GT_access='GENERATION_ONLY',inverse_and_score_not_run=True,
   observations_completed=len(costs),additive=True))
 return 0 if status=='COMPLETE' else 1

if __name__=='__main__':raise SystemExit(main())
