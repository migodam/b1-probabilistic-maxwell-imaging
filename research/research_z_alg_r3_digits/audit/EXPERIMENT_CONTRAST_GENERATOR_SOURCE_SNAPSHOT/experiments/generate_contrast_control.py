"""Post-GT two-case controlled contrast sensitivity; no main replacement."""
import json,time,uuid,shutil
from pathlib import Path
import numpy as np
from scipy.linalg import solve
ROOT=Path(__file__).resolve().parents[1]
SOURCES=['obs_a63df7f20566','obs_42c48e795766']
def dump(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2)+'\n')
def main():
 solverfreeze=ROOT/'audit/FULLJ_DIAGNOSTIC_CONTRAST_SOURCE_SNAPSHOT/FREEZE.json'
 if not solverfreeze.exists():raise RuntimeError('Solver/config freeze required before new Y')
 out=ROOT/'public_data/contrast_control';truth=ROOT/'private_truth/contrast_control'
 if out.exists() and list(out.glob('*.npz')):raise RuntimeError('Refusing overwrite diagnostic observations')
 paths=['experiments/generate_contrast_control.py','contrast_control.py','audit/FULLJ_DIAGNOSTIC_CONTRAST_CONFIG.json','src/solver.py','src/physics.py','src/multiscale.py','src/__init__.py'];snap=ROOT/'audit/EXPERIMENT_CONTRAST_GENERATOR_SOURCE_SNAPSHOT';files=[]
 for name in paths:
  p=ROOT/name;q=snap/name;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q);files.append(dict(path=name,snapshot=str(q.relative_to(ROOT)),bytes=p.stat().st_size))
 dump(ROOT/'audit/EXPERIMENT_CONTRAST_GENERATOR_FREEZE.json',dict(frozen_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),before_new_Y=True,post_main_GT=True,selected_source_ids=SOURCES,selection='two public high-misfit cases, postGT diagnostic; not a new blind benchmark',formula='m_low=1+1.3*(m_GT-min(m_GT))/(max(m_GT)-min(m_GT)) on ALL pixels',files=files,scope='physical regime sensitivity only; cannot replace original8-18 targets'))
 ops=np.load(ROOT/'public_data/operators.npz');G,H,U0=ops['G'],ops['H'],ops['U0'];mapping=json.loads((ROOT/'truth_after_freeze/holdout/mapping_seedless.json').read_text());lookup={r['observation_id']:r for r in mapping['records']};records=[];total=0;started=time.perf_counter();out.mkdir(parents=True,exist_ok=True);truth.mkdir(parents=True,exist_ok=True)
 for source in SOURCES:
  old=np.load(ROOT/'truth_after_freeze/holdout'/f'{source}.npz');high=old['m'];bg=float(high.min());peak=float(high.max());m=1+1.3*(high-bg)/(peak-bg);ident='contrast_'+uuid.uuid4().hex[:12];t=time.perf_counter();A=np.eye(len(m),dtype=complex)-m[:,None]*G;b=m[:,None]*U0;j=solve(A,b,assume_a='gen',check_finite=False);Y=H@j;elapsed=time.perf_counter()-t;total+=elapsed
  np.savez_compressed(out/f'{ident}.npz',Y=Y,sigma=0.,observation_id=np.array(ident));np.savez_compressed(truth/f'{ident}.npz',m=m,mask=old['mask'],Y_clean=Y)
  records.append(dict(observation_id=ident,digit=lookup[source]['digit'],source_observation_id=source,source_material_min=bg,source_material_max=peak,transformed_material_min=float(m.min()),transformed_material_max=float(m.max()),generation_seconds=elapsed,state_relative_residual=float(np.linalg.norm(A@j-b)/np.linalg.norm(b))))
 dump(truth/'mapping_PRIVATE.json',dict(scope='POSTGT physical contrast sensitivity, same2morphologies notmainrecovery',records=records))
 receipt=dict(scope='POSTGT controlled physical-regime sensitivity; notmainreplacement/notnewblindholdout',count=2,public_ids=[r['observation_id'] for r in records],additional_PDE_LU=2,additional_complex_RHS=24,generation_seconds=total,total_preparation_seconds=time.perf_counter()-started,public_payload_keys=['Y','sigma','observation_id'],same_operators_geometry_frequency=True,noise=0,mask_preserved=True,solver_truth_access='NONE_NEW_ANONYMOUS_Y_ONLY')
 dump(ROOT/'audit/EXPERIMENT_GENERATION_CONTRAST_CONTROL.json',receipt);print(json.dumps(receipt))
if __name__=='__main__':main()
