"""Public postseal endpoint metadata only; never imports physics or reads GT."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,csv,time,sys
import numpy as np
HERE=Path(__file__).absolute().parent;ROOT=HERE.parents[2];REPO=ROOT.parents[1]
FIELDS=('projected_gradient_inf','raw_gradient_inf','exact_lower_active','exact_upper_active','near_lower_1e_minus8','near_upper_1e_minus8','wrong_sign_exact_lower','wrong_sign_exact_upper','wrong_sign_exact_total')
SAVED=('R11_ACTUAL_PAID_INCUMBENT.npz','R11_LAST_PAID_STATE.npz')
def read(p):return json.loads(Path(p).read_text())
def resolve(p):
 p=Path(str(p).replace('\\','/'));p=p if p.is_absolute() else REPO/p
 if {'sealed_truth','data_private','private'} & {x.lower() for x in p.parts}:raise ValueError('Private path forbidden')
 return p
def load_material_gradient(p):
 with np.load(p,allow_pickle=False) as z:
  m=z['m'].copy() if 'm' in z.files else None;g=z['gradient'].copy() if 'gradient' in z.files else None
 return m,g
def strict_match(cap,m,g):
 if m is None:return False,'NO_SAVED_MATERIAL'
 if cap.shape!=m.shape or not np.array_equal(cap,m):return False,'MATERIAL_NOT_ELEMENTWISE_EXACT'
 if g is None:return False,'EXACT_MATERIAL_BUT_NO_SAVED_GRADIENT'
 if g.shape!=cap.shape or np.iscomplexobj(g) or not np.isfinite(g).all():return False,'INVALID_SAVED_REAL_GRADIENT'
 return True,'EXACT_MATERIAL_AND_SAVED_GRADIENT'
def endpoint(directory):
 directory=resolve(directory);result=dict(output_dir=str(directory),available=False,reason=None,gradient_source=None,candidate_match_records=[],gradient_units='Unscaled originalPhi gradient including original prior, not1e6 optimizer units',new_forward=0,new_adjoint=0,GT_access='NONE',**{k:None for k in FIELDS})
 capfile=directory/'R11_CAP_PREDICTION.npz'
 if not capfile.is_file():result['reason']='NO_SEALED_CAP_MATERIAL';return result
 cap,_=load_material_gradient(capfile)
 if cap is None or cap.shape!=(576,) or np.iscomplexobj(cap) or not np.isfinite(cap).all():result['reason']='INVALID_CAP_REAL_MATERIAL';return result
 picked=None
 for name in SAVED:
  p=directory/name;item=dict(path=str(p),exists=p.is_file(),material_exact=False,gradient_exists=False,reason='SAVED_RECORD_MISSING')
  if p.is_file():
   try:
    m,g=load_material_gradient(p);ok,why=strict_match(cap,m,g);item.update(material_exact=m is not None and m.shape==cap.shape and np.array_equal(m,cap),gradient_exists=g is not None,reason=why)
    if ok and picked is None:picked=(p,g)
   except (OSError,ValueError,KeyError) as exc:item['reason']='INVALID_SAVED_RECORD: '+repr(exc)
  result['candidate_match_records'].append(item)
 if picked is None:result['reason']='NO_EXACT_CAP_MATERIAL_WITH_VALID_SAVED_GRADIENT';return result
 p,g=picked;lo,hi=.02,20.;lower=cap==lo;upper=cap==hi;wrong_low=lower&(g<0);wrong_high=upper&(g>0)
 result.update(available=True,reason='EXACT_SAVED_ENDPOINT_MATCH',gradient_source=str(p),projected_gradient_inf=float(np.max(np.abs(cap-np.clip(cap-g,lo,hi)))),raw_gradient_inf=float(np.max(np.abs(g))),exact_lower_active=int(lower.sum()),exact_upper_active=int(upper.sum()),near_lower_1e_minus8=int(np.count_nonzero(np.abs(cap-lo)<=1e-8)),near_upper_1e_minus8=int(np.count_nonzero(np.abs(cap-hi)<=1e-8)),wrong_sign_exact_lower=int(wrong_low.sum()),wrong_sign_exact_upper=int(wrong_high.sum()),wrong_sign_exact_total=int(wrong_low.sum()+wrong_high.sum()))
 return result
def gate(p,expected):
 x=read(p);status=x.get('status')
 accepted=(status=='ALL96_PUBLIC_DEPLOYMENTS_SEALED_BEFORE48_GT' and expected==288 and x.get('primary_tasks_terminal')==288 and x.get('truth_access')=='NONE') or (status=='PASS_ALL288_TERMINAL' and expected==288 and x.get('tasks')==288 and x.get('truth_opened') is False)
 if not accepted:raise RuntimeError('Requires completed all288 public seal/gate metadata; no endpoint arrays read')
 return x
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--manifest',type=Path,required=True);ap.add_argument('--public-seal',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();wall=time.perf_counter();tasks=read(a.manifest)['tasks']
 if len(tasks)!=288 or len({t['task_id'] for t in tasks})!=288 or len({t['case_id'] for t in tasks})!=48:raise ValueError('Original288 scientific identities/48materials required')
 seal=gate(a.public_seal,288)
 # The production metadata-only recovery helper preserves scientific IDs while
 # selecting actual valid retry directories; no interrupted directory fallback.
 sys.path.insert(0,str(ROOT/'experiments'))
 from twod_recovery import active_task
 a.output.mkdir(parents=True,exist_ok=False);rows=[]
 for original in tasks:
  t=active_task(original);d=resolve(t['output_dir']);r=read(d/'RECEIPT.json') if (d/'RECEIPT.json').is_file() else {}
  item=dict(task_id=t['task_id'],case_id=t['case_id'],method=t.get('material_readout',{}).get('label',t['arm']),engineering_attempt_id=t.get('engineering_attempt_id'),scientific_original_output=original['output_dir'],receipt_status=r.get('status'))
  if r.get('cap_prediction_status')!='SEALED_PAID_INCUMBENT':item.update(available=False,reason='NO_PUBLIC_SEALED_CAP_RECEIPT',**{k:None for k in FIELDS})
  else:item.update(endpoint(d))
  rows.append(item)
 (a.output/'ENDPOINT_BOOKKEEPING.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with (a.output/'ENDPOINT_BOOKKEEPING.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in rows)
 receipt=dict(status='POST_ALL288_PUBLIC_SEAL_ENDPOINT_BOOKKEEPING_COMPLETE',created_utc=datetime.now(timezone.utc).isoformat(),source_public_seal=str(a.public_seal),seal_status=seal['status'],tasks=288,independent_materials=48,available=sum(r['available'] for r in rows),unavailable=sum(not r['available'] for r in rows),new_forward=0,new_adjoint=0,GT_access='NONE',selected_candidate_changed=False,CPU_seconds_near_write=time.process_time(),wall_seconds=time.perf_counter()-wall,scope='Posthoc numerical endpoint description only; no optimization or convergence/success gate')
 (a.output/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
