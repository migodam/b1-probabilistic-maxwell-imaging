#!/usr/bin/env python3
"""ZERO-PHYSICS closed-star norm audit. No solve/inverse/LU or running-chain IO."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
from pathlib import Path
import sys,json,time,traceback,resource
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,write
from p3r4.decision_refresh import decision_from_phi_interval
from p3r4.decision_refresh_tight import (inverse_norm_upper_bounds,neumann_constants,
    phi_interval_from_saved_residual)

JOB='decision_bound_tightening_r1'
RESULT=OUT/'AUDIT/decision_bound_tightening'
STAR=OUT/'results/II/decision_refresh_star_v1'
STAR_JOB=OUT/'II/jobs/decision_refresh_star_v1'
PARENTS=tuple(f'R4DEV_{i:05d}' for i in range(8))


def stats(values):
 a=np.asarray(values,float)
 return dict(n=len(a),minimum=float(a.min()),median=float(np.median(a)),maximum=float(a.max())) if len(a) else dict(n=0)


def decision(interval,base_phi,logu):
 return 'UNKNOWN' if not interval['available'] else decision_from_phi_interval(interval['phi_lower'],interval['phi_upper'],base_phi,logu)


def append(path,row):
 with Path(path).open('a') as f:f.write(json.dumps(row,allow_nan=False)+'\n');f.flush()


def audit():
 if RESULT.exists():raise FileExistsError('audit directory exists; preserve prior development evidence')
 closed=json.loads((STAR_JOB/'CLOSED.json').read_text())
 if closed['status']!='COMPLETE':raise ValueError('only CLOSED star is allowlisted')
 RESULT.mkdir(parents=True)
 book=Book('II',JOB,max_wall=900,analysis_only=True)
 status='COMPLETE';error=None;started=time.perf_counter();allrows=[];summaries=[];missing=[]
 oldmodule=ROOT/'src/p3r4/decision_refresh.py';original_stat=(oldmodule.stat().st_size,oldmodule.stat().st_mtime_ns)
 try:
  write(RESULT/'PROTOCOL.json',dict(scope='ZERO_PHYSICS_CLOSED_STAR_ONLY',source=str(STAR),parents=PARENTS,
    old_module_default_unchanged=True,running_chain_read=False,GT_read=False,
    candidates='only saved star proposal chi, no new/scaled/interpolated physical candidate',
    norms='min(Frobenius,sqrt(induced1*inducedInfinity)); eta0=||I-BL0||F; eta=eta0+b*||diagDchi G||F',
    residuals='same saved Richardson history and GMRES terminal residual scalars; no new approximate-state iteration',
    Full_reference='already-paid Full endpoint predictions from CLOSED star job; no queries',
    missing='record MISSING, never invert or acquire missing columns',
    floating_scope='CONDITIONAL_EXACT_ARITHMETIC_FLOAT64_NOT_DIRECTED_INTERVAL_CERTIFIED',
    no_new_complex_RHS=True,no_new_LU=True,no_linear_solve=True,wall_cap_seconds=900,
    development_revision=1,positive_chain_claim=False))
  with np.load(OUT/'inputs/background/weak.npz',allow_pickle=False) as z:G=z['G'].copy();S=z['S'].copy()
  cfg=json.loads((OUT/'configs/TARGET.json').read_text());sigma=cfg['noise']['sigma'];SF=float(np.linalg.norm(S))
  for parent in PARENTS:
   folder=STAR/parent;inverse_path=folder/'PAID_R0.npz';cache_path=OUT/'results/contexts/cache'/f'{parent}.npz'
   if not inverse_path.exists() or not cache_path.exists():
    item=dict(parent=parent,status='MISSING',inverse=str(inverse_path),cache=str(cache_path));missing.append(item);write(RESULT/f'{parent}_MISSING.json',item);continue
   tick=time.perf_counter()
   with np.load(inverse_path,allow_pickle=False) as z:B=z['R0'].copy();chi0=z['base_chi'].copy();q0=z['base_q'].copy();stored_eta0=float(z['left_inverse_defect_F'])
   with np.load(cache_path,allow_pickle=False) as z:
    if not np.array_equal(z['chi'],chi0) or not np.array_equal(z['q'],q0):raise ValueError('closed inverse/base cache mismatch')
    saved_base_prediction=z['prediction'].copy()
   with np.load(folder/'FROZEN_PROPOSALS.npz',allow_pickle=False) as z:
    chis=z['chi'].copy();qs=z['q'].copy();base_phi=float(z['base_Phi']);logus=z['logu'].copy()
   with np.load(OUT/'inputs/online'/f'{parent}.npz',allow_pickle=False) as z:
    if set(z.files)!={'observation','target_id','sigma'} or str(z['target_id'])!=cfg['target_id'] or float(z['sigma'])!=sigma:raise ValueError('online observation keys or target mismatch')
    y=z['observation'].copy()
   records=[json.loads(line) for line in (folder/'PAIRED_ROWS.jsonl').read_text().splitlines()]
   if len(records)!=len(chis):raise ValueError('closed candidate count mismatch')
   L0=np.eye(len(chi0))-chi0[:,None]*G
   common_tick=time.perf_counter();eta0=float(np.linalg.norm(np.eye(len(chi0))-B@L0));common_seconds=time.perf_counter()-common_tick
   frobenius=float(np.linalg.norm(B));bounds=inverse_norm_upper_bounds(B,known_frobenius=frobenius)
   b_tight=bounds['selected_upper'];parentrows=[];norm_start=time.perf_counter()
   for record in records:
    i=record['candidate'];chi=chis[i]
    if float(record['logu'])!=float(logus[i]):raise ValueError('fixed uniform mismatch')
    delta=float(np.linalg.norm((chi-chi0)[:,None]*G))
    constants={'F':neumann_constants(eta0,frobenius,delta),'min_F_1inf':neumann_constants(eta0,b_tight,delta)}
    endpoint=STAR_JOB/'endpoints'/f"{record['full_endpoint']}.npz"
    with np.load(endpoint,allow_pickle=False) as z:
     if not np.array_equal(z['chi'],chi) or not np.array_equal(z['q'],qs[i]):raise ValueError('Full endpoint candidate mismatch')
     prediction=z['prediction'].copy()
    difference=prediction-y;phi=float(np.vdot(difference,difference).real/sigma**2)
    Full_accept=bool(logus[i]<min(0.,base_phi-phi))
    for arm in ('richardson','gmres'):
     original=record[arm]
     history=original['refinement_history'] if arm=='richardson' else [dict(iteration=0,rho_F=original['rho_F'],residual_white_norm=float(np.sqrt(2*original['approximate_phi'])))]
     variants={}
     for label,c in constants.items():
      entries=[];first=None
      for h in history:
       interval=phi_interval_from_saved_residual(h['residual_white_norm'],h['rho_F'],sigma,SF,c)
       answer=decision(interval,base_phi,float(logus[i]));contains=(not interval['available']) or interval['phi_lower']<=phi<=interval['phi_upper']
       item=dict(iteration=int(h['iteration']),decision=answer,contains_Full=bool(contains),**interval)
       entries.append(item)
       if first is None and answer!='UNKNOWN':first=item
      selected=entries[-1] if first is None else first
      resolved=selected['decision']!='UNKNOWN'
      variants[label]=dict(eta_upper=c['eta_upper'],inverse_upper=c['inverse_upper'],valid_window=c['sufficient_condition_valid'],
       first_resolved_iteration=selected['iteration'] if resolved else None,decision_before_Full=selected['decision'],
       fallback=not resolved,wrong_resolved_decision=bool(resolved and (selected['decision']=='ACCEPT')!=Full_accept),
       interval_misses_on_saved_history=sum(not z['contains_Full'] for z in entries),saved_history_count=len(entries),
       first_interval_width=None if entries[0]['phi_upper'] is None else entries[0]['phi_upper']-entries[0]['phi_lower'],
       selected_interval_width=None if selected['phi_upper'] is None else selected['phi_upper']-selected['phi_lower'],history=entries)
     row=dict(parent=parent,candidate=i,arm=arm,Full_Phi=phi,Full_accept=Full_accept,logu=float(logus[i]),
       delta_operator_F=delta,source_Full=str(endpoint),original_record_fullPhi_difference=phi-record['Full_Phi'],
       F=variants['F'],min_F_1inf=variants['min_F_1inf'])
     append(RESULT/'REPLAY_ROWS.jsonl',row);allrows.append(row);parentrows.append(row)
   norm_seconds=time.perf_counter()-norm_start
   summary=dict(parent=parent,status='COMPLETE',saved_candidates=len(records),norm_bounds=bounds,
      eta0_recomputed=eta0,eta0_saved=stored_eta0,eta0_difference=eta0-stored_eta0,
      common_B_L0_matrix_product_wall_seconds=common_seconds,
      safe_perturbation_F_threshold_F=(1-eta0)/frobenius,safe_perturbation_F_threshold_tight=(1-eta0)/b_tight,
      candidate_replay_wall_seconds=norm_seconds,total_parent_wall_seconds=time.perf_counter()-tick,
      cache_bytes=B.nbytes,source_inverse=str(inverse_path),base_observation_Phi_difference=float(np.linalg.norm(saved_base_prediction-y)**2/sigma**2-base_phi))
   for arm in ('richardson','gmres'):
    rows=[r for r in parentrows if r['arm']==arm];summary[arm]={}
    for label in ('F','min_F_1inf'):
     v=[r[label] for r in rows]
     summary[arm][label]=dict(valid_windows=sum(z['valid_window'] for z in v),fallbacks=sum(z['fallback'] for z in v),
       wrong_resolved_decisions=sum(z['wrong_resolved_decision'] for z in v),interval_misses=sum(z['interval_misses_on_saved_history'] for z in v),
       resolved_iteration=stats([z['first_resolved_iteration'] for z in v if z['first_resolved_iteration'] is not None]),
       first_interval_width=stats([z['first_interval_width'] for z in v if z['first_interval_width'] is not None]))
   summaries.append(summary);write(RESULT/f'{parent}_SUMMARY.json',summary)
   book.counts['saved_inverse_arrays']+=1;book.counts['saved_candidates']+=len(records);book.counts['added_norm_reduction_entries']+=bounds['profile']['row_column_sum_entries']
   print(json.dumps(dict(parent=parent,tightening_factor=bounds['tightening_factor'],richardson=summary['richardson'])),flush=True)
  aggregate=dict(status='COMPLETE_DEVELOPMENT_AUDIT_NOT_CHAIN_RESULT',parents=summaries,missing=missing,
    independent_source_parents=len(summaries),saved_proposals=sum(s['saved_candidates'] for s in summaries),
    norm_tightening_factors=stats([s['norm_bounds']['tightening_factor'] for s in summaries]),
    added_norm_setup_wall_seconds=sum(s['norm_bounds']['profile']['wall_seconds'] for s in summaries),
    common_eta0_matrix_product_wall_seconds=sum(s['common_B_L0_matrix_product_wall_seconds'] for s in summaries),
    new_Maxwell_states=0,new_complex_RHS=0,new_LU=0,new_inverse_columns=0,new_physical_evaluations=0,
    all_linear_solve_calls=0,GT_read=False,running_chain_read_or_modified=False,source_closed_receipts_modified=False,
    old_default_module_unchanged=(oldmodule.stat().st_size,oldmodule.stat().st_mtime_ns)==original_stat,
    optional_revision='src/p3r4/decision_refresh_tight.py; default old class unchanged',
    norm_scope='CONDITIONAL_EXACT_ARITHMETIC_FLOAT64_NOT_DIRECTED_INTERVAL_CERTIFIED',
    positive_chain_speedup='NOT_TESTED_NOT_CLAIMED',current_potential_interval_handling='unchanged; chain validity remains caller responsibility')
  for arm in ('richardson','gmres'):
   rows=[z for z in allrows if z['arm']==arm];aggregate[arm]={}
   for label in ('F','min_F_1inf'):
    v=[z[label] for z in rows]
    aggregate[arm][label]=dict(valid_windows=sum(z['valid_window'] for z in v),fallbacks=sum(z['fallback'] for z in v),
      wrong_resolved_decisions=sum(z['wrong_resolved_decision'] for z in v),interval_misses=sum(z['interval_misses_on_saved_history'] for z in v),
      resolved_iteration=stats([z['first_resolved_iteration'] for z in v if z['first_resolved_iteration'] is not None]),
      first_interval_width=stats([z['first_interval_width'] for z in v if z['first_interval_width'] is not None]))
   ratios=[z['F']['first_interval_width']/z['min_F_1inf']['first_interval_width'] for z in rows if z['F']['first_interval_width'] is not None and z['min_F_1inf']['first_interval_width'] is not None and z['min_F_1inf']['first_interval_width']>0]
   aggregate[arm]['first_width_improvement_ratio']=stats(ratios)
   aggregate[arm]['earlier_resolution_candidates']=sum(z['min_F_1inf']['first_resolved_iteration'] is not None and z['F']['first_resolved_iteration'] is not None and z['min_F_1inf']['first_resolved_iteration']<z['F']['first_resolved_iteration'] for z in rows)
  write(RESULT/'RESULTS.json',aggregate)
 except BaseException as exc:
  status='ERROR';error=repr(exc);write(RESULT/'ERROR.json',dict(error=error,traceback=traceback.format_exc()));raise
 finally:
  book.close(status,error)
  write(RESULT/'FINAL_JOB_COST.json',dict(status=status,error=error,elapsed_seconds=time.perf_counter()-started,
      new_Maxwell_states=0,new_complex_RHS=0,new_LU=0,receipt=str(book.jobdir/'CLOSED.json')))

if __name__=='__main__':audit()
