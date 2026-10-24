#!/usr/bin/env python3
"""Registered pre-outcome applicability grid; prepare-only absent parent release."""
from pathlib import Path
import sys,json,os,argparse,subprocess,time,csv
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'tools')]
import numpy as np
import run_r4_interval_chain_global_gmres as common
from p3r4.core import OUT,Book,write
from p3r4.evaluation import checkpoint_identity

engine=common.engine;EXPERIMENT='interval_chain_boundary_v1';ARMS=common.ARMS
PARENTS=tuple(f'R4DEV_{i:05d}' for i in range(56,60));BETAS=(.005,.01,.02);SEEDS=(20261025181,20261025182,20261025183)
CASES=tuple(dict(parent=p,parent_index=i,beta=beta,seed=seed,steps=256,case_id=f'b{beta:g}_s{seed}_{p}') for beta in BETAS for seed in SEEDS for i,p in enumerate(PARENTS))
CAPS=dict(full=73728,richardson=75424,global_gmres=75424,global_gmres_accuracy=75424)
ROOT_RESULT=OUT/'results/II'/EXPERIMENT;AUDIT=OUT/'AUDIT/interval_boundary';PREPARED=AUDIT/'FROZEN_PROTOCOL.json';MARKER=OUT/'results/II/interval_boundary_pair'
engine.EXPERIMENT=EXPERIMENT;engine.PARENTS=PARENTS;engine.RESULT=ROOT_RESULT;engine.CAP=CAPS;engine.STEPS=256
engine.LEDGER_STAGE='II';engine.CASE_DESIGNS=CASES;engine.REUSE_ORIGINAL_R0_PER_PARENT=True;engine.MATERIAL_MOVE_AUDIT=True
engine.CONTEXT_DIR=OUT/'results/contexts/dev';engine.CACHE_DIR=OUT/'results/contexts/cache'
engine.FROZEN_DESIGN_EXTENSION=dict(engine.FROZEN_DESIGN_EXTENSION,experiment=EXPERIMENT,case_designs=CASES,beta=BETAS,seed=SEEDS,source_pCN_beta=BETAS,grid_order='beta,seed,parent ascending',
 cold_R0='576RHS peroriginalparent perapproxarm once,2304RHS perarm; shared immutable across9samebasecases; never rebase',
 prospective_expansion=True,scope='failureboundary,notnewposteriorreference; no benefit requirement awayfrombeta.002',full_receives_same_original_anchor_cache_permissions=True)

def protocol():
 return dict(status='FROZEN_PREPARED_NOT_RELEASED',experiment=EXPERIMENT,campaign_id='P3LR4-20261025',ledger_stage='II',parents=PARENTS,
  betas=BETAS,seeds=SEEDS,steps_per_case=256,cases=CASES,arms=ARMS,transitions_per_arm=9216,
  target_config=json.loads((OUT/'configs/TARGET.json').read_text()),target_identity=checkpoint_identity(OUT/'configs/TARGET.json'),
  combined_complex_RHS_cap=300000,arm_complex_RHS_caps=CAPS,whole_arm_wall_seconds=1200,RSS_cap_bytes=3*1024**3,
  original_R0_setup_RHS_perapproxarm=2304,original_R0_lifetime='one immutable original perparent;9radii/seedcases share; fresh endpointbasis percase, no crosschi basis/rebase',
  random_stream='SeedSequence([seed,parent_index]) reset eachradius/arm;66normal thenuniform;cases start originalanchor independently',
  shared_rights='allFull/approx have samepaid originalanchor/cache/receiverdata;eachapprox pays ownR0;R0GSmatrixapply0newsolve,actualnewcomplexLS1RHS reservedfirst',
  strict_fallback='failed inverse bound or unresolved MH aftermax4 calls nativeFull8sources; ifbudget exhausted STOP without decision or relaxation',
  accuracy_arm='bothendpoints whiteerrorupper<=.05 andPhiwidth<=.05 beforeMH,thenmax4orFull',
  partial_policy='fixednontransferable armcaps; retain allattempts/failures/rejections/fields/proofs; markpartialandunrun cases,never quotarescue/rerun',
  independence='4originalDEVparents,9settings/seedreplicates are NOT36independentparents',
  coordinator_wall='external wholecoordinator inclallarm starts,sealing,state-field publication,postclose comparison/report;NOT an additional childwall debit',
  source_modules=['tools/run_r4_interval_chain_boundary.py','tools/run_r4_interval_chain_fast.py','tools/run_r4_interval_chain_global_gmres.py'],
  release_required=['matchingreleased:true','prepared_protocol_identity','queue_dependency_closed_receipt'],
  no_GT=True,posterior_reference='NONE',float64_certification=False,no_forced_speedup=True,new_physics_in_preparation=0)

def check_release(path):
 r=json.loads(path.read_text())
 if r.get('released') is not True or r.get('experiment')!=EXPERIMENT or r.get('campaign_id')!='P3LR4-20261025':raise PermissionError('matching explicit parent release required')
 p=json.loads(PREPARED.read_text())
 if r.get('prepared_protocol_identity')!=checkpoint_identity(PREPARED) or p['target_config']!=json.loads((OUT/'configs/TARGET.json').read_text()):raise ValueError('frozen protocol/target changed')
 dependency=Path(r.get('queue_dependency_closed_receipt','')).resolve()
 if not dependency.is_file() or dependency.name!='CLOSED.json' or OUT not in dependency.parents:raise ValueError('parentVIcontroller must close before released physical scheduling')
 if json.loads(dependency.read_text()).get('status') not in ('COMPLETE','ERROR','FAILED','STOPPED'):raise ValueError('dependency is not closed')
 return r

def rows(path):return [json.loads(s) for s in path.read_text().splitlines()] if path.exists() else []
def postclose():
 if (ROOT_RESULT/'COMPARISON.json').exists():raise FileExistsError('preserve first audit')
 b=Book('II',EXPERIMENT+'_postclose',max_wall=600,analysis_only=True);status='COMPLETE';error=None
 try:
  timing=json.loads((ROOT_RESULT/'ARM_PROCESS_WALL.json').read_text());gates={};table=[];case_rows=[];misses=[];proof_count=0;actual_groups=0
  # Actual reservations, including a failed/last partial action, are retained.
  import fcntl
  with (OUT/'budget.lock').open('a') as lock:
   fcntl.flock(lock,fcntl.LOCK_SH)
   with (OUT/'RESOURCE_LEDGER.jsonl').open('rb') as f:f.seek(0,2);prefix=f.tell()
   fcntl.flock(lock,fcntl.LOCK_UN)
  ledger={a:dict(RHS=0,small=0,physics=0.,cases={}) for a in ARMS};jobmap={EXPERIMENT+'_'+a:a for a in ARMS}
  with (OUT/'RESOURCE_LEDGER.jsonl').open('rb') as f:
   while f.tell()<prefix:
    line=f.readline()
    if not line.endswith(b'\n'):break
    if EXPERIMENT.encode() not in line:continue
    event=json.loads(line);arm=jobmap.get(event['job'])
    if arm is None:continue
    cost=event.get('RHS',0);isLS=event.get('role')=='samecache_global_GMRES_reduced_LS';ledger[arm]['RHS']+=cost;ledger[arm]['small']+=cost if isLS else 0;ledger[arm]['physics']+=event.get('physics_wall_seconds',0.)
    case_id=event.get('context',{}).get('case_id')
    if case_id is not None:
     entry=ledger[arm]['cases'].setdefault(case_id,dict(RHS=0,small=0,states=0));entry['RHS']+=cost;entry['small']+=cost if isLS else 0;entry['states']+=event.get('states',0)
  for arm in ARMS:
   path=ROOT_RESULT/arm/'ARM_COMPLETE.json'
   if not path.exists():table.append(dict(arm=arm,status='NOT_RUN',RHS=0,states=0,LU=0));continue
   g=json.loads(path.read_text());gates[arm]=g;c=json.loads(Path(g['receipt']).read_text())
   if c['status']!=g['status']:raise ValueError('closed arm status mismatch')
   if gates.get('full') and g['input_identities']!=gates['full']['input_identities']:raise ValueError('original inputs differ')
   for identity in g['all_array_seals']:
    if checkpoint_identity(identity['path'])!=identity:raise ValueError('sealed array changed')
   small=0;setup=0;setupwall=0.;fallbacks=0;accepted=0;completed=0;ref_total=0;proposal_values=[];move_values=[]
   for case in CASES:
    d=ROOT_RESULT/arm/case['case_id'];records=rows(d/'CHAIN_REFERENCES.jsonl');summary=json.loads((d/'SUMMARY.json').read_text()) if (d/'SUMMARY.json').exists() else None
    completed+=len(records);ref_total+=len(records);accepted+=sum(z['accepted'] for z in records);fallbacks+=sum(z['Full_calls'] for z in records) if arm!='full' else 0
    props=[z['proposal_material_RMS'] for z in records];moves=[z['material_move_RMS'] for z in records];proposal_values+=props;move_values+=moves
    if summary:setup+=summary['setup'].get('RHS',0)
    setupmeta=summary['setup'] if summary is not None else json.loads((d/'SETUP.json').read_text()) if (d/'SETUP.json').exists() else None
    if setupmeta:setupwall+=setupmeta.get('wall_including_R0_boundcopy_and_persistence',0.)
    for p in rows(d/'GLOBAL_GMRES_REDUCED_PROOFS.jsonl'):
     if p['event']=='global_GMRES_refined' and not p['solution_cache_hit'] and p['iteration']>0:small+=1
    coverage='COMPLETE' if len(records)==256 and (d/'CHAIN.npz').exists() else 'PARTIAL' if d.exists() else 'NOT_RUN'
    cost=ledger[arm]['cases'].get(case['case_id'],dict(RHS=0,small=0,states=0))
    case_rows.append(dict(arm=arm,case_id=case['case_id'],parent=case['parent'],beta=case['beta'],seed=case['seed'],coverage=coverage,completed_steps=len(records),
     accepts=sum(z['accepted'] for z in records),Full_fallbacks=sum(z['Full_calls'] for z in records) if arm!='full' else 0,
     case_RHS=cost['RHS'],reduced_LS_RHS=cost['small'],Full_attempts=cost['states'],setup_RHS=cost['RHS']-cost['small']-8*cost['states'],
     recorded_setup_seconds=None if setupmeta is None else setupmeta.get('wall_including_R0_boundcopy_and_persistence',0.),R0_reuse_hit=None if setupmeta is None else setupmeta.get('shared_original_R0_cache_hit'),
     mean_proposal_material_RMS=float(np.mean(props)) if props else None,mean_actual_material_move_RMS=float(np.mean(moves)) if moves else None,
     max_proposal_material_RMS=max(props) if props else None,terminal_state_RMS_from_originalanchor=records[-1]['state_material_RMS_from_anchor'] if records else None,
     case_wholewall_seconds=None if summary is None else summary['wholeparentwall']))
    goldpath=ROOT_RESULT/'full'/case['case_id']/'CHAIN.npz'
    if arm=='full' or not goldpath.exists():continue
    with np.load(goldpath,allow_pickle=False) as z:gold={k:z[k] for k in z.files}
    goldrefs=rows(goldpath.parent/'CHAIN_REFERENCES.jsonl');goldphi={'initial':float(gold['current_phi_before'][0])};goldphi.update({f'c{i:04d}':float(v) for i,v in enumerate(gold['candidate_phi'])})
    for i,r in enumerate(records):
     if r['logu']!=goldrefs[i]['logu'] or r['accepted']!=bool(gold['accept'][i]) or not np.array_equal(r['x'],gold['x_path'][i+1]):misses.append(dict(arm=arm,case=case['case_id'],step=i,kind='random/flag/path'))
     with np.load(d/'attempts'/f'c{i:04d}.npz',allow_pickle=False) as a,np.load(goldpath.parent/'attempts'/f'c{i:04d}.npz',allow_pickle=False) as a0:
      if any(not np.array_equal(a[k],a0[k]) for k in ('x','q','chi','logu')):misses.append(dict(arm=arm,case=case['case_id'],step=i,kind='candidate'))
     if arm=='global_gmres_accuracy' and not (r['history'][-1]['current_accuracy'] and r['history'][-1]['candidate_accuracy']):misses.append(dict(arm=arm,case=case['case_id'],step=i,kind='accuracy'))
    for p in rows(d/'ENDPOINT_PROOFS.jsonl'):
     proof_count+=1;lo,hi=p['Phi_interval'];lo=0. if lo is None else lo;hi=float('inf') if hi is None else hi
     if p['endpoint_id'] in goldphi and not lo<=goldphi[p['endpoint_id']]<=hi:misses.append(dict(arm=arm,case=case['case_id'],step=p['step'],kind='potentialproof'))
   proof_small=small;small=ledger[arm]['small'];primary=g['costs']['RHS']-small;inverse=primary-8*g['costs']['states']
   if ledger[arm]['RHS']!=g['costs']['RHS']:raise AssertionError('closed receipt/ledger reservations differ')
   if inverse<0 or inverse>2304 or inverse%576:raise AssertionError('partial/cold inverse reconciliation')
   table.append(dict(arm=arm,status=g['status'],RHS=g['costs']['RHS'],primary_grid_RHS=primary,reduced_complex_LS_RHS=small,reduced_successful_published_proofs=proof_small,cold_inverse_RHS=inverse,
    states=g['costs']['states'],LU=g['costs']['LU'],completed_steps=completed,planned_steps=9216,accepts=accepted,completed_step_Full_fallbacks=fallbacks,all_Full_attempts=g['costs']['states'],
    whole_process_seconds=timing[arm]['seconds'],controlled_seconds=c['controlled_wall_seconds'],physical_operation_seconds=ledger[arm]['physics'],recorded_setup_including_warm_copy_seconds=setupwall,CPU_seconds=c['profile']['CPU_seconds'],RSS_peak_bytes=c['profile']['RSS_peak_bytes'],
    mean_proposal_material_RMS=float(np.mean(proposal_values)) if proposal_values else None,mean_actual_material_move_RMS=float(np.mean(move_values)) if move_values else None))
  total=sum(t['RHS'] for t in table)
  if total>300000:raise AssertionError('combined cap')
  complete=len(gates)==4 and all(g['status']=='COMPLETE' for g in gates.values())
  result=dict(status='NO_GO_NUMERICAL_MISMATCH' if misses else 'PASS_SAME_KERNEL_GRID_NUMERICS' if complete else 'PARTIAL_CAP_OR_ERROR_NOT_FULL_GRID_PASS',
   misses=misses,all_potential_proofs_checked=proof_count,costs=table,cases=case_rows,combined_complex_RHS=total,
   independent_parents=4,settings_are_not_independent_parents=True,posterior_reference='NONE',float64_certified=False,new_postclose_RHS=0)
  write(ROOT_RESULT/'COMPARISON.json',result)
  with (ROOT_RESULT/'PER_CASE.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(case_rows[0]));w.writeheader();w.writerows(case_rows)
  report=['# 前瞻较大步长四臂：实际封账附录','',result['status'],f'组合 RHS {total}/300000；数值不一致 {len(misses)}。完整/部分/未运行范围见 PER_CASE.csv。',
   '各臂独立进程及固定不可转移子预算。重复半径/种子不是独立父对象；所有尝试、拒绝、失败/未解析状态、最终场和 proof 都保留。包络失败必须 Full 或停止，不作近似判决。',
   '共享 R0 冷拥有费用每近似臂最多2304RHS，跨相同原始父锚点9配置复用，不是免费构造、重基或跨候选 basis。Full 同样拥有原在线锚点与接收器缓存。',
   '报告普通Krylov/固定精度/decision实际费用，不预设beta.002之外仍有收益；这里没有后验参考、独立样本、MCSE或定向float64证书。',
   '完整外部协调器时间（含postclose、全部状态/场数组与proof发布）由最后 COORDINATOR_COMPLETE.json 记录，不在已计费的子臂时间上再次加账。']
  report += ['|臂|状态|完成步|接受|已完成步Full回退|全Full尝试|总RHS|冷R0RHS|小LSRHS|整进程秒|设置记录秒|平均proposalRMS|平均实际moveRMS|',
   '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
  for r in table:
   if r['status']=='NOT_RUN':report.append(f"|{r['arm']}|NOT_RUN||||||||||||")
   else:report.append('|'+ '|'.join(str(r[k]) for k in ('arm','status','completed_steps','accepts','completed_step_Full_fallbacks','all_Full_attempts','RHS','cold_inverse_RHS','reduced_complex_LS_RHS','whole_process_seconds','recorded_setup_including_warm_copy_seconds','mean_proposal_material_RMS','mean_actual_material_move_RMS'))+'|')
  (ROOT_RESULT/'REPORT_ZH.md').write_text('\n\n'.join(report)+'\n')
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)

def pair(release):
 check_release(release)
 if MARKER.exists() or ROOT_RESULT.exists():raise FileExistsError('no replay')
 MARKER.mkdir(parents=True);start=time.perf_counter();env=dict(os.environ,R4_ARRAY_CODEC='stored');times={};records=[]
 for arm in ARMS:
  tick=time.perf_counter();p=subprocess.run([sys.executable,'-u',str(Path(__file__).resolve()),arm,'--release',str(release.resolve())],cwd=ROOT,env=env)
  times[arm]=dict(seconds=time.perf_counter()-tick,returncode=p.returncode);records.append(dict(arm=arm,returncode=p.returncode))
  gate=ROOT_RESULT/arm/'ARM_COMPLETE.json'
  if p.returncode and (not gate.exists() or json.loads(gate.read_text())['status']!='STOPPED'):break
 write(ROOT_RESULT/'ARM_PROCESS_WALL.json',times)
 tick=time.perf_counter();p=subprocess.run([sys.executable,'-u',str(Path(__file__).resolve()),'compare'],cwd=ROOT,env=env);times['postprocessing']=dict(seconds=time.perf_counter()-tick,returncode=p.returncode)
 payload=dict(status='CLOSED_WITH_POSTPROCESSING' if p.returncode==0 else 'POSTPROCESSING_FAILED',arm_processes=records,process_wall=times,
  coordinator_whole_wall_seconds=time.perf_counter()-start,scope='allarmstart/coldsetup/saves/sealing+postclose+reports included; childwallnotdoubledebited')
 write(ROOT_RESULT/'COORDINATOR_COMPLETE.json',payload);write(MARKER/'COMPLETE.json',payload)
 if p.returncode:raise SystemExit(p.returncode)

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','pair','compare',*ARMS));p.add_argument('--release',type=Path);a=p.parse_args()
 if a.action=='prepare':
  if PREPARED.exists():raise FileExistsError('preserve preparation')
  write(PREPARED,protocol());write(OUT/'configs/INTERVAL_CHAIN_BOUNDARY_RELEASE_TEMPLATE.json',dict(released=False,experiment=EXPERIMENT,campaign_id='P3LR4-20261025',prepared_protocol_identity=checkpoint_identity(PREPARED),
   queue_dependency_closed_receipt='<parentVIcontrollerCLOSED.json>',source=str(Path(__file__).resolve()),scope='inert; main schedules after VI controller closes'))
  print('Prepared boundary only;physics0.');return
 if a.action=='compare':postclose();return
 if a.release is None:p.error('parent release required; no autoqueue or launch')
 check_release(a.release)
 if a.action=='pair':pair(a.release)
 else:engine.run(a.action,a.release)

if __name__=='__main__':main()
