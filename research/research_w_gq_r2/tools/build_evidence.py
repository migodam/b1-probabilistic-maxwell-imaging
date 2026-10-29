"""Assemble saved evidence; no physics, no hidden truth-selected method changes."""
from pathlib import Path
import json,csv,time,platform,sys,numpy as np
from r2.accounting import BudgetBook,plain
ROOT=Path(__file__).resolve().parents[1]
def read(p):return json.loads((ROOT/p).read_text())
def dump(p,x):(ROOT/p).write_text(json.dumps(x,default=plain,indent=2,allow_nan=False)+'\n')
def table(p,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (ROOT/p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v,default=plain) if isinstance(v,(list,dict)) else v for k,v in r.items()})
def cost():
    b=BudgetBook();rows=b.records();caps=b.auth;ph={p:b.snapshot(p) for p in caps['caps']}
    jobs=[json.loads(p.read_text()) for p in sorted((ROOT/'results/jobs').glob('*/RECEIPT.json'))]
    allcpu=sum(j['process_cpu_seconds'] for j in jobs);counts=b.snapshot();peak=max(j.get('polled_peak_RSS_bytes',0) for j in jobs)
    checks=[r for r in rows if r.get('label')=='r2_checked_linear_action']
    result=dict(counted=counts,phases=ph,phase_caps=caps['caps'],all_caps=caps['total'],RSS_cap_bytes=caps['rss_bytes'],peak_job_RSS_bytes=peak,jobs=jobs,completed_wrapper_process_CPU_seconds=allcpu,total_wall_seconds_serial_jobs=sum(j['wall_seconds'] for j in jobs),failed_jobs=[j['job'] for j in jobs if j['exit_code']],accounting='wrapper process CPU includes child CPU and monitoring once. Reservation totals preserved; actual successful linear actions separately logged.',unmetered_gap='initial PDF extraction/source inspection, tool-service network latency, interactive source authoring and final readback are outside wrapper CPU; do not call ledger complete human/research elapsed effort.',old_bank_paid_historical=dict(full_fe=2950,rhs=11800,CPU_seconds=9288.02304,not_new_cost=True),cold_Bayesian='NOT_RUN: reliable pre-full-table bound not obtained; new cold Bayes quality/cost not claimed',publication='LOCAL_ONLY_NO_AUTOMATIC_PUSH',physics_closed=read('PHYSICS_CLOSED.json'))
    assert all(ph[p]['rhs']<=caps['caps'][p]['rhs'] and ph[p]['full_fe']<=caps['caps'][p]['full_fe'] and ph[p]['process_cpu_seconds']<=caps['caps'][p]['process_cpu_seconds'] for p in ph)
    assert peak<=caps['rss_bytes']
    dump('PHYSICS_COST.json',result);return result
def main():
    cost()
    e1=[dict(experiment='resolved_microstructure',status='NUMERICALLY_SUPPORTED',**r) for r in read('results/micro/RESULTS.json')['rows']]
    e1 += [dict(experiment='symmetry',status='NUMERICALLY_SUPPORTED_LOCAL_ONLY',**r) for r in read('results/symmetry/RESULTS.json')]
    with (ROOT/'results/imaging/TARGET_RESULTS.csv').open() as f:e1 += [dict(experiment='finite_online',status='MODEL_MISMATCH_DIAGNOSTIC',**r) for r in csv.DictReader(f)]
    first=read('results/continuous_fwi/ONLINE_OUTPUTS.json');third=read('results/continuous_fwi_third/ONLINE_OUTPUTS.json');allstarts=first+third;chosen=[]
    truth=read('results/continuous_fwi/TRUTH_OFFLINE.json')['target']
    for m in sorted(set(r['method'] for r in allstarts)):
        candidates=[dict(r,run='original') for r in first if r['method']==m]+[dict(r,run='missing_orientation_common_start') for r in third if r['method']==m]
        c=min(candidates,key=lambda r:r['best_objective']);t=c['best_target'];row=dict(method=m,selected_by='minimum saved full objective, never truth',selected_run=c['run'],selected_start=c['start'],position_error=float(np.linalg.norm(np.asarray(t['pose'])-truth['pose'])),radius_error=abs(t['radius']-truth['radius']),eps_real_error=abs(t['contrast']-truth['contrast']),loss_error=abs(t['loss']-truth['loss']),shape_correct=t['shape']==truth['shape'],objective=c['best_objective'],stop=c['stop'],predicted=t)
        chosen.append(row);e1.append(dict(experiment='continuous_combined_measurement_selected',status='OPTIMIZATION_NOT_REACHED',**row))
    dump('results/CONTINUOUS_COMBINED_SELECTED.json',chosen);table('E1_RESULTS.csv',e1)
    e2=[dict(experiment='A0_profile',**r) for r in read('results/compensation/CURVES.json')]
    e2 += [dict(experiment='reoptimized_nested_profile',**r) for r in read('results/compensation_extended/RESULTS.json')['curves']]
    e2.append(dict(experiment='A2_final_profile',**read('results/compensation_final/RESULTS.json')))
    e2.append(dict(experiment='A2_third_refinement',**read('results/final_pair_refinement/RESULTS.json')))
    with (ROOT/'results/e2_micro_reference/CASE_RESULTS.csv').open() as f:e2 += [dict(experiment='offline_physical_finite_reference',**r) for r in csv.DictReader(f)]
    table('E2_RESULTS.csv',e2)
    # Explicitly save both data-only Schur matrices; the original JSON omits Kw.
    d=np.load(ROOT/'results/compensation/BASE.npz');n=d['eps'].size;ws=np.r_[np.full(n,.5),np.full(n,.05)];ts=np.array([.2,.2,2.,.08,.04]);data={};summ=[]
    for name,A,B in [('A0',d['A'][1].reshape(4,16,-1)[:2].reshape(32,-1),d['B'][1].reshape(4,16,-1)[:2].reshape(32,-1)),('A1',d['A'][1],d['B'][1]),('A2',np.concatenate(d['A']),np.concatenate(d['B']))]:
        A=np.vstack((A.real,A.imag))/.01*ws;B=np.vstack((B.real,B.imag))/.01*ts
        def projector(X):
            u,s,v=np.linalg.svd(X,full_matrices=False);r=s>np.finfo(float).eps*max(X.shape)*s[0];return u[:,r],s
        UA,sa=projector(A);UB,sb=projector(B);RA=A-UB@(UB.T@A);RB=B-UA@(UA.T@B)
        kw=RA.T@RA;kx=RB.T@RB;data[name+'_K_W_given_X_data']=kw;data[name+'_K_X_given_W_data']=kx
        data[name+'_A']=A;data[name+'_B']=B
        summ.append(dict(contract=name,target_eigenvalues=np.linalg.eigvalsh(kx),wall_nonzero_gram_eigenvalues=np.linalg.eigvalsh(RA@RA.T),machine_rank_scope='floating matrix only; all singularvalues retained in original INFORMATION; no uniform grid/derivative bound'))
    np.savez_compressed(ROOT/'results/BIDIRECTIONAL_SCHUR_SUPPLEMENT.npz',**data);dump('results/BIDIRECTIONAL_SCHUR_SUPPLEMENT.json',summ)
    e3=[]
    for p,run in [('results/e3_return_operator/PER_RANK_ERRORS.json','original'),('results/e3_return_transpose_repair/PER_RANK_ERRORS.json','transpose_repair')]:
        for r in read(p):e3.append(dict(experiment='return_representation',run=run,wall=r['wall_id'],target=r['target_id'],method=r['method'],rank=r['rank'],status=r['status'],**r.get('metrics',{})))
    with (ROOT/'results/e3_logic/offline_scores/TARGET_RESULTS.csv').open() as f:e3 += [dict(experiment='oldbank_consensus_replay',**r) for r in csv.DictReader(f)]
    e3.append(dict(experiment='cold_Bayesian',status='NOT_RUN',reason='No affordable reliable consensus/contrast tail bound before all queries'))
    table('E3_RESULTS.csv',e3)
    with (ROOT/'QUERY_LEDGER.jsonl').open('w') as f:
        for p in sorted((ROOT/'results/e3_logic/online').glob('*/QUERY_TRACE.json')):
            for r in json.loads(p.read_text()):f.write(json.dumps(dict(case=p.parent.name,identity='REPLAY_NOT_NEW_PDE_SPEED_EVIDENCE',**r))+'\n')
    ledger=list(csv.DictReader((ROOT/'research/delegated/evidence/EVIDENCE_REUSE_MAP.csv').open()))
    for phase,path,claim in [('E1','results/micro/RESULTS.json','finite resolved microstructure equivalence'),('E1','results/imaging/RESULTS.json','unknown macro finite imaging failure'),('E2','results/final_pair_refinement/RESULTS.json','new finite A2 compensated witness'),('E2','results/e2_micro_reference/SUMMARY.json','target/wall risk separation finite class'),('E3','results/e3_return_transpose_repair/SUMMARY.json','exact self + return representation'),('E3','results/e3_logic/offline_scores/SUMMARY.json','conditional consensus does not certify stopping')]:ledger.append(dict(phase=phase,path=path,claim=claim,evidence='NEW_R2_SAVED_EXECUTION'))
    table('EVIDENCE_LEDGER.csv',ledger)
    gates=dict(overall='PHYSICAL_MECHANISMS_SUPPORTED_ONLINE_AND_NET_SPEEDUP_NOT_ESTABLISHED',A_wrong_wall_target_task='NUMERICALLY_SUPPORTED_FINITE_TASK_CLASS',B_unknown_wall_online='NOT_ESTABLISHED; finite macro model mismatch and capped continuous optimization',C_nonlinear_compensation='NUMERICALLY_SUPPORTED_DISCRETE_FINITE_WITNESS, not continuum exact equality',D_target_vs_wall='NUMERICALLY_SUPPORTED_4x12_PROPER_FINITE_PRIOR',E_return_only='NUMERICALLY_SUPPORTED_SAME_DISCRETE_TARGETS; rank8 maxwhite .000204769',F_target_consensus='NOT_ESTABLISHED_AS_UNIFORM_ONLINE_BOUND; retain macro masses and within-group tails',G_speed='NET_SPEEDUP_NOT_ESTABLISHED',H_attribution='homogenization, symmetry, IFT, Schur/T-matrix and complex symmetric projection are standard tools; new finite mechanism evidence only',I_limits=['2D nonmagnetic scalar TM','static calibrated instrument','finite prior and finite task family','continuous target FWI not converged','no uniform continuous PDE guarantee','no broad unknown support success'],J_paper='Develop constrained finite target resolution through unresolved materials; no independent faster Bayesian algorithm claim',physics=read('PHYSICS_CLOSED.json'),gates_unchanged=read('CONTRACT.json')['gates'],new_locked_evaluation=False,new_NN=False,new_3D=False,automatic_publication=False)
    dump('GO_NO_GO.json',gates)
    print(json.dumps(dict(E1_rows=len(e1),E2_rows=len(e2),E3_rows=len(e3),chosen_continuous=chosen),default=plain))
if __name__=='__main__':main()
