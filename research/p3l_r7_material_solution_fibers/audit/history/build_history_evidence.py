from pathlib import Path
import subprocess,json,csv
root=Path.cwd(); out=root/'research/p3l_r7_material_solution_fibers/audit/history'
stages=[
('F','b1-research-f','research_f_intrinsic_transport'),('F3','b1-research-f3','research_f3_quotient_exploration'),('P1','b1-p1-material-target','p1_material_target_20261018'),('P2','b1-p2-material-fiber','p2_material_fiber_20261019'),('P3L','b1-p3l-local-posterior','p3l_local_posterior_imaging'),('R2','b1-p3l-r2-risk-posterior','p3l_r2_risk_posterior_exploration'),('R4','b1-p3l-r4-conditional-transport','p3l_r4_conditional_transport'),('R5','b1-p3l-r5-parallel-ab','p3l_r5_parallel_ab'),('R5B','b1-p3l-r5-b-confirm','p3l_r5_b_independent_confirmation'),('R6','b1-p3l-r6-task-optimal','p3l_r6_task_optimal_maxwell')]
def git(*args,cwd=root):return subprocess.check_output(['git',*args],cwd=cwd)
ident=[]
for stage,wt,d in stages:
 w=root.parent/wt; head=git('rev-parse','HEAD',cwd=w).decode().strip(); f=f'research/{d}/RESEARCH_REPORT_ZH.md'
 row={'stage':stage,'worktree':str(w),'head':head,'branch':git('branch','--show-current',cwd=w).decode().strip(),'report':f,'report_blob':git('rev-parse',f'{head}:{f}').decode().strip(),'current_bytes_equal_fixed_head':git('show',f'{head}:{f}')==(root/f).read_bytes(),'historical_worktree_bytes_equal_fixed_head':git('show',f'{head}:{f}')==(w/f).read_bytes(),'report_last_commit':git('log','-1','--format=%H',head,'--',f).decode().strip()};ident.append(row)
(out/'FIXED_IDENTITIES.json').write_text(json.dumps({'current_head':git('rev-parse','HEAD').decode().strip(),'R7T_scientific_commit':git('rev-parse','d1bdc8aceae4cd6a843f873b61b019d5eee3cd1d').decode().strip(),'stages':ident},ensure_ascii=False,indent=2)+'\n')
# All statements are transcriptions/limited summaries; no raw GT or physics arrays read.
rows=[]
def add(i,s,file,a,b,status,fact,limit):
 c=next(x['head'] for x in ident if x['stage']==s)
 p=f'research/{next(x[2] for x in stages if x[0]==s)}/{file}'
 fixed=git('show',f'{c}:{p}').decode(); lines=fixed.splitlines()
 assert 1<=a<=b<=len(lines),(p,a,b,len(lines))
 assert (root/p).read_bytes()==git('show',f'{c}:{p}'),p
 rows.append(dict(id=i,stage=s,commit=c,file=p,line_start=a,line_end=b,status=status,fact=fact,limit=limit,excerpt='\n'.join(f'{n}: {lines[n-1]}' for n in range(a,b+1))))
add('F01','F','RESEARCH_REPORT_ZH.md',43,54,'QUALIFIED_GEOMETRY_MISSING_VECTOR','weak/ring material second-order normal fractions .0032863%/14.77154%; old Maxwell vector missing','material fractions are not response fractions; full vector attribution unavailable')
add('F02','F','RESEARCH_REPORT_ZH.md',80,85,'FAIL','ring/weak stationary gate failed despite residuals inside noise ball','do not retroactively relax stationary threshold')
add('F03','F','RESEARCH_REPORT_ZH.md',120,137,'NOT_RUN_NO_GO','F2 four arms not run; learned Explorer no-go; exact fiber unconfirmed','not a proof no finite uncertainty exists')
add('F301','F3','RESEARCH_REPORT_ZH.md',21,39,'VERIFIED_IDENTITY','weight-scale quotient preserves material image; induced rho prior nonuniform','labels/collisions remain; not global smooth one-to-one atlas')
add('F302','F3','RESEARCH_REPORT_ZH.md',43,54,'PASS_NEW_GATE','quotient anchors pass newly frozen projected/KKT criteria','old eight-dimensional F Gate2 FAIL remains; likelihood anchors not posterior MAP')
add('F303','F3','RESEARCH_REPORT_ZH.md',64,82,'CANDIDATE_ONLY','weak native C and D task-scale endpoints competitive; ring not','three nearby numerical endpoints are not posterior modes')
add('F304','F3','RESEARCH_REPORT_ZH.md',102,120,'NO_GO_HIGHER_ORDER_INCREMENT','eight-step MH unqualified; same-RHS native found 2 competitive weak candidates versus D 1','finite endpoint discovery not posterior accuracy')
add('P101','P1','RESEARCH_REPORT_ZH.md',7,16,'MATERIAL_READY_PHYSICS_NOT_RUN','material stage generated 192 materials; real physics initially NOT_RUN','must distinguish later separately authorized physics phase')
add('P102','P1','PHYSICS_QUALIFICATION_REPORT_ZH.md',7,14,'A_PASS_B_FAIL','A residual22.2923 inside23.7854; B309.7873 fails; one feasible anchor','A not certified MAP; no alternative worlds or posterior probabilities')
add('P103','P1','PHYSICS_QUALIFICATION_REPORT_ZH.md',32,50,'LIMITED_QUALIFICATION','both endpoint Dg ranks66; 26 states/976 RHS; mechanism comparisons NOT_RUN','B failure is optimization at frozen cap; not absence of legal material')
add('P201','P2','RESEARCH_REPORT_ZH.md',3,15,'REDIRECT','15/24 feasible anchors; 69 task-scale alternatives overall; heldout29; best-of4 oracle MAE gain1.35%','mixed archive/budgets; no matched-method victory or posterior mass')
add('P202','P2','RESEARCH_REPORT_ZH.md',80,88,'NO_GO_CURRENT_LEARNING','global methods no stable extra branches; short MH mixing insufficient; Gaussian/Flow no best-of4 gain','not proof of no remote branch or learnability impossibility')
add('P203','P2','RESEARCH_REPORT_ZH.md',96,117,'NO_GO_A17','A17 fallback17/48 and setup-amortized cost ~11.1x optimized Full; spectral consistency standard','E evaluator uses offline oracle error checks; not online controller speed proof')
add('P204','P2','RESEARCH_REPORT_ZH.md',119,139,'ISOLATION_DEVIATION_CLOSED','dataset and inversion shared OS process; 2236 charged states/46728RHS','code-path separation does not satisfy planned process isolation')
add('P3L01','P3L','RESEARCH_REPORT_ZH.md',3,7,'REDIRECT','128-parent approximate mean MSE .0615021 versus point .0110309 (+457.5%)','does not refute correct Bayes mean optimality')
add('P3L02','P3L','RESEARCH_REPORT_ZH.md',31,47,'UNQUALIFIED','all11 Full references unqualified;44chains155648 production states;SMC128 one ancestor','retained states not ESS; stable Phi not material mixing')
add('P3L03','P3L','RESEARCH_REPORT_ZH.md',93,105,'NO_GO_COST','decoder-aware surrogate reduces response error but not material MCSE, wall increases','standard exact-Hamiltonian correction not new theorem; bracket not exact fiber proof')
add('R201','R2','RESEARCH_REPORT_ZH.md',5,7,'SUPERVISED_RISK_GO','fresh256 CNN MSE .007571318 versus .010131196,25.27% reduction, adjusted interval22.34–27.99%','39 worse objects retained; learned maps not posterior samples or all noise-compatible')
add('R202','R2','RESEARCH_REPORT_ZH.md',58,68,'PARTIAL_MECHANISM','sensitivity vs random decoder envelope .70% improvement;backprojection independent increment unresolved','static frame/channel contribution not dynamic Maxwell theorem')
add('R203','R2','RESEARCH_REPORT_ZH.md',78,94,'FINITE_NOISE_CANDIDATES_FULL_UNQUALIFIED','1934/2048 legal selected endpoints are noise/likelihood competitive at .02/.04 steps; trajectory risk fails','counts endpoints not parents/modes; weak projector noninvolutivity not exact-kernel nonintegrability')
add('R204','R2','CAMPAIGN_CLOSED.json',1,12,'CLOSED','5956328 RHS322717 attempts/LU;13906.677 controlled CPUwall','final closure differs from report per-stage wall snapshot; old balances not reused')
add('R401','R4','RESEARCH_REPORT_ZH.md',33,46,'RISK_GO','512+256 risk gains vs R2 13.10%/11.36%;vs strongerCNN5.36%/4.90%','Spline/Gaussian fusion difference intervals cross0')
add('R402','R4','RESEARCH_REPORT_ZH.md',66,91,'RELATIVE_DENSITY_GO_ABSOLUTE_FAIL','Full66Spline NLL77.50922/77.27755 versus83.07743/83.44856;CFMenergy~21% gain','absolute response calibration severe failure;decoder-only CFM missing; density tolerance FAIL;8 solves NOT_RUN')
add('R403','R4','RESEARCH_REPORT_ZH.md',105,124,'LIMITED_COMPUTE_GO_PATH_NO_GO','8192-step same-path replay71.397→24.153seconds,zero decision difference;111 candidate archive','short beta=.002 kernel;no mixing gain;foursegment2 final supported versus94 direct')
add('R501','R5','RESEARCH_REPORT_ZH.md',3,10,'A_REDIRECT_B_CONDITIONAL_GO','A regular-region numerical moments improve;B VAL retention92.96%;fresh256+256 NOT_RUN','development evidence cannot substitute fresh independent confirmation')
add('R502','R5','RESEARCH_REPORT_ZH.md',64,81,'LIMITED_A','local response moments accurate; QR-RTO regular success but host0/8','coverage, branch normalization,TV and20%equalaccuracycost not established;Full66NOT_RUN')
add('R503','R5','RESEARCH_REPORT_ZH.md',118,142,'MODULE_COMPRESSION_CPU_FAIL','provider width1288.45us versusK32332.26us;pipeline349.73versus640.91us','shared anchor cost remains;CPU/MPS3e-6FAIL retained;new confirmation not yet run')
add('R5B01','R5B','RESEARCH_REPORT_ZH.md',4,12,'INDEPENDENT_RISK_GO','new512 frozen objects;student3.11%lowerMSE vs StrongCNN;Q64 .006854734','retention80.70% point is not >=80% qualification')
add('R5B02','R5B','RESEARCH_REPORT_ZH.md',34,46,'HOLD','retention adjusted interval60.28–113.12%;80%contrast lowerbound -5.43036e-5','no ensemble rescue;standard squareloss decomposition explains94.85%seed-gap but no new theorem')
add('R5B03','R5B','RESEARCH_REPORT_ZH.md',50,64,'CLOSED_LIMITED_SCOPE','1999504RHS56734states/LU;shared anchor dominates system;fullposterior not qualified','143596deployments only teacherCPU recovery;not total economic break-even')
add('R601','R6','RESEARCH_REPORT_ZH.md',23,51,'A_REDIRECT','stratified localQ1/Q2 moments and copula-controlled Maxwell counterexample','not Full66;no global tail/mode bound;no taskaccuracy cost separation;4D not expanded')
add('R602','R6','RESEARCH_REPORT_ZH.md',55,69,'B_NO_GO_CURRENT_OBJECTIVES','fourN1024objectivefamilies do not beat historicalN32768deployment;moment regression not trustworthy covariance','does not show all probability learning impossible')
add('R603','R6','RESEARCH_REPORT_ZH.md',75,92,'ORIGINAL_FROZEN_PASS','512newparents in2noise cohorts;lowcost24+64 q .005335107/.009555500 vs full .005463559/.009534329;wall34.551%/32.571%less','sameendpoint q independent risk gain low-noise only;two noise cohorts not pooled')
add('R604','R6','RESEARCH_REPORT_ZH.md',98,104,'STRENGTHENED_CONTROL_HOLD','Full-q low-noise.005315940;candidate1%NIuppercontrast+1.26745e-5;newwindow20%savinglowerbound16.82%','original PASS preserved;strongerNI and across-window20% not established;Direct competitive')
add('R605','R6','RESEARCH_REPORT_ZH.md',110,130,'CLOSED_REDIRECT','6860984RHS210355states/LU;33trainingstages;GO-TAP REDIRECT','noFull66qualification/calibration/layout/largermeshguarantee;R5B80%HOLD unchanged')
with (out/'HISTORICAL_EVIDENCE_TABLE.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
(out/'HISTORICAL_EVIDENCE_TABLE.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n')
# Read machine summaries only, never truth files or raw experiment arrays.
snap={}
for stage,wt,d in stages:
 for name in ['PHYSICS_DECISION.json','GATE_DECISION.json','GATE_SUMMARY.json','CAMPAIGN_CLOSED.json','RUN_CLOSED.json','RUN_FINAL_CLOSED.json','RESOURCE_SUMMARY.json','RESOURCE_SUMMARY_FINAL.json','RESULTS_SUMMARY.json']:
  p=root/'research'/d/name
  if p.exists():snap[str(p.relative_to(root))]=json.loads(p.read_text())
(out/'MACHINE_SUMMARY_READBACK.json').write_text(json.dumps(snap,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'identities':len(ident),'evidence_rows':len(rows),'machine_summaries':len(snap),'all_report_bytes_match':all(x['current_bytes_equal_fixed_head'] and x['historical_worktree_bytes_equal_fixed_head'] for x in ident)},indent=2))
