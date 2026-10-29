from core import *
from invert_digits import Problem
import collections,csv,datetime

def main():
 logs=[P/'data/COST_LEDGER.jsonl',P/'agent_a/COST_LEDGER.jsonl',P/'agent_b/KERNEL_COST_LEDGER.jsonl',P/'agent_b/SUPPLEMENT_COST_LEDGER.jsonl',P/'agent_c/path_results/COST_LEDGER.jsonl',P/'agent_c/API_COST_LEDGER.jsonl']+sorted((P/'results').glob('*/*/*/COST_LEDGER.jsonl'))
 target=P/'COST_LEDGER.jsonl';tot=collections.Counter();by=collections.Counter();luorder=collections.Counter();seq=0;groups={}
 with target.open('w') as out:
  for path in logs:
   group=path.relative_to(P).parts[0];gg=groups.setdefault(group,dict(RHS=0,LU=0,physical_charge_records=0))
   for li,line in enumerate(path.read_text().splitlines(),1):
    a=json.loads(line);rhs=int(a.get('RHS',0));lu=int(a.get('LU',0))
    if rhs==0 and lu==0:continue
    d=a.get('RHS_by_order',{})
    if not d:
     order=a.get('order',a.get('factor_order',None));d={str(order):rhs}
    if len(d)!=1:raise RuntimeError(('ambiguousorder',path,li,d))
    k=next(iter(d));by[k]+=rhs;luorder[k]+=lu;tot['RHS']+=rhs;tot['LU']+=lu;gg['RHS']+=rhs;gg['LU']+=lu;gg['physical_charge_records']+=1;seq+=1
    out.write(json.dumps(dict(sequence=seq,utc=a.get('utc'),raw_source=str(path.relative_to(P)),raw_line=li,charge_RHS=rhs,charge_LU=lu,order=k,event=a.get('event','DIRECT_EVALUATION'),context=a.get('context'),label=a.get('label'),scope='oneactualpaidattempt;summarysnapshotsnotcharges'))+'\n')
 costs=dict(status='FAIL_CAP_EXCEEDED',RHS_conservative_all_orders=tot['RHS'],LU_all_orders=tot['LU'],RHS_by_order=dict(by),LU_by_order=dict(luorder),group_totals=groups,charged_attempt_records=seq,ledger_sources=[str(p.relative_to(P)) for p in logs],new_RHS576=by.get('576',0),new_LU576=luorder.get('576',0),caps=dict(RHS=400000,LU=12000,RHS576=1024),RHS_cap_pass=tot['RHS']<=400000,LU_cap_pass=False,LU_excess=tot['LU']-12000,history_host_setup=dict(RHS576=264,LU576=1,new_this_campaign=False),notes=['Do not reinterpretallorderLUasLU256.','ConservativecrossorderRHSnotFLOP-equivalent.','Originalperjobledgersappendonly;rootnormalizationgeneratedfromexplicitchargeeventsonly.','No physicsafterSTOP,alladditionalanalysiszeroRHS.','RootglobalLUcapnotenforcedinmultiplejobdrivers;branchcapsaloneinsufficient.'])
 write(P/'COST_SUMMARY.json',costs)
 s=load_setup();phase=[]
 for ds in ['N','P']:
  problem=Problem(ds,s)
  for i in range(9):
   oid=f'{ds}_{i:03d}';g=dict(np.load(P/'data/truth'/f'{oid}.npz'));z=dict(np.load(P/'results'/ds/'FULL256'/oid/'PREDICTION.npz'));meta=json.loads((P/'results'/ds/'FULL256'/oid/'PREDICTION_SEALED.json').read_text());y=np.load(P/'data/online'/f'{oid}.npz')['y']
   xg=np.r_[g['occupancy'][s['ROI']],0.,0.] if ds=='N' else g['xgt']
   gtphi=phi(g['Fmean'],y);gtr=problem.reg(xg)[0];rr=problem.reg(z['x'])[0]
   history=json.loads((P/'results'/ds/'FULL256'/oid/'HISTORY.json').read_text())
   phase.append(dict(dataset=ds,digit=i+1,GT_Phi=gtphi,GT_regularizer=gtr,GT_total_objective=gtphi+gtr,reconstruction_Phi=meta['Phi'],reconstruction_regularizer=rr,reconstruction_total_objective=meta['Phi']+rr,objective_at_reconstruction_below_GT=bool(meta['Phi']+rr<gtphi+gtr),status=meta['status'],accepted_steps=sum(h['accepted'] for h in history),last_relative_objective_change=(history[-1]['objective']-history[-1]['new_objective'])/max(1.,abs(history[-1]['objective'])),phase='offlineGTdiagnosticonly,notonlinechoice'))
 with (P/'audit/PRIOR_AND_CONVERGENCE_DIAGNOSTICS.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(phase[0]));w.writeheader();w.writerows(phase)
 write(P/'audit/FINAL_RESOURCE_RECONCILIATION.json',dict(costs=costs,methods_sealed=72,GT_scoring_after_seal=True,physics_stopped=True,original_budget_not_rewritten=True))
 print(json.dumps({k:v for k,v in costs.items() if k not in ['RHS_by_order','LU_by_order','ledger_sources']},indent=2))
if __name__=='__main__':main()
