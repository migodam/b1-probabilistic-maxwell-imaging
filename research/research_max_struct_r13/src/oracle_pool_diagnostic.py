"""OFFLINE_DIAGNOSTIC: stream scored pools to separate search and selection.

GT-selected candidate indices are never deployed or called algorithm results.
No truth NPZ/solver is opened; inputs must already have authorized scoring.
Candidate counts differ and all corresponding generation/audit costs are paid.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import argparse,csv,json,time
from pathlib import Path
from seal_results import read,validate_seal


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--campaign',action='append',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    cpu,wall=time.process_time(),time.perf_counter();root=a.root.resolve();sources=[]
    for name in a.campaign:
        folder=root/'analysis'/name;receipt=read(folder/'SCORING_RECEIPT.json');validate_seal(read(receipt['authorized_seal']),True);sources.append((name,folder))
    from score_results import table
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);output=[];processed=0
    yes=lambda value:str(value).lower()=='true'
    for campaign,folder in sources:
        primary={r['task_id']:r for r in read(folder/'SCORES.json') if r['view']=='primary'};pools={tid:dict(legal_paid_count=0,noise_feasible_paid_count=0,best=None,best_feasible=None) for tid in primary}
        with (folder/'ALL_PAID_AND_OVERRUN_CANDIDATES.csv').open(newline='') as f:
            for row in csv.DictReader(f):
                processed+=1
                if not(yes(row.get('eligible')) and yes(row.get('legal'))):continue
                tid=row['task_id'];pool=pools[tid];pool['legal_paid_count']+=1
                feasible=yes(row.get('noise_feasible'))
                if feasible:pool['noise_feasible_paid_count']+=1
                if not row.get('NRMSE'):continue
                item={k:row.get(k) for k in ('candidate','phase','NRMSE','RMSE','IoU','Phi','D','R','CPU_seconds','residual_over_delta')}
                item['NRMSE']=float(item['NRMSE'])
                if pool['best'] is None or item['NRMSE']<pool['best']['NRMSE']:pool['best']=item
                if feasible:
                    if pool['best_feasible'] is None or item['NRMSE']<pool['best_feasible']['NRMSE']:pool['best_feasible']=item
        for tid,selected in primary.items():
            pool=pools[tid]
            for label,key in [('all_legal_paid','best'),('legal_paid_discrepancy_feasible','best_feasible')]:
                best=pool[key];nrmse=selected.get('NRMSE')
                output.append(dict(campaign=campaign,task_id=tid,case_id=selected['case_id'],method=selected['method'],physics_kind=selected['physics_kind'],cohort=selected['cohort'],view='OFFLINE_GT_ORACLE_NOT_DEPLOYABLE',pool=label,available=best is not None,legal_paid_candidates=pool['legal_paid_count'],discrepancy_feasible_paid_candidates=pool['noise_feasible_paid_count'],primary_NRMSE=nrmse,oracle_NRMSE=best['NRMSE'] if best else None,primary_minus_oracle_NRMSE=nrmse-best['NRMSE'] if best and nrmse is not None else None,oracle_candidate=best,interpretation='Diagnostic search coverage versus frozen physical selection only; GT selects this row. Not a deployable success, not an equal-candidate-count comparison, and not evidence of free search.'))
    table(out,'ORACLE_POOL_DIAGNOSTIC',output)
    receipt=dict(status='OFFLINE_DIAGNOSTIC_NOT_DEPLOYABLE',campaigns=a.campaign,candidate_rows_streamed=processed,diagnostic_rows=len(output),new_forward=0,new_optimization=0,raw_truth_NPZ_read=False,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
