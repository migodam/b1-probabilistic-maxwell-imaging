"""Produce report tables from completed experiment receipts; no solver replay."""
from pathlib import Path
from collections import defaultdict
import csv,json,time,shutil
import numpy as np
from .budget import ROOT,atomic_json

def read(p):return json.loads(Path(p).read_text())
def table(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,keys);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()})

def main():
    ledger=[json.loads(z) for z in (ROOT/'COST_LEDGER.jsonl').read_text().splitlines()]
    by=defaultdict(lambda:dict(LU=0,RHS=0,RHS576=0,kernel_wall_seconds=0.,reservations=0))
    for row in ledger:
        if row['event']=='RESERVE':
            v=by[row['group']]
            for k in ('LU','RHS','RHS576'):v[k]+=row['request'][k]
            v['reservations']+=1
        elif row['event']=='COMPLETE':by[row['group']]['kernel_wall_seconds']+=row['wall_seconds']
    stage=defaultdict(lambda:dict(LU=0,RHS=0,RHS576=0,kernel_wall_seconds=0.))
    for group,v in by.items():
        k='P1' if group.startswith('P1/') else 'N2' if group.startswith('N2/') else 'N1' if group.startswith('N/') else 'small_controls'
        for field in stage[k]:stage[k][field]+=v[field]
    state=read(ROOT/'GLOBAL_BUDGET_STATE.json')
    summary=dict(actual=state['used'],limits=state['limits'],stopped=state['stop'],stage_totals=dict(stage),groups=dict(by),source='onefinalaggregationoforiginalRESERVE/COMPLETErecords; no new physical/numerical verification',historical_host_setup_RHS576=264,historical_R9_original_anchor_RHS=6136,historical_R9_original_anchor_LU=175,history_recharged_as_new=False,physics_end=read(ROOT/'audit/PHYSICS_STOP.json')['utc'],new_GPU_MPS_training=0,cloud_purchases=0,scope='physical LU and complex RHS of all orders; real normal Cholesky/SVD/QR are separate CPU operations, not additional Maxwell RHS')
    atomic_json(ROOT/'RESOURCE_SUMMARY.json',summary)
    atomic_json(ROOT/'GLOBAL_BUDGET_AUDIT.json',dict(status='WITHIN_CAPS',final_actual=state['used'],limits=state['limits'],atomic_prepayment=True,failed_attempts_retained=True,original_R11_failure_preserved=True,source='finalexperimentcostaggregation; earliercontrolreceipt in audit/ATOMIC_BUDGET_CONTROL.json',no_new_validation_campaign=True))
    table(ROOT/'COST_BY_PHYSICAL_GROUP.csv',[dict(group=k,**v) for k,v in by.items()])
    p1=read(ROOT/'agent_c/p1_score/SUMMARY.json')['datasets']['P1']
    methods=[]
    for arm,v in p1.items():
        records=[read(p) for p in (ROOT/'results/P1').glob(f'*/{arm}/SEALED.json')]
        methods.append(dict(dataset='P1_DEV4',method=arm,parents=4,GT_coverage_eps075=v['coverage_eps075'],mean_distinct_qualified_including_anchor=v['qualified_unique_mean'],mean_oracle_MSE=v['nearest_MSE_mean'],mean_online_selected_MSE=v['selected_MSE_mean'],new_method_RHS=sum(r['new_RHS'] for r in records),new_method_LU=sum(r['new_LU'] for r in records),mean_standalone_RHS=float(np.mean([r['standalone_RHS'] for r in records])),mean_wall_seconds=float(np.mean([r['wall_seconds'] for r in records])),qualification='EXPOSED_DEV; matchedcapsnotidenticalactualcost'))
    n2=read(ROOT/'agent_c_max/n2_score/SUMMARY.json')['datasets']['N']
    nr=list(csv.DictReader((ROOT/'figures/n2_online_and_clouds/N2_SELECTED_STRUCTURE.csv').open()))
    for arm,v in n2.items():
        records=[read(p) for p in (ROOT/'results/N2').glob(f'*/{arm}/SEALED.json')];rows=[r for r in nr if r['method']==arm]
        methods.append(dict(dataset='N2_DESIGNED9',method=arm,parents=9,GT_coverage_eps075=v['coverage_eps075'],mean_distinct_qualified_including_anchor=v['qualified_unique_mean'],mean_oracle_MSE=None,mean_online_selected_MSE=v['selected_MSE_mean'],mean_occupancy_IoU=float(np.mean([float(r['occupancy_IoU']) for r in rows])),mean_material_IoU=float(np.mean([float(r['material_IoU']) for r in rows])),selected_noise_pass=sum(r['selected_noise_pass'].lower()=='true' for r in rows),new_method_RHS=sum(r['new_RHS'] for r in records),new_method_LU=sum(r['new_LU'] for r in records),mean_standalone_RHS=float(np.mean([r['new_RHS']+40 for r in records])),mean_wall_seconds=float(np.mean([r['wall_seconds'] for r in records])),qualification='DESIGNED_NATIVE258_NOT_P1; nofreshpopulationclaim'))
    table(ROOT/'METHOD_COMPARISON.csv',methods)
    shutil.copy2(ROOT/'figures/n2_online_and_clouds/N2_SELECTED_STRUCTURE.csv',ROOT/'DIGIT_PER_NUMBER_RESULTS.csv')
    groups=read(ROOT/'agent_c/p1_score/GROUP_SCORES.json')
    table(ROOT/'MATERIAL_COVERAGE_COMPARISON.csv',groups)
    table(ROOT/'DEPTH_ABLATION.csv',read(ROOT/'agent_c/p1_score/LAYER_SCORES.json'))
    table(ROOT/'SUBSPACE_EVOLUTION.csv',read(ROOT/'agent_c/p1_diagnostics/J_PROJECTOR_EVOLUTION.json'))
    table(ROOT/'NONLINEAR_SECANT_RESULTS.csv',read(ROOT/'agent_c/p1_diagnostics/FINITE_SECANTS.json'))
    shutil.copy2(ROOT/'figures/n2_online_and_clouds/DIGITS_123456789_RECONSTRUCTION.png',ROOT/'DIGITS_123456789_RECONSTRUCTION.png')
    result=dict(id='P3L_R12_BF',scope='4exposedP1DEV+9frozendesignedNdigits, notindependentpopulationconfirmation',P1=p1,N2=n2,resource=state['used'],GO_HIERARCHY='NO_GO_FOR_TESTED_MECHANISM',GO_BILATERAL_NEW_INFORMATION='NO:exactFullJfactorization',GO_FIXEDPOINT_FIBER='NOT_ESTABLISHED_FOR_P1',GO_DIGIT_TOPOLOGY='NO_TWO_HOLE_DIGIT8_IN_42_STORED_CANDIDATES',GO_NEW_COVERAGE='NO_GO',GO_FRESH_CONFIRMATION='NOT_RUN_GATE_FAILED',GO_TAP='HOLD',FULL66_POSTERIOR='UNQUALIFIED_NOT_TARGET_BLOCKER',H3='NOT_RUN_NO_GATE',shared_material_TSOM='NOT_RUN; R11GalerkinfailuredoesnotrefuteTSOM',probability_accuracy='notrequiredforfinitequalifiedsets; nofrequencyaspopstposterior',formula_atlas_count=167)
    atomic_json(ROOT/'RESULTS_SUMMARY.json',result)
    print(json.dumps(dict(resource=state['used'],stage=stage,methods=len(methods)),indent=2))
if __name__=='__main__':main()
