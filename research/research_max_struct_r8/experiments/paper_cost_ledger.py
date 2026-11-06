"""Aggregate paid costs without counting reused checkpoints twice.

Measured-CPU gaps are retained as None; observed process CPU snapshots are
explicit lower bounds for interrupted development attempts.
"""
import json
from pathlib import Path
from research.research_max_struct_r8.experiments.paper_benchmark import ROOT,OUT,write

def read(p):return json.loads(p.read_text())
def seconds(token):
    parts=token.split(':');return sum(float(v)*60**i for i,v in enumerate(reversed(parts)))

def run():
    rows=[]
    for p in OUT.glob('*/*/RECEIPT.json'):
        if not (p.parent/'PREDICTION_SEAL.json').exists():continue
        r=read(p);rows.append(dict(path=str(p.relative_to(ROOT)),kind='32grid_final_comparative_and_executed',CPU_seconds=r['cpu_seconds'],wall_seconds=r['wall_seconds'],LU=r['physical_cost'].get('LU',0),complex_RHS=r['physical_cost'].get('actual_complex_RHS',0)))
    for group in ['original64','strong_corrector']:
        for p in (OUT/group).glob('*/*/RECEIPT.json'):
            if not (p.parent/'PREDICTION_SEAL.json').exists():continue
            r=read(p);reuse=r.get('reused_auxiliary_cpu_seconds',0);old=r.get('auxiliary',{}).get('physical_cost',{}) if reuse else {}
            rows.append(dict(path=str(p.relative_to(ROOT)),kind=group,comparative_algorithm_CPU_seconds=r['cpu_seconds'],CPU_seconds=r['cpu_seconds']-reuse,wall_seconds=r['wall_seconds']-r.get('reused_auxiliary_wall_seconds',0),reused_auxiliary_CPU_seconds=reuse,LU=r['physical_cost'].get('LU',0)-old.get('LU',0),complex_RHS=r['physical_cost'].get('actual_complex_RHS',0)-old.get('actual_complex_RHS',0),comparative_algorithm_complex_RHS=r['physical_cost'].get('actual_complex_RHS',0),RHS_note='reusedauxRHS subtracted fromexecutedrow; interruptedattempt partialRHS still lowerbound/notrecorded'))
    for p in OUT.glob('*/DATA_RECEIPT.json'):
        r=read(p);rows.append(dict(path=str(p.relative_to(ROOT)),kind='shared150_forward_data_setup',CPU_seconds=r['cpu_seconds'],wall_seconds=r['wall_seconds'],LU=0,complex_RHS=r['forward']['physical_complex_RHS']))
    for p in (ROOT/'delegated/paper_baselines/forward64_setup_attempt').glob('*/DATA_RECEIPT.json'):
        r=read(p);rows.append(dict(path=str(p.relative_to(ROOT)),kind='development64_forward_setup',CPU_seconds=r['cpu_seconds'],wall_seconds=r['wall_seconds'],LU=0,complex_RHS=r['forward']['physical_complex_RHS']))
    for p in (OUT/'original64').glob('*/*/FAILURE.json'):
        r=read(p);rows.append(dict(path=str(p.relative_to(ROOT)),kind='failed_development_attempt',CPU_seconds=r.get('cpu_seconds'),wall_seconds=r['wall_seconds'],physical_cost=r.get('physical_cost'),CPU_status='NOT_RECORDED if None; do not infer CPU from wall'))
    for p in (OUT/'original64').glob('*/*/INTERRUPTED*.json'):
        r=read(p);tokens=r['observed_process_elapsed_CPU_RSS'].split();rows.append(dict(path=str(p.relative_to(ROOT)),kind='interrupted_engineering_attempt',CPU_seconds=seconds(tokens[1]),wall_seconds=seconds(tokens[0]),CPU_status='observed processCPU lowerbound beforeownworker SIGTERM; not exact exitCPU',RSS_KB=int(tokens[2])))
    grid=OUT/'OFFLINE_GRID_ADEQUACY.json'
    if grid.exists():
        for r in read(grid)['rows']:rows.append(dict(kind='offline_paid_grid_adequacy',side=r['side'],CPU_seconds=r['cost']['cpu_seconds'],wall_seconds=r['cost']['wall_seconds'],LU=0,complex_RHS=r['cost']['physical_complex_RHS']))
    for kind in ['NN_INITIALIZATION','NN_AUG_INITIALIZATION']:
        f=OUT/'austria_eps4_400MHz'/kind/'RECEIPT.json'
        if f.exists():
            r=read(f);rows.append(dict(kind='shared_frozen_NN_features_weights_load_inference',initialization=kind,path=str(f.relative_to(ROOT)),CPU_seconds=r['total_cpu_seconds'],wall_seconds=r['total_wall_seconds'],LU=0,complex_RHS=0,exclusions='Python/PyTorch imports and parenttraining costs'))
    write(OUT/'COST_LEDGER.json',dict(status='ACTUAL_COSTS_WITH_EXPLICIT_GAPS',rows=rows,known_CPU_seconds_sum=sum(r['CPU_seconds'] for r in rows if r.get('CPU_seconds') is not None),CPU_missing_rows=[r['path'] for r in rows if r.get('CPU_seconds') is None],warning='Do not sum comparative_algorithm_CPU when reused checkpoints already appear in interruptedattempt processCPU. Longcorrector fees incremental fromearly80. CPUlowerbounds/gaps and sharedoperatorassembly exclusions prevent claiming exactgrandtotal.',external_training_cost_receipt='RESULTS/training/physics/RECEIPT.json (parent ledger)'))

if __name__=='__main__':run()
