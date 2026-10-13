"""Measured zero-physics CDF probes and actual public observation replay."""
from pathlib import Path
import json,time
import numpy as np
from src.b1.e0 import dump,read_observation,ledger_snapshot
from src.b1.continuous import FrozenDensity

def run(root):
    root=Path(root);stage=root/'research/b1_next/e0';out=stage/'results'
    summary=json.loads((out/'SUMMARY.json').read_text());stats={};start=time.perf_counter();cpu=time.process_time()
    before=ledger_snapshot(root);polynomial=np.load(root/'results/refinement/FORWARD_POLYNOMIAL.npz');rows=[]
    for name,case in summary['cases'].items():
        density=FrozenDensity(polynomial['coeff7'],polynomial['bounds'],read_observation(root,name),stats,'poly7')
        for source in ['legacy_original','legacy_refined']:
            for h,old in enumerate(case[source]['conditional_material']):
                norm=np.exp(case['result']['log_evidence'][h]-density.shifts[h])
                for axis in range(2):
                    for probability,cut in zip([.05,.5,.95],old['quantiles_05_50_95'][axis]):
                        z=(cut-density.lo[axis])/density.width[axis];box=np.array([[0.,1.],[0.,1.]])
                        box[axis,1]=z
                        value=float(density.rectangle(h,box,97,'direct_historical_CDF_probe')[0]/norm)
                        rows.append(dict(case=name,source=source,H=h,coordinate=axis,saved_claimed_quantile=probability,
                                         cutoff=cut,continuous_target_CDF=value,difference=value-probability,
                                         method='direct truncated 97-point tensor integral, not grid interpolation',
                                         target='frozen polynomial7 continuous posterior',certified=False))
    probe_wall=time.perf_counter()-start;probe_cpu=time.process_time()-cpu
    dump(out/'CONTINUOUS_CDF_PROBES.json',dict(rows=rows,new_physics=0,polynomial_counts=stats,wall_seconds=probe_wall,cpu_seconds=probe_cpu))
    t=time.perf_counter();c=time.process_time()
    from src.b1.cli_infer import run as imaging
    config=json.loads((root/'configs/protocol.json').read_text());cli=out/'cli_gain_replay'
    imaging(root,config,root/'results/cases/gain_marginalized/posterior.npz',cli,'continuous')
    receipt=json.loads((cli/'RUN_RECEIPT.json').read_text());obs=read_observation(root,'gain_marginalized');npz=np.load(cli/'posterior.npz')
    error=float(np.max(np.abs(np.array(receipt['result']['topology_probabilities'])-summary['cases']['gain_marginalized']['result']['topology_probabilities'])))
    agreement={k:bool(np.array_equal(npz[k],v)) for k,v in obs.items()}
    if not all(agreement.values()) or error>1e-12:raise RuntimeError('Actual public gain replay changed target')
    cli_wall=time.perf_counter()-t;cli_cpu=time.process_time()-c
    dump(out/'PUBLIC_INTERFACE_REPLAY.json',dict(all_metadata_equal=all(agreement.values()),metadata_equal=agreement,
          gain_sd=float(receipt['measurement_target']['gain_sd']),max_probability_difference=error,
          material_summary_status=receipt['result']['material_summary_status'],new_physics=0,
          polynomial_counts=receipt['reference']['polynomial_response_counts'],wall_seconds_including_images=cli_wall,cpu_seconds=cli_cpu))
    for phase,wall,used_cpu,counts in [('continuous_CDF_probes',probe_wall,probe_cpu,stats),
          ('actual_public_gain_replay_and_images',cli_wall,cli_cpu,receipt['reference']['polynomial_response_counts'])]:
        row=dict(phase=phase,status='PASS',wall_seconds=wall,cpu_seconds=used_cpu,new_physics_evaluations=0,
                 polynomial_response_evaluations=counts['polynomial_response_evaluations'])
        with (out/'JOB_LEDGER.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
    after=ledger_snapshot(root)
    if before!=after:raise RuntimeError('historical physics ledger changed')
    resources=json.loads((out/'RESOURCE_SUMMARY.json').read_text());resources['supplemental_CDF_and_public_replay']=dict(
         wall_seconds=probe_wall+cli_wall,cpu_seconds=probe_cpu+cli_cpu,
         polynomial_response_evaluations=stats['polynomial_response_evaluations']+receipt['reference']['polynomial_response_counts']['polynomial_response_evaluations'],new_physics=0)
    resources['agent_report_attempts']=dict(native_capacity_failures=3,monetary_and_model_service_time='unmeasured',
                                          fallback='parent Codex writes report; no provider/model change',new_physics=0)
    dump(out/'RESOURCE_SUMMARY.json',resources);summary['resources']=resources;summary['public_interface_replay']='PUBLIC_INTERFACE_REPLAY.json';dump(out/'SUMMARY.json',summary)
    with (out/'FAILURES.jsonl').open('a') as f:
        f.write(json.dumps(dict(phase='native_report_delegation',attempts=3,error='Selected model is at capacity',
                               recovered_by='parent Codex; no model/provider change',new_physics=0,model_service_cost='unmeasured'))+'\n')
    print(json.dumps(dict(probes=len(rows),metadata_replay_max_probability_difference=error,new_physics=0,
                         wall_seconds=time.perf_counter()-start)),flush=True)

if __name__=='__main__':run(Path.cwd())
