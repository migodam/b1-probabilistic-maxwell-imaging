"""Read-only experiment accounting, with explicit hardware/measurement gaps."""
from pathlib import Path
import json,datetime
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'RESULTS/vector_extension';D=ROOT/'delegated/3d_extension'
def read(p):return json.loads(Path(p).read_text())
def duration(s):return sum(float(x)*60**i for i,x in enumerate(reversed(s.split(':'))))

def main():
    tasks=read(ROOT/'configs/vector_extension_all29.json')['tasks'];rows=[]
    for t in tasks:
        r=read(OUT/t['case_id']/t['arm']/'RECEIPT.json');rows.append(dict(case_id=t['case_id'],arm=t['arm'],platform=r.get('platform','Mac_AppleSilicon'),
            status=r['status'],CPU_seconds=r['phase2_incremental_CPU_seconds'],wall_seconds=r['phase2_incremental_wall_seconds'],counts=r.get('counts',{})))
    grouped={}
    for r in rows:
        g=grouped.setdefault(r['platform'],dict(completed_scientific_tasks=0,whole_worker_CPU_seconds=0,physical_LU=0,physical_source_RHS=0,physical_adjoint_RHS=0))
        g['completed_scientific_tasks']+=1;g['whole_worker_CPU_seconds']+=r['CPU_seconds']
        for counter in ('physical_LU','physical_source_RHS','physical_adjoint_RHS'):g[counter]+=r['counts'].get(counter,0)
    abort=read(D/'LOCAL_REMOTE_MIGRATION_ABORT.json');snap={}
    for line in abort['process_snapshot'].splitlines():
        tokens=line.split();snap[int(tokens[0])]=dict(observed_CPU_lower_bound_seconds=duration(tokens[3]),elapsed=tokens[2])
    old=read(D/'LOCAL_AMENDMENT_COORDINATOR_BEFORE_MIGRATION.json');mie=read(OUT/'validation_refinement16/RECEIPT.json');final=read(OUT/'FINALIZATION_COST.json');remote=read(D/'xinan_resolved_queue/PARENT_COST.json');transport=read(D/'REMOTE_TRANSPORT.json')
    local_aborted=snap[46361]['observed_CPU_lower_bound_seconds']+snap[46362]['observed_CPU_lower_bound_seconds']
    timing={line.split()[0]:float(line.split()[1]) for line in transport.get('return_transport_timing','').splitlines() if line.split() and line.split()[0] in ('real','user','sys')}
    return_transport_cpu=timing.get('user',0)+timing.get('sys',0)
    render=read(OUT/'RENDER_COST.json') if (OUT/'RENDER_COST.json').exists() else {'CPU_seconds':0,'scope':'notrun'}
    mac_extra=old['original_parent_CPU_observed_lower_bound_seconds']+snap[39368]['observed_CPU_lower_bound_seconds']+snap[46345]['observed_CPU_lower_bound_seconds']+local_aborted+mie['CPU_seconds']+final['CPU_seconds']+return_transport_cpu+render['CPU_seconds']
    remote_extra=remote['parent_CPU_seconds']+transport.get('remote_result_packaging_CPU_seconds',0)
    costs=dict(status='COMPLETED_WITH_EXPLICIT_MEASUREMENT_LIMITS',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        task_counts=dict(primary_Mac=24,resolution_control_XINAN=5,aborted_Mac_resource_handoff=2,Mie_validation_LU=3),
        scientific_workers_by_hardware=grouped,
        Mac_incremental_ancillary=dict(original24_parent_CPU_lower_bound=old['original_parent_CPU_observed_lower_bound_seconds'],
            original24_parent_additional_CPU_bound=old['original_parent_CPU_unobserved_upper_bound_extra_seconds'],
            source='ps sampled old parent until exit; includes original8/10 Mie, originalgeneration, imports, queue; do not add their nested receipts again',
            amendment_metadata_coordinator_CPU_lower_bound=snap[39368]['observed_CPU_lower_bound_seconds'],
            local_control_parent_CPU_lower_bound=snap[46345]['observed_CPU_lower_bound_seconds'],
            local_control_parent_scope='includes one necessary fine10/coarse8 resolvedhollow generation, imports, abortedqueue; generation not repeated remotely',
            aborted_workers_CPU_lower_bound=local_aborted,aborted_available_physics_lowerbounds=abort['partial_counters_lowerbounds'],
            Mie16_CPU_seconds=mie['CPU_seconds'],offline_seal_score_figures_CPU_seconds=final['CPU_seconds'],postscore_render_CPU_seconds=render['CPU_seconds'],return_transport_CPU_seconds=return_transport_cpu),
        XINAN_ancillary=dict(queue_parent_CPU_seconds=remote['parent_CPU_seconds'],packaging_CPU_seconds=transport.get('remote_result_packaging_CPU_seconds')),
        retained_total_CPU_lower_bounds_by_hardware=dict(Mac_AppleSilicon=grouped['Mac_AppleSilicon']['whole_worker_CPU_seconds']+mac_extra,
            XINAN_Windows_i9=grouped['XINAN_Windows_i9']['whole_worker_CPU_seconds']+remote_extra),
        transport=transport,
        missing_or_bounded_costs=['Original already-running parent CPU has an explicit polling bound, not exact postmortem measurement.',
            'Killed Mac control/metadata processes use immediately pre-signal observed CPU lower bounds; terminal unobserved CPU is not falsely zero.',
            'Small metadata ps/upload/remote-extract helpers were not separately CPU-metered; transport record gives wall times. These are unmeasured overhead, not scientific zero-cost operations.'],
        comparison_scope='Within same object and hardware only. Crossplatform CPU sum is consumption bookkeeping, not a pooled speedup claim. Original24 Mac predictions/configs preserved; five main resolvedhollow arms all XINAN.',
        validation_generation_scope='9 newlypaid forwardLU outside inversions: Mie8/10/16=3; twooriginalnewcasesfine/coarse=4; resolvedhollowfine/coarse=2. Legacyfourcubesfine10 not rerun.',
        phaseI_modified=False,GT_policy='all29 public predictions sealed before offline scoring; no GTonline',rows=rows)
    (OUT/'PHASE2_COMPLETE.json').write_text(json.dumps(costs,indent=2)+'\n');return costs

if __name__=='__main__':print(json.dumps(main(),indent=2))
