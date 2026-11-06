"""Metadata-only independent physical-refinement ledger; old9 untouched."""
import time
T=time.perf_counter();C=time.process_time()
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'RESULTS/original_vector3d/physical_refinement';D=ROOT/'delegated/3d_extension/physical_refinement'
def read(p):return json.loads(Path(p).read_text())
def source(p):
    try:return str(Path(p).relative_to(ROOT))
    except ValueError:return str(Path(p).relative_to(ROOT.parents[1]))
def main():
    rows=[];whole=read(OUT/'RUN_PROCESS_COST.json');assert whole['status']=='COMPLETE'
    def add(p,platform,cpu,wall,scope,operation,counts,additive=True,included_in=None):
        rows.append(dict(artifact=source(p),platform=platform,CPU_seconds=cpu,wall_seconds=wall,additive=additive,included_in=included_in,scope=scope,operation=operation,counts=counts))
    add(OUT/'RUN_PROCESS_COST.json','XINAN',whole['CPU_seconds'],whole['wall_seconds'],whole['scope'],whole['operation'],whole['counts'])
    for name in ('air','GlobalDual_stage2','AplusGlobal_stage2'):
        p=OUT/name/'RECEIPT.json';r=read(p);add(p,'XINAN',r['CPU_seconds'],r['wall_seconds'],'Nested perstart interval, uniqueinitial included; initialfun cachehit not a second solve','ONE_PHYSICAL_START',r['counts'],False,'RUN_PROCESS_COST.json')
    ev=read(OUT/'PAID_EVALUATIONS.json');we=[x for x in ev if x['owner'].startswith('gradient_wiring_')]
    add(OUT/'GRADIENT_WIRING.json','XINAN',sum(x['CPU_seconds'] for x in we),sum(x['wall_seconds'] for x in we),'Nested air analytic gradient plus one qualified FD forward batch','GRADIENT_WIRING',dict(primal_RHS=sum(x['actual_primal_RHS'] for x in we),adjoint_RHS=sum(x['actual_adjoint_RHS'] for x in we)),False,'RUN_PROCESS_COST.json')
    for p in [OUT/'REMOTE_FETCH_PROCESS_COST.json',OUT/'LOCAL_REPORT_PROCESS_COST.json',OUT/'LOCAL_FIGURE_REVISION_ATTEMPT1_COST.json',OUT/'LOCAL_FIGURE_REVISION_COST.json',OUT/'LOCAL_FINAL_REPORT_COST.json']:
        r=read(p);add(p,r['platform'],r['CPU_seconds'],r['wall_seconds'],r['scope'],r['operation'],r['counts'])
    r=read(ROOT.parents[1]/'EXTRACTION_COST.json');add(ROOT.parents[1]/'EXTRACTION_COST.json','XINAN',r['CPU_seconds'],r['wall_seconds'],'Remote isolated public-code/data archive extraction only','PUBLIC_BUNDLE_EXTRACTION',{})
    for p in [D/'TRANSPORT_CONTROLLER_COST.json',D/'FETCH_CONTROLLER_COST.json']:
        r=read(p)
        for key in ('local_controller_CPU_seconds','local_transport_child_CPU_seconds'):
            add(p,'Mac',r[key],None,'Own '+key+' only; remote worker/helper CPU excluded','CONNECTED_PUBLIC_TRANSPORT',{})
    add(D/'LAUNCH.json','Mac',None,None,'Syntax-only AST checks and first connection CUDA availability/affinity probe before transport, whole CPU unmeasured; no scientific RHS','UNMEASURED_PREFLIGHT',dict(new_physical_RHS=0))
    add(D/'LOCAL_ENGINEERING_FAILURES.json','Mac',None,.348026792,'First failed offline guard plus dependent ledger call; no truth read, no new physical RHS','FAILED_LOCAL_PATH_PORTABILITY',dict(new_physical_RHS=0))
    for f in read(D/'LOCAL_ENGINEERING_FAILURES.json').get('attempts',[])[1:]:
        add(D/'LOCAL_ENGINEERING_FAILURES.json','Mac',f.get('CPU_seconds'),f.get('wall_seconds'),f['reason'],f['operation'],dict(new_physical_RHS=0))
    history=read(OUT/'COST_COMPILER_HISTORY.json') if (OUT/'COST_COMPILER_HISTORY.json').exists() else []
    history.append(dict(invocation=len(history)+1,CPU_seconds=time.process_time()-C,wall_seconds=time.perf_counter()-T,source='Measured metadata compiler pre-export CPU'))
    (OUT/'COST_COMPILER_HISTORY.json').write_text(json.dumps(history,indent=2)+'\n')
    for h in history:add(OUT/'COST_COMPILER_HISTORY.json','Mac',h['CPU_seconds'],h['wall_seconds'],'Independent metadata compiler invocation '+str(h['invocation'])+'; excludes physics/score/transport','COMPILE_INDEPENDENT_LEDGER',dict(new_physical_RHS=0))
    by={p:sum(x['CPU_seconds'] for x in rows if x['additive'] and x['platform']==p and x['CPU_seconds'] is not None) for p in ('Mac','XINAN')}
    result=dict(status='COMPLETE_RECORDED_WORK_KNOWN_CPU_LOWER_BOUND',classification='PUBLIC_DEVELOPMENT_SHORT_PHYSICAL_REFINEMENT',rows=rows,by_platform=by,known_CPU_seconds=sum(by.values()),
        actual_full_primal_RHS=whole['counts'].get('primal_RHS',0),actual_full_adjoint_RHS=whole['counts'].get('adjoint_RHS',0),
        nonadditive_note='Perarm/wiring/setup/evaluation/cache views are already inside single GPU process. Wall intervals are not campaign makespan and are not summed. No LU performed.',
        unknowns=['Unmeasured connection-probe CPU and pre-export process/output tails remain unknown, not zero; no claim of total GPU elapsed from host CPU/wall.'],
        old_nine_canonical_modified=False)
    (OUT/'ACCOUNTING_FINAL.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],by_platform=by,known_CPU_seconds=result['known_CPU_seconds'],primal=result['actual_full_primal_RHS'],adjoint=result['actual_full_adjoint_RHS'])))
if __name__=='__main__':main()
