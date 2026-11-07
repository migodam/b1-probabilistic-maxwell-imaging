"""Own Windows R10 metadata snapshot only; never loads NPZ or GT."""
import time
T=time.perf_counter()
from pathlib import Path
import json,os
import psutil
ROOT=Path(__file__).resolve().parents[1]


def load(path):return json.loads(path.read_text())


def main():
    snapshots=[];parents=[];active=[]
    for folder in sorted((ROOT/'delegated/vector').glob('windows_*_queue*')):
        rows=folder/'PROCESS_RECEIPTS.json'
        if rows.exists():
            for row in load(rows):snapshots.append(dict(source=str(rows.relative_to(ROOT)),**row))
        single=folder/'PROCESS_COST.json'
        if single.exists():snapshots.append(dict(source=str(single.relative_to(ROOT)),label=folder.name,**load(single)))
        p=folder/'CONTROLLER_COST.json'
        if p.exists():parents.append(dict(source=str(p.relative_to(ROOT)),**load(p)))
        marker=folder/'ACTIVE.json'
        if marker.exists() and not p.exists():
            info=load(marker);info.setdefault('label',folder.name)
            try:
                process=psutil.Process(info['pid']);command=' '.join(process.cmdline())
                if 'MAX_STRUCT_R10_VECTOR_20261107' not in command or 'vector_run.py' not in command:continue
                cpu=process.cpu_times();memory=process.memory_info()
                active.append(dict(source=str(marker.relative_to(ROOT)),**info,CPU_seconds=cpu.user+cpu.system,
                    wall_seconds=time.time()-process.create_time(),RSS_bytes=memory.rss,peak_RSS_bytes=memory.peak_wset,
                    status='RUNNING_CUMULATIVE_LOWER_BOUND_NOT_ADDITIVE_WITH_FINAL'))
            except psutil.NoSuchProcess:pass
    endpoints=[];seals=[]
    for folder in sorted((ROOT/'RESULTS/vector').glob('*')):
        receipt=folder/'RECEIPT.json'
        if receipt.exists() and folder.name.startswith(('v04_','v08_')):
            d=load(receipt);corrector=d.get('physical_corrector',{})
            endpoints.append(dict(task_id=folder.name,CPU_seconds=d.get('CPU_seconds'),wall_seconds=d.get('wall_seconds'),
                Phi=d.get('physical_Phi'),PG_inf=d.get('PG_inf'),counts=d.get('counts'),stages=d.get('stages'),
                physical_corrector={k:v for k,v in corrector.items() if k!='events'},
                CUDA_event_stream_interval_ms=d.get('CUDA_event_stream_interval_ms'),
                GPU_peak_reserved_bytes_to_date=d.get('GPU_peak_reserved_bytes_to_date'),
                SEALED=(folder/'SEALED.json').exists(),GT_score='NOT_RUN'))
            if (folder/'SEALED.json').exists():seals.append(folder.name)
    process_headers=[]
    for name in ('generation_windows_v04_attempt1_missing_dependency','generation_windows_v04','generation_gpu_v08','case_v04_windows_process','case_v08_gpu_process','case_v08_cuda_dense_process'):
        p=ROOT/'RESULTS/vector'/name/'PROCESS_COST.json'
        if p.exists():
            row=load(p);row.pop('actual_solve_records',None)
            process_headers.append(dict(source=str(p.relative_to(ROOT)),**row))
    current=psutil.Process();cpu=current.cpu_times()
    receipt=dict(status='READ_ONLY_OWN_R10_METADATA_SNAPSHOT',completed_process_rows=snapshots,controller_rows=parents,
        active_process_snapshots=active,process_cost_headers=process_headers,endpoints=endpoints,sealed=seals,
        helper_CPU_seconds=cpu.user+cpu.system,helper_wall_seconds=time.perf_counter()-T,platform='XINAN',
        GT_arrays_read=False,new_solvers=0,snapshot_CPU_semantics='Active counters replace earlier same-process counters, never sum snapshots/final together')
    out=ROOT/'delegated/vector/REMOTE_STATUS_SNAPSHOT.json';out.write_text(json.dumps(receipt,indent=2)+'\n')
    with (ROOT/'delegated/vector/STATUS_HELPERS.jsonl').open('a') as handle:handle.write(json.dumps(dict(pid=os.getpid(),CPU_seconds=cpu.user+cpu.system,wall_seconds=receipt['helper_wall_seconds'],platform='XINAN',scope='Own one metadata helper process, no scientific child CPU'))+'\n')
    print(json.dumps(dict(sealed=seals,active=[dict(pid=p['pid'],label=p['label'],CPU_seconds=p['CPU_seconds'],wall_seconds=p['wall_seconds']) for p in active],helper_CPU_seconds=receipt['helper_CPU_seconds'])))


if __name__=='__main__':main()
