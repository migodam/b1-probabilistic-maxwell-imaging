"""One read-only Windows processor-information snapshot; no configuration edits."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,subprocess,json,time,platform

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--phase',required=True);a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    start=datetime.now(timezone.utc).isoformat();cpu=time.process_time()
    command=('Get-CimInstance Win32_PerfFormattedData_Counters_ProcessorInformation | '
        'Select-Object Name,PercentProcessorTime,ProcessorFrequency,PercentProcessorPerformance,'
        'PercentofMaximumFrequency,PercentPerformanceLimit,PerformanceLimitFlags | ConvertTo-Json -Compress')
    result=subprocess.run(['powershell','-NoProfile','-Command',command],capture_output=True,text=True)
    raw=json.loads(result.stdout.lstrip('\ufeff')) if result.returncode==0 and result.stdout.strip() else None
    value=dict(status='READ_ONLY_SNAPSHOT' if raw is not None else 'SNAPSHOT_QUERY_FAILED',phase=a.phase,
        started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),returncode=result.returncode,
        counter_source='Win32_PerfFormattedData_Counters_ProcessorInformation raw named fields',
        values=raw,stderr=result.stderr,scope='Cost environment; not a scientific algorithm intervention',
        power_frequency_affinity_settings_changed=False,own_process_CPU_seconds=time.process_time()-cpu,
        query_helper_CPU_seconds=None,query_helper_CPU_scope='PowerShell child CPU not separately measured here; transport child/runtime are nested',platform=platform.platform())
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.open('x').write(json.dumps(value,indent=2)+'\n')
    print(json.dumps(dict(status=value['status'],phase=a.phase,path=str(a.output),UTC=value['finished_UTC'])),flush=True)

if __name__=='__main__':main()
