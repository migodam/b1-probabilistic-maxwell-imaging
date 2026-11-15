"""Small offline-review transport plus separate full derived-analysis ZIP.

The original full-state campaign archives are separate and unchanged. This
utility runs after the authorized scorer completes; it never reads truth.
"""
import argparse,json,time,zipfile
from pathlib import Path


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--campaign',required=True);a=p.parse_args();r=a.root.resolve();campaign=a.campaign
    analysis=r/'analysis'/campaign;receipt=json.loads((analysis/'RECEIPT.json').read_text())
    if receipt['status']!='AUTHORIZED_EXTENSION_OFFLINE_SCORING_COMPLETE':raise ValueError('Completed offline analysis required')
    cpu=time.process_time();wall=time.perf_counter();dest=r/'archives';dest.mkdir(exist_ok=True)
    excluded={'NATIVE_STAGE_HISTORY.json','NATIVE_STAGE_HISTORY.csv','PAID_CANDIDATE_HISTORY.json','PAID_CANDIDATE_HISTORY.csv'}
    compact=[p for p in analysis.rglob('*') if p.is_file() and p.name not in excluded]
    for directory in sorted((r/'RESULTS'/campaign).iterdir()):
        compact.extend(p for p in directory.glob('*.json') if p.name in ('RECEIPT.json','ACTUAL_CONFIG.json','WRAPPER_PROVENANCE.json'))
    compact.extend(p for p in (r/'queues'/campaign).glob('*.json') if p.name in ('RECEIPT.json','FROZEN_MANIFEST.json','RESOURCE_PROBE.json'))
    compact.extend(p for p in (r/'audit').glob('*independent*SEAL*.json'))
    compact.extend(p for p in (r/'audit').glob('*INDEPENDENT*.json'))
    config=r/'configs'/(campaign+'.json');compact.append(config)
    files=list(dict.fromkeys(compact));paths=[]
    for label,items in [('REVIEW',files),('ANALYSIS_FULL',[p for p in analysis.rglob('*') if p.is_file()])]:
        out=dest/(campaign+'_'+label+'.zip');before_cpu=time.process_time();before_wall=time.perf_counter()
        with zipfile.ZipFile(out,'x',zipfile.ZIP_DEFLATED,compresslevel=1) as z:
            for p in sorted(items):
                z.write(p,p.relative_to(r).as_posix(),compress_type=zipfile.ZIP_STORED if p.suffix.lower() in ('.npz','.png','.pdf') else zipfile.ZIP_DEFLATED)
        paths.append(dict(name=out.name,bytes=out.stat().st_size,members=len(items),
            CPU_seconds=time.process_time()-before_cpu,wall_seconds=time.perf_counter()-before_wall))
    result=dict(status='REVIEW_AND_FULL_DERIVED_ANALYSIS_ARCHIVED',campaign=campaign,assets=paths,
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,new_physics=0,
        truth_access='NONE',selected_material_transport='No material array changes; full original field archive remains authoritative')
    (dest/(campaign+'_REVIEW_EXPORT_RECEIPT.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
