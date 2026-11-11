"""Compact public scoring handoff; preserve full raw trajectories remotely.

Copies already sealed materials/scalar metadata and four frozen prefixes.
No GT, physics/model imports, optimizers, or new forward evaluations.
"""
from pathlib import Path
import json,zipfile,io,time
import numpy as np
ROOT=Path(__file__).resolve().parent

def main():
    archive=ROOT/'R12_SELECTED_PUBLIC_OUTPUTS_20261111.zip'
    if archive.exists():raise FileExistsError(archive)
    rows=[];started=time.perf_counter()
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for stage in ('dev6','confirm24'):
            sealed=ROOT/('SEALED_'+stage.upper()+'_PUBLIC_COHORT.json')
            cohort=json.loads(sealed.read_text())
            if cohort['status']!='COMPLETE_PUBLIC_COHORT_SEALED_BEFORE_GT':raise ValueError('Public cohort is not sealed')
            z.write(sealed,sealed.name)
            for row in cohort['rows']:
                task=ROOT/'RESULTS'/stage/row['task_id'];dest=Path(stage)/row['task_id']
                # Float64 material bytes are copied without reoptimization/reselection.
                with np.load(row['prediction'],allow_pickle=False) as data:m=data['m'].copy()
                buf=io.BytesIO();np.savez_compressed(buf,m=m);z.writestr(str(dest/'primary.npz'),buf.getvalue())
                secondary=task/'NOISE_FEASIBLE_MIN_PRIOR.npz'
                if secondary.exists():
                    with np.load(secondary,allow_pickle=False) as data:m=data['m'].copy()
                    buf=io.BytesIO();np.savez_compressed(buf,m=m);z.writestr(str(dest/'secondary.npz'),buf.getvalue())
                keep=('RECEIPT.json','R11_CAP_PREDICTION.json','R11_CPU_CHECKPOINTS.json','R11_CHECKPOINT_MATERIALS.npz',
                    'paid_prefix.npz','NOISE_FEASIBLE_MIN_PRIOR.json','NOISE_FEASIBLE_MIN_PRIOR_TIMEPREFIX.json',
                    'NOISE_FEASIBLE_MIN_PRIOR_TIMEPREFIX.npz','CONFIG.json','FROZEN_CONFIG.json')
                for name in keep:
                    if (task/name).exists():z.write(task/name,str(dest/name))
                metadata=dict(row,compact_primary=str(dest/'primary.npz'),
                    compact_secondary=str(dest/'secondary.npz') if secondary.exists() else None,
                    original_complete_remote_directory=str(task))
                if row['method'] in ('dr','multistart','nn_direct'):
                    history=json.loads((task/'history.json').read_text())
                    eligible=[r for r in history if r['eligible']]
                    metadata['primary_paid_audit']=min(eligible,key=lambda r:r['Phi'])
                rows.append(metadata)
        for folder in ('queues','engineering'):
            for path in (ROOT/folder).rglob('*'):
                if path.is_file() and path.suffix in ('.json','.jsonl','.py') and 'CHECKPOINTS' not in path.parts:
                    z.write(path,str(path.relative_to(ROOT)))
        for name in ('ROOT_CORE_FREEZE.json','SECONDARY_SELECTION_FREEZE.json'):
            z.write(ROOT/name,name)
        manifest=dict(status='SEALED_PUBLIC_SELECTED_OUTPUTS',rows=rows,predictions=len(rows),GT='NONE',new_FW=0,
            scope='Selected full-precision materials, original receipts and four frozen CPU prefixes; complete raw m/J/Q/field and all paid trajectories retained unchanged on remote.',
            export_process_CPU_seconds_near_write=time.process_time(),wall_seconds=time.perf_counter()-started)
        z.writestr('COMPACT_PUBLIC_MANIFEST.json',json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(status=manifest['status'],archive=str(archive),bytes=archive.stat().st_size,
        predictions=len(rows),CPU_seconds=time.process_time(),GT='NONE',new_FW=0)),flush=True)

if __name__=='__main__':main()
