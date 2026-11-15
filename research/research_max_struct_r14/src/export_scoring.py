"""Material-only transport of sealed predictions; full field originals retained.

No selection, truth, optimizer, forward model, numerical cast or field recompute.
"""
import argparse,json,shutil,time,zipfile
from pathlib import Path


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--campaign',required=True);a=p.parse_args()
    import numpy as np
    r=a.root;q=r/'queues'/a.campaign;terminal=json.loads((q/'RECEIPT.json').read_text())
    if terminal.get('running') or terminal.get('pending_tasks',0):raise RuntimeError('Only terminal predictions are exported')
    cpu=time.process_time();wall=time.perf_counter();dest=r/'exports'/(a.campaign+'_scoring');dest.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((q/'FROZEN_MANIFEST.json').read_text());files=[]
    for task in manifest['tasks']:
        folder=r/'RESULTS'/a.campaign/task['task_id'];target=dest/'RESULTS'/a.campaign/task['task_id'];target.mkdir(parents=True)
        for name in ['RECEIPT.json','ACTUAL_CONFIG.json','HISTORY.json','NATIVE_HISTORY.json','WRAPPER_PROVENANCE.json','PREFIXES.json','PROJECTION_HISTORY.json']:
            source=folder/name
            if source.exists():shutil.copy2(source,target/name)
        for name in ('initial','raw','corrected','primary','secondary'):
            source=folder/(name+'.npz')
            if not source.exists():continue
            with np.load(source,allow_pickle=False) as z:
                material=z['m'];np.savez_compressed(target/(name+'.npz'),m=material)
            metadata=folder/(name+'.json')
            if metadata.exists():shutil.copy2(metadata,target/metadata.name)
            files.append(dict(original=str(source.relative_to(r)),original_file_bytes=source.stat().st_size,
                              exported=str((target/(name+'.npz')).relative_to(dest)),retained_keys=['m'],
                              material_shape=list(material.shape),dtype=str(material.dtype),value_transform='NONE'))
        for name in ('CANDIDATE_POOL.npz','CANDIDATES.npz'):
            if (folder/name).exists():shutil.copy2(folder/name,target/name)
    (dest/'queues'/a.campaign).mkdir(parents=True)
    for name in ('RECEIPT.json','FROZEN_MANIFEST.json','RESOURCE_PROBE.json'):
        if (q/name).exists():shutil.copy2(q/name,dest/'queues'/a.campaign/name)
    (dest/'configs').mkdir();shutil.copy2(q/'FROZEN_MANIFEST.json',dest/'configs'/(a.campaign+'.json'))
    receipt=dict(status='DERIVED_MATERIAL_ONLY_NO_NUMERICAL_CHANGE',campaign=a.campaign,files=files,
                 tasks=len(manifest['tasks']),full_originals_retained=str(r/'RESULTS'/a.campaign),
                 CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
                 truth_access='NONE',new_physical_LU=0,new_complex_RHS=0,scope='offline scoring transport; not full field archive')
    (dest/'OFFLINE_EXPORT_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    path=r/'archives'/(a.campaign+'_scoring_only.zip');zcpu=time.process_time();zwall=time.perf_counter()
    with zipfile.ZipFile(path,'x',zipfile.ZIP_DEFLATED,compresslevel=1) as z:
        for f in sorted(dest.rglob('*')):
            if f.is_file():z.write(f,str(f.relative_to(dest)),compress_type=zipfile.ZIP_STORED if f.suffix=='.npz' else zipfile.ZIP_DEFLATED)
    receipt.update(zip_CPU_seconds=time.process_time()-zcpu,zip_wall_seconds=time.perf_counter()-zwall,archive_bytes=path.stat().st_size)
    (r/'archives'/(a.campaign+'_SCORING_EXPORT_RECEIPT.json')).write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='files'}))


if __name__=='__main__':main()
