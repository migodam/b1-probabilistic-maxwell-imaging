import time,resource,subprocess,zipfile,json
from pathlib import Path
T=time.perf_counter();C=time.process_time();HERE=Path(__file__).parent;ROOT=HERE.parents[2];REPO=ROOT.parents[1]
KEY='/Users/migodam/.ssh/heatcryst_xinan_ed25519';HOST='jayzh@100.121.97.39';REMOTE='D:/AI/MAX_STRUCT_R9_PHYSREF_20261106'
try:
    subprocess.run(['scp','-i',KEY,str(HERE/'fetch_remote.py'),HOST+':'+REMOTE+'/fetch_remote.py'],check=True)
    subprocess.run(['ssh','-i',KEY,HOST,'D:/python/python.exe '+REMOTE+'/fetch_remote.py'],check=True)
    subprocess.run(['scp','-i',KEY,HOST+':'+REMOTE+'/COMPLETED_PUBLIC_OUTPUTS.zip',str(HERE/'COMPLETED_PUBLIC_OUTPUTS.zip')],check=True)
    with zipfile.ZipFile(HERE/'COMPLETED_PUBLIC_OUTPUTS.zip') as z:
        assert all('sealed_truth' not in p and not p.startswith('/') and '..' not in Path(p).parts for p in z.namelist())
        assert not (ROOT/'RESULTS/original_vector3d/physical_refinement').exists()
        z.extractall(REPO)
    status='COMPLETE'
finally:
    c=resource.getrusage(resource.RUSAGE_CHILDREN)
    (HERE/'FETCH_CONTROLLER_COST.json').write_text(json.dumps(dict(status=locals().get('status','FAILED'),platform='Mac',local_controller_CPU_seconds=time.process_time()-C,local_transport_child_CPU_seconds=c.ru_utime+c.ru_stime,wall_seconds=time.perf_counter()-T,GT_access=False),indent=2))
