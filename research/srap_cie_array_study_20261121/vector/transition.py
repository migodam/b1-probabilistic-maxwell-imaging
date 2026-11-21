"""Stop only this branch's two v1 processes and preserve paid artifacts."""
from pathlib import Path
import psutil,time,json,shutil
root=Path(__file__).resolve().parent;arc=root/'development_ring';arc.mkdir(exist_ok=True)
proc=[]
for p in psutil.process_iter(['pid','name']):
    try:
        if p.name().lower()!='python.exe' or Path(p.cwd()).resolve()!=root:continue
        args=p.cmdline();scripts={Path(v).name for v in args[1:]}
        if not scripts.intersection({'invert.py','finish_remote.py'}):continue
        cpu=p.cpu_times();proc.append(dict(pid=p.pid,cmd=args,wall_since_start=time.time()-p.create_time(),cpu=cpu.user+cpu.system))
        p.terminate();p.wait(10)
    except (psutil.NoSuchProcess,psutil.AccessDenied):continue
for f in ['predictions','data']:
    target=arc/f;target.mkdir(exist_ok=True)
    for p in (root/f).glob('plane_*'):shutil.move(str(p),target/p.name)
for fn in ['physics.py','generate.py','invert.py','projection.py','FREEZE.json','data/generation.json']:
    a=root/fn;b=arc/fn;b.parent.mkdir(exist_ok=True);shutil.copyfile(a,b)
partial=[a.name for a in (arc/'predictions').glob('*') if not (a/'SEALED.json').exists()]
(arc/'transition.json').write_text(json.dumps(dict(proc=proc,time_unix=time.time(),partial_ring=partial,reason='Coplanar single ring is development only; same hemisphere retained.'),indent=2))
print(proc,partial)
