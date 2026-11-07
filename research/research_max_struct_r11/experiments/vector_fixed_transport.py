"""Finite R11 public fixed-state CUDA deployment, preserves R10 files."""
import time,subprocess,resource,sys,json,zipfile,argparse
from pathlib import Path
T=time.perf_counter();C=time.process_time();ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
p=argparse.ArgumentParser();p.add_argument('--tag',default='fixed_state');p.add_argument('--basis',default='research/research_max_struct_r10/data_public/vector/operator_spectral/GIN_TOP32_CAPTURE.npz');args=p.parse_args()
OUT=ROOT/('delegated/vector/'+args.tag+'_transport');OUT.mkdir(parents=True,exist_ok=False)
REMOTE='D:/AI/MAX_STRUCT_R10_VECTOR_20261107';HOST='jayzh@100.121.97.39';KEY='/Users/migodam/.ssh/heatcryst_xinan_ed25519'
SSH=['ssh','-i',KEY,'-o','BatchMode=yes','-o','ServerAliveInterval=15','-o','ServerAliveCountMax=3',HOST]
paths=[ROOT/'src/vector_denominator.py',ROOT/'experiments/vector_fixed_state.py',ROOT/'delegated/vector/PUBLIC_V04_FIXED_STATE.npz']
steps=[]
def run(command,label):
 start=time.perf_counter()
 with (OUT/(label+'.log')).open('w') as log:
  process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT);status=process.wait()
 steps.append(dict(operation=label,exit=status,wall_seconds=time.perf_counter()-start))
 if status:raise RuntimeError(label+' failed: '+str(status))
 print(json.dumps(steps[-1]),flush=True)
try:
 with zipfile.ZipFile(OUT/'SOURCE.zip','w',zipfile.ZIP_DEFLATED) as z:
  for path in paths:z.write(path,str(path.relative_to(REPO)))
  z.write(REPO/args.basis,'research/research_max_struct_r11/delegated/vector/'+args.tag+'_BASIS.npz')
 run(['scp','-i',KEY,str(OUT/'SOURCE.zip'),HOST+':'+REMOTE+'/R11_'+args.tag+'_SOURCE.zip'],'upload')
 code="import zipfile;zipfile.ZipFile('"+REMOTE+"/R11_"+args.tag+"_SOURCE.zip').extractall('"+REMOTE+"')"
 run(SSH+['D:/python/python.exe -c "'+code+'"'],'extract')
 cmd='D:/python/python.exe -u '+REMOTE+'/research/research_max_struct_r11/experiments/vector_fixed_state.py --basis '+REMOTE+'/research/research_max_struct_r11/delegated/vector/'+args.tag+'_BASIS.npz --state '+REMOTE+'/research/research_max_struct_r11/delegated/vector/PUBLIC_V04_FIXED_STATE.npz --output delegated/vector/'+args.tag+'_cuda --device cuda'
 run(SSH+[cmd],'run')
 code="import zipfile;from pathlib import Path;r=Path('"+REMOTE+"');z=zipfile.ZipFile(r/'R11_"+args.tag+"_RESULT.zip','w',zipfile.ZIP_DEFLATED);[z.write(p,str(p.relative_to(r))) for p in (r/'research/research_max_struct_r11/delegated/vector/"+args.tag+"_cuda').rglob('*') if p.is_file()];z.close()"
 run(SSH+['D:/python/python.exe -c "'+code+'"'],'pack')
 run(['scp','-i',KEY,HOST+':'+REMOTE+'/R11_'+args.tag+'_RESULT.zip',str(OUT/'RESULT.zip')],'download')
 with zipfile.ZipFile(OUT/'RESULT.zip') as z:z.extractall(REPO)
 status='COMPLETE'
except BaseException as exc:status='FAILED';print(repr(exc),flush=True)
finally:
 child=resource.getrusage(resource.RUSAGE_CHILDREN)
 (OUT/'TRANSPORT_COST.json').write_text(json.dumps(dict(status=status,platform='Mac',local_controller_CPU_seconds=time.process_time()-C,local_transport_child_CPU_seconds=child.ru_utime+child.ru_stime,wall_seconds=time.perf_counter()-T,steps=steps,scope='Local connected SSH/SCP only; remote scientific CPU excluded'),indent=2)+'\n')
if status!='COMPLETE':raise SystemExit(1)
