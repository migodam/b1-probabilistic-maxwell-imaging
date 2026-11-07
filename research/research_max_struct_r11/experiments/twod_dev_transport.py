"""Finite explicitly authorized R11 exposed DEV deployment; no truth access."""
import json,subprocess,time,resource,tarfile,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
OUT=ROOT/'delegated/twod/dev_transport';OUT.mkdir(parents=True,exist_ok=True)
HOST='jayzh@100.121.97.39';KEY='/Users/migodam/.ssh/heatcryst_xinan_ed25519';REMOTE='D:/AI/MAX_STRUCT_R11_DEV_20261107'
def run(cmd,label):
 t,c=time.perf_counter(),time.process_time();before=resource.getrusage(resource.RUSAGE_CHILDREN)
 p=subprocess.run(cmd,capture_output=True,text=True,encoding='utf-8',errors='replace');after=resource.getrusage(resource.RUSAGE_CHILDREN)
 with (OUT/'TRANSPORT_COSTS.jsonl').open('a') as f:f.write(json.dumps(dict(label=label,returncode=p.returncode,wall_seconds=time.perf_counter()-t,controller_CPU_seconds=time.process_time()-c,child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime))+'\n')
 if p.returncode:raise RuntimeError(p.stderr+'\n'+p.stdout)
 return p.stdout

def ssh(command,label):return run(['ssh','-i',KEY,'-o','BatchMode=yes','-o','ServerAliveInterval=30',HOST,command],label)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['prepare','deploy','probe','replay','lagged','fetch']);args=ap.parse_args()
 if args.mode=='prepare':
  cfg=json.loads((REPO/'research/research_max_struct_r10/configs/twod_solver.json').read_text());cfg.update(run_id='r11_exposed_dev_480CPU',whole_child_CPU_budget_seconds=480,CPU_budget_checkpoints_seconds=[60,120,240,480],workers=2,freeze_authority='Root direct authorized exposed DEV; no GT selection',experiment_scope='EXPOSED_DEV_NOT_CONFIRMATION')
  (ROOT/'configs/twod_dev_solver.json').write_text(json.dumps(cfg,indent=2)+'\n')
  tasks=[]
  for case,label in [('r10_twod_0001','full_replay'),('r10_twod_0001','lagged'),('r10_twod_0004','lagged'),('r10_twod_0022','lagged'),('r10_conditional_0000','lagged'),('r10_conditional_0003','lagged')]:
   kind='twod' if case.startswith('r10_twod') else 'conditional_confirmation'
   ident=case+'__'+label+'__480CPU'
   tasks.append(dict(task_id=ident,case_id=case,arm='a_plus',initialization='uniform',initial_material=.1,dpcs_rule='noise_ball_99',public_files=['research/research_max_struct_r10/data_public/'+kind+'/operators.npz','research/research_max_struct_r10/data_public/'+kind+'/'+case+'.npz'],output_dir='research/research_max_struct_r11/RESULTS/exposed_dev/'+ident,experiment_scope='EXPOSED_DEV_NOT_CONFIRMATION',material_readout=dict(module_file='research/research_max_struct_r11/src/material_readouts.py',callable='full_material_readout' if label=='full_replay' else 'lagged_material_readout',prepare='prepare_full' if label=='full_replay' else 'prepare_lagged',label=label)))
  for kind,selected in [('replay',tasks[:1]),('lagged',tasks[1:])]:
   (ROOT/('configs/twod_dev_'+kind+'_tasks.json')).write_text(json.dumps(dict(status='ROOT_FROZEN_EXPOSED_DEV',tasks=selected),indent=2)+'\n')
  files=[]
  for study in ['research_max_struct_r2','research_max_struct_r4','research_max_struct_r5','research_max_struct_r6','research_max_struct_r7','research_max_struct_r8','research_fbg_r3_experiments','research_z_alg_r4_diversity','research_z_alg_r4_basin']:
   files.extend((REPO/'research'/study).rglob('*.py'))
  files += [REPO/'research/research_max_struct_r10/src/twod_adapter.py',REPO/'research/research_max_struct_r10/src/lowrank_denominator.py']
  files += list((ROOT/'src').glob('*.py'))+[ROOT/'experiments/twod_run.py',ROOT/'experiments/twod_queue.py']+list((ROOT/'configs').glob('twod_dev*.json'))
  for task in tasks:files += [REPO/p for p in task['public_files']]
  with tarfile.open(OUT/'inverse_only.tgz','w:gz') as tar:
   for p in sorted(set(files)):
    if any(x in str(p) for x in ['sealed_truth','data_private']):raise ValueError('private')
    tar.add(p,arcname=str(p.relative_to(REPO)),recursive=False)
  (OUT/'PACKAGE_RECEIPT.json').write_text(json.dumps(dict(status='PREPARED',truth_access='NONE',files=[dict(path=str(p.relative_to(REPO)),bytes=p.stat().st_size) for p in sorted(set(files))],bytes=(OUT/'inverse_only.tgz').stat().st_size),indent=2)+'\n')
  print('PREPARED',len(files));return
 if args.mode=='probe':
  code="import json,psutil; print(json.dumps(dict(load=psutil.cpu_percent(interval=2,percpu=True),python=[dict(pid=p.pid,affinity=p.cpu_affinity()) for p in psutil.process_iter() if 'python' in p.name().lower()])))"
  s=ssh('D:/python/python.exe -c "'+code+'"','probe');(OUT/'P8_P10_PROBE.json').write_text(s);print(s);return
 if args.mode=='deploy':
  ssh('powershell -NoProfile -Command "New-Item -ItemType Directory -Force '+REMOTE+' | Out-Null"','mkdir')
  run(['scp','-i',KEY,str(OUT/'inverse_only.tgz'),HOST+':'+REMOTE+'/inverse_only.tgz'],'upload')
  code="import tarfile,shutil,pathlib; tarfile.open('"+REMOTE+"/inverse_only.tgz').extractall('"+REMOTE+"'); pathlib.Path('"+REMOTE+"/vendor').mkdir(exist_ok=True); shutil.copy2('D:/AI/MAX_STRUCT_R10_2D_20261107/vendor/threadpoolctl.py','"+REMOTE+"/vendor/threadpoolctl.py')"
  print(ssh('D:/python/python.exe -c "'+code+'"','extract'));return
 if args.mode in ('replay','lagged'):
  study=REMOTE+'/research/research_max_struct_r11';kind=args.mode
  command='D:/python/python.exe '+study+'/experiments/twod_queue.py --execute --receipt-dir '+study+'/delegated/twod/dev_'+kind+'_queue --vendor-path '+REMOTE+'/vendor --manifest '+study+'/configs/twod_dev_'+kind+'_tasks.json --config '+study+'/configs/twod_dev_solver.json --cpus '+('8' if kind=='replay' else '8 10')
  s=ssh(command,'run_'+kind);(OUT/(kind+'_stdout.log')).write_text(s);print(s[-10000:]);return
 if args.mode=='fetch':
  code="import tarfile,pathlib; r=pathlib.Path('"+REMOTE+"'); t=tarfile.open(str(r/'r11_dev_outputs.tgz'),'w:gz'); [t.add(r/p,arcname=p) for p in ['research/research_max_struct_r11/RESULTS','research/research_max_struct_r11/delegated/twod'] if (r/p).exists()]; t.close()"
  ssh('D:/python/python.exe -c "'+code+'"','package_outputs')
  run(['scp','-i',KEY,HOST+':'+REMOTE+'/r11_dev_outputs.tgz',str(OUT/'r11_dev_outputs.tgz')],'fetch_outputs')
  with tarfile.open(OUT/'r11_dev_outputs.tgz') as tar:tar.extractall(REPO)
  print('FETCHED')
if __name__=='__main__':main()
