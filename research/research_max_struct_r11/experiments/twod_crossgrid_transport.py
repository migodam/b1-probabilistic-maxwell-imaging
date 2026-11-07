"""Public-only supplementary64 controller, gated on primary288 completion.
No private transport or scoring; frozen inputs and production solver unchanged.
"""
import argparse,json,tarfile,time,subprocess,resource,os
from pathlib import Path
import twod_dev_transport as t
from twod_postprocess_cost import attach_external_cpu
D=t.ROOT/'delegated/twod/crossgrid8_paired';P=t.ROOT/'delegated/diagnostics/crossgrid8_protocol'
R=t.REMOTE+'/research/research_max_struct_r11';Q=R+'/delegated/twod/crossgrid8_paired/queue'
def dump(name,x):
 D.mkdir(parents=True,exist_ok=True);(D/name).write_text(json.dumps(x,indent=2)+'\n')
def remote(code,label):return t.ssh('D:/python/python.exe -c "'+code+'"',label)
def prepare():
 archive=P/'PUBLIC_INVERSE_INPUTS_V2.tgz'
 with tarfile.open(archive) as f:
  names=f.getnames()
  if any('sealed_truth' in n or 'data_private' in n or Path(n).is_absolute() or '..' in Path(n).parts for n in names):raise ValueError('Unsafe/nonpublic package')
  if len(names)!=19:raise ValueError('Expected frozen17 NPZ plus2 contracts')
 t.run(['scp','-i',t.KEY,str(archive),t.HOST+':'+t.REMOTE+'/crossgrid8_public_v2.tgz'],'crossgrid64_upload_public')
 code="import tarfile,time,json; c=time.process_time(); f=tarfile.open('"+t.REMOTE+"/crossgrid8_public_v2.tgz'); f.extractall('"+t.REMOTE+"'); print(json.dumps(dict(CPU_seconds=time.process_time()-c)))"
 cost=json.loads(remote(code,'crossgrid64_extract_public'))
 dump('PREPARED.json',dict(status='PUBLIC_INPUTS_UPLOADED_NO_INVERSION',UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),bytes=archive.stat().st_size,files=names,remote_extraction_CPU_seconds=cost['CPU_seconds'],truth_access='NONE',start_gate='All288 designated science terminal, primary public seals retained, original workers naturally exited'))
def gate():
 primary=t.ROOT/'delegated/twod/independent48'
 status=json.loads((t.ROOT/'delegated/twod/dev_transport/INDEPENDENT48_STATUS.json').read_text())
 seals=[primary/p for p in ('secondary/public_seal/RECEIPT.json','secondary/global_appendix/public_seal/RECEIPT.json','actual_cpu_fairness/public_seal/RECEIPT.json')]
 if status['completed']!=288 or status['running'] or status['status'] not in ('COMPLETED','COMPLETED_WITH_ERRORS') or not all(p.exists() for p in seals):return False
 syncs=[json.loads(p.read_text()) for p in primary.glob('SYNC_288_*.json')]
 if not any(s.get('all288_terminal') for s in syncs):return False
 code="import psutil,json; print(json.dumps([dict(pid=p.pid,args=p.cmdline()) for p in psutil.process_iter() if 'python' in p.name().lower() and any(a.endswith(('twod_run.py','twod_queue.py')) for a in p.cmdline())]))"
 live=json.loads(remote(code,'crossgrid64_primary_worker_exit_gate'));dump('PRIMARY_RELEASE_PROBE.json',dict(active_original_workers=live,UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),truth_access='NONE'))
 return not live
def execute():
 if not (D/'PREPARED.json').exists():raise RuntimeError('Prepare public archive first')
 dump('ACTIVE_CONTROLLER.json',dict(pid=os.getpid(),UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),mode='WAIT_FOR_PRIMARY288_THEN_RUN64_NO_GT'))
 while not gate():time.sleep(60)
 # The scientific queue is finite and refuses existing output folders.
 if (D/'DISPATCH_STARTED.json').exists():raise RuntimeError('Already dispatched; inspect existing queue, never duplicate')
 dump('DISPATCH_STARTED.json',dict(UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),tasks=64,materials=8,cpus=[0,2,8,10],truth_access='NONE',primary288_gate='PASSED'))
 command='D:/python/python.exe '+R+'/experiments/twod_queue.py --execute --receipt-dir '+Q+' --vendor-path '+t.REMOTE+'/vendor --manifest '+R+'/delegated/diagnostics/crossgrid8_protocol/TASKS_V2.json --config '+R+'/delegated/diagnostics/crossgrid8_protocol/SOLVER_CONFIG_V2.json --cpus 0 2 8 10 > '+R+'/delegated/twod/crossgrid8_paired/controller_stdout.log 2> '+R+'/delegated/twod/crossgrid8_paired/controller_stderr.log'
 remote("import pathlib; pathlib.Path('"+R+"/delegated/twod/crossgrid8_paired').mkdir(parents=True,exist_ok=True)",'crossgrid64_controller_mkdir')
 t.ssh(command,'crossgrid64_run_finite_queue')
 inventory=json.loads(remote("import json,pathlib; print(pathlib.Path('"+Q+"/RECEIPT.json').read_text())",'crossgrid64_terminal_inventory'))
 if inventory['status'] not in ('COMPLETED','COMPLETED_WITH_ERRORS') or inventory['running'] or len(inventory['jobs'])!=64:raise RuntimeError('Nonterminal64; retain and inspect')
 items=['research/research_max_struct_r11/delegated/twod/crossgrid8_paired','research/research_max_struct_r11/RESULTS/crossgrid8_paired']
 code="import pathlib,tarfile,time,json; c=time.process_time(); root=pathlib.Path('"+t.REMOTE+"'); f=tarfile.open(str(root/'crossgrid64_public_outputs.tgz'),'w:gz'); [f.add(root/p,arcname=p) for p in "+repr(items)+"]; f.close(); print(json.dumps(dict(CPU_seconds=time.process_time()-c)))"
 packed=json.loads(remote(code,'crossgrid64_package_public_outputs'))
 archive=D/'PUBLIC_OUTPUTS.tgz';t.run(['scp','-i',t.KEY,t.HOST+':'+t.REMOTE+'/crossgrid64_public_outputs.tgz',str(archive)],'crossgrid64_fetch_public_outputs')
 c=time.process_time()
 with tarfile.open(archive) as f:
  if any('sealed_truth' in p for p in f.getnames()):raise RuntimeError('Private archive member')
  f.extractall(t.REPO,filter='data')
 extract=time.process_time()-c
 runtime=json.loads((t.ROOT/'configs/twod_local_runtime.json').read_text())['python'];cost=D/'PUBLIC_SEAL_PROCESS_COST.json';out=D/'public_seal';before=resource.getrusage(resource.RUSAGE_CHILDREN);c=time.process_time();wall=time.perf_counter()
 process=subprocess.run([runtime,str(P/'seal_paired_grid8.py'),'--queue',str(D/'queue/RECEIPT.json'),'--output',str(out)],capture_output=True,text=True);after=resource.getrusage(resource.RUSAGE_CHILDREN)
 dump(cost.name,dict(returncode=process.returncode,whole_child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,controller_CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-wall,stdout=process.stdout,stderr=process.stderr,truth_access='NONE'))
 if process.returncode:raise RuntimeError('Public seal failed: '+process.stderr)
 attach_external_cpu(out,cost)
 dump('COSTS.json',dict(status='ALL64_TERMINAL_PUBLIC_SEALED_GT_HELD',scientific_terminal=64,scientific_failed=sum(j['status']!='FINISHED' for j in inventory['jobs']),independent_materials=8,child_GetProcessTimes_CPU_seconds=sum(j['child_process_CPU_seconds'] for j in inventory['jobs']),queue_controller_CPU_seconds=inventory['controller_CPU_seconds'],remote_packaging_CPU_seconds=packed['CPU_seconds'],local_extraction_CPU_seconds=extract,source_transport_ledger=str(D/'TRANSPORT_COSTS.jsonl'),public_seal_external_cost=str(cost),GT_scoring='NOT_RUN',truth_access='NONE'))
 print('ALL64_TERMINAL_PUBLIC_SEALED_GT_HELD',flush=True)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['prepare','follow']);a=ap.parse_args();D.mkdir(parents=True,exist_ok=True);t.OUT=D
 if a.mode=='prepare':prepare()
 else:execute()
if __name__=='__main__':main()
