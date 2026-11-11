"""Finite four-P-core execution; same frozen paired tasks, no GT or retries."""
import os,sys,time,json,argparse,subprocess,importlib.util
from pathlib import Path
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
def main():
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--legacy-repo',type=Path,required=True);a=p.parse_args()
 import psutil
 root=a.root.resolve();legacy=a.legacy_repo.resolve();q=root/'queues/main';q.mkdir(parents=True,exist_ok=False);(q/'logs').mkdir()
 spec=importlib.util.spec_from_file_location('existing_external_cpu',legacy/'research/research_max_struct_r11/experiments/twod_queue.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
 tasks=json.loads((root/'configs/TASKS.json').read_text())['tasks'];cores=[0,2,4,6]
 probe=dict(selected_cores=cores,CPU_percent=psutil.cpu_percent(interval=1,percpu=True),other_python=[dict(pid=p.pid,command=p.cmdline(),affinity=p.cpu_affinity()) for p in psutil.process_iter() if 'python' in p.name().lower() and p.pid!=os.getpid()])
 (q/'RESOURCE_PROBE.json').write_text(json.dumps(probe,indent=2))
 env=os.environ.copy();env['PYTHONPATH']=os.pathsep.join([str(legacy),str(legacy/'vendor')]);pending=list(tasks);running={};jobs=[];start=time.perf_counter();events=(q/'EVENTS.jsonl').open('x')
 def save(status):
  d=dict(status=status,total=len(tasks),pending=len(pending),running=[dict(task_id=x['task']['task_id'],pid=p.pid,cpu=x['cpu']) for p,x in running.items()],jobs=jobs,wall_seconds=time.perf_counter()-start,controller_CPU_seconds=time.process_time(),finished_child_CPU_seconds=sum(x['external'].get('child_process_CPU_seconds') or 0 for x in jobs))
  temp=q/'RECEIPT.tmp';temp.write_text(json.dumps(d,indent=2));temp.replace(q/'RECEIPT.json')
 while pending or running:
  for core in cores:
   if not pending:break
   if core in [x['cpu'] for x in running.values()]:continue
   task=pending.pop(0);argv=[str(root/'src/run_paired.py'),'--root',str(root),'--legacy-repo',str(legacy),'--task-id',task['task_id'],'--cpu',str(core)]
   boot='import os,sys,runpy,psutil;psutil.Process().cpu_affinity(['+str(core)+']);psutil.Process().nice(psutil.NORMAL_PRIORITY_CLASS);sys.argv='+repr(argv)+';runpy.run_path('+repr(argv[0])+',run_name="__main__")'
   command=[sys.executable,'-c',boot];log=(q/'logs'/(task['task_id']+'.log')).open('x');p=subprocess.Popen(command,env=env,cwd=root,stdout=log,stderr=subprocess.STDOUT)
   running[p]=dict(task=task,cpu=core,log=log,start=time.perf_counter(),command=command);events.write(json.dumps(dict(event='START',task_id=task['task_id'],pid=p.pid,cpu=core))+'\n');events.flush()
  for p,x in list(running.items()):
   code=p.poll()
   if code is None:continue
   x['log'].close();receipt=root/'RESULTS'/x['task']['task_id']/'RECEIPT.json';child=json.loads(receipt.read_text()) if receipt.exists() else None
   row=dict(task_id=x['task']['task_id'],case_id=x['task']['case_id'],method=x['task']['method'],cpu=x['cpu'],returncode=code,external=old.external_cpu(p),wall_seconds=time.perf_counter()-x['start'],command=x['command'],child_status=child['status'] if child else None)
   jobs.append(row);del running[p];events.write(json.dumps(dict(event='END',**row))+'\n');events.flush();print(json.dumps({k:row[k] for k in ('task_id','returncode','external','child_status')}),flush=True)
  save('RUNNING')
  if running:time.sleep(1)
 save('COMPLETED' if all(j['returncode']==0 for j in jobs) else 'COMPLETED_WITH_ERRORS');events.close()
if __name__=='__main__':main()
