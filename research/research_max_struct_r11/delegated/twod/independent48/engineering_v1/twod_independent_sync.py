"""Sync public queue metadata/new terminal task files; sealed truth only after all288 terminal."""
import argparse,json,tarfile,time,subprocess,sys,resource
from pathlib import Path
import twod_dev_transport as t
D=t.ROOT/'delegated/twod/independent48'
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--finished',action='store_true');ap.add_argument('--include-sealed-truth',action='store_true');a=ap.parse_args();D.mkdir(parents=True,exist_ok=True)
 code="import pathlib,json; r=pathlib.Path('"+t.REMOTE+"'); qs=[r/'research/research_max_struct_r11/delegated/twod/independent48/queue/RECEIPT.json',r/'research/research_max_struct_r11/delegated/twod/independent48/interleaved_queue/RECEIPT.json']; q=[json.loads(p.read_text()) for p in qs if p.exists()]; print(json.dumps(dict(queues=[dict(status=x['status'],pending=x.get('pending'),running=x.get('running'),jobs=x['jobs']) for x in q])))"
 snapshot=json.loads(t.ssh('D:/python/python.exe -c "'+code+'"','independent48_sync_inventory'))
 jobs=[j for q in snapshot['queues'] for j in q['jobs']];ids=[j['task_id'] for j in jobs]
 if len(ids)!=len(set(ids)):raise ValueError('Duplicate jobs; audit before sync')
 all_terminal=len(jobs)==288 and all(not q['running'] for q in snapshot['queues']) and snapshot['queues'][-1]['status'] in ('COMPLETED','COMPLETED_WITH_ERRORS')
 if a.include_sealed_truth and not all_terminal:raise RuntimeError('Sealed truth byte copy waits for all288 terminal')
 held_for_root_review=bool(a.include_sealed_truth and (not (D/'ROOT_POSTSEAL_SCORER_APPROVED.json').exists() or not (D/'secondary/public_seal/RECEIPT.json').exists()))
 if held_for_root_review:a.include_sealed_truth=False  # Root specifically holds even the private package pending review.
 items=['research/research_max_struct_r11/delegated/twod/independent48','research/research_max_struct_r11/configs/twod_independent48_interleaved_remaining_tasks.json']
 new=[]
 if a.finished:
  for j in jobs:
   if j['status']=='LAUNCH_FAILED':continue
   p=j['output_dir']
   if not (t.REPO/p).exists():new.append(p)
  items+=new
 if a.include_sealed_truth:items.append('research/research_max_struct_r11/sealed_truth/independent48')
 tag=str(len(jobs))+'_'+str(time.time_ns());archive='independent48_sync_'+tag+'.tgz'
 code="import pathlib,tarfile,time,json,psutil; psutil.Process().cpu_affinity([8,10]); c=time.process_time(); r=pathlib.Path('"+t.REMOTE+"'); f=tarfile.open(str(r/'"+archive+"'),'w:gz'); [f.add(r/p,arcname=p) for p in "+repr(items)+" if (r/p).exists()]; f.close(); print(json.dumps(dict(packaging_CPU_seconds=time.process_time()-c,bytes=(r/'"+archive+"').stat().st_size)))"
 packed=json.loads(t.ssh('D:/python/python.exe -c "'+code+'"','independent48_package_sync'))
 target=t.OUT/archive;t.run(['scp','-i',t.KEY,t.HOST+':'+t.REMOTE+'/'+archive,str(target)],'independent48_fetch_sync');cpu=time.process_time()
 with tarfile.open(target) as f:f.extractall(t.REPO,filter='data')
 extraction_cpu=time.process_time()-cpu
 if all_terminal and (a.finished or (D/'secondary/public_seal/RECEIPT.json').exists()):
  public_seal=D/'secondary/public_seal'
  if not (public_seal/'RECEIPT.json').exists():
   science_python=json.loads((t.ROOT/'configs/twod_local_runtime.json').read_text())['python']  # Preserve venv symlink entry, never resolve().
   before=resource.getrusage(resource.RUSAGE_CHILDREN);control_cpu=time.process_time();seal_wall=time.perf_counter()
   process=subprocess.run([science_python,str(t.ROOT/'experiments/twod_secondary_seal.py'),'--manifest',str(t.ROOT/'configs/twod_independent48_all_tasks.json'),'--queue',str(D/'queue/RECEIPT.json'),'--queue',str(D/'interleaved_queue/RECEIPT.json')],capture_output=True,text=True)
   after=resource.getrusage(resource.RUSAGE_CHILDREN)
   (D/'secondary/SEAL_PROCESS_COST.json').write_text(json.dumps(dict(returncode=process.returncode,whole_child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,controller_CPU_seconds=time.process_time()-control_cpu,wall_seconds=time.perf_counter()-seal_wall,scope='Disjoint local wait4/RUSAGE_CHILDREN delta for one NumPy public-only seal child; child internal timer is nested',truth_access='NONE',stdout=process.stdout,stderr=process.stderr),indent=2)+'\n')
   if process.returncode:raise RuntimeError('Secondary public seal failed; GT stays closed: '+process.stderr)
   print(process.stdout,flush=True)
  elif json.loads((public_seal/'RECEIPT.json').read_text()).get('status')!='ALL96_PUBLIC_DEPLOYMENTS_SEALED_BEFORE48_GT':raise RuntimeError('Invalid partial secondary seal; retain and audit before GT')
 receipt=dict(status='SYNCED',terminal_tasks=len(jobs),all288_terminal=all_terminal,new_task_directories=new,remote_packaging_CPU_seconds=packed['packaging_CPU_seconds'],bytes=packed['bytes'],local_extraction_CPU_seconds=extraction_cpu,sealed_truth_access='SEALED_BYTE_COPY_NO_NUMERIC_GT_ACCESS' if a.include_sealed_truth else 'NONE',GT_scoring='NOT_RUN',sealed_truth_byte_copy_held_for_root_review=held_for_root_review,archive=str(target))
 (D/('SYNC_'+tag+'.json')).write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
