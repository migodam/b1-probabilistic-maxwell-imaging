"""Fetch only newly sealed exposed DEV task directories, never GT or live arrays."""
import json,tarfile,time
from pathlib import Path
import twod_dev_transport as t

def main():
 code="import pathlib,json; r=pathlib.Path('"+t.REMOTE+"'); qs=[r/'research/research_max_struct_r11/delegated/twod/dev_replay_queue/RECEIPT.json',r/'research/research_max_struct_r11/delegated/twod/dev_lagged_queue/RECEIPT.json']; print(json.dumps([j['output_dir'] for q in qs if q.exists() for j in json.loads(q.read_text())['jobs'] if j['status']=='FINISHED']))"
 paths=json.loads(t.ssh('D:/python/python.exe -c "'+code+'"','completed_inventory'))
 new=[p for p in paths if not (t.REPO/p).exists()]
 if not new:print('NO_NEW_SEALED_TASKS');return
 archive='completed_new_'+str(len(paths))+'.tgz'
 code="import tarfile,pathlib,time,json; c=time.process_time(); r=pathlib.Path('"+t.REMOTE+"'); t=tarfile.open(str(r/'"+archive+"'),'w:gz'); [t.add(r/p,arcname=p) for p in "+repr(new)+"]; t.close(); print(json.dumps(dict(packaging_CPU_seconds=time.process_time()-c)))"
 s=t.ssh('D:/python/python.exe -c "'+code+'"','package_new_completed')
 (t.OUT/('package_'+archive+'.json')).write_text(s)
 t.run(['scp','-i',t.KEY,t.HOST+':'+t.REMOTE+'/'+archive,str(t.OUT/archive)],'fetch_new_completed')
 started=time.process_time()
 with tarfile.open(t.OUT/archive) as a:a.extractall(t.REPO,filter='data')
 (t.OUT/('extract_'+archive+'.json')).write_text(json.dumps(dict(CPU_seconds=time.process_time()-started,new_tasks=new),indent=2)+'\n')
 print(json.dumps(dict(synced=new,all_completed=len(paths))))
if __name__=='__main__':main()
