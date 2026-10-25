"""Delivery only: four bounded-memory streams, closed markers, name/size audit."""
import argparse
import concurrent.futures as cf
import fcntl
import json
import os
from pathlib import Path
import subprocess
import tempfile
import threading
import uuid
import time

REPO = 'migodam/b1-probabilistic-maxwell-imaging'
MAX_METADATA = 16 * 1024 * 1024


def write(path, obj):
    tmp = path.with_name(path.name + '.pending')
    with tmp.open('w') as f:
        json.dump(obj, f, indent=2)
        f.write('\n'); f.flush(); os.fsync(f.fileno())
    os.replace(tmp, path)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bundle', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--tag', default='research-p3lr4-20261025')
    p.add_argument('--max-wall', type=float, default=43200)
    p.add_argument('--max-cost', type=float, default=86400, help='remaining aggregate seconds: coordinator plus upload streams')
    a = p.parse_args()
    if not 0 < a.max_wall <= 43200:
        p.error('--max-wall must be positive and at most 43200 seconds')
    if not 0 < a.max_cost <= 86400:
        p.error('--max-cost must be positive and at most 86400 seconds')
    a.bundle = a.bundle.resolve(); a.out = a.out.resolve()
    if a.out == a.bundle or a.bundle in a.out.parents:
        p.error('receipts must be outside immutable bundle')
    a.out.mkdir(parents=True, exist_ok=True)
    receipt_root=a.out
    run_id=time.strftime('%Y%m%dT%H%M%S')+'-'+uuid.uuid4().hex[:8]
    a.out=a.out/'runs'/run_id
    a.out.mkdir(parents=True)
    (a.out/'jobs').mkdir(); (a.out/'logs').mkdir()
    lock = (receipt_root/'RUN.lock').open('a+')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    start = time.monotonic(); stop = threading.Event()
    deadline = start + a.max_wall
    budget_lock=threading.Lock(); active={}; stream_spent=0.0
    prior={}; previous_states=[]
    for path in (receipt_root/'PROGRESS.json',receipt_root/'UPLOAD_RECEIPT.json'):
        if path.exists():
            old=json.loads(path.read_text()); previous_states.append(old)
            for key in ('cumulative_wall_seconds','cumulative_aggregate_cost_seconds'):
                prior[key]=max(prior.get(key,0.0),old.get(key,0.0))
    if previous_states:
        latest=max(previous_states,key=lambda x:x.get('checkpoint_epoch',0))
        if latest.get('status')=='UPLOADING':
            gap=max(0,time.time()-latest.get('checkpoint_epoch',start))
            prior['cumulative_wall_seconds']=prior.get('cumulative_wall_seconds',0)+gap
            prior['cumulative_aggregate_cost_seconds']=prior.get('cumulative_aggregate_cost_seconds',0)+5*gap
    # --max-wall / --max-cost are invocation allowances. Previously recorded usage
    # is also subtracted from the absolute 12h / 24h ceilings on restart.
    wall_allowance=min(a.max_wall,43200-prior.get('cumulative_wall_seconds',0))
    cost_allowance=min(a.max_cost,86400-prior.get('cumulative_aggregate_cost_seconds',0))
    deadline=start+wall_allowance
    def cost():
        with budget_lock:
            now=time.monotonic()
            return now-start+stream_spent+sum(now-t for t in active.values())
    def guard():
        if time.monotonic()>=deadline-5 or cost()>=cost_allowance-15:
            stop.set(); raise TimeoutError('delivery wall or aggregate cost cap (early safety margin)')
    jobs = {}; assets = {}; claimed = {}; futures = {}
    status = 'ERROR'; error = None; release_id = None
    def elapsed(): return time.monotonic() - start
    def remaining():
        guard()
        value = deadline-time.monotonic()
        if value <= 0: raise TimeoutError('delivery wall cap')
        return value
    def api(path, paginate=False):
        cmd = ['gh', 'api', path]
        if paginate: cmd += ['--paginate', '--slurp']
        with tempfile.TemporaryFile() as out, tempfile.TemporaryFile() as err:
            r = subprocess.run(cmd, stdout=out, stderr=err, timeout=min(90, remaining(), max(.1,(cost_allowance-cost()-15)/5)))
            out.seek(0, 2); size = out.tell()
            err.seek(0, 2); n = err.tell(); err.seek(max(0, n-4000)); msg=err.read().decode(errors='replace')
            if r.returncode: raise RuntimeError('gh metadata read failed: '+msg)
            if size > MAX_METADATA: raise RuntimeError('remote metadata cap exceeded')
            out.seek(0); return json.load(out)
    def remote():
        pages=api('repos/'+REPO+'/releases/'+str(release_id)+'/assets?per_page=100', True)
        found={}
        for page in pages:
            for v in page:
                if v['name'] in found: raise ValueError('duplicate remote asset name: '+v['name'])
                found[v['name']]={'bytes':v['size'], 'state':v['state'], 'id':v['id']}
        return found
    def check(found, name, size):
        v=found.get(name)
        if v is None: return False
        if v['bytes'] != size: raise ValueError('remote size mismatch; no clobber: '+name)
        if v['state'] != 'uploaded': raise ValueError('remote asset state '+v['state']+'; no clobber: '+name)
        return True
    def upload(name, size):
        began=time.monotonic(); worker=threading.current_thread().name
        row={'file':name,'bytes':size,'worker':worker,'attempts':[], 'status':'ERROR'}
        try:
            if check(remote(),name,size): row['status']='ALREADY_UPLOADED_NAME_SIZE_MATCH'; return row
            for attempt in range(1,4):
                if stop.is_set(): raise RuntimeError('coordinator stopped after another failure')
                remaining()
                log=a.out/'logs'/(name+'.attempt'+str(attempt)+'.log')
                t=time.monotonic(); code=None; uncertainty=None
                with log.open('ab') as f:
                    try:
                        guard()
                        child=subprocess.Popen(['gh','release','upload',a.tag,str(a.bundle/name),'--repo',REPO], stdout=f, stderr=f)
                        with budget_lock: active[worker]=t
                        attempt_deadline=min(deadline-5,t+1200)
                        try:
                            while child.poll() is None:
                                guard()
                                if time.monotonic()>=attempt_deadline: raise TimeoutError('1200s attempt limit')
                                try: child.wait(timeout=.2)
                                except subprocess.TimeoutExpired: pass
                            code=child.returncode
                        except BaseException as exc:
                            child.terminate()
                            try: child.wait(timeout=2)
                            except subprocess.TimeoutExpired: child.kill(); child.wait()
                            uncertainty=repr(exc)+'; child terminated and reaped'
                        finally:
                            with budget_lock:
                                nonlocal stream_spent
                                stream_spent += time.monotonic()-active.pop(worker)
                    except subprocess.TimeoutExpired:
                        uncertainty='upload timeout; child terminated and reaped'
                    f.flush(); os.fsync(f.fileno())
                record={'attempt':attempt,'returncode':code,'uncertainty':uncertainty,'elapsed_seconds':time.monotonic()-t,'log':str(log)}
                row['attempts'].append(record)
                # Verification must succeed before any retry. Unknown remote state fails closed.
                found=remote()
                if check(found,name,size):
                    row['status']='UPLOADED_NAME_SIZE_VERIFIED' if code==0 else 'UPLOADED_VERIFIED_AFTER_UNCERTAIN_RETURN'
                    return row
                if code==0: raise RuntimeError('upload returned success but remote asset absent: '+name)
                if attempt<3:
                    if stop.wait(min(30,remaining())): raise RuntimeError('coordinator stopped')
            raise RuntimeError('upload failed after 3 attempts: '+name)
        except BaseException as exc:
            row['error']=repr(exc); raise
        finally:
            row['elapsed_seconds']=time.monotonic()-began
            row['upload_stream_seconds']=sum(v['elapsed_seconds'] for v in row['attempts'])
            write(a.out/'jobs'/(name+'.json'),row)
    def summary(state):
        workers={}
        for job in jobs.values():
            worker=workers.setdefault(job['worker'],{'job_elapsed_seconds':0.0,'upload_stream_seconds':0.0,'jobs':0,'attempts':0})
            worker['job_elapsed_seconds']+=job['elapsed_seconds']
            worker['upload_stream_seconds']+=job['upload_stream_seconds']
            worker['jobs']+=1; worker['attempts']+=len(job['attempts'])
        return {'status':state,'tag':a.tag,'repo':REPO,'bundle':str(a.bundle),'max_parallel_upload_streams':4,'checkpoint_epoch':time.time(),
                'wall_seconds':elapsed(),'aggregate_cost_seconds':cost(),
                'cumulative_wall_seconds':prior.get('cumulative_wall_seconds',0)+elapsed(),
                'cumulative_aggregate_cost_seconds':prior.get('cumulative_aggregate_cost_seconds',0)+cost(),
                'wall_allowance_seconds':wall_allowance,'cost_allowance_seconds':cost_allowance,'run_id':run_id,
                'sum_job_elapsed_seconds':sum(v['elapsed_seconds'] for v in jobs.values()),
                'sum_upload_stream_seconds':sum(v['upload_stream_seconds'] for v in jobs.values()),
                'workers':workers,'assets':assets,'claimed':claimed,'jobs':jobs,'error':error,'hashes_used':False,'no_clobber':True,
                'RHS':0,'training':0,'scope':'delivery only; does not alter scientific qualification'}
    pool=cf.ThreadPoolExecutor(max_workers=4,thread_name_prefix='delivery')
    try:
        release_pages=api('repos/'+REPO+'/releases?per_page=100', True)
        matches=[v for page in release_pages for v in page if v.get('tag_name')==a.tag]
        if len(matches)!=1: raise ValueError('expected one authenticated release for exact tag; got '+str(len(matches)))
        release=matches[0]
        if release.get('draft') is not True: raise ValueError('existing Release is not DRAFT')
        release_id=release['id']
        while True:
            remaining()
            for future in list(futures):
                if future.done():
                    name=futures.pop(future)
                    try: row=future.result()
                    except BaseException:
                        row=json.loads((a.out/'jobs'/(name+'.json')).read_text()); jobs[name]=row
                        raise
                    jobs[name]=row; assets[name]={'bytes':row['bytes'],'status':row['status']}
            markers={}
            for marker in sorted((a.bundle/'closed_parts').glob('*.json')):
                v=json.loads(marker.read_text()); name=v['file']; size=v['bytes']
                if not isinstance(name,str) or Path(name).name!=name or '.writing' in name or '#' in name or '?' in name:
                    raise ValueError('unsafe closed part name')
                if type(size) is not int or size<0: raise ValueError('invalid closed part bytes')
                part=a.bundle/name
                if part.is_symlink() or not part.is_file() or part.stat().st_size != size:
                    raise ValueError('closed local part size/type mismatch: '+name)
                if name in markers and markers[name]!=size: raise ValueError('conflicting marker')
                markers[name]=size
                if name in claimed and claimed[name]!=size: raise ValueError('closed marker changed: '+name)
            receipt=a.bundle/'RECEIPT.json'; expected=None
            if receipt.exists():
                done=json.loads(receipt.read_text())
                if done['status'] != 'COMPLETE': raise ValueError('packager closed unsuccessfully; preserve partial assets')
                expected={v['file']:v['bytes'] for v in done['parts']}
                if len(expected)!=len(done['parts']) or expected!=markers:
                    raise ValueError('COMPLETE pack receipt and closed marker map differ')
            for name,size in markers.items():
                if len(futures)>=4: break
                if name not in claimed:
                    claimed[name]=size; futures[pool.submit(upload,name,size)]=name
            write(a.out/'PROGRESS.json',summary('UPLOADING'))
            write(receipt_root/'PROGRESS.json',summary('UPLOADING'))
            if expected is not None and set(assets)==set(expected) and not futures:
                final_release=api('repos/'+REPO+'/releases/'+str(release_id))
                if not final_release.get('draft') or final_release['id']!=release_id or final_release.get('tag_name')!=a.tag:
                    raise ValueError('Release identity/draft state changed')
                audit=remote()
                for name,size in expected.items():
                    if not check(audit,name,size): raise ValueError('final remote asset missing: '+name)
                write(a.out/'REMOTE_ASSET_AUDIT.json',{'status':'NAME_SIZE_STATE_VERIFIED','expected':expected,'remote_assets':audit,'wall_seconds':elapsed()})
                status='ALL_PARTS_UPLOADED_NAME_SIZE_VERIFIED'; break
            stop.wait(min(.5,remaining()))
    except BaseException as exc:
        error=repr(exc); stop.set()
    finally:
        # Active children finish or hit the global deadline; no additional uploads are submitted.
        pool.shutdown(wait=True, cancel_futures=True)
        for name in claimed:
            jobpath=a.out/'jobs'/(name+'.json')
            if jobpath.exists(): jobs[name]=json.loads(jobpath.read_text())
        for name,row in jobs.items():
            if row['status']!='ERROR': assets[name]={'bytes':row['bytes'],'status':row['status']}
        write(a.out/'UPLOAD_RECEIPT.json',summary(status))
        write(receipt_root/'UPLOAD_RECEIPT.json',summary(status))
    if error: raise SystemExit(error)

if __name__=='__main__': main()
