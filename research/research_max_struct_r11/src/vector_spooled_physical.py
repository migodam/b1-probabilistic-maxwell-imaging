"""Thin CPU-only storage subclass: historical J/Q views in one temporary map.

All forward/gradient/cache decisions remain in the qualified R10 objective.
The returned records keep the same arrays' shapes, dtype and numerical values.
"""
import os,time,tempfile,json
from pathlib import Path
import numpy as np
from research.research_max_struct_r10.src.vector_physical import PhysicalObjective


class SpooledPhysicalObjective(PhysicalObjective):
    def __init__(self,model,Y,*,upper,capacity,directory,owner,verify_first=False):
        if model.device.type!='cpu':raise ValueError('Spool optimization is CPU-only')
        wall,cpu=time.perf_counter(),time.process_time()
        super().__init__(model,Y,upper=upper)
        directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
        handle=tempfile.NamedTemporaryFile(prefix=owner+'_'+str(os.getpid())+'_',suffix='.dat',dir=directory,delete=False)
        self.spool_path=Path(handle.name);handle.close();self.capacity=int(capacity)
        self.spool=np.memmap(self.spool_path,mode='w+',dtype=np.complex128,shape=(self.capacity,2,model.M,self.P))
        self.views={};self.verify_first=verify_first;self.roundtrip_verified=False;self.archive_verified=False
        self.storage=dict(setup_CPU_seconds=time.process_time()-cpu,setup_wall_seconds=time.perf_counter()-wall,
            write_CPU_seconds=0.,write_wall_seconds=0.,flush_CPU_seconds=0.,flush_wall_seconds=0.,
            bytes_logically_written=0,record_writes=0,cache_hit_rewrites_skipped=0,
            capacity_states=self.capacity,capacity_bytes=int(self.spool.nbytes),single_mmap_file=True,
            shape=[self.capacity,2,model.M,self.P],dtype='complex128',new_physical_RHS=0)

    def evaluate(self,*args,**kwargs):
        record=super().evaluate(*args,**kwargs)
        wall,cpu=time.perf_counter(),time.process_time();slot=int(record['id'][5:])
        if not 0<=slot<self.capacity:raise RuntimeError('Frozen800request+initial/stage spool bound exceeded')
        views=self.views.get(slot)
        if views is None:
            views=(self.spool[slot,0],self.spool[slot,1]);self.views[slot]=views
        if record['J'] is views[0] and record['Q'] is views[1]:
            self.storage['cache_hit_rewrites_skipped']+=1
        else:
            for key in ('J','Q'):
                if record[key].shape!=(self.model.M,self.P) or record[key].dtype!=np.complex128:
                    raise ValueError('J/Q spool must preserve exact3N x60 complex128 arrays')
            views[0][...]=record['J'];views[1][...]=record['Q']
            if self.verify_first and not self.roundtrip_verified:
                if not np.array_equal(views[0],record['J']) or not np.array_equal(views[1],record['Q']):
                    raise RuntimeError('Normal first-state spool roundtrip differs')
                self.roundtrip_verified=True
            record['J'],record['Q']=views
            self.storage['bytes_logically_written']+=2*self.model.M*self.P*np.dtype(np.complex128).itemsize
            self.storage['record_writes']+=1
            fw,fc=time.perf_counter(),time.process_time();self.spool.flush()
            self.storage['flush_CPU_seconds']+=time.process_time()-fc;self.storage['flush_wall_seconds']+=time.perf_counter()-fw
        dc,dw=time.process_time()-cpu,time.perf_counter()-wall
        self.storage['write_CPU_seconds']+=dc;self.storage['write_wall_seconds']+=dw
        self.events[-1].update(storage_spool_CPU_seconds=dc,storage_spool_wall_seconds=dw,
            storage_scope='Post-physics exact J/Q storage; included in outer completion CPU/wall; do not add twice')
        return record

    def verify_initial_archive(self,path,record):
        if self.verify_first and not self.archive_verified:
            wall,cpu=time.perf_counter(),time.process_time()
            with np.load(path,allow_pickle=False) as z:
                if not np.array_equal(z['J'],record['J']) or not np.array_equal(z['Q'],record['Q']):
                    raise RuntimeError('Normal first-state archived J/Q differs from spool view')
            self.archive_verified=True
            self.storage['first_archive_check_CPU_seconds']=time.process_time()-cpu
            self.storage['first_archive_check_wall_seconds']=time.perf_counter()-wall
            (Path(path).parent/'MATERIAL_SPOOL_NORMAL_FIRST_ARCHIVE_CHECK.json').write_text(json.dumps(dict(
                status='EXACT_NORMAL_ARCHIVE_ROUNDTRIP',roundtrip_verified=self.roundtrip_verified,archive_verified=True,
                new_physical_RHS=0,scope='Already-paid normal initial state; exact J/Q view and saved-array check, not a new forward probe'),indent=2)+'\n')

    def finish_spool(self):
        """Called only after all original arrays, receipts and seals are saved."""
        wall,cpu=time.perf_counter(),time.process_time();self.spool.flush()
        allocated=int(self.spool_path.stat().st_blocks*512)
        for record in self.records.values():
            record.pop('J',None);record.pop('Q',None)
        self.pending_current.clear();self.records.clear();self.views.clear()
        self.spool._mmap.close();self.spool=None;self.spool_path.unlink()
        result=dict(self.storage,used_states=len({e['state_id'] for e in self.events}),
            allocated_temporary_disk_bytes_at_close=allocated,temporary_path=str(self.spool_path),temporary_file_deleted=True,
            close_CPU_seconds=time.process_time()-cpu,close_wall_seconds=time.perf_counter()-wall,
            first_roundtrip_verified=self.roundtrip_verified,first_archive_verified=self.archive_verified,
            semantics='No mathematical/precision/cache-hit/step/budget change; all required archived arrays preserved')
        return result
