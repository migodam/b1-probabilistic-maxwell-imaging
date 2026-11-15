"""R7 dispatch: inherited algorithms, native FFT actions, observed paid budgets."""
import time
import numpy as np
from pathlib import Path
from research.research_max_struct_r6.src import benchmark_solver as r6
from research.research_max_struct_r8.src.cie_fft import run_cie_fft

ARMS=('cie_single','cie_dual','global_dual','a_material','a_plus','a_plus_sourcewise','physical_lbfgs')
SOURCE_PATHS=dict(r6.SOURCE_PATHS,a_plus='R7 Mobius/MM or projected A readout; frozen config names method',physical_lbfgs='Original Maxwell L-BFGS-B from same initial material',FFT='R7 native Fourier action over R6 zero-padded disk Green')
difference=r6.difference
refine_observed=r6.refine_observed

def arm_config(arm,common,Y,E,continuation_stages=None):
    base='a_material' if arm in ('a_plus','a_plus_sourcewise') else 'global_dual' if arm=='physical_lbfgs' else arm
    cfg=r6.arm_config(base,common,Y,E,continuation_stages)
    cfg['arm']='a_plus' if arm=='a_plus_sourcewise' else arm
    if arm=='a_plus_sourcewise':cfg['sourcewise_cg']=True
    return cfg

def run_auxiliary(model,Y,m0,arm,config):
    if arm=='physical_lbfgs':
        t,c=time.perf_counter(),time.process_time()
        before=dict(model.cost);state=model.state(m0)
        value=model.physical_score(m0,Y)
        return dict(status='DIRECT_PHYSICS_NO_AUXILIARY',history=[],
            candidates=[dict(m=m0.copy(),objective=value,stage=-1,reason='common_initial_material')],
            cost=dict(),execution_receipt=dict(status='DIRECT_PHYSICS_NO_AUXILIARY',wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,physical_model_delta=difference(dict(model.cost),before)))
    cfg=dict(config)
    if arm in ('a_plus','a_plus_sourcewise'):
        from research.research_max_struct_r7.src.aplus_readout import mobius_material_readout,projected_material_readout
        mode=cfg.get('aplus_mode','mobius')
        cfg['material_readout']=mobius_material_readout if mode=='mobius' else projected_material_readout
    return run_cie_fft(model,Y,m0,cfg)

class PaidStateObserver(r6.PaidStateObserver):
    """Store incumbent material every ten already paid physical evaluations.

    Float32 material snapshots are only offline budget grading outputs; all
    online optimization, original Phi and final predictions remain float64.
    """
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.budget_rows=[];self.budget_material=[];self.budget_calls=0
        self.budget_cpu_offset=0.

    def record(self,*args,**kwargs):
        row=super().record(*args,**kwargs)
        self.budget_calls+=1
        if self.budget_calls==1 or self.budget_calls%10==0 or self.phase=='auxiliary':
            self._budget_store()
        return row

    def _budget_store(self):
        self.budget_rows.append((self.budget_cpu_offset+time.process_time()-self.cpu,
            self.model.cost.get('LU',0),self.model.cost.get('actual_complex_RHS',0),self.best['objective']))
        self.budget_material.append(self.best['m'].astype(np.float32))

    def save_budget_states(self,out):
        self._budget_store()
        np.savez_compressed(Path(out)/'BUDGET_STATES.npz',
            columns=np.array(['CPU_seconds_including_process_prefix','LU','complex_RHS','original_Phi']),
            costs=np.array(self.budget_rows),m=np.array(self.budget_material))
