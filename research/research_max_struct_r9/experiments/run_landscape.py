"""Execute one public R9 continuation task; no truth or scoring access."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[k]='1'
import argparse,json,sys,time,traceback
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))

def main():
    wall,cpu=time.perf_counter(),time.process_time()
    ap=argparse.ArgumentParser();ap.add_argument('--manifest',type=Path,required=True);ap.add_argument('--task-id',required=True)
    args=ap.parse_args(); manifest=json.loads(args.manifest.read_text());task=next(t for t in manifest['tasks'] if t['task_id']==args.task_id)
    out=Path(task['output_dir']);out=out if out.is_absolute() else REPO/out;out.mkdir(parents=True,exist_ok=False)
    receipt=dict(status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),task=task,
        config=manifest['landscape_config'],git_base='efaad956f',truth_access='NONE',physical_convention='exp(-i omega t), outgoing Hankel1, real eps_r-1, fixed mu',
        prior='Inherited identical smooth+TV',material_class='576 independent real isotropic cells',noise_scope='public noise; original normalized Phi is not a new whitened likelihood')
    def save(): (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    save()
    try:
        import numpy as np
        from threadpoolctl import threadpool_limits
        from research.research_max_struct_r8.experiments.run_case import load_public
        from research.research_fbg_r3_experiments.src.physics import MaxwellModel
        from research.research_max_struct_r9.src.landscape_solver import LandscapeRun
        with threadpool_limits(limits=1):
            public=load_public(task['public_files']);model=MaxwellModel(public);Y=np.asarray(public['Y'],complex)
            m0=np.full(model.n,float(task.get('initial_material',.1)))
            if task.get('initialization_file'):
                with np.load(task['initialization_file']) as d:m0=np.asarray(d['m_init'],float)
            run=LandscapeRun(model,Y,out)
            receipt.update(pre_solver_CPU_seconds=run.cpu-cpu,
                           pre_solver_wall_seconds=run.wall-wall)
            result=run.run(m0,task['algorithm'],manifest['landscape_config'])
            for name in ('history','budget'):
                (out/(name.upper()+'.json')).write_text(json.dumps(result.pop(name),indent=2)+'\n')
            receipt.update(result)
    except Exception as exc:
        receipt.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc());raise
    finally:
        receipt.update(process_CPU_seconds=time.process_time()-cpu,process_wall_seconds=time.perf_counter()-wall,
                       ended_utc=datetime.now(timezone.utc).isoformat());save()

if __name__=='__main__':main()
