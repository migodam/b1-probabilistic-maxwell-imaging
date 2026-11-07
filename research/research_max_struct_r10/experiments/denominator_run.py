"""Full historical DEV inversion with one changed material readout.

No truth inputs. Current/basis/schedule/physical corrector inherited verbatim.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,json,time,argparse,importlib.util
from types import SimpleNamespace
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--historical-repo',required=True);ap.add_argument('--case',required=True)
    ap.add_argument('--basis',choices=['twofold','green_equal_rank','green_strong'],required=True)
    ap.add_argument('--projection',choices=['galerkin','right'],default='galerkin')
    ap.add_argument('--retain-self',action='store_true')
    ap.add_argument('--public-files',nargs='+',help='Portable public observation/operator paths; no truth')
    ap.add_argument('--output-dir',help='Fresh output directory for reproduction')
    args=ap.parse_args();old=Path(args.historical_repo)
    from research.research_max_struct_r10.src.lowrank_denominator import prepare_bases,LowRankDenominator,matrix_material_readout
    from research.research_max_struct_r8.src import benchmark_solver as solver
    from threadpoolctl import threadpool_limits
    import numpy as np
    task0=json.loads((old/'research/research_max_struct_r9/RESULTS/development'/f'{args.case}__global_dual__fixed_001'/'PUBLIC_TASK.json').read_text())
    if args.public_files:task0['public_files']=args.public_files
    label=f'{args.case}__{args.basis}__{args.projection}'
    if args.retain_self:label+='__retain_self'
    out=Path(args.output_dir).absolute() if args.output_dir else ROOT/'RESULTS/denominator_dev'/label
    inp=out.parent/(out.name+'__inputs');inp.mkdir(parents=True,exist_ok=False)
    cfg=json.loads((REPO/'research/research_max_struct_r9/configs/inherited_solver.json').read_text())
    cfg['run_id']='r10_projected_green_denominator_dev';cfg['candidate']=dict(basis=args.basis,projection=args.projection,
        threshold=.1,inner_maxiter=8,retain_self=args.retain_self,frozen_before_complete_trajectory=True,truth_access='NONE',development_objects='historically exposed')
    (inp/'CONFIG.json').write_text(json.dumps(cfg,indent=2))
    task=dict(task0,task_id=label,arm='a_plus',config=str(inp/'CONFIG.json'),output_dir=str(out))
    (inp/'MANIFEST.json').write_text(json.dumps(dict(tasks=[task]),indent=2))
    def dispatch(model,Y,m0,arm,config):
        started=time.process_time();bases,info=prepare_bases(model.G,model.S,.1)
        cls=LowRankDenominator
        if args.retain_self:
            from research.research_max_struct_r10.src.self_retained_denominator import SelfRetainedDenominator
            cls=SelfRetainedDenominator
        op=cls(model.G,bases[args.basis],args.projection)
        (out/'DENOMINATOR_SETUP.json').write_text(json.dumps(dict(**info,
            projection=args.projection,basis=args.basis,complete_setup_cpu_seconds=time.process_time()-started,
            history_state_term_material='LEGACY_ONEHOP_DIAGNOSTIC_ONLY; actual candidate loss in material_final_block_loss'),indent=2))
        np.savez_compressed(out/'DENOMINATOR_BASIS.npz',U=op.U,B=op.B)
        supplied=dict(config,material_readout=matrix_material_readout,matrix_denominator=op,matrix_inner_maxiter=8)
        return solver.run_cie_fft(model,Y,m0,supplied)
    solver.run_auxiliary=dispatch
    spec=importlib.util.spec_from_file_location('inherited_r8',REPO/'research/research_max_struct_r8/experiments/run_case.py')
    runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
    with threadpool_limits(limits=1):
        receipt,path=runner.execute(SimpleNamespace(manifest=str(inp/'MANIFEST.json'),public_list=None,task_id=label,
            config=str(inp/'CONFIG.json'),output_dir=str(out)))
    receipt.update(study='R10',candidate_label=label,truth_access='NONE',scoring='NOT_RUN',
                   inherited_git_base_metadata=receipt.pop('git_base',None),baseline='R9 same fixed.001 DPCS GlobalDual current and inherited original physical corrector')
    runner.dump(out/'RECEIPT.json',receipt)
    print(json.dumps(dict(task=label,status=receipt['status'],CPU_seconds=receipt.get('cpu_seconds'),Phi=receipt.get('physical_objective'))),flush=True)


if __name__=='__main__':main()
