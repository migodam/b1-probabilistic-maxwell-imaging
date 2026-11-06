"""Two exposed DEV controls: old one-hop metric, same eight-step inner LBFGS.

Question frozen before launch: can replacing the A+ inner optimizer alone
explain the matrix-denominator rescue? No fresh-data tuning or truth input.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import sys,json,argparse,importlib.util
from types import SimpleNamespace
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--historical-repo',required=True)
    ap.add_argument('--case',required=True)
    ap.add_argument('--public-files',nargs='+',help='Portable public observation/operator paths; no truth')
    ap.add_argument('--output-dir',help='Fresh reproduction output')
    args=ap.parse_args();old=Path(args.historical_repo)
    from research.research_max_struct_r8.src import benchmark_solver as solver
    from research.research_max_struct_r4.src.cie_onehop import onehop_material_readout
    from threadpoolctl import threadpool_limits
    task0=json.loads((old/'research/research_max_struct_r9/RESULTS/development'/f'{args.case}__global_dual__fixed_001'/'PUBLIC_TASK.json').read_text())
    if args.public_files:task0['public_files']=args.public_files
    label=f'{args.case}__onehop_lbfgs8_corrected_dispatch'
    out=Path(args.output_dir).absolute() if args.output_dir else ROOT/'RESULTS/onehop_solver_control'/label
    inp=out.parent/(out.name+'__inputs');inp.mkdir(parents=True,exist_ok=False)
    cfg=json.loads((REPO/'research/research_max_struct_r9/configs/inherited_solver.json').read_text())
    cfg.update(run_id='r10_onehop_solver_only_dev',material_inner_maxiter=8,
               material_inner_ftol=1e-12,material_inner_gtol=1e-8,
               cache_material_value_gradient=True)
    # R8 constructs the arm schedule from common first; only aplus_options
    # overrides that constructed schedule. Top-level keys alone are not used.
    cfg.setdefault('aplus_options',{}).update(material_inner_maxiter=8,
        material_inner_ftol=1e-12,material_inner_gtol=1e-8,cache_material_value_gradient=True)
    cfg['candidate']=dict(metric='unchanged original onehop',inner='bounded LBFGS8',
        seed='same ordinary CIE seed or old point by original onehop objective',
        frozen_before_complete_trajectory=True,truth_access='NONE',
        scientific_question='Does inner optimizer alone explain the matrix denominator development rescue?',
        scope='two historical exposed objects; excluded from fresh confirmation')
    (inp/'CONFIG.json').write_text(json.dumps(cfg,indent=2))
    task=dict(task0,task_id=label,arm='a_plus',config=str(inp/'CONFIG.json'),output_dir=str(out))
    (inp/'MANIFEST.json').write_text(json.dumps(dict(tasks=[task]),indent=2))
    def dispatch(model,Y,m0,arm,config):
        return solver.run_cie_fft(model,Y,m0,dict(config,material_readout=onehop_material_readout))
    solver.run_auxiliary=dispatch
    spec=importlib.util.spec_from_file_location('r8_onehop_control',REPO/'research/research_max_struct_r8/experiments/run_case.py')
    runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
    with threadpool_limits(limits=1):
        receipt,path=runner.execute(SimpleNamespace(manifest=str(inp/'MANIFEST.json'),public_list=None,
            task_id=label,config=str(inp/'CONFIG.json'),output_dir=str(out)))
    receipt.update(study='R10',candidate_label=label,truth_access='NONE',scoring='NOT_RUN',
                   inherited_git_base_metadata=receipt.pop('git_base',None),
                   baseline='same inherited global current/basis/corrector, original onehop metric, eight-step LBFGS')
    runner.dump(out/'RECEIPT.json',receipt)
    print(json.dumps(dict(task=label,status=receipt['status'],CPU_seconds=receipt.get('cpu_seconds'),
                         Phi=receipt.get('physical_objective'))),flush=True)


if __name__=='__main__':main()
