"""Thin path router for existing entry points; explicit fresh output required.

No scheduler, installs, numerical changes, or automatic multi-job execution.
Run from any cwd using the complete repository and existing CUDA environment.
"""
import argparse,importlib,sys,os,json,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
PREFIX='research.research_max_struct_r9.experiments.'
MODES={
 'benchmark':('benchmark_original_cuda',None,'CUDA_WHOLE_RHS_BENCHMARK.json'),
 'generate':('vector_original_cuda','generate','public/observations.npz'),
 'qualify':('vector_original_cuda','qualify-discretization','DISCRETIZATION_QUALIFICATION.json'),
 'nine':('run_original_cm_nine',None,'cohorts'),
 'score':('score_original_cm_nine',None,'NINE_OFFLINE_SCORES.json'),
 'plot':('plot_original_cm_nine',None,'NINE_PLOT_COST.json'),
 'cost':('finalize_original_cm_accounting',None,'ACCOUNTING_FINAL.json'),
 'igt-prepare':('diagnose_original_igt','prepare','igt_source59/prepare_PROCESS_COST.json'),
 'igt-solve':('diagnose_original_igt','solve','igt_source59/solve_PROCESS_COST.json')}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run-dir',required=True,type=Path);p.add_argument('mode',choices=MODES);a=p.parse_args()
    out=a.run_dir.expanduser().resolve();name,sub,guard=MODES[a.mode]
    if (out/guard).exists():raise SystemExit('Refusing to overwrite existing sealed/recorded output: '+str(out/guard))
    if a.mode=='generate' and (out/'generate_CUDA_COST.json').exists():
        receipt=json.loads((out/'generate_CUDA_COST.json').read_text())
        if receipt.get('status')=='COMPLETE':raise SystemExit('Completed generation receipt exists; use a fresh directory')
        history=out/'resumed_generation_attempts';index=0
        while (history/('attempt%03d'%index)).exists():index+=1
        archive=history/('attempt%03d'%index);archive.mkdir(parents=True)
        for item in ['generate_CUDA_COST.json','generation_CUDA_progress.jsonl','GENERATION_SETUP_CUDA.json']:
            if (out/item).exists():shutil.move(str(out/item),str(archive/item))
    if a.mode=='qualify' and (out/'qualify-discretization_CUDA_COST.json').exists():
        raise SystemExit('Prior reference attempt is recorded; refusing silent overwrite/repetition')
    if a.mode in ('benchmark','generate','qualify','nine','igt-prepare','igt-solve'):
        if os.name!='nt':raise SystemExit('Recorded numerical runner requires Windows psutil affinity/peak_wset. No automatic platform/backend substitution.')
        import psutil
        if not set([12,14]).issubset(set(psutil.Process().cpu_affinity())):raise SystemExit('Recorded P-core affinity12/14 unavailable; explicit host-only adaptation required before running.')
    module=importlib.import_module(PREFIX+name)
    if name=='diagnose_original_igt':module.PARENT=out;module.OUT=out/'igt_source59'
    else:module.OUT=out
    if a.mode=='plot':module.FIG=out/'FIGURES'
    sys.argv=[name]+([sub] if sub else [])
    module.main()

if __name__=='__main__':main()
