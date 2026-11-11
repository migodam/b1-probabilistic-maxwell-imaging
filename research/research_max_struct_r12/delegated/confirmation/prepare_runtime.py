"""Rebase public-only task paths on the execution host; no new data/inversion."""
from pathlib import Path
import argparse,json,sys,os
ROOT=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--legacy-repo',type=Path,required=True)
    ap.add_argument('--stage',choices=('dev6','confirm24','all'),default='dev6')
    ap.add_argument('--search-runner',type=Path)
    ap.add_argument('--weights',type=Path)
    a=ap.parse_args();legacy=a.legacy_repo.resolve()
    stages=('dev6','confirm24') if a.stage=='all' else (a.stage,)
    package_search=ROOT.parent/'research/research_max_struct_r12/experiments/run_confirmation_search.py'
    local_search=ROOT.parents[1]/'experiments/run_confirmation_search.py'
    search=(a.search_runner or (package_search if package_search.is_file() else local_search)).resolve()
    package_weights=ROOT/'assets/physics_weights.pt'
    weights=(a.weights or (package_weights if package_weights.is_file() else legacy/'research/research_max_struct_r8/RESULTS/training/physics/weights.pt')).resolve()
    if not search.is_file() or not weights.is_file():raise FileNotFoundError('Qualified search runner and frozen physics weights required')
    if not (legacy/'research/research_max_struct_r8/experiments/run_case.py').is_file():raise FileNotFoundError(legacy)
    results=[]
    for stage in stages:
      public=ROOT/'data_public'/stage
      source=json.loads((public/'manifest.json').read_text())
      tasks=[];baseline=[]
      config=ROOT/('DEV6_BASELINE_SOLVER.json' if stage=='dev6' else 'CONFIRM24_BASELINE_SOLVER.json')
      if not config.is_file():raise FileNotFoundError(config)
      for index,record in enumerate(source['records']):
        case=record['case_id']
        files=[public/'operators.npz',public/(case+'.npz')]
        if not all(p.is_file() for p in files):raise FileNotFoundError(files)
        methods=('cie_single','cie_dual','global_dual','a_plus','fullG19','dr','multistart','nn_direct')
        offset=index%8
        for method in methods[offset:]+methods[:offset]:
            ident=case+'__'+method+'__240CPU'
            task=dict(task_id=ident,case_id=case,arm='a_plus' if method=='fullG19' else method,
                algorithm_label=method,dpcs_rule='noise_ball_99',initialization='uniform',initial_material=.1,
                public_files=list(map(str,files)),input_base=str(legacy),config=str(config),
                output_dir=str(ROOT/'RESULTS'/stage/ident),experiment_scope='R12_'+stage.upper(),
                executor_kind='search' if method in ('dr','multistart','nn_direct') else 'baseline')
            if method=='fullG19':task['material_readout']=dict(module_file=str(legacy/'research/research_max_struct_r11/src/material_readouts.py'),
                callable='full_material_readout',prepare='prepare_full',label='fullG19')
            if task['executor_kind']=='search':
                task.update(runner=str(search),public_dir=str(public),method=method,
                    weights=str(weights) if method=='nn_direct' else None)
            else:
                task['runner']=str(ROOT/'run_baseline.py');baseline.append(task)
            tasks.append(task)
      destinations=[ROOT/('RUNTIME_'+stage.upper()+'_BASELINE_TASKS.json'),ROOT/('RUNTIME_'+stage.upper()+'_EIGHT_METHOD_TASKS.json')]
      for dest in destinations:
        if dest.exists():raise FileExistsError(dest)
      common=dict(status='PUBLIC_RUNTIME_PATHS_NO_GT',stage=stage,legacy_repo=str(legacy),
            dispatch='Finite FIFO; case-index cyclic method order, fixed before inversion',truth_access='NONE')
      destinations[0].write_text(json.dumps(dict(common,tasks=baseline),indent=2)+'\n')
      destinations[1].write_text(json.dumps(dict(common,tasks=tasks),indent=2)+'\n')
      results.append(dict(stage=stage,baseline_manifest=str(destinations[0]),eight_method_manifest=str(destinations[1]),
        tasks=len(tasks),baseline_tasks=len(baseline),config=str(config)))
    print(json.dumps(dict(status='READY_NOT_LAUNCHED',stages=results,
        required_environment=dict(R12_LEGACY_REPO=str(legacy)),python=sys.executable,vendor_path=str(legacy/'vendor'),
        search_runner=str(search),physics_weights=str(weights),truth_access='NONE')),flush=True)

if __name__=='__main__':main()
