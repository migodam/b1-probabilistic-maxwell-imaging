"""Rebase recorded public tasks onto an extracted R9 checkout; no solves/GT.

Original execution manifests remain immutable evidence. This writes new task
manifests with fresh output directories and bundled public-file references.
"""
import argparse,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cohort',choices=['development','confirmation'],required=True)
    ap.add_argument('--destination',required=True,type=Path)
    args=ap.parse_args();dest=args.destination.resolve();dest.mkdir(parents=True,exist_ok=False)
    names=(['baseline_development','baseline_historical_fixed','baseline_cost_controls','landscape_development']
           if args.cohort=='development' else ['confirmation_baselines','confirmation_landscape','confirmation_single_noise'])
    for name in names:
        data=json.loads((ROOT/'configs'/(name+'.json')).read_text())
        for t in data['tasks']:
            t['output_dir']=str(dest/'outputs'/t['task_id'])
            if args.cohort=='development':
                t['public_files']=[str(ROOT/'data_public/development'/Path(f).name) for f in t['public_files']]
            else:
                t['public_files']=[str(ROOT/('data_public/operators' if 'operators' in Path(f).parts else 'data_public/confirmation')/Path(f).name) for f in t['public_files']]
            if 'arm' in t:t['config']=str(ROOT/'configs/inherited_solver.json')
        (dest/(name+'.json')).write_text(json.dumps(data,indent=2)+'\n')
    print('Prepared '+str(len(names))+' public manifests; choose this checkout as --r8-repo for baseline runner. No solver or GT access.')

if __name__=='__main__':main()
