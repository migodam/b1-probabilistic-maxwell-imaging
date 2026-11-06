"""Prepare portable R10 public manifests with fresh outputs; no inversion/GT."""
from pathlib import Path
import argparse,json

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--name',required=True,help='New reproduction folder name, e.g. repeat01')
    args=ap.parse_args()
    if Path(args.name).name!=args.name:raise ValueError('Use a folder name without path separators')
    destination=ROOT/'reproduction'/args.name
    destination.mkdir(parents=True,exist_ok=False)
    for filename in ('twod_tasks.json','denominator_confirmation_tasks.json'):
        manifest=json.loads((ROOT/'configs'/filename).read_text())
        manifest['reproduction_status']='PREPARED_NOT_EXECUTED'
        manifest['input_base']='.'
        for task in manifest['tasks']:
            task['input_base']='.'
            task['output_dir']=str((ROOT/'RESULTS/twod'/args.name/task['task_id']).relative_to(REPO))
        (destination/filename).write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(status='PREPARED_NOT_EXECUTED',destination=str(destination),
                         truth_access='NONE',forward_calls=0)))


if __name__=='__main__':main()
