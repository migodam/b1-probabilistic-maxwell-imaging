"""Package frozen raw evidence without rerunning, hashing or comparing experiments."""
import json,zipfile,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[1]

def main():
    begun=time.perf_counter();out=ROOT/'release';out.mkdir(exist_ok=True)
    files=[p for p in ROOT.rglob('*') if p.is_file() and not any(x in p.relative_to(ROOT).parts for x in ('private_truth','__pycache__','release','.venv'))]
    history=ROOT.parent/'research_z_alg_r3_digits'
    support=[]
    for name in ('src/__init__.py','src/physics.py','src/solver.py','src/multiscale.py','experiments/generate_data.py','experiments/score_predictions.py','public_data/operators.npz'):
        support.append(history/name)
    for split,case in [('dev_v2','dev_2a5c221387b1'),('holdout','obs_c3604cf3555d'),('holdout','obs_42c48e795766'),('holdout','obs_a63df7f20566')]:
        support.extend([history/'truth_after_freeze'/split/(case+'.npz'),history/'truth_after_freeze'/split/'mapping_seedless.json',history/'public_data'/('dev' if split=='dev_v2' else split)/(case+'.npz')])
        for method in ('fullj_native','fullj_reciprocal','globaldual'):support.append(history/'results'/split/method/(case+'.json'))
    allfiles=list(dict.fromkeys(files+support))
    full=out/'R4_BASIN_FULL_REPRODUCIBILITY.zip'
    with zipfile.ZipFile(full,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5) as z:
        for p in allfiles:z.write(p,p.relative_to(REPO).as_posix())
    reading=[p for p in files if (p.suffix.lower() in ('.md','.html','.png','.pdf','.csv') and 'published_results' not in p.parts and 'results' not in p.parts) or p.name in ('EXPERIMENT_RESULTS.json','RESULTS_SUMMARY.json','COST_TOTALS.json')]
    light=out/'R4_BASIN_READING_KIT.zip'
    with zipfile.ZipFile(light,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5) as z:
        for p in reading:z.write(p,p.relative_to(REPO).as_posix())
    manifest=dict(files=[dict(path=p.relative_to(REPO).as_posix(),bytes=p.stat().st_size) for p in allfiles],
       full_archive=dict(name=full.name,bytes=full.stat().st_size,members=len(allfiles)),
       reading_kit=dict(name=light.name,bytes=light.stat().st_size,members=len(reading)),
       package_seconds=time.perf_counter()-begun,
       privacy='private_truth and generation entropy excluded; only already-disclosed seedless materials',
       validation='No hashes, repeated content comparisons or new scientific experiments; raw evidence preserved by packaging.')
    (out/'BUNDLE_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:v for k,v in manifest.items() if k!='files'},indent=2),flush=True)

if __name__=='__main__':main()
