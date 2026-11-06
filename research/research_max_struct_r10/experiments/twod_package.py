"""Make an inverse-only portable bundle; never read fresh truth or generator labels."""
from pathlib import Path
import argparse,json,tarfile,time,io
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]


def main():
    cpu,wall=time.process_time(),time.perf_counter()
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir',type=Path,default=ROOT/'delegated/twod/transport')
    ap.add_argument('--manifest',type=Path,default=ROOT/'configs/twod_tasks.json')
    ap.add_argument('--extra-file',type=Path,action='append',default=[])
    ap.add_argument('--public-dir',type=Path,default=ROOT/'data_public/twod')
    args=ap.parse_args();out=args.output_dir
    if not out.absolute().is_relative_to(ROOT/'delegated/twod'):raise ValueError('Output outside owned tree')
    out.mkdir(parents=True,exist_ok=False)
    import threadpoolctl
    inherited=['research_max_struct_r2','research_max_struct_r4','research_max_struct_r5',
        'research_max_struct_r6','research_max_struct_r7','research_max_struct_r8',
        'research_fbg_r3_experiments','research_z_alg_r4_diversity','research_z_alg_r4_basin']
    files=[]
    for name in inherited:files.extend((REPO/'research'/name).rglob('*.py'))
    files.extend([ROOT/'src/twod_adapter.py',ROOT/'experiments/twod_run.py',
        ROOT/'experiments/twod_queue.py',ROOT/'configs/twod_solver.json',args.manifest.absolute()])
    files.extend(p.absolute() for p in args.extra_file)
    files.extend(args.public_dir.absolute().glob('*.npz'))
    files.append(args.public_dir.absolute()/'manifest.json')
    archive=out/'inverse_only.tgz';records=[]
    with tarfile.open(archive,'w:gz') as tar:
        for p in sorted(set(files)):
            rel=p.relative_to(REPO)
            if any(x in rel.parts for x in ('sealed_truth','data_private','private_truth')):
                raise ValueError('Private input in package')
            if p.name in ('twod_generator.json','twod_materials.py','twod_generate.py','twod_conditional_generator.json','twod_conditional_generate.py'):
                raise ValueError('Fresh generator configuration/code is excluded from inverse bundle')
            tar.add(p,arcname=str(rel),recursive=False)
            records.append(dict(path=str(rel),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns))
        module=Path(threadpoolctl.__file__);tar.add(module,arcname='vendor/threadpoolctl.py')
    receipt=dict(status='INVERSE_ONLY_BUNDLE_PREPARED',archive=str(archive),bytes=archive.stat().st_size,
        files=records,truth_access='NONE',fresh_generator_code_included=False,
        threadpoolctl=dict(version=threadpoolctl.__version__,source=str(module),action='copied existing pure Python module; no install'),
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        SHA_checks='NONE',inverse_runs=0)
    (out/'PACKAGE_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('status','archive','bytes','CPU_seconds')}))


if __name__=='__main__':main()
