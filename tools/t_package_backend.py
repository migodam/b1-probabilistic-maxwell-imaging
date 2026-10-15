"""Copy the original reusable backend and paid input subset; zero solves."""
from pathlib import Path
import shutil,json,platform,sys
import numpy as np,scipy
from research_t.target import OUT,A30,RELEASE,CASES
from research_t.budget import dump

def main():
    dst=OUT/'T_PHYSICS_HANDOFF/vendor_backend/a30_root';dst.mkdir(parents=True,exist_ok=True)
    paths=[]
    for name in ['__init__.py','physics.py','representations.py']:
        source=A30/'src/a30'/name;target=dst/'src/a30'/name
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target);paths.append(str(target.relative_to(OUT)))
    source=A30/'vendor/frozen_a29/src';target=dst/'vendor/frozen_a29/src'
    shutil.copytree(source,target,dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    for p in sorted(target.rglob('*.py')):paths.append(str(p.relative_to(OUT)))
    for bg in ['weak','ring']:
        target=dst/'inputs/a29_release/inputs/a25/background'/f'{bg}.npz';target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(OUT/'T_PHYSICS_HANDOFF/backend_assets/background'/f'{bg}.npz',target)
    for case in CASES:
        target=dst/'inputs/a29_release/inputs/a28/results/a28/exact_states'/f'{case}.npz'
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(RELEASE/'inputs/a28/results/a28/exact_states'/f'{case}.npz',target)
        target=dst/'inputs/a29_release/results/a29/pools'/f'{case}.npz';target.parent.mkdir(parents=True,exist_ok=True)
        with np.load(RELEASE/'results/a29/pools'/f'{case}.npz',allow_pickle=False) as z:
            np.savez_compressed(target,observation=z['observation'],sigma=z['sigma'],n3=z['n3'])
    for name in ['LICENSE','LICENSE.md','LICENSE.txt','COPYING']:
        for label,root in [('a30',A30),('frozen_a29',A30/'vendor/frozen_a29')]:
            if (root/name).exists():
                target=dst/'source_licenses'/label/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(root/name,target)
    dump(OUT/'T_PHYSICS_HANDOFF/vendor_backend/SOURCE_PROVENANCE.json',dict(
        upstream_a30_commit='717d3d26a290248bd2b9aeb55f3d03694b407b84',upstream_a29_commit='5d5e31fcfc56ac8ef9d64c46b4b7bbdedef96bb4',
        operation='exact source copies plus explicit already-paid input subset, no hashing/no solving',
        source_files=paths,pool_subset_fields=['observation','sigma','n3'],
        cached_states='8 exact original exposed N3 caches; no old truth/teacher bank copied',
        known_background='sanitized required physical arrays and metadata, no teacher J/step',
        historical_setup_cost='already paid, not quantified here; not free new training data',
        license_status='existing license files preserved if present; no new grant inferred',new_full_state_attempts=0,new_full_solver_RHS=0))
    dump(OUT/'EXECUTION_ENVIRONMENT.json',dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,
                                             architecture=platform.machine(),platform=platform.platform(),
                                             GPU_training=0,BLAS_threads=1))
if __name__=='__main__':main()
