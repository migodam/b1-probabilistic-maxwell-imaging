"""Collect completed artifacts without reading private arrays in the solver."""
from pathlib import Path
import zipfile,subprocess,sys
root=Path(__file__).resolve().parent
for p in ['COHORT_SEAL.json','objects.csv','scoring_cost.json']:
    if not (root/p).is_file():raise FileNotFoundError(p)
subprocess.run([sys.executable,str(root/"qualify_finer.py")],check=True)
with zipfile.ZipFile(root/'results_bundle.zip','w',zipfile.ZIP_DEFLATED) as z:
    for folder in ['predictions','data','sealed_truth','development_ring']:
        for p in sorted((root/folder).rglob('*')):
            if p.is_file():z.write(p,p.relative_to(root))
    for p in ['COHORT_SEAL.json','objects.csv','scoring_cost.json']:z.write(root/p,p)
