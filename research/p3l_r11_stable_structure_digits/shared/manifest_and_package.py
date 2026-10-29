"""Content identities and compact lossless evidence package, no physics."""
from core import *
import subprocess,zipfile,datetime,platform

def main():
 ignore={'packaging','__pycache__','.DS_Store','PACKAGE.log'};items=[]
 for p in sorted(P.rglob('*')):
  if not p.is_file() or any(a in ignore for a in p.relative_to(P).parts) or p.name in ['SOURCE_MANIFEST.json','PUBLICATION.md','PUBLIC_READBACK.json','WEB_READBACK.json']:continue
  items.append(dict(path=str(p.relative_to(REPO)),bytes=p.stat().st_size,blake2b_256=blake(p),role='RAW_LOSSLESS_EVIDENCE' if p.suffix in ['.npz','.jsonl'] else 'CODE_CONFIG_REPORT_FIGURE_OR_COMPACT_RESULT'))
 sources=dict(id='R11_STABLE_DIGITS_SOURCE_AND_EVIDENCE_V1',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),base_commit='d6120d4bc0f9d8699b9e788da7abee345752b8f6',historical_R4_source_commit='092da74e320de38c086d07f24a1d98ca4892e34a',source_identity='exactBLAKE2bfileidentities;releasepinsGitcommitwithoutself-referentialmanifesthash',immutable_history=True,files=items,files_total=len(items),bytes_total=sum(x['bytes'] for x in items),external_reused_inputs=[dict(path=str(R9/'results/stable_roi/SETUP.npz'),bytes=(P/'inputs/FIXED_HOST_SETUP.npz').stat().st_size,blake2b_256=blake(P/'inputs/FIXED_HOST_SETUP.npz'),historical_physics_RHS576=264,new_RHS=0),dict(path=str(R4/'src/p1_material/prior.py'),blake2b_256=blake(P/'inputs/PRIOR_SOURCE.py'))],resource_status='FAIL_CAP_EXCEEDED',python=sys.version,platform=platform.platform(),cpu='AppleM4',MPS_GPU_hours=0,actual_workers='three native gpt-6.1-sol ultra, Main currentCodexproviderunchanged')
 write(P/'SOURCE_MANIFEST.json',sources)
 out=P/'packaging';out.mkdir(exist_ok=True)
 include=[]
 for p in P.rglob('*'):
  if p.is_file() and not any(a in ignore for a in p.relative_to(P).parts):include.append(p)
 at=REPO/'research/maxwell_formula_atlas'
 include += list(at.glob('*.md'))+list(at.glob('*.html'))+[at/'FORMULA_REGISTRY.json',at/'SOURCE_FORMULA_INDEX.csv']
 paths=sorted(set(str(p.relative_to(REPO)) for p in include))
 (out/'FILES.txt').write_text('\n'.join(paths)+'\n')
 archive=out/'R11_stable_structure_digits_complete.tar.zst'
 with archive.open('wb') as f:
  tar=subprocess.Popen(['tar','-cf','-','-T',str(out/'FILES.txt')],cwd=REPO,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  zstd=subprocess.run(['zstd','-3','-T1','-c'],stdin=tar.stdout,stdout=f,stderr=subprocess.PIPE)
  tar.stdout.close();err=tar.stderr.read();code=tar.wait()
  if code or zstd.returncode:raise RuntimeError((code,err.decode(),zstd.stderr.decode()))
 reading=out/'R11_reading_kit.zip'
 with zipfile.ZipFile(reading,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as f:
  for p in include:
   if p.suffix in ['.md','.html','.png','.csv','.py'] and 'results/' not in str(p.relative_to(REPO)) or p.name in ['EXPERIMENT_CONTRACT.json','EXPERIMENT_RESULTS.json','COST_SUMMARY.json','RESULTS_SUMMARY.json','SOURCE_MANIFEST.json']:
    f.write(p,str(p.relative_to(REPO)))
 assets=[dict(name=p.name,bytes=p.stat().st_size,blake2b_256=blake(p)) for p in [archive,reading]]
 write(out/'R11_RELEASE_MANIFEST.json',dict(release='p3lr11-stable-digits-20261030',created_utc=time.time(),scientific_base=sources['base_commit'],science_version='release target Git commit, recorded at publication',archive_paths='repository-relative, nooverwritinghistory',file_count=len(paths),assets=assets,resource_status='FAIL_CAP_EXCEEDED',new_RHS=363240,new_LU=28727,unqualified_claims=['full66posterior','newmaterialfibers','universalwholeimagingspeedup','strongTAPnovelty']))
 print(json.dumps(dict(files=len(paths),assets=assets),indent=2))
if __name__=='__main__':main()
