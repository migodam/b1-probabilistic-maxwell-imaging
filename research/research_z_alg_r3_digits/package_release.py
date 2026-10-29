#!/usr/bin/env python3
"""Package existing frozen evidence. Never runs physics or reads private truth."""
from pathlib import Path
import argparse, json, zipfile, datetime

ROOT=Path(__file__).resolve().parent
EXCLUDE={'private_truth','__pycache__','.git','release','node_modules'}
READ_EXT={'.md','.html','.png','.pdf','.csv','.json'}

def files():
    out=[]
    for p in sorted(ROOT.rglob('*')):
        rel=p.relative_to(ROOT)
        if not p.is_file() or any(x in EXCLUDE for x in rel.parts):continue
        if p.suffix in {'.zip','.pyc'} or p.name=='.DS_Store':continue
        out.append(p)
    return out

def reading(p):
    rel=p.relative_to(ROOT)
    if p.suffix not in READ_EXT:return False
    if len(rel.parts)==1:return True
    return rel.parts[0] in {'figures','theory'} or str(rel) in {
        'audit/FINAL_SCIENTIFIC_REVIEW.md','audit/FINAL_SCIENTIFIC_REVIEW.html',
        'audit/SCALING_REVIEW.md','audit/SCALING_REVIEW.html'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--build',action='store_true');a=ap.parse_args()
    chosen=files();manifest={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'scope':'Existing scientific files; private seeds, private truth and caches excluded. No physics rerun.',
        'files':[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size} for p in chosen]}
    print(json.dumps({'files':len(chosen),'bytes':sum(x['bytes'] for x in manifest['files']),
                      'build':a.build,'private_truth_included':False}))
    if not a.build:return
    dest=ROOT/'release';dest.mkdir(exist_ok=True)
    for name,items in [('FULL_REPRODUCIBILITY.zip',chosen),('Reading_Kit.zip',[p for p in chosen if reading(p)])]:
        target=dest/name
        if target.exists():raise FileExistsError('Preserve prior archive: '+str(target))
        with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            for p in items:z.write(p,arcname=str(Path('research_z_alg_r3_digits')/p.relative_to(ROOT)))
            z.writestr('research_z_alg_r3_digits/BUNDLE_MANIFEST.json',json.dumps(manifest,indent=2))
        if target.stat().st_size>=100_000_000:raise RuntimeError('Release asset exceeds100MB')
        with zipfile.ZipFile(target) as z:
            if any('/private_truth/' in x for x in z.namelist()):raise RuntimeError('Private path')
            if z.testzip() is not None:raise RuntimeError('ZIP integrity failure')
        print(json.dumps({'archive':name,'bytes':target.stat().st_size,'members':len(items)+1}))
    (dest/'BUNDLE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))

if __name__=='__main__':main()
