"""Existing-repository evidence assets; no physical replay or file hashes."""
from pathlib import Path
import sys,json,zipfile
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from material_fiber.core import OUT,write

def main():
    release='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p2-material-fiber-20261019'
    selected={'OLD_B'}|set(json.loads((OUT/'A/ANCHORS.json').read_text()).values())|set(json.loads((OUT/'C/BEST_POINTS.json').read_text()).values())|{'C_00008'}
    dest=ROOT/'publication_local';dest.mkdir(exist_ok=True);assets=[]
    def raw(p):
        rel=p.relative_to(OUT)
        return len(rel.parts)>1 and rel.parts[0] in 'ABCDE' and rel.parts[1] in ['endpoints','attempts','evaluations']
    for kind,name in [('reports','P2_Material_Fiber_Code_Reports_20261019.zip'),('raw','P2_Material_Fiber_Raw_Maxwell_20261019.zip')]:
        archive=dest/name;n=0
        with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=3,allowZip64=True) as z:
            for p in sorted(OUT.rglob('*')):
                if not p.is_file() or p.name=='executor.lock' or p.name.endswith('.log') or p.name in ['PUBLICATION_RECEIPT.json','WEB_READBACK.json','EVIDENCE_RELEASE.json']:continue
                if '__pycache__' in p.parts:continue
                if 'cache' in p.relative_to(OUT).parts and p.stem not in selected:continue
                if (kind=='raw')!=raw(p):continue
                z.write(p,p.relative_to(ROOT));n+=1
            if kind=='reports':
                for folder in ['src/material_fiber','src/p1_material','src/a17x']:
                    for p in (ROOT/folder).glob('*.py'):z.write(p,p.relative_to(ROOT));n+=1
                for p in (ROOT/'tools').glob('*material_fiber*.py'):z.write(p,p.relative_to(ROOT));n+=1
                for name2 in ['publish_lxz.py','publish_research_f.py']:
                    p=ROOT/'tools'/name2;z.write(p,p.relative_to(ROOT));n+=1
                for p in (ROOT/'src/lxz').glob('*.py'):z.write(p,p.relative_to(ROOT));n+=1
        assets.append(dict(file=archive.name,bytes=archive.stat().st_size,files=n))
    write(OUT/'EVIDENCE_RELEASE.json',dict(release=release,assets=assets,all_2012_native_raw_endpoint_arrays_included=True,all_paid_attempts_included=True,raw_Gaussian_uniform_generation_and_noise_included=True,selected_cache_keys=sorted(selected),other_LU_cache='retained locally, derivable computational cache; not raw material/prediction loss',no_hash_checks=True,new_repo=False,publication_independent_of_scientific_verdict=True))
    print(json.dumps(dict(assets=assets),ensure_ascii=False))
if __name__=='__main__':main()
