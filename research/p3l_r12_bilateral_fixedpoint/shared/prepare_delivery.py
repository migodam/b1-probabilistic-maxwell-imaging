"""Build one reading/reuse ZIP and manifests; no scientific checks or replay."""
from pathlib import Path
import json,hashlib,zipfile,time
from .budget import ROOT,atomic_json
def main():
    files=[];code=[]
    for f in ROOT.rglob('*'):
        if not f.is_file():continue
        rel=f.relative_to(ROOT)
        if any(x in rel.parts for x in ('packaging','__pycache__','states')) or f.name.endswith('.tmp') or f.name=='GLOBAL_BUDGET.lock':continue
        files.append(f)
        if f.suffix=='.py':code.append(dict(path=str(rel),bytes=f.stat().st_size,blake2b=hashlib.blake2b(f.read_bytes(),digest_size=32).hexdigest()))
    raw=json.loads((ROOT/'packaging/RAW_ARCHIVE_CREATION.json').read_text())
    atomic_json(ROOT/'SOURCE_MANIFEST.json',dict(created_utc=time.time(),base_commit='63f8707c3ad06a66ca9b992ece118b2e0e3100fb',code_identity=code,original_P1_prior='inputs/PRIOR_SOURCE.py + PRIOR_CONFIG.json + PRIOR_POINTS.npz',original_fixed_host_cache='inputs/FIXED_HOST_SETUP.npz',reused_exposed_development='inputs/reuse/P1',frozen_N_truth_and_anchor='inputs/reuse/N',source_changes='P1 initialimplementation, serialization-onlyresume, nativeN1boundaryfailure, explicitpreGTN2revision; defaultoldbehaviorretained',data_predictions='rootCLOUD.npz identities in each SEALED.json',raw_assets=[dict(name=r['name'],bytes=r['bytes'],relative_source_paths=r['source_paths']) for r in raw],no_new_SCIENTIFIC_verification=True))
    target=ROOT/'packaging/R12_reading_and_reuse.zip'
    with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for f in sorted(set(files+[ROOT/'SOURCE_MANIFEST.json'])):
            z.write(f,str(Path('research')/ROOT.name/f.relative_to(ROOT)))
    assets=[dict(name=r['name'],bytes=r['bytes'],type='losslessrawphysicalarrays',restore='extractatB1repositoryroot') for r in raw]+[dict(name=target.name,bytes=target.stat().st_size,type='reports+code+compactclouds+reuseinputs; NOTallrawJstates',restore='extractatB1repositoryroot')]
    atomic_json(ROOT/'packaging/R12_RELEASE_MANIFEST.json',dict(id='P3L_R12_BF',created_utc=time.time(),source_commit='seePUBLICATION.md; frozenbeforeupload',raw_physical_results_unchanged_since_stop=True,assets=assets,scientific_scope='DEV_ONLY+DESIGNED_DIGITS; hierarchicalmethodNO_GO; nofullposterior',actual_cost=json.loads((ROOT/'GLOBAL_BUDGET_STATE.json').read_text())['used'],reading_kit_not_full_evidence=True,no_duplicate_array_comparison_requested=True))
    print(json.dumps(dict(reading_files=len(files)+1,reading_bytes=target.stat().st_size,all_asset_bytes=sum(r['bytes'] for r in assets)),indent=2))
if __name__=='__main__':main()
