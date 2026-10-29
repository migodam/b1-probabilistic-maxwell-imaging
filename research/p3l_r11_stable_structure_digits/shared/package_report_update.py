"""Report-only amendment package; never alters frozen physics arrays."""
from pathlib import Path
import json,time,hashlib,zipfile,subprocess
P=Path(__file__).resolve().parents[1];REPO=P.parents[1]
files=[]
for pat in ['*.md','*.html']:files+=list(P.glob(pat))+list((P/'GPT_PRO_HANDOFF').glob(pat))
files += [P/f for f in ['SOURCE_MANIFEST.json','EXPERIMENT_RESULTS.json','RESULTS_SUMMARY.json','DIGIT_PER_NUMBER_RESULTS.csv','shared/score_and_plot.py','shared/build_site.py','shared/package_report_update.py','figures/DECODER_3_5_VISUAL_AUDIT.png','audit/DECODER_VISUAL_QUALIFICATION.json','audit/REPORT_UPDATE_DATA_IMMUTABILITY.json','audit/FINAL_QA_RESOLUTION.json','audit/INITIAL_RELEASE_NOTES.md','audit/REPORT_UPDATE_PACKAGING_FAILURE.json']]
files=sorted(set(files));out=P/'packaging/R11_report_qualification_update.zip'
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in files:z.write(f,str(f.relative_to(REPO)))
a=dict(id='R11_REPORT_ONLY_QUALIFICATION_UPDATE',created_utc=time.time(),report_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True,cwd=REPO).strip(),immutable_execution_data_commit='7d4d1fa2930afacdd8f683c3025e5538b4f3be8a',initial_release_assets_unchanged=True,changes='P3visualrepresentationfailthisfit,per-digitChinesevisualexplanation,clarifymetadata;noarray/physicalsolver/config/likelihoodgate/materialMSEchange',new_RHS=0,new_LU=0,assets=[dict(name=out.name,bytes=out.stat().st_size,blake2b_256=hashlib.blake2b(out.read_bytes(),digest_size=32).hexdigest())],files=[dict(path=str(f.relative_to(REPO)),bytes=f.stat().st_size,blake2b_256=hashlib.blake2b(f.read_bytes(),digest_size=32).hexdigest()) for f in files])
(P/'packaging/R11_REPORT_UPDATE_MANIFEST.json').write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(dict(bytes=out.stat().st_size,files=len(files),CRCbad=zipfile.ZipFile(out).testzip(),blake2b_256=a['assets'][0]['blake2b_256'])))
