"""Archive only the closed R4 project; preserve source, raw arrays and failures."""
from pathlib import Path
import json,time,zipfile,csv,sys,os,resource

ROOT=Path(__file__).resolve().parents[1]

def main():
    closed=ROOT/'RUN_CLOSED.json'
    if not closed.exists():raise RuntimeError('Close scientific execution before final archive')
    delivery=ROOT/'delivery';delivery.mkdir(exist_ok=True)
    target=delivery/'RESEARCH_W_GSDA_R4_20261030.zip'
    if target.exists():raise FileExistsError('Do not overwrite a completed delivery; choose a new explicit archive version')
    paths=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file():continue
        relative=p.relative_to(ROOT)
        if relative.parts[0]=='delivery' or '__pycache__' in relative.parts:continue
        if p.name in ['.DS_Store','.fem.lock','.ledger.lock'] or p.suffix in ['.pyc','.pyo']:continue
        paths.append((p,relative))
    manifest=[dict(path=str(rel),bytes=p.stat().st_size) for p,rel in paths]
    (delivery/'ARCHIVE_MANIFEST.json').write_text(json.dumps(dict(files=manifest,file_count=len(paths),uncompressed_bytes=sum(x['bytes'] for x in manifest),hash_checks='NOT_PERFORMED_NOT_REQUESTED',scope='Only new R4; all original scientific evidence remains in original directories'),ensure_ascii=False,indent=2)+'\n')
    started=time.time()
    with zipfile.ZipFile(target,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=5,allowZip64=True) as archive:
        for p,rel in paths:archive.write(p,Path('research_w_gsda_r4')/rel)
        archive.write(delivery/'ARCHIVE_MANIFEST.json','research_w_gsda_r4/delivery/ARCHIVE_MANIFEST.json')
    # Read central-directory sizes and required entries, not gratuitous hash campaigns.
    with zipfile.ZipFile(target) as archive:
        infos={x.filename:x for x in archive.infolist()}
        required=['START_HERE.md','WEB_CHAT_HANDOFF.md','REPORT_ZH.md','QUESTIONS_FOR_WEB_CHAT.md','GO_NO_GO.json','DIGITS_1_TO_9_CONTRAST.png','DIGITS_1_TO_9_METRICS.csv','PHYSICS_COST.json','RUN_CLOSED.json','tools/deploy.py','src/gsda/currents.py']
        required += [f'results/{method}/c{i:03d}/endpoint.npz' for method in ['B0_FULL_FWI','B4_GSDA','FULL_GN_TV_SUPPLEMENT'] for i in range(1,10)]
        missing=[x for x in required if 'research_w_gsda_r4/'+x not in infos]
        if missing:raise RuntimeError('Archive missing required actual artifacts: '+repr(missing))
        mismatched=[str(rel) for p,rel in paths if infos['research_w_gsda_r4/'+str(rel)].file_size!=p.stat().st_size]
        if mismatched:raise RuntimeError('Archive size mismatch: '+repr(mismatched))
    receipt=dict(status='CREATED_REQUIRED_CONTENTS_PRESENT',archive=str(target),archive_bytes=target.stat().st_size,file_count=len(paths)+1,uncompressed_bytes=sum(x['bytes'] for x in manifest),wall_seconds=time.time()-started,packaging_process_CPU_seconds=time.process_time(),peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),packaging_CPU_not_in_physics_or_inference_cost=True,required_entries_checked=len(required),required_actual_digit_endpoints=27,no_hash_campaign=True,not_published=True)
    (delivery/'DELIVERY_RECEIPT.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
    web=delivery/'RESEARCH_W_GSDA_R4_WEB_CHAT.zip'
    omitted=[];webpaths=[]
    for p,rel in paths:
        binary=p.suffix in ['.npz','.npy']
        if binary and ('public_reference' in rel.parts or 'preconditioner' in rel.parts):
            omitted.append(dict(path=str(rel),bytes=p.stat().st_size,reason='Bulk propagation cache retained in full archive; metadata and cost kept here'))
        else:webpaths.append((p,rel))
    webmanifest=delivery/'WEB_CHAT_SCOPE.json'
    webmanifest.write_text(json.dumps(dict(status='ANALYSIS_PACKAGE_WITH_ALL_SELECTED_ENDPOINTS_AND_COSTS',full_archive=target.name,omitted_bulk_arrays=omitted,raw_measurements_and_selected_materials_preserved=True,all_27_primary_digit_endpoints_preserved=True,full_jacobian_snapshots_preserved=True,spectral_and_full_absorption_validation_arrays_preserved=True,not_all_cached_operators_in_this_smaller_package=True),ensure_ascii=False,indent=2)+'\n')
    with zipfile.ZipFile(web,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=5,allowZip64=True) as archive:
        for p,rel in webpaths:archive.write(p,Path('research_w_gsda_r4')/rel)
        archive.write(webmanifest,'research_w_gsda_r4/delivery/WEB_CHAT_SCOPE.json')
        archive.write(delivery/'DELIVERY_RECEIPT.json','research_w_gsda_r4/delivery/FULL_ARCHIVE_RECEIPT.json')
    with zipfile.ZipFile(web) as archive:
        names=set(archive.namelist())
        for name in required:
            if 'research_w_gsda_r4/'+name not in names:raise RuntimeError('Web package missing primary artifact '+name)
    receipt['web_chat_archive']=str(web);receipt['web_chat_archive_bytes']=web.stat().st_size
    receipt['web_chat_omitted_bulk_cache_files']=len(omitted)
    receipt['packaging_process_CPU_seconds']=time.process_time()
    receipt['wall_seconds']=time.time()-started
    (delivery/'DELIVERY_RECEIPT.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(receipt,ensure_ascii=False))

if __name__=='__main__':main()
