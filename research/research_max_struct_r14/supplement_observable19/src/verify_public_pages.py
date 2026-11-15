#!/usr/bin/env python3
"""Anonymous supplement Pages readback; transport checks, no science or hashes.

Adapted from R14/src/verify_public_pages.py: identity encoding, 45-second
anonymous requests, ranged 64-byte media checks and retained failure rows.
Unlike the parent, pages and public sources must match current local bytes.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path
import struct
import time
from urllib.parse import quote, urlsplit
from urllib.request import Request, urlopen

ROOT=Path(__file__).resolve().parents[1]
BASE_URL='https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r14/supplement_observable19/'
BLOCKED={'sealed_truth','private_operational','results','predictions','chat_readback','credentials','secrets'}
DOCS=['README.md','RESEARCH_REPORT_ZH.md','THEORY_CLARIFICATION_ZH.md','EXPERIMENT_CONTRACT.md',
      'REPRODUCE.md','COST_LEDGER.md','FAILURES.md','CHAT_HANDOFF.md','DOWNLOADS.json']
CSVS=['SUMMARY.csv','ALL_VIEWS.csv','CORRECTOR_EFFECTS.csv','RECOVERY_GAIN_LOSS.csv',
      'NORMALIZED_GAIN_LOSS.csv','PAIRED_DIFFERENCES.csv','NORMALIZED_PAIRED_DIFFERENCES.csv',
      'COST_LEDGER.csv','QUEUE_COSTS.csv','ISSUES.csv','SELECTOR_AUDIT.csv']


def safe_local(path,root):
    rel=Path(path)
    if rel.is_absolute() or any(x in ('.','..') or x.lower() in BLOCKED for x in rel.parts):
        raise ValueError('Nonpublic or outside path: '+str(path))
    local=root/rel
    if any((root/Path(*rel.parts[:i])).is_symlink() for i in range(1,len(rel.parts)+1)):
        raise ValueError('Symlink source excluded: '+str(path))
    if not local.resolve().is_relative_to(root):raise ValueError('Path outside supplement')
    return local


def verify(row):
    cpu,wall=time.process_time(),time.perf_counter()
    result={k:v for k,v in row.items() if k!='local'}
    try:
        local=row['local']
        if not local.is_file():raise FileNotFoundError(str(local))
        request=Request(row['url'],headers={'User-Agent':'MAX-STRUCT-supplement-anonymous-Pages-audit',
            'Accept-Encoding':'identity',**({'Range':'bytes=0-63'} if row['kind']=='media' else {})})
        with urlopen(request,timeout=45) as response:
            body=response.read(64) if row['kind']=='media' else response.read()
            result.update(status=response.status,final_url=response.geturl(),content_type=response.headers.get('Content-Type'),
                          bytes_read=len(body),header_length=response.headers.get('Content-Length'))
        if row['kind']=='media':
            with local.open('rb') as source:expected=source.read(64)
            result['first64_bytes_match_local']=body==expected
            if local.suffix.lower()=='.png' and len(body)>=24 and body.startswith(b'\x89PNG\r\n\x1a\n'):
                result['PNG_dimensions']=list(struct.unpack('>II',body[16:24]))
            if local.suffix.lower()=='.pdf':result['PDF_signature']=body.startswith(b'%PDF-')
            result['verified']=result['status'] in (200,206) and body==expected
        else:
            expected=local.read_bytes()
            result.update(local_bytes=len(expected),full_content_bytes_match_local=body==expected)
            result['verified']=result['status']==200 and body==expected
        if not result['verified']:result['error']='HTTP status or current-local-byte comparison failed'
    except Exception as exc:
        result.update(verified=False,error=repr(exc))
    result.update(wall_seconds=time.perf_counter()-wall)
    # Whole-process CPU is measured once in main; per-thread CPU would overlap.
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--base-url',default=BASE_URL)
    parser.add_argument('--workers',type=int,choices=range(1,9),default=8)
    args=parser.parse_args();cpu,wall=time.process_time(),time.perf_counter();root=args.root.resolve()
    results=[];rows={};counts={};inventory_path=root/'audit/site_implementation/BUILD_RECEIPT.json'
    try:
        parsed=urlsplit(args.base_url)
        if parsed.scheme!='https' or not parsed.netloc or parsed.username or parsed.password or parsed.query or parsed.fragment:
            raise ValueError('Base URL must be an anonymous HTTPS directory URL')
        base=args.base_url.rstrip('/')+'/'
        inventory=json.loads(inventory_path.read_text(encoding='utf8'))
        def add(path,kind):
            path=str(path).replace('\\','/')
            local=safe_local(path,root)
            rows.setdefault(path,dict(path=path,kind=kind,url=base+quote(path,safe='/._-'),local=local))
        for path in inventory['generated_pages']:add(path,'page')
        media=set(inventory['inventory']['images'])|set(inventory['inventory']['standalone_pdfs'])
        for path in list(media):
            # The builder links a same-stem PDF without making a second card.
            pdf=str(Path(path).with_suffix('.pdf'))
            if safe_local(pdf,root).is_file():media.add(pdf)
        for path in sorted(media):add(path,'media')
        for path in DOCS:add(path,'source')
        for name in CSVS:add('analysis/observable19_v1/'+name,'source')
        for path in ['analysis/fixed_state_summary/SUMMARY.json','analysis/fixed_state_summary/SUMMARY.md',
                     'analysis/fixed_state_summary/DIAGNOSTIC_COMPACT.csv']:
            add(path,'source')
        for directory,suffix in [('src','.py'),('configs','.json')]:
            folder=root/directory
            if folder.is_dir() and not folder.is_symlink():
                for path in sorted(folder.iterdir()):
                    if path.suffix==suffix and path.is_file():add(str(path.relative_to(root)),'source')
        counts=dict(generated_pages=len(inventory['generated_pages']),PNG_items=sum(Path(p).suffix.lower()=='.png' for p in media),
                    PDF_items=sum(Path(p).suffix.lower()=='.pdf' for p in media))
        expected=dict(generated_pages=9,PNG_items=13,PDF_items=3)
        for key,value in expected.items():
            if counts[key]!=value:results.append(dict(path='BUILD_RECEIPT.json inventory:'+key,kind='inventory',verified=False,
                                                      expected=value,actual=counts[key],error='Actual inventory count differs from requested final publication set'))
        with ThreadPoolExecutor(max_workers=args.workers) as executor:results+=list(executor.map(verify,rows.values()))
    except Exception as exc:
        results.append(dict(path=str(inventory_path.relative_to(root)),kind='setup',verified=False,error=repr(exc)))
    failures=[row for row in results if not row['verified']]
    receipt=dict(status='ALL_REQUESTED_SUPPLEMENT_PAGES_AND_MEDIA_READ_BACK' if not failures else 'PUBLIC_ACCESS_FAILURES_RETAINED',
                 created_utc=datetime.now(timezone.utc).isoformat(),base_url=args.base_url,build_receipt=str(inventory_path),
                 items=len(results),counts=counts,results=results,failures=failures,authentication='NONE',
                 comparison='Media first64 bytes; pages and public sources full bytes against local snapshot; no inherited parent report markers',
                 CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,physics=0,scoring=0,
                 truth_arrays_loaded=False,hash_operations='NONE',
                 visual_QA='Separate parent review; transport matches do not establish full media visual review or scientific validity')
    out=root/'audit/publication';out.mkdir(parents=True,exist_ok=True)
    (out/'PUBLIC_PAGES_READBACK.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print(json.dumps({k:receipt[k] for k in ['status','items','counts','CPU_seconds','wall_seconds']},indent=2))
    return 1 if failures else 0


if __name__=='__main__':raise SystemExit(main())
