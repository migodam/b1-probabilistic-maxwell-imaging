"""Anonymous access check for actually built R15 Pages, without science or GT.

Uses the existing site build inventory. Reads complete small HTML and prefixes
of figures/tables; does not redownload released multi-GB scientific archives.
"""
from __future__ import annotations
import argparse,json,time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request,urlopen


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--base',default='https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r15/')
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();root=a.root.resolve();cpu=time.process_time();wall=time.perf_counter()
    build=json.loads((root/'SITE_BUILD_RECEIPT.json').read_text())
    paths=set(build['generated_pages'])
    for group in build['gallery_inventory']:
        for name in group['figures']:
            paths.add(name)
            pdf=Path(name).with_suffix('.pdf').as_posix()
            if (root/pdf).is_file():paths.add(pdf)
    paths.update(['CHAT_HANDOFF.md','RESEARCH_REPORT_ZH.md','THEORY_DERIVATIONS.md',
                  'PRO_PACKET_CONTENTS.md','DOWNLOADS.json','src/twofold_solver.py',
                  'analysis/twofold_exposed24_v1/PRIMARY_MATRIX24.csv',
                  'analysis/twofold_exposed24_v1/ALL_VIEWS.csv',
                  'analysis/twofold_exposed24_v1/RECOVERY_GAIN_LOSS.csv'])
    def one(name):
        local=root/name;limit=1_000_000 if local.suffix=='.html' else 256
        row=dict(path=name,url=a.base+quote(name,safe='/_.-'),local_bytes=local.stat().st_size,authentication='NONE')
        try:
            req=Request(row['url'],headers={'User-Agent':'MAX-STRUCT-R15-publication-audit','Accept-Encoding':'identity'})
            with urlopen(req,timeout=45) as response:
                data=response.read(limit);status=response.status;headers={k.lower():v for k,v in response.headers.items()}
            with local.open('rb') as f:prefix=f.read(limit)
            row.update(HTTP_status=status,read_bytes=len(data),prefix_matches_local=data==prefix,
                       Content_Length=headers.get('content-length'),content_type=headers.get('content-type'))
            row['verified']=status==200 and data==prefix
        except Exception as e:row.update(verified=False,error=repr(e))
        return row
    with ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(one,sorted(paths)))
    failures=[r for r in rows if not r['verified']]
    receipt=dict(status='ANONYMOUS_PAGES_READBACK_COMPLETE' if not failures else 'ACCESS_OR_DEPLOYMENT_MISMATCH_RETAINED',
                 created_utc=datetime.now(timezone.utc).isoformat(),base=a.base,files=len(rows),rows=rows,failures=failures,
                 CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
                 scope='Full small HTML and original first256 bytes of figures/PDF/CSV/code, not full arrays or scientific revalidation',
                 visual_check='SEPARATE_PARENT_BROWSER_REVIEW',physics=0,GT_access='NONE')
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open('x') as f:json.dump(receipt,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({k:receipt[k] for k in ('status','files','CPU_seconds','wall_seconds')}))
    return 1 if failures else 0


if __name__=='__main__':raise SystemExit(main())
