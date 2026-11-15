"""Anonymous readback of actual generated pages and gallery image/PDF headers.

Uses the existing site inventory, without scientific arrays or hashes. Header
checks establish transport/access, not visual interpretation or physics validity.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path
import struct
import time
from urllib.parse import quote
from urllib.request import Request, urlopen


def verify(row):
    result={k:v for k,v in row.items() if k!="local"}
    try:
        request=Request(row["url"],headers={"User-Agent":"MAX-STRUCT-public-Pages-audit",
            "Accept-Encoding":"identity",**({"Range":"bytes=0-63"} if row["kind"]=="media" else {})})
        with urlopen(request,timeout=45) as response:
            body=response.read(64) if row["kind"]=="media" else response.read()
            result.update(status=response.status,content_type=response.headers.get("Content-Type"),
                bytes_read=len(body),header_length=response.headers.get("Content-Length"))
        if row["kind"]=="media":
            with row["local"].open("rb") as source:local=source.read(64)
            result["prefix_matches_original"]=body==local
            if row["path"].endswith(".png") and body.startswith(b"\x89PNG\r\n\x1a\n"):
                result["PNG_dimensions"]=list(struct.unpack(">II",body[16:24]))
            if row["path"].endswith(".pdf"):result["PDF_signature"]=body.startswith(b"%PDF-")
            result["verified"]=result["status"] in (200,206) and body==local
        else:
            text=body.decode("utf-8",errors="replace")
            result["verified"]=result["status"]==200 and ("<html" in text.lower() or not row["path"].endswith(".html"))
            if row["path"]=="index.html":
                result["report_markers"]={marker:marker in text for marker in ("MAX-STRUCT R14","r14_new_0021","已核验发布附件")}
                result["verified"]=result["verified"] and all(result["report_markers"].values())
    except Exception as error:
        result.update(verified=False,error=repr(error))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument("--base-url",default="https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r14/")
    parser.add_argument("--label",default="PUBLIC_PAGES_READBACK")
    args=parser.parse_args();root=args.root.resolve();cpu,wall=time.process_time(),time.perf_counter()
    inventory=json.loads((root/"SITE_BUILD_RECEIPT.json").read_text())
    rows=[]
    def add(path,kind):rows.append({"path":path,"kind":kind,"url":args.base_url+quote(path,safe="/._-"),"local":root/path})
    for path in inventory["generated_pages"]:add(path,"page")
    for group in inventory["gallery_inventory"]:
        for path in group["figures"]:add(path,"media")
    for path in ("CHAT_HANDOFF.md","PRO_HANDOFF.md","RESEARCH_REPORT_ZH.md","DOWNLOADS.json",
        "src/dualmetric_core.py","configs/independent24_v1.json"):add(path,"source")
    with ThreadPoolExecutor(max_workers=8) as executor:results=list(executor.map(verify,rows))
    failures=[r for r in results if not r["verified"]]
    receipt=dict(status="ALL_REQUESTED_PUBLIC_PAGES_AND_MEDIA_READ_BACK" if not failures else "PUBLIC_ACCESS_FAILURES_RETAINED",
        created_utc=datetime.now(timezone.utc).isoformat(),base_url=args.base_url,items=len(results),
        generated_pages=len(inventory["generated_pages"]),gallery_items=sum(len(g["figures"]) for g in inventory["gallery_inventory"]),
        results=results,failures=failures,authentication="NONE",physics=0,scoring=0,hash_operations="NONE",
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        visual_QA="Separate parent inspection; header matching does not establish full visual review")
    out=root/"audit/publication";out.mkdir(exist_ok=True,parents=True)
    (out/(args.label+".json")).write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({k:receipt[k] for k in ("status","items","generated_pages","gallery_items","CPU_seconds","wall_seconds")}))
    return 1 if failures else 0


if __name__=="__main__":raise SystemExit(main())
