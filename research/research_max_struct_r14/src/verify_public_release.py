"""Read an already published release anonymously and write actual download links.

Presentation/publication audit only. No credentials, scientific arrays, scoring,
hashes or optimizer are used. Large assets are opened with a small Range request.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import time
from urllib.request import Request, urlopen


def read_request(url, limit=None):
    headers={"User-Agent":"MAX-STRUCT-public-access-audit", "Accept-Encoding":"identity"}
    if limit is not None: headers["Range"]=f"bytes=0-{limit-1}"
    with urlopen(Request(url, headers=headers), timeout=60) as response:
        data=response.read() if limit is None else response.read(limit)
        return dict(status=response.status, headers=dict(response.headers), data=data)


def verify_one(asset):
    row={"name":asset["name"],"url":asset["browser_download_url"],"bytes":asset["size"],
         "authentication":"NONE","verified":False}
    try:
        result=read_request(row["url"],64)
        headers={k.lower():v for k,v in result["headers"].items()}
        content_range=headers.get("content-range")
        match=re.search(r"/(\d+)$",content_range or "")
        actual=int(match.group(1)) if match else int(headers.get("content-length","-1"))
        prefix=result["data"]
        zip_ok=prefix.startswith(b"PK\x03\x04") if row["name"].lower().endswith(".zip") else True
        row.update(http_status=result["status"],content_range=content_range,
                   header_total_bytes=actual,bytes_read=len(prefix),zip_header_ok=zip_ok,
                   verified=result["status"] in (200,206) and actual==row["bytes"] and zip_ok)
    except Exception as error:
        row["error"]=repr(error)
    return row


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument("--repo",default="migodam/b1-probabilistic-maxwell-imaging")
    parser.add_argument("--tag",default="research-max-struct-r14-dualmetric-20261115")
    args=parser.parse_args(); root=args.root.resolve()
    cpu,wall=time.process_time(),time.perf_counter()
    api_url=f"https://api.github.com/repos/{args.repo}/releases/tags/{args.tag}"
    metadata=json.loads(read_request(api_url)["data"])
    if metadata.get("draft") or metadata.get("tag_name")!=args.tag:
        raise RuntimeError("Require actual public release metadata for this tag")
    with ThreadPoolExecutor(max_workers=4) as executor:
        assets=list(executor.map(verify_one,metadata["assets"]))
    failures=[row for row in assets if not row["verified"]]
    receipt=dict(status="ANONYMOUS_READBACK_COMPLETE" if not failures else "ACCESS_FAILURES_RETAINED",
        created_utc=datetime.now(timezone.utc).isoformat(),release_id=metadata["id"],
        release_url=metadata["html_url"],tag=args.tag,public_metadata_url=api_url,
        asset_count=len(assets),bytes=sum(row["bytes"] for row in assets),assets=assets,
        authentication="NONE",failures=failures,new_physics=0,GT_access="NONE",hash_operations="NONE",
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        scope="Anonymous metadata plus each asset header/size and ZIP signature; not a full redownload or scientific revalidation")
    out=root/"audit/publication";out.mkdir(parents=True,exist_ok=True)
    (out/"PUBLIC_ASSET_READBACK.json").write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+"\n")
    if not failures:
        downloads={k:receipt[k] for k in ("created_utc","release_id","release_url","tag","asset_count","bytes","assets","scope")}
        downloads["status"]="PUBLISHED_AND_ANONYMOUSLY_READ_BACK"
        (root/"DOWNLOADS.json").write_text(json.dumps(downloads,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({k:receipt[k] for k in ("status","asset_count","bytes","CPU_seconds","wall_seconds")}))
    return 1 if failures else 0


if __name__=="__main__": raise SystemExit(main())
