"""Anonymous public Release/HEAD checks, deliberately no file digests."""
from pathlib import Path
import json,urllib.request
ROOT=Path(__file__).resolve().parents[1]
def main():
    prefix='https://api.github.com/repos/migodam/b1-probabilistic-maxwell-imaging/releases/tags/'
    tag='research-w-bgga-20261028'
    request=urllib.request.Request(prefix+tag,headers={'User-Agent':'Research-W-BGGA-anonymous-readback'})
    with urllib.request.urlopen(request,timeout=30) as r: release=json.load(r)
    expected=json.loads((ROOT/'PUBLICATION_ASSETS.json').read_text())['assets']
    rows=[]
    for local in expected:
        asset=next(a for a in release['assets'] if a['name']==local['asset'])
        req=urllib.request.Request(asset['browser_download_url'],method='HEAD',
            headers={'User-Agent':'Research-W-BGGA-anonymous-readback'})
        with urllib.request.urlopen(req,timeout=30) as r:
            length=r.headers.get('Content-Length')
            row=dict(name=asset['name'],url=asset['browser_download_url'],anonymous=True,
                status=r.status,API_bytes=asset['size'],HEAD_bytes=None if length is None else int(length),
                expected_bytes=local['bytes'])
        row['passed']=row['status']==200 and row['API_bytes']==row['expected_bytes'] and row['HEAD_bytes'] in (None,row['expected_bytes'])
        rows.append(row)
    result=dict(status='PASS' if all(r['passed'] for r in rows) else 'PARTIAL',tag=tag,
        release_url=release['html_url'],target_commitish=release['target_commitish'],assets=rows,
        no_file_digest_checks=True,new_physics=False)
    (ROOT/'ASSET_READBACK.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
if __name__=='__main__':main()
