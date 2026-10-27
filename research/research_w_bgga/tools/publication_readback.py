"""Anonymous HTTP readback only; no credentials, physics or file digests."""
from pathlib import Path
import json,subprocess,time,urllib.request,urllib.error
ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
def main():
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip()
    prefix='https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/'+commit+'/research/research_w_bgga/'
    urls=dict(raw_report=prefix+'RESEARCH_W_BGGA_REPORT_ZH.md',
        resource_plan=prefix+'RESOURCE_REQUEST.md',
        repository='https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/'+commit+'/research/research_w_bgga',
        page='https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_w_bgga/')
    rows=[]
    for kind,url in urls.items():
        row=dict(kind=kind,url=url,anonymous=True)
        for attempt in range(6 if kind=='page' else 1):
            try:
                request=urllib.request.Request(url,headers={'User-Agent':'Research-W-BGGA-public-readback'})
                with urllib.request.urlopen(request,timeout=20) as response:
                    data=response.read();status=response.status
                content=data.decode('utf-8',errors='replace')
                if kind in ('raw_report','page') and '266,326' not in content:raise ValueError('current BGGA report marker missing')
                if kind=='raw_report' and data!=(ROOT/'RESEARCH_W_BGGA_REPORT_ZH.md').read_bytes():raise ValueError('raw report differs from local source')
                row.update(status=status,bytes=len(data),passed=status==200,attempts=attempt+1);break
            except Exception as exc:
                row.update(passed=False,error=repr(exc),attempts=attempt+1)
                if kind=='page' and attempt<5:time.sleep(10)
        rows.append(row)
    result=dict(source_commit=commit,status='PASS' if all(x['passed'] for x in rows) else 'PARTIAL',
        checks=rows,no_sha256_check=True,scientific_gates_unchanged=True)
    (ROOT/'PUBLICATION_READBACK.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False))
if __name__=='__main__':main()
