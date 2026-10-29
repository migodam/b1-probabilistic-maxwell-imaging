"""Anonymous streaming publication verification; no extra local raw-asset copies."""
from pathlib import Path
import argparse,json,urllib.request,hashlib,time,concurrent.futures
P=Path(__file__).resolve().parents[1];BASE='https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/';RAW='https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/'
def get(url,expected=None):
 start=time.monotonic();h=hashlib.blake2b(digest_size=32);n=0
 with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'R9-anonymous-evidence-readback'}),timeout=120) as r:
  code=r.status;ctype=r.headers.get('Content-Type');declared=r.headers.get('Content-Length')
  for b in iter(lambda:r.read(8<<20),b''):n+=len(b);h.update(b)
 rec=dict(url=url,status=code,bytes=n,blake2b_256=h.hexdigest(),content_type=ctype,declared_bytes=declared,wall_seconds=time.monotonic()-start,anonymous=True)
 if expected:rec['matches_expected']=n==expected['archive_bytes'] and h.hexdigest()==expected['blake2b_256']
 return rec
def assets():
 m=json.loads((P/'EVIDENCE_MANIFEST.json').read_text());out=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
  jobs={pool.submit(get,a['url'],a):a['file'] for a in m['assets']}
  for j in concurrent.futures.as_completed(jobs):
   name=jobs[j]
   try:rec=dict(file=name,**j.result())
   except Exception as e:rec=dict(file=name,error=repr(e),matches_expected=False)
   out.append(rec);(P/'packaging/PUBLIC_ASSET_READBACK.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(rec),flush=True)
 return out
def pages():
 paths=['p3l_r9_stabilized_subspace/','p3l_r9_stabilized_subspace/EVIDENCE_PACKAGE.html','p3l_r9_stabilized_subspace/GO_NO_GO.html','p3l_r9_stabilized_subspace/GPT_PRO_HANDOFF/START_HERE.html','maxwell_formula_atlas/','maxwell_formula_atlas/FORMULA_ATLAS.html','maxwell_formula_atlas/08_TASK_COVERAGE_AND_RANDOMIZED_OBJECTIVES.html']
 out=[]
 for x in paths:out.append(get(BASE+x))
 for x in ['p3l_r9_stabilized_subspace/RESEARCH_REPORT_ZH.md','p3l_r9_stabilized_subspace/EVIDENCE_MANIFEST.json','maxwell_formula_atlas/FORMULA_ATLAS.md']:
  rec=get(RAW+x);f=P.parents[1]/'research'/x;h=hashlib.blake2b(f.read_bytes(),digest_size=32).hexdigest();rec['matches_local']=rec['bytes']==f.stat().st_size and rec['blake2b_256']==h;out.append(rec)
 (P/'packaging/PUBLIC_PAGE_READBACK.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True);return out
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('part',choices=['assets','pages']);v=a.parse_args();globals()[v.part]()
