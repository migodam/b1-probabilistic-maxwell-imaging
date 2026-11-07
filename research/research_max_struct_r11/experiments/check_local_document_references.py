"""Check local Markdown links in R11 public-facing root docs and absorptive handoff."""
from pathlib import Path
from urllib.parse import unquote,urlsplit
import json,re,time
ROOT=Path(__file__).absolute().parents[1]
def main():
 docs=sorted(ROOT.glob('*.md'))+[ROOT/'delegated/diagnostics/absorptive_dev/FACTUAL_SUMMARY.md']
 rows=[]
 for doc in docs:
  for match in re.finditer(r'\[[^\]]*\]\((<[^>]*>|[^\s)]+)(?:\s+"[^"]*")?\)',doc.read_text()):
   raw=match.group(1).strip('<>');url=urlsplit(raw)
   if url.scheme or url.netloc or not url.path:continue
   p=Path(unquote(url.path));target=p if p.is_absolute() else doc.parent/p
   rows.append(dict(document=str(doc.relative_to(ROOT)),line=doc.read_text()[:match.start()].count('\n')+1,target=raw,resolved=str(target.absolute()),exists=target.exists()))
 sources=json.loads((ROOT/'FIGURES/absorptive_dev/MATERIAL_SOURCE_TABLE.json').read_text())
 raw_missing=[r for r in sources if not r['raw_material_path'] or not Path(r['raw_material_path']).exists()]
 missing=[r for r in rows if not r['exists']]
 result=dict(status='PASS' if not missing and not raw_missing else 'MISSING_REFERENCES_REPORTED',scope='Local Markdown link file/directory existence only, root R11 docs plus absorptive factual handoff; no remote URL or anchor validation',documents=len(docs),local_links=len(rows),missing=missing,raw_material_paths_checked=len(sources),missing_raw_material_paths=raw_missing,rows=rows,new_solver=0,new_FW=0,GT_arrays_read=0,CPU_seconds_near_write=time.process_time(),CPU_scope='Own process_time lower bound, no external OS measurement')
 (ROOT/'delegated/diagnostics/absorptive_dev/LOCAL_DOCUMENT_REFERENCE_CHECK.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k!='rows'}))
if __name__=='__main__':main()
