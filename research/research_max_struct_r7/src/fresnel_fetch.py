"""Download the four public 2001 dielectric TM files from their official page.
Does not solve challenges, use credentials, or compute checksums.
"""
from pathlib import Path
import re,html,json,urllib.request
PAGE='https://iopscience.iop.org/article/10.1088/0266-5611/17/6/301/data'
def main():
 out=Path(__file__).resolve().parents[1]/'data_public'/'fresnel';out.mkdir(parents=True,exist_ok=True)
 text=urllib.request.urlopen(PAGE,timeout=30).read().decode();rows=[]
 for name in ['dielTM_dec4f.exp','dielTM_dec8f.exp','twodielTM_4f.exp','twodielTM_8f.exp']:
  links=re.findall(r'href=["\']([^"\']+)["\']',text)
  url=next(html.unescape(x) for x in links if '/'+name+'?' in x)
  content=urllib.request.urlopen(url,timeout=30).read()
  if not content.startswith(b'#'):raise RuntimeError('Expected original ASCII measurement header')
  (out/name).write_bytes(content);rows.append(dict(file=name,bytes=len(content),official_page=PAGE,url_path=url.split('?')[0]))
 (out/'DOWNLOAD_RECEIPT.json').write_text(json.dumps(rows,indent=2)+'\n')
if __name__=='__main__':main()
