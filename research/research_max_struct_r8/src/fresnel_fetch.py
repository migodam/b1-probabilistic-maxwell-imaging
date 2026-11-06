"""Acquire original official Fresnel files without redistributing measurements.

An HTTP 200 challenge is explicitly an access failure, never data. Signed
publisher links are regenerated from the canonical supplementary page.
"""
from pathlib import Path
import urllib.request,urllib.error,re,html,json,time,argparse,zipfile,io
ROOT=Path(__file__).resolve().parents[1]
CAMPAIGNS={'2001':('https://iopscience.iop.org/article/10.1088/0266-5611/17/6/301/data',['dielTM_dec8f.exp','twodielTM_8f.exp']),
 '2005':('https://iopscience.iop.org/article/10.1088/0266-5611/21/6/S09/data',['FoamTwinDielTM.exp'])}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--campaign',choices=CAMPAIGNS,default='2005');args=ap.parse_args()
 page,names=CAMPAIGNS[args.campaign];start=time.perf_counter();rows=[];out=ROOT/'data_local/fresnel';out.mkdir(parents=True,exist_ok=True)
 try:
  text=urllib.request.urlopen(page,timeout=30).read().decode()
  if 'captcha' in text.lower():raise RuntimeError('OFFICIAL_ACCESS_BLOCKED_CAPTCHA; no measurements acquired')
  links=[html.unescape(x) for x in re.findall(r'href=["\']([^"\']+)["\']',text)]
  archive=None;archive_url=None
  if args.campaign=='2005':
   zlinks=[x for x in links if '/supp.zip?' in x or x.endswith('/supp.zip')]
   if zlinks:
    archive_url=zlinks[0];archive=zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(archive_url,timeout=30).read()))
  for name in names:
   if archive is not None:
    members=[x for x in archive.namelist() if Path(x).name==name]
    if len(members)!=1:raise RuntimeError('Original file missing/ambiguous in official supp.zip: '+name)
    content=archive.read(members[0])
    if not content.startswith(b'#'):raise RuntimeError('Original ASCII header absent in official archive: '+name)
    (out/name).write_bytes(content);rows.append(dict(file=name,status='DOWNLOADED_FROM_OFFICIAL_ARCHIVE',bytes=len(content),official_page=page,url_path=archive_url.split('?')[0],archive_member=members[0]));continue
   matches=[x for x in links if '/'+name+'?' in x or x.endswith('/'+name)]
   if not matches:raise RuntimeError('Original file link absent from official page: '+name)
   content=urllib.request.urlopen(matches[0],timeout=30).read()
   if not content.startswith(b'#'):raise RuntimeError('Original ASCII header absent: '+name)
   (out/name).write_bytes(content);rows.append(dict(file=name,status='DOWNLOADED',bytes=len(content),official_page=page,url_path=matches[0].split('?')[0]))
 except Exception as exc:rows.append(dict(status='FAILED',official_page=page,error=str(exc)))
 (ROOT/'delegated/fresnel'/f'DOWNLOAD_{args.campaign}_RECEIPT.json').write_text(json.dumps(dict(attempts=rows,wall_seconds=time.perf_counter()-start),indent=2)+'\n')
 print(json.dumps(rows))
if __name__=='__main__':main()
