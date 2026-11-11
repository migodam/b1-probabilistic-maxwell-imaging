"""Restore selected R12 logical archives from published transport parts.

Only restores files; it never reruns experiments. The output directory is explicit
so callers can select an external disk. Requires zstd only when --decode is used.
"""
import argparse,json,subprocess,urllib.request
from pathlib import Path

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--manifest',type=Path,default=Path(__file__).with_name('PUBLICATION_PARTS_MANIFEST.json'))
 p.add_argument('--asset',action='append',required=True,help='Exact logical archive name; may be repeated')
 p.add_argument('--output-dir',type=Path,required=True)
 p.add_argument('--decode',action='store_true',help='Decode .zst into ZIP; retain the verified logical wrapper')
 a=p.parse_args();m=json.loads(a.manifest.read_text());selected={x['name']:x for x in m['logical_assets']}
 a.output_dir.mkdir(parents=True,exist_ok=True)
 for name in a.asset:
  x=selected[name];dest=a.output_dir/name
  if not dest.exists() or dest.stat().st_size!=x['bytes']:
   partial=dest.with_name(dest.name+'.download')
   with partial.open('wb') as out:
    for part in x['parts']:
     count=0
     with urllib.request.urlopen(part['url']) as response:
      while block:=response.read(8*1024*1024):out.write(block);count+=len(block)
     if count!=part['bytes']:raise RuntimeError(f'Incomplete part: {part["name"]}')
   if partial.stat().st_size!=x['bytes']:raise RuntimeError(f'Incomplete archive: {name}')
   partial.replace(dest)
  print(dest,flush=True)
  if a.decode and name.endswith('.zst'):
   target=dest.with_suffix('')
   subprocess.run(['zstd','-d','-f',str(dest),'-o',str(target)],check=True)
   print(target,flush=True)
if __name__=='__main__':main()
