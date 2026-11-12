"""Check local HTML/Markdown targets among allowlisted publication documents.

No network access, file changes, publishing or scientific reads. Optional strict
Git target membership; existing links outside R13 are repository dependencies.
"""
import argparse,json,re,time
from pathlib import Path
from urllib.parse import unquote,urlsplit
from html.parser import HTMLParser
class Links(HTMLParser):
 def __init__(self):super().__init__();self.targets=[]
 def handle_starttag(self,tag,attrs):
  self.targets +=[v for k,v in attrs if k in ['href','src'] and v]
def main():
 c,w=time.process_time(),time.perf_counter();p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--allowlist',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--strict-allowlist',action='store_true');a=p.parse_args();root=a.root.resolve();allowed=set(json.loads(a.allowlist.read_text())['git_allowlist']);issues=[];external=[];checked=0;membership=[]
 for name in sorted(allowed):
  f=root/name
  if f.suffix not in ['.md','.html'] or not f.is_file():continue
  text=f.read_text(errors='replace')
  if f.suffix=='.html':parser=Links();parser.feed(text);targets=parser.targets
  else:targets=[x[1] or x[2] for x in re.findall(r'(!?)\[[^\]]*\]\((?:<([^>]+)>|([^\s)]+))(?:\s+"[^"]*")?\)',text)]
  for target in targets:
   target=target.strip();url=urlsplit(target)
   if url.scheme or url.netloc or not url.path:continue
   path=unquote(url.path)
   if path.startswith('/'):
    issues.append(dict(source=name,target=target,reason='absolute local/root path not portable'));continue
   dest=(f.parent/path).resolve();checked+=1
   if not dest.exists():issues.append(dict(source=name,target=target,reason='missing local target'));continue
   if not dest.is_relative_to(root):external.append(dict(source=name,target=target,reason='existing external repository dependency'));continue
   rel=dest.relative_to(root).as_posix()
   if dest.is_file() and rel not in allowed:membership.append(dict(source=name,target=target,reason='target exists but is not Git allowlisted'))
 report=dict(status='PASS' if not issues and not membership else 'LOCAL_TARGETS_EXIST_BUT_PUBLICATION_MEMBERSHIP_PENDING' if not issues and not a.strict_allowlist else 'ISSUES_RETAINED',checked_local_links=checked,broken_or_absolute=issues,not_allowlisted=membership,existing_external_repository_dependencies=external,strict_allowlist=a.strict_allowlist,CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w,host='Apple_local',network_checks=False,science_or_GT_read=False)
 a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['status','checked_local_links','CPU_seconds']}));print(json.dumps(dict(broken_or_absolute=len(issues),not_allowlisted=len(membership))))
if __name__=='__main__':main()
