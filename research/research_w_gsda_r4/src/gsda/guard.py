"""Online-only file-access guard: no truth, glyph bank, or scoring reference."""
from pathlib import Path
import sys,os,json,time
class OnlineGuard:
 def __init__(self,root):
  self.root=Path(root).resolve();self.count=0;self.denied=[];self.cache={};self.hook_wall_seconds=0.
  self.blocked=tuple(str(self.root/x) for x in ('inputs/private','evaluation','results/absorption_controls','results/spectral','results/independent_validation'))
 def hook(self,event,args):
  if event!='open' or not args:return
  x=args[0]
  if not isinstance(x,(str,bytes,os.PathLike)):return
  started=time.perf_counter();key=os.fsdecode(x)
  try:
   p=self.cache.get(key)
   if p is None:p=str(Path(key).resolve());self.cache[key]=p
  except Exception:return
  self.hook_wall_seconds+=time.perf_counter()-started;self.count+=1
  if any(p==b or p.startswith(b+os.sep) for b in self.blocked) or 'DIGIT_PATTERNS_1_TO_9' in p or p.endswith(('/digit_patterns.py','/gsda/generate.py','/gsda/score.py','/tools/generate_inputs.py','/tools/score_all.py')):
   self.denied.append(p);raise PermissionError('Online reconstruction cannot read private truth or glyph bank: '+p)
 def install(self):sys.addaudithook(self.hook);return self
 def save(self,path):
  Path(path).write_text(json.dumps(dict(time_unix=time.time(),online_opens=self.count,denied=self.denied,unique_resolved_paths=len(self.cache),measured_hook_wall_seconds=self.hook_wall_seconds,cache_assumption='online program does not mutate symlink targets',scope='Python open audit; no OS isolation, not a proof against malicious code'),indent=2)+'\n')
