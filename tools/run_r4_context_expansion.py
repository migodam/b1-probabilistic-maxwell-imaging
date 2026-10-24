"""One-off chunked context expansion under the common physical lock/ledger."""
from pathlib import Path
import subprocess,sys,json
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r4_conditional_transport'
def run(args,label):
 p=subprocess.run([sys.executable,*args],cwd=ROOT)
 print(json.dumps(dict(task=label,returncode=p.returncode)),flush=True)
 if p.returncode:raise SystemExit(p.returncode)
def main():
 a=json.loads((OUT/'configs/CONTEXT_EXPANSION_01.json').read_text())
 if a['status']!='RELEASED_DEVELOPMENT_ONLY':raise ValueError('expansiongate')
 for group,n in [('train',a['train']),('validation',a['validation'])]:
  for start in range(0,n,a['chunk']):
   # A separately frozen physical decision experiment may run between chunks.
   # It owns its own Book and is never hidden inside context costs.
   priority=OUT/'configs/INTERVAL_CHAIN_RELEASE.json'
   if priority.exists() and not (OUT/'results/II/interval_chain_pair/PAIR_COMPLETE.json').exists():
    run(['tools/run_r4_interval_chain.py','--release',str(priority)],'registered_interval_chain_pair')
   run(['tools/run_r4_contexts.py',group,'--start',str(start),'--count',str(min(a['chunk'],n-start))],f'contexts_{group}_{start}')
 (OUT/'results/CONTEXT_EXPANSION_COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',train=a['train'],validation=a['validation'],GT_read=False),indent=2)+'\n')
if __name__=='__main__':main()
