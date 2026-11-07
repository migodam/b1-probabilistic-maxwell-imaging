"""Wait for main public seals, then zero-FW endpoint metadata with external CPU."""
import json,time,subprocess,sys,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'delegated/twod/independent48'
def main():
 out=D/'endpoint_bookkeeping'
 if out.exists():raise RuntimeError('Preserve existing endpoint attempt; do not overwrite')
 (D/'ENDPOINT_WAITER.json').write_text(json.dumps(dict(pid=os.getpid(),truth_access='NONE',status='WAITING_FOR_ALL_MAIN_PUBLIC_SEALS'))+'\n')
 seals=[D/p for p in ('secondary/public_seal/RECEIPT.json','secondary/global_appendix/public_seal/RECEIPT.json','actual_cpu_fairness/public_seal/RECEIPT.json')]
 while not all(p.exists() for p in seals):time.sleep(60)
 subprocess.run([sys.executable,str(ROOT/'experiments/twod_postprocess_run.py'),'--kind','endpoint-bookkeeping','--','--manifest',str(ROOT/'configs/twod_independent48_all_tasks.json'),'--public-seal',str(seals[0]),'--output',str(out)],check=True)
 subprocess.run([sys.executable,str(ROOT/'experiments/twod_independent_costs.py')],check=True)
 print('PUBLIC_ONLY_ENDPOINT_BOOKKEEPING_COMPLETE_NO_GT',flush=True)
if __name__=='__main__':main()
