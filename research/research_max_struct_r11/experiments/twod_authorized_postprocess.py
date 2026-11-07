"""Execute root-authorized main48-only postseal pipeline; no supplementary GT."""
import json,sys,subprocess,time
from pathlib import Path
import twod_dev_transport as t
from twod_recovery import queue_paths
D=t.ROOT/'delegated/twod/independent48'
RUNTIME=json.loads((t.ROOT/'configs/twod_local_runtime.json').read_text())['python']
def run(label,args):
 p=subprocess.run(args,capture_output=True,text=True,errors='replace');(D/(label+'_stdout.log')).write_text(p.stdout);(D/(label+'_stderr.log')).write_text(p.stderr)
 print(json.dumps(dict(stage=label,returncode=p.returncode,stdout_tail=p.stdout[-2000:],stderr_tail=p.stderr[-2000:])),flush=True)
 if p.returncode:raise RuntimeError(label+' failed; original attempt preserved')
def wrap(kind,args):return [RUNTIME,str(t.ROOT/'experiments/twod_postprocess_run.py'),'--kind',kind,'--',*map(str,args)]
def main():
 auth=t.ROOT/'audit/INDEPENDENT48_ROOT_POSTSEAL_AUTHORIZATION.json';a=json.loads(auth.read_text())
 if a['status']!='ROOT_EXPLICIT_POSTSEAL_AUTHORIZATION' or a['scientific_predictions_terminal']!=288:raise PermissionError('Explicit root main48 authorization required')
 (D/'ROOT_POSTSEAL_SCORER_APPROVED.json').write_text(json.dumps(dict(status=a['status'],source=str(auth),scope='MAIN48_ONLY; no supplementary paired-grid or3D GT',created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())),indent=2)+'\n')
 if not (t.ROOT/'sealed_truth/independent48/generator_metadata.json').exists():run('authorized_private48_copy',[RUNTIME,str(t.ROOT/'experiments/twod_independent_sync.py'),'--finished','--include-sealed-truth'])
 args=['--manifest',t.ROOT/'configs/twod_independent48_all_tasks.json','--output',D/'scored','--authorize-postseal','ROOT_EXPLICIT_POSTSEAL_AUTHORIZATION','--authorization-note',str(auth)]
 for q in queue_paths():args+=['--queue',q]
 run('primary_score',wrap('primary',args));run('costs_after_primary',[RUNTIME,str(t.ROOT/'experiments/twod_independent_costs.py')])
 run('secondary_reindex',wrap('secondary',['--public-seal',D/'secondary/public_seal','--primary-scores',D/'scored','--output',D/'secondary/scored','--authorize-postscore','ROOT_APPROVED_SECONDARY_POSTSCORE']))
 run('all48_plots',wrap('plots',['--scores',D/'scored','--manifest',t.ROOT/'configs/twod_independent48_all_tasks.json','--truth-dir',t.ROOT/'sealed_truth/independent48','--output',D/'plots','--language','en','--authorize-postscore','ROOT_APPROVED_POSTSCORE_PLOTS']))
 run('costs_after_all_postprocess',[RUNTIME,str(t.ROOT/'experiments/twod_independent_costs.py')]);print('MAIN48_AUTHORIZED_POSTPROCESS_COMPLETE',flush=True)
if __name__=='__main__':main()
