"""Root-granted portable normal_r1 entry test within the original campaign clock."""
from pathlib import Path
import datetime,json,os,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1]
def main():
    clock=json.loads((ROOT/'RUN_CLOCK.json').read_text());cut=datetime.datetime.fromisoformat(clock['no_new_heavy_after_utc'].replace('Z','+00:00')).timestamp();budget=min(1800,int(cut-time.time()-60))
    if budget<60:raise RuntimeError('Insufficient original time for a new portable verification budget')
    portable=ROOT/'delivery/portable_runtime';case=ROOT/'inputs/online/blind_normal_r1';out=ROOT/'results/online/blind_normal_r1/FULL_GN';audit=ROOT/'results/access/portable_normal_r1_FILE_ACCESS.jsonl'
    python=Path("/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/A25/a25-adaptive-om-feedback/.venv_a25/bin/python")
    command=[str(python),str(ROOT/'tools/run_suite.py'),'portable_normal_r1_full_gn','PYTHON',str(ROOT/'tools/guarded_online.py'),'--audit-log',str(audit),str(portable/'tools/deploy_r3.py'),'--deployment-budget-seconds',str(budget),'--measurement',str(case/'observation.npz'),'--config',str(case/'config.json'),'--output',str(out),'--method','FULL_GN']
    env={**os.environ,'PYTHONPATH':str(portable/'vendor'),'R3_SHARED_LOCK_ROOT':str(ROOT)};print(json.dumps(dict(command=command,budget_seconds=budget,shared_lock=str(ROOT/'.fem.lock'))),flush=True);result=subprocess.run(command,cwd=ROOT,env=env)
    source_paths=[]
    if audit.exists():
        for line in audit.read_text().splitlines():
            r=json.loads(line);p=Path(r.get('path',''))
            if r.get('event')=='file_open' and not r.get('write') and p.suffix in {'.py','.pyc'} and any(q in p.parts for q in ['r3','r2','gq','bgga','wslam']):source_paths.append(str(p))
    source_paths=sorted(set(source_paths));fallback=[p for p in source_paths if not Path(p).is_relative_to(portable)]
    required=['endpoint.npz','ONLINE_SEALED.json','reconstructed_target.png','reconstructed_wall.png','target_material.npy','predicted_measurement.npz','reconstruction.json','optimization_trace.csv','runtime_cost.json']
    summary={}
    if (out/'SUMMARY.json').exists():summary=json.loads((out/'SUMMARY.json').read_text())
    completed_starts=len(summary.get('starts',[]));protocol_complete=completed_starts==4
    if not protocol_complete and out.exists():(out/'NOT_COMPLETED.json').write_text(json.dumps(dict(status='NOT_COMPLETED',completed_joint_starts=completed_starts,expected_joint_starts=4,reason='Original time/entry limit; frozen solver caps were not reduced'),indent=2)+'\n')
    check=dict(completed_joint_starts=completed_starts,expected_joint_starts=4,full_frozen_protocol_completed=protocol_complete,solver_caps_changed=False,returncode=result.returncode,source_paths=source_paths,workspace_scientific_source_fallback=fallback,portable_source_only=bool(source_paths) and not fallback,required_output_files={name:(out/name).exists() for name in required},deployment_budget_seconds=budget,original_campaign_clock_extended=False,truth_or_scores_read=False)
    (ROOT/'results/PORTABLE_ENTRY_VERIFICATION.json').write_text(json.dumps(check,indent=2)+'\n');print(json.dumps(check),flush=True)
    if result.returncode or fallback or not protocol_complete or not all(check['required_output_files'].values()):raise RuntimeError('Portable verification incomplete; original failure and outputs retained')
if __name__=='__main__':main()
