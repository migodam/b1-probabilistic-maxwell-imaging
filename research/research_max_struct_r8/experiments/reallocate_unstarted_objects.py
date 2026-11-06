"""Move only wholly unstarted objects to the now-idle local CPU.

The numerical runner and budgets are unchanged. The existing remote finite
dispatcher records explicit relocation acknowledgements instead of executing
those objects. Credentials are command-line inputs and are never archived.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,subprocess,time
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--host',required=True);ap.add_argument('--identity',required=True);a=ap.parse_args()
    candidates=[f'r8_test_{i:04d}' for i in [33,34,35,37,38,39,41]]
    code=r'''from pathlib import Path
from datetime import datetime,timezone
import json,psutil,time
primary=Path('D:/AI/MAX_STRUCT_R8_20261106/research/research_max_struct_r8')
anchor=Path('D:/AI/MAX_STRUCT_R8_ANCHOR_20261106/research/research_max_struct_r8')
candidates=CANDIDATES
seen={p.name.split('__')[0] for p in (primary/'RESULTS/warm_primary').iterdir() if p.is_dir()}
active=set()
for p in psutil.process_iter(['cmdline']):
 c=p.info['cmdline'] or []
 if '--task-id' in c:active.add(c[c.index('--task-id')+1].split('__')[0])
eligible=sorted(set(candidates)-seen-active)
if set(eligible)!=set(candidates):raise RuntimeError('Some proposed objects have started; do not split their hardware: '+str(set(candidates)-set(eligible)))
record=dict(status='WHOLE_OBJECT_RELOCATION_BEFORE_ANY_SOLVE',created_utc=datetime.now(timezone.utc).isoformat(),objects=eligible,source='XINAN',destination='Mac',truth_access='NONE',algorithm_change='NONE; dispatch-only acknowledgement',reason='Mac original285 and other workstreams complete; reuse idle CPU for entire unstarted objects; no outcome selection')
block="""    control=ROOT/'audit/RELOCATE_PUBLIC_OBJECTS.json'
    if control.exists():
        relocation=json.loads(control.read_text())
        row=choose_task(json.loads(resolve(args.manifest or args.public_list).read_text()),args.task_id)
        if row['case_id'] in relocation['objects']:
            out=resolve(args.output_dir or row['output_dir']);out.mkdir(parents=True,exist_ok=False)
            rec=dict(status='RELOCATED_BEFORE_EXECUTION',task_id=args.task_id,case_id=row['case_id'],arm=row['arm'],truth_access='NONE',physical_calls=0,LU=0,complex_RHS=0,CPU_seconds=None,reason='Entire previously unstarted object delegated to Mac; actual scientific receipt is produced there. Interpreter/dispatch overhead not independently metered.')
            dump(out/'RECEIPT.json',rec)
            print(json.dumps(rec),flush=True)
            return 0
"""
for root in (primary,anchor):
 audit=root/'audit';audit.mkdir(exist_ok=True,parents=True)
 path=root/'experiments/run_case.py';old=path.read_text()
 if "control=ROOT/'audit/RELOCATE_PUBLIC_OBJECTS.json'" in old:raise RuntimeError('Already patched')
 (audit/('run_case_before_reallocation_'+('anchor' if root==anchor else 'primary')+'.py')).write_text(old)
 new=old.replace('    receipt,out=execute(args)',block+'    receipt,out=execute(args)')
 if new==old:raise RuntimeError('Runner entrypoint not found')
 temp=path.with_suffix('.py.reallocate');temp.write_text(new);temp.replace(path)
 (audit/'RELOCATE_PUBLIC_OBJECTS.json').write_text(json.dumps(record,indent=2))
# An already launched process would have created its RUNNING receipt before
# numerical imports. Check again after the dispatcher switch; no process killed.
time.sleep(1)
conflicts=[]
for p in (primary/'RESULTS/warm_primary').glob('*/RECEIPT.json'):
 if p.parent.name.split('__')[0] not in eligible:continue
 d=json.loads(p.read_text())
 if d['status']!='RELOCATED_BEFORE_EXECUTION':conflicts.append(d['task_id'])
if conflicts:raise RuntimeError('Relocation race, local work not launched: '+str(conflicts))
print(json.dumps(record))
'''.replace('CANDIDATES',repr(candidates))
    started=time.perf_counter();p=subprocess.run(['ssh','-i',a.identity,'-o','BatchMode=yes','-o','ConnectTimeout=10',a.host,'D:/python/python.exe','-B','-'],input=code,text=True,capture_output=True,timeout=60)
    if p.returncode:raise RuntimeError(p.stdout+'\n'+p.stderr)
    rec=json.loads(p.stdout);objects=set(rec['objects']);tasks=[]
    for manifest in ['warm_all_manifest.json','warm_anchor_confirmation_manifest.json']:
        tasks.extend(t for t in json.loads((ROOT/'configs'/manifest).read_text())['tasks'] if t['case_id'] in objects)
    # Every selected local output must be genuinely absent; a partial remote
    # import or previous local solve would invalidate whole-object reassignment.
    if any((REPO/t['output_dir']).exists() for t in tasks):raise RuntimeError('A relocated output already exists locally')
    path=ROOT/'configs/warm_reallocated_manifest.json';path.write_text(json.dumps(dict(config='research/research_max_struct_r8/configs/solver.json',tasks=tasks),indent=2))
    rec.update(tasks=len(tasks),primary_tasks=sum('/warm_primary/' in t['output_dir'] for t in tasks),anchor_tasks=sum('/warm_anchor_confirmation/' in t['output_dir'] for t in tasks),control_wall_seconds=time.perf_counter()-started,cpu_seconds=None,original_assignment_preserved=True)
    (ROOT/'audit/HARDWARE_REALLOCATION.json').write_text(json.dumps(rec,indent=2));(ROOT/'audit/warm_reallocated_slots.json').write_text(json.dumps(dict(slots=4)))
    print(json.dumps(rec))

if __name__=='__main__':main()
