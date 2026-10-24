"""Explicit metadata-only repair after generator COMPLETE, no inference replay."""
from pathlib import Path
import json,sys,subprocess
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,Book,write
from p3r4.evaluation import checkpoint_identity

def main():
    root=OUT/'results/postfit_dev_v2';cfg=OUT/'configs'
    if json.loads((root/'GENERATOR_COMPLETE.json').read_text())['status']!='COMPLETE':raise ValueError('closed generator required')
    if (root/'paid_probability').exists():raise ValueError('metadata repair only before ANY physical phase/folder')
    original=cfg/'POSTFIT_DEV_V2_PAID_RELEASE.json';bad=json.loads(original.read_text())
    good=checkpoint_identity(root/'FINAL_PREINFERENCE_SPEC.json')
    if bad.get('preinference_spec')!=good or bad['final_spec']['bytes'] is not None:raise ValueError('unexpected release failure')
    log=OUT/'logs/postfit_execution.log';receipt=OUT/'I/jobs/postfit_dev_v2_score/CLOSED.json'
    elapsed=(log.stat().st_mtime_ns-receipt.stat().st_mtime_ns)/1e9
    if not 0<=elapsed<15 or 'paid release must bind final preinference spec' not in log.read_text():raise ValueError('reconcile actual pre-Book failure first')
    book=Book('I','postfit_release_metadata_reconciliation',analysis_only=True,max_wall=60)
    book._update('FAILED_COMMAND_OVERHEAD_UPPER_BOUND',controlled=elapsed,
        physical_actions=0,scope='after last successful scoring receipt to pre-Book phase guard failure; conservative bound, no double counting child inference',source_log=str(log))
    recovery=dict(reason='parent coordinator wrote correct final identity under wrong metadata key; strict guard rejected before Book/physics',
        original_release=checkpoint_identity(original),new_physics_RHS=0,new_inference_calls=0,
        failed_guard_wall_upper_bound=elapsed,all_original_sources_and_samples_unchanged=True)
    write(OUT/'AUDIT/POSTFIT_RELEASE_BINDING_RECOVERY.json',recovery);book.close()
    for phase,key in [('paid_probability','PAID'),('corrected_MH','MH')]:
        p=cfg/f'POSTFIT_DEV_V2_{key}_RELEASE_BINDING_FIXED.json'
        if p.exists():raise FileExistsError('no automatic physical replay')
        r=json.loads((cfg/f'POSTFIT_DEV_V2_{key}_RELEASE_TEMPLATE.json').read_text())
        r.update(released=True,final_spec=good,metadata_recovery=recovery)
        p.write_text(json.dumps(r,indent=2)+'\n')
        q=subprocess.run([sys.executable,'-u','tools/r4_postfit_dev.py',phase,'--release',str(p)],cwd=ROOT)
        if q.returncode:raise SystemExit(q.returncode)
    write(OUT/'results/POSTFIT_DEV_V2_EXECUTION_COMPLETE.json',dict(status='COMPLETE',generator_NOT_replayed=True,
        scope='development diagnostics only; no qualified posterior reference'))

if __name__=='__main__':main()
