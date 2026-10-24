"""Parent-released fixed DEV diagnostics; no training or automatic replay."""
from pathlib import Path
import json,subprocess,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.evaluation import checkpoint_identity
OUT=ROOT/'research/p3l_r4_conditional_transport'

def call(phase,path):
    r=subprocess.run([sys.executable,'-u','tools/r4_postfit_dev.py',phase,'--release',str(path)],cwd=ROOT)
    if r.returncode:raise SystemExit(r.returncode)

def main():
    # This wrapper is the parent's concrete release, not an adaptive search.
    cfg=OUT/'configs';name='POSTFIT_DEV_V2_GENERATOR_RELEASE.json'
    path=cfg/name
    if path.exists():raise FileExistsError('no replay of released sequence')
    result=json.loads((OUT/'results/physics_transport_paired_v2/CONTROLLER_RESULT.json').read_text())
    if result['status']!='COMPLETE':raise ValueError('all actual paired jobs closed first')
    r=json.loads((cfg/'POSTFIT_DEV_V2_GENERATOR_RELEASE_TEMPLATE.json').read_text())
    r['released']=True
    r['MH_validation_binding']=dict(method_id='raw_spline__physics',validation_only=True,before_DEV_inference=True,
        reason='Among the two raw512 physics-refined models on identical VAL2048/fixed32 Full roster, rawSpline has smaller frozen composite430.8349 versus613.8660. Anchor VAL256 is not pooled. This chooses a diagnostic proposal, not a qualified posterior; all3 family replication gates FAILED and remain failed.')
    path.write_text(json.dumps(r,indent=2)+'\n')
    call('generator',path)
    identity=checkpoint_identity(OUT/'results/postfit_dev_v2/FINAL_PREINFERENCE_SPEC.json')
    for phase,template,key in [('paid_probability','PAID','PAID'),('corrected_MH','MH','MH')]:
        src=cfg/f'POSTFIT_DEV_V2_{template}_RELEASE_TEMPLATE.json'
        p=cfg/f'POSTFIT_DEV_V2_{key}_RELEASE.json'
        if p.exists():raise FileExistsError(p)
        q=json.loads(src.read_text());q['released']=True
        q['final_spec']=identity
        p.write_text(json.dumps(q,indent=2)+'\n')
        call(phase,p)
    (OUT/'results/POSTFIT_DEV_V2_EXECUTION_COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',
        scope='exposed development calibration and valid finite-chain diagnostic only',new_locked_GT=False))+'\n')

if __name__=='__main__':main()
