"""Serial physical release after first probability diagnostics; no GT input."""
from pathlib import Path
import json,sys,subprocess,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r4_conditional_transport'
def call(args):
    p=subprocess.run([sys.executable,*args],cwd=ROOT)
    if p.returncode:raise SystemExit(p.returncode)

def main():
    release=OUT/'configs/GEOMETRY_KERNEL_RELEASE.json'
    if release.exists():raise FileExistsError('physical release already registered')
    release.write_text(json.dumps(dict(released=True,experiment='frozen_geometry_probability_v1',campaign_id='P3LR4-20261025',
        scope='finite-run exact-prior-reversible kernel comparison, not qualifiedFULLreference',
        new_RHS_cap=262144,controlled_wall_cap=1800,earlier_teacher_acquisition_RHS=459712,
        kernels=['isotropic','weak','bracket_enriched','random_enriched'],no_GT=True),indent=2)+'\n')
    request=OUT/'configs/REPRESENTATIVE_DEV_REOPT_CNN.json'
    request.write_text(json.dumps(dict(request_id='DEV_REOPT_CNN',estimator_id='local_cnn_s2026102501',
        predictions_npz=str(OUT/'results/development_risk_first/online/local_cnn_s2026102501.npz'),
        mu_layout='complex_N576',parents=[f'R4DEV_{k:05d}' for k in range(16)],
        outputs=['prediction_raw','prediction_clipped'],frozen_estimator=True,no_GT=True,
        state_cap=16,RHS_cap=256,gradient_mode='source_adjoint8',
        estimator_scope='reoptimized R2 function class, risk action notposterior',
        online_gate=str(OUT/'results/development_risk_first/online/ONLINE_COMPLETE.json')),indent=2)+'\n')
    call(['tools/run_r4_representatives.py','--prepare','--manifest',str(request)])
    while not (OUT/'results/probability_physics_first/COMPLETE.json').exists():time.sleep(10)
    call(['tools/r4_geometry_probability.py','run','--prepared',str(OUT/'AUDIT/geometry_probability/PREPARED_V1_VALIDATED'),
          '--release',str(release),'--max-wall','1800','--out',str(OUT/'results/III/frozen_geometry_probability_v1')])
    call(['tools/run_r4_representatives.py','--run',
        '--prepared-request',str(OUT/'results/II/representative_preparation/DEV_REOPT_CNN/FROZEN_REQUEST.json'),
        '--job','representatives_DEV_REOPT_CNN','--max-wall','1200'])
    (OUT/'results/GEOMETRY_AND_REPRESENTATIVES_COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',no_GT=True),indent=2)+'\n')
if __name__=='__main__':main()
