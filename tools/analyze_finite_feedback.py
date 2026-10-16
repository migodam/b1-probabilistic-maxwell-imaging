"""Postprocess saved paid arrays only. No new backend import or solver call."""
from pathlib import Path
import sys,json,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from finite_feedback.common import OUT,dump,mass,pack,rms
from lxz.local2d import MaterialPrior

def main():
    start=time.perf_counter();m=json.loads((OUT/'METRICS.json').read_text());rows=[]
    for key in ['250033_weak','250034_ring']:
        geo=np.load(OUT/f'results/{key}_CHART.npz');J=np.load(OUT/f'results/{key}_JACOBIAN.npz')
        base=np.load(OUT/f'results/endpoints/{key}_base.npz');a=np.load(OUT/f'inputs/{key}.npz');sigma=float(a['sigma'])
        for t in [.25,.5,.75,1.]:
            ep=np.load(OUT/f'results/endpoints/{key}_weak_t{t:g}.npz');h=ep['chi']-base['chi'];hlin=t*geo['D']@geo['v']
            actual=pack(ep['prediction']-base['prediction'])*np.sqrt(2)/sigma
            affine=J['J_white_mass']@mass(h);tangent=geo['A']@geo['v']*t
            rows.append(dict(anchor=key,t=t,material_RMS=rms(h),linear_chart_material_RMS=rms(hlin),
                material_secant_tangent_cosine=float(np.vdot(h,hlin).real/np.linalg.norm(h)/np.linalg.norm(hlin)),
                material_decoder_remainder_RMS=rms(h-hlin),frozen_chart_response=float(np.linalg.norm(tangent)),
                actual_material_increment_J_response=float(np.linalg.norm(affine)),actual_full_response=float(np.linalg.norm(actual)),
                frozen_chart_response_relative_error=float(np.linalg.norm(tangent-actual)/np.linalg.norm(actual)),
                full_to_frozen_chart_amplitude_ratio=float(np.linalg.norm(actual)/np.linalg.norm(tangent)),
                Maxwell_only_linearization_error=float(np.linalg.norm(affine-actual)/np.linalg.norm(actual))))
    old=ROOT.parents[1]/'A17/N3_APPLICATION_R1/raw/p250065_weak_a1_snr20_shared.npz'
    oldJ=np.load(old)['J_full'];newJ=np.load(OUT/'results/N3_JACOBIAN.npz')['J_white_mass']
    metrics=dict(new_physics=0,new_RHS=0,role='posthoc decomposition of already frozen numerical diagnostics',
        chart_vs_Maxwell=rows,N3_J_vs_independent_historical_relative_error=float(np.linalg.norm(newJ-oldJ)/np.linalg.norm(oldJ)),
        rational_worst_relative_error=max(r['rational_relative_error'] for r in m['feedback']),
        first_feedback_vs_receiver_J_worst_relative_error=max(r['first_term_J_consistency'] for r in m['feedback']),
        max_saved_backward_residual=max(r['backward_residual'] for r in json.loads((OUT/'results/ENDPOINTS.json').read_text())),
        new_endpoint_arrays=len(list((OUT/'results/endpoints').glob('*.npz'))),
        prior_or_noise_changed=False,wall_seconds=time.perf_counter()-start)
    dump(OUT/'POSTHOC_CONSISTENCY.json',metrics)
    print(json.dumps(metrics))

if __name__=='__main__':main()
