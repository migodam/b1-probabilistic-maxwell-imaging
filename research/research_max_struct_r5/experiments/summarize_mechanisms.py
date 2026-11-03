"""Post-seal summaries only: no Maxwell evaluation or online selection."""
from pathlib import Path
import json, time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def main():
    wall,cpu=time.perf_counter(),time.process_time()
    scores=json.loads((ROOT/'RESULTS/OFFLINE_TRAJECTORY_SCORES.json').read_text())['rows']
    rows=[]; endpoints=[]; counter=[]
    for score in scores:
        if score['source']!='sealed_replay_checkpoint':continue
        path=ROOT/'RESULTS/trajectory_diagnostics'/score['run_id']/score['checkpoint']
        d=json.loads(path.with_suffix('.json').read_text());p=d['physical'];g=d['gradients']['plus']
        row=dict(run_id=score['run_id'],case=score['case'],arm=score['arm'],stage=score['stage'],iteration=score['iteration'],
            cross=d['cross_value'],actual_deltaPhi=p['actual_deltaPhi'],delta_nrmse=score['delta_nrmse'],
            transport_step_cos=d['frozen_LS']['plus']['transport']['cosine_actual_step'],
            connection_to_exact_state_norm=g['connection']['norm']/g['exact_state']['norm'],
            Phi_after=p['Phi_after'],relative_residual=p['relative_residual_after'],NRMSE_after=score['after']['amplitude_normalized_nrmse'],
            projected_gradient_inf=p['box_projected_gradient_inf'],active_upper=p['active_upper_before'],active_lower=p['active_lower_before'])
        rows.append(row)
        if score['iteration']==500:
            z=np.load(path.with_suffix('.npz'));upper=z['m_before']>=20-1e-10
            row=dict(row,endpoint_geometry=d['endpoint_geometry'],
                upper_physical_negative_gradient_inward=int(np.sum(upper&(z['g_physical']>0))),
                upper_aux_negative_gradient_outward=int(np.sum(upper&((z['g_exact_state_plus']+z['g_prior'])<0))))
            endpoints.append(row)
            for kind,cf in d['counterfactual_readouts'].items():
                offline=score['counterfactuals'][kind]
                counter.append(dict(run_id=score['run_id'],case=score['case'],owner_arm=score['arm'],stage=score['stage'],
                    kind=kind,deltaPhi=cf['deltaPhi'],delta_nrmse=offline['delta_nrmse'],
                    cosine_to_negative_physical_gradient=cf['step']['cosine_negative_physical_gradient']))
    summary=dict(classification='OFFLINE_DIAGNOSTIC_POST_SEAL',independent_material_objects=4,actual_steps=len(rows),
        positive_cross_count=sum(r['cross']>0 for r in rows),cross_min=min(r['cross'] for r in rows),cross_max=max(r['cross'] for r in rows),
        transport_step_cos_min=min(r['transport_step_cos'] for r in rows),transport_step_cos_max=max(r['transport_step_cos'] for r in rows),
        actual_Phi_decrease_and_material_worsening_count=sum(r['actual_deltaPhi']<0 and r['delta_nrmse']>0 for r in rows),
        actual_Phi_increase_and_material_improvement_count=sum(r['actual_deltaPhi']>0 and r['delta_nrmse']<0 for r in rows),
        counterfactual_Phi_decrease_and_material_worsening_count=sum(r['deltaPhi']<0 and r['delta_nrmse']>0 for r in counter),
        counterfactuals_count=len(counter),rows=rows,endpoints=endpoints,counterfactuals=counter,
        warning='Checkpoints and counterfactuals are correlated observations, never independent statistical samples. Endpoints are pre-readout m_before, not corrected final outputs.')
    (ROOT/'RESULTS/MECHANISM_SUMMARY.json').write_text(json.dumps(summary,indent=2))
    (ROOT/'audit/MECHANISM_SUMMARY_RECEIPT.json').write_text(json.dumps(dict(status='COMPLETED',
        truth_access='READ_ONLY_PREVIOUSLY_SEALED_OFFLINE_SCORES',new_physical_LU=0,new_complex_RHS=0,
        wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu),indent=2))
    print(json.dumps({k:v for k,v in summary.items() if k not in ('rows','endpoints','counterfactuals')}))

if __name__=='__main__':main()
