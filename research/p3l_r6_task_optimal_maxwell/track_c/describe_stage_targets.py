"""Saved CDEV64 only: Native/Probit target mismatch and warm-start allocation.

No truth/Target/Jacobian solve/import, no new optimization. Native likelihood is
read at its actual endpoint. Clipped property penalties refer ONLY to Original
Probit initialization; the first Probit likelihood is read from its paid trace.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
import csv,json
from pathlib import Path
import numpy as np
from scipy.special import ndtri
TRACK=Path(__file__).resolve().parent
CONFIGS=('native24_probit64','native48_probit64','native96_probit64')


def main():
    destination=TRACK/'stage_target_mismatch'
    if destination.exists():raise FileExistsError('new descriptive analysis; preserve prior outputs')
    destination.mkdir()
    parents=[r['parent'] for r in json.loads((TRACK.parent/'configs/parents_CDEV.json').read_text())['parents']]
    if len(parents)!=64:raise ValueError('all64 inclusive required')
    rows=[];finals={};reference={}
    for parent in parents:
        with np.load(TRACK/'contexts/CDEV'/parent/'native96_probit64/CONTEXT.npz',allow_pickle=False) as z:
            reference[parent]={k:z[k].copy() for k in ('q0','chi0','x')}
    for config in CONFIGS:
        for parent in parents:
            folder=TRACK/'contexts/CDEV'/parent/config
            report=json.loads((folder/'REPORT.json').read_text());native=report['original_native_report'];probit=report['original_probit_report']
            with np.load(folder/'CONTEXT.npz',allow_pickle=False) as z:
                q,chi,nq,nchi,x=(z[k].copy() for k in ('q0','chi0','native_q','native_chi','x'))
            u=nq[64:];initial_prop_x=ndtri(np.clip(u,1e-4,1-1e-4))
            first=None
            for line in (folder/'original_optimizer/probit/MODE_HISTORY.jsonl').read_text().splitlines():
                event=json.loads(line)
                if event['event']=='EVALUATED':first=event;break
            if first is None:raise ValueError('paid initial Probit endpoint required')
            accepted=sum(json.loads(line).get('phase')=='probit' for line in (folder/'path/PATH.jsonl').read_text().splitlines())
            native_rhs=8*native['fresh_forward_states']+32*native['receiver_J32_calls']
            probit_rhs=report['path_resource_cost']['RHS']-native_rhs
            probit_fresh=probit['fresh_state_attempts']
            derivative_rhs=probit_rhs-8*probit_fresh
            if derivative_rhs<0 or derivative_rhs%32:raise ValueError('Original Probit receiver32 accounting mismatch')
            ref=reference[parent]
            row=dict(parent=parent,config=config,native_cap=report['config']['native_cap'],
                native_Phi=native['phi'],native_objective=native['objective'],native_field_norm=float(np.linalg.norm(nq[:64])),
                native_field_prior_energy=.5*float(nq[:64]@nq[:64]),native_u_real=float(u[0]),native_u_imag=float(u[1]),
                native_real_at_bound=bool(u[0]<=1e-6 or u[0]>=1-1e-6),native_imag_at_bound=bool(u[1]<=1e-6 or u[1]>=1-1e-6),
                native_any_property_at_bound=bool(np.any((u<=1e-6)|(u>=1-1e-6))),
                Probit_initial_clipped_property_x_norm_squared=float(initial_prop_x@initial_prop_x),
                Probit_initial_clipped_property_prior_energy=.5*float(initial_prop_x@initial_prop_x),
                first_paid_Probit_Phi=float(first['Phi']),first_paid_Probit_transformed_energy=float(first['transformed_energy']),
                initialization_clip_only=True,native_RHS=native_rhs,native_fresh_states=native['fresh_forward_states'],
                native_wall_seconds=native['wall_seconds'],Probit_RHS=probit_rhs,Probit_forward_RHS=8*probit_fresh,
                Probit_fresh_forward_attempts=probit_fresh,Probit_receiver_J32_calls=derivative_rhs//32,Probit_accepted_steps=accepted,
                Probit_optimizer_nfev=probit['optimizer'].get('nfev'),Probit_optimizer_status=probit['optimizer']['status'],
                final_Phi=None,
                final_field_norm=float(np.linalg.norm(q[:64])),final_property_x_norm_squared=float(x[64:]@x[64:]),
                final_Probit_transformed_energy=None,
                final_vs_Full_field_latent_l2=float(np.linalg.norm(q[:64]-ref['q0'][:64])),
                final_vs_Full_property_l2=float(np.linalg.norm(q[64:]-ref['q0'][64:])),
                final_vs_Full_material_complex_RMS=float(np.linalg.norm(chi-ref['chi0'])/24),
                final_vs_Full_material_complex_MAE=float(np.mean(abs(chi-ref['chi0']))),
                Native_to_final_field_latent_l2=float(np.linalg.norm(q[:64]-nq[:64])),
                Native_to_final_material_complex_RMS=float(np.linalg.norm(chi-nchi)/24))
            # Read actual best mode energy/likelihood rather than treating the
            # terminal SciPy cost (possibly a different endpoint) as Phi.
            best=json.loads((folder/'original_optimizer/probit/BEST_MODE.json').read_text())
            row['final_Probit_transformed_energy']=best['transformed_energy']
            row['final_Phi']=best['transformed_energy']-.5*float(x@x)
            rows.append(row)
    with (destination/'PARENT_STAGE_TARGETS.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    summaries=[]
    for config in CONFIGS:
        selected=[r for r in rows if r['config']==config]
        summary=dict(config=config,parents=len(selected),property_bound_parent_rate=float(np.mean([r['native_any_property_at_bound'] for r in selected])),
            real_property_bound_rate=float(np.mean([r['native_real_at_bound'] for r in selected])),
            imag_property_bound_rate=float(np.mean([r['native_imag_at_bound'] for r in selected])))
        for key in ('native_Phi','native_field_norm','Probit_initial_clipped_property_x_norm_squared','first_paid_Probit_Phi',
                    'first_paid_Probit_transformed_energy','native_RHS','native_wall_seconds','Probit_RHS','Probit_fresh_forward_attempts',
                    'Probit_accepted_steps','final_Phi','final_field_norm','final_vs_Full_field_latent_l2','final_vs_Full_property_l2',
                    'final_vs_Full_material_complex_RMS','final_vs_Full_material_complex_MAE','Native_to_final_material_complex_RMS'):
            values=np.array([r[key] for r in selected]);summary[key]=dict(mean=float(values.mean()),median=float(np.median(values)),min=float(values.min()),max=float(values.max()))
        summaries.append(summary)
    result=dict(status='SAVED_CDEV64_DESCRIPTIVE_STAGE_TARGET_ANALYSIS',parents=64,configs=list(CONFIGS),summaries=summaries,
        native_target='Phi(q;y)+0.5||z64||^2 with uniform u2 in native box',Probit_target='Phi(from_Probit(x);y)+0.5||x66||^2',
        property_penalty_convention='||ndtri(clip(u,1e-4,1-1e-4))||^2; clipping is Original initialization only; not a prior/target modification',
        first_Probit_likelihood_convention='actual first paid EVALUATED trace; never reuse unclipped NativePhi as clipped-point Phi',
        final_likelihood_convention='BEST_MODE transformed_energy minus0.5||saved x||^2; SciPy terminalcost is not Phi',
        material_difference_scope='candidate final anchor minus saved Full anchor, no GT accuracy claim',new_RHS=0,new_optimizer_steps=0,GT_read=False,D_GT_read=False,
        interpretation='descriptive target-mismatch/warm-start allocation evidence only; no causal proof and no change to either MAP law')
    (destination/'SUMMARY.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
