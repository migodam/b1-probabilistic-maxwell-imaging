"""Saved-array response decomposition; no Target, LU solve, GT or training."""
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_key]='1'
import csv
import json
import sys
from pathlib import Path
import numpy as np
from scipy.special import ndtr
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import OUT, Book, pack, write, save
from p1_material.prior import MaterialPrior, saturation, saturation_prime

def quantiles(values):
    a=np.asarray(values,dtype=float)
    return dict(zip(('min','p10','p25','median','p75','p90','p99','max'),
                    np.quantile(a,[0,.1,.25,.5,.75,.9,.99,1]).tolist()))

def main():
    dest=OUT/'B/full_proposal_curvature_audit'
    if dest.exists(): raise FileExistsError('No overwrite: '+str(dest))
    book=Book('B','proposal_curvature_audit',analysis_only=True,max_wall=600)
    status='COMPLETE';error=None
    try:
        dest.mkdir(parents=True,exist_ok=False)
        job=OUT/'B/jobs/B_full_nonlinear_DEV030_r0'
        protocol=json.loads((job/'PROTOCOL.json').read_text())
        anchor_path=Path(protocol['anchor'])
        with np.load(anchor_path,allow_pickle=False) as z:
            anchor={k:z[k].copy() for k in ('q','chi','prediction','field','receiver_T')}
        with np.load(job/'FROZEN_PROPOSAL.npz',allow_pickle=False) as z:
            J=z['J'].copy();x0=z['x0'].copy();frozen_q=z['native_anchor_q'].copy()
        cfg=json.loads((OUT/'configs/TARGET.json').read_text())
        points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
        prior=MaterialPrior(cfg,points);sigma=cfg['noise']['sigma']
        with np.load(OUT/f"inputs/online/{protocol['parent']}.npz",allow_pickle=False) as z:
            y=z['observation'].copy()
        q0=anchor['q'];chi0=anchor['chi'];prediction0=anchor['prediction']
        if not np.array_equal(q0,frozen_q): raise ValueError('Frozen J native anchor differs from paid anchor')
        C=np.einsum('rp,ps->rsp',anchor['receiver_T'],anchor['field']).reshape(256,576)
        D=prior.derivative(q0)
        def white(a): return pack(a)*np.sqrt(2)/sigma
        J_check=np.column_stack([white((C@D[:,k]).reshape(32,8)) for k in range(66)])
        j_error=float(np.linalg.norm(J-J_check)/max(np.linalg.norm(J),1e-30))
        if j_error>1e-10: raise ValueError('Saved J does not match native-anchor material tangent')
        r0=white(prediction0-y);phi0=.5*float(r0@r0)
        native_field=cfg['field_mean']+prior.A@q0[:64]
        fraction0=saturation(native_field/prior.temperature)*prior.mask
        df0=saturation_prime(native_field/prior.temperature)/prior.temperature*prior.mask
        delta0=complex(prior.real_interval[0]+np.diff(prior.real_interval)[0]*q0[64],
                       prior.imag_interval[0]+np.diff(prior.imag_interval)[0]*q0[65])
        init_q=x0.copy();init_q[64:]=ndtr(init_q[64:])
        files=sorted((job/'endpoints').glob('*.npz'))
        rows=[];vectors=[]
        for path in files:
            with np.load(path,allow_pickle=False) as z:
                q=z['q'];chi=z['chi'];prediction=z['prediction']
            dq=q-q0;dchi=chi-chi0
            exact=white(prediction-prediction0)
            material=white((C@dchi).reshape(32,8))
            latent=J@dq
            maxwell=exact-material;decoder=material-latent;total=exact-latent
            field=cfg['field_mean']+prior.A@q[:64]
            fraction=saturation(field/prior.temperature)*prior.mask
            ddelta=complex(np.diff(prior.real_interval)[0]*dq[64],np.diff(prior.imag_interval)[0]*dq[65])
            saturation_material=delta0*(fraction-fraction0-df0*(prior.A@dq[:64]))
            interaction_material=ddelta*(fraction-fraction0)
            saturation_response=white((C@saturation_material).reshape(32,8))
            interaction_response=white((C@interaction_material).reshape(32,8))
            e2=float(total@total);m2=float(maxwell@maxwell);d2=float(decoder@decoder)
            cross=2*float(maxwell@decoder);denom=max(e2,1e-30)
            exact_norm=float(np.linalg.norm(exact))
            def phi_change(response): return .5*float((r0+response)@(r0+response))-phi0
            def regimes(f): return np.where(f/prior.temperature<=-1,0,np.where(f/prior.temperature>=1,2,1))
            oldreg,newreg=regimes(native_field)[prior.mask],regimes(field)[prior.mask]
            row=dict(endpoint=path.stem,norm_dz=float(np.linalg.norm(dq[:64])),norm_du=float(np.linalg.norm(dq[64:])),
                material_secant_norm=float(np.linalg.norm(dchi)),exact_response_norm=exact_norm,
                latent_response_norm=float(np.linalg.norm(latent)),material_response_norm=float(np.linalg.norm(material)),
                total_error_norm=float(np.sqrt(e2)),maxwell_remainder_norm=float(np.sqrt(m2)),decoder_secant_error_norm=float(np.sqrt(d2)),
                total_error_relative_exact=float(np.sqrt(e2)/max(exact_norm,1e-30)),
                maxwell_squared_share=m2/denom,decoder_squared_share=d2/denom,cross_squared_share=cross/denom,
                error_cross_cosine=float((maxwell@decoder)/max(np.sqrt(m2*d2),1e-30)),
                cancellation_fraction=float(-cross/max(m2+d2,1e-30)),
                delta_phi_exact=phi_change(exact),delta_phi_latent=phi_change(latent),delta_phi_material=phi_change(material),
                phi_exact=phi0+phi_change(exact),saturation_response_norm=float(np.linalg.norm(saturation_response)),
                property_saturation_interaction_norm=float(np.linalg.norm(interaction_response)),
                decoder_material_error_norm=float(np.linalg.norm(dchi-D@dq)),
                saturation_regime_changed=int(np.sum(oldreg!=newreg)),entered_plateau=int(np.sum((oldreg==1)&(newreg!=1))),
                left_plateau=int(np.sum((oldreg!=1)&(newreg==1))),fraction_delta_rms=float(np.sqrt(np.mean((fraction-fraction0)[prior.mask]**2))),
                decomposition_closure=float(np.linalg.norm(total-maxwell-decoder)),
                decoder_split_closure=float(np.linalg.norm(decoder-saturation_response-interaction_response)))
            rows.append(row);vectors.append(np.stack((exact,material,latent)))
        if not rows: raise ValueError('No paid endpoints')
        with (dest/'ENDPOINTS.csv').open('x',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        save(dest/'RESPONSE_VECTORS.npz',endpoint=np.array([r['endpoint'] for r in rows]),
             vectors=np.asarray(vectors),vector_axis=np.array(['exact_delta_F_white','material_C_delta_chi_white','latent_J_delta_q']),
             anchor_whitened_residual=r0)
        summed={k:sum(r[k]**2 for r in rows) for k in ('total_error_norm','maxwell_remainder_norm','decoder_secant_error_norm')}
        total_sq=summed['total_error_norm'];maxwell_sq=summed['maxwell_remainder_norm'];decoder_sq=summed['decoder_secant_error_norm']
        cross_sq=total_sq-maxwell_sq-decoder_sq
        result=dict(status='COMPLETE_ARRAY_DECOMPOSITION_NOT_POSTERIOR_QUALIFICATION',n_endpoint=len(rows),
          source_job=str(job),anchor=str(anchor_path),physics_calls=0,new_RHS=0,new_LU=0,GT_read=False,training_updates=0,
          packing='source-major Re32 then Im32; sqrt(2)/sigma; Phi=.5||white(F-y)||²',
          anchor_phi=phi0,saved_J_relative_reconstruction_error=j_error,
          initialization=dict(native_anchor_q=q0.tolist(),transformed_initialization_q=init_q.tolist(),
            q_difference_norm=float(np.linalg.norm(init_q-q0)),native_u_boundary=bool(np.any((q0[64:]==0)|(q0[64:]==1))),
            inward_clip_flag=protocol['coordinate_initialization_inward_only'],
            note='All response deltas and J use native paid anchor, not transformed/clipped initialization. No physics evaluated for initialization.'),
          decomposition='exact-latent=(exact-material)+(material-latent); material=C_anchor delta_chi; C[r,s,p]=T[r,p] field[p,s]',
          decoder_split='delta0*(fraction-fraction0-fraction_derivative0*A*delta_z) + delta_property*(fraction-fraction0); includes saturation curvature plus property/fraction bilinear interaction',
          pooled_squared_error=dict(total=total_sq,maxwell=maxwell_sq,decoder=decoder_sq,cross=cross_sq,
              maxwell_share=maxwell_sq/max(total_sq,1e-30),decoder_share=decoder_sq/max(total_sq,1e-30),cross_share=cross_sq/max(total_sq,1e-30)),
          quantiles={k:quantiles([r[k] for r in rows]) for k in rows[0] if k!='endpoint'},
          limitations=['4100 includes chain initial evaluations and proposed points; not4100 independent posterior samples.',
             'Exact response means saved paid FULL endpoint prediction minus paid native-anchor prediction; no new solve.',
             'Material response is anchor tangent in unconstrained complex chi; its remainder includes nonlinear Maxwell propagation, not decoder curvature.',
             'Squared component shares can exceed1; cross term records cancellation or reinforcement. Norm shares alone are not additive causal attribution.',
             'All diagnostics are relative to one fixed native anchor and existing proposal sizes, not infinitesimal derivative tests or alternate proposal trials.',
             'No posterior mixing, novelty or GO decision follows from this mechanism decomposition.'])
        write(dest/'RESULT.json',result)
        book.counts.update(saved_endpoints_read=len(rows),array_C_material_matvec=3*len(rows),array_J_latent_matvec=len(rows),array_D_material_matvec=len(rows),zero_physics_analysis=1)
        lines=['# FULL proposal response decomposition','',f'Paid endpoints: {len(rows)}. New physics/RHS/LU/training: 0. No ground truth read.',
          '',f'Native-anchor Phi: {phi0:.8g}; frozen J reconstruction relative error: {j_error:.3g}.',
          f'Native versus transformed initialization q difference: {np.linalg.norm(init_q-q0):.6g}. All deltas use native anchor.',
          '',f'Pooled squared response-error shares: Maxwell {maxwell_sq/total_sq:.6g}; decoder {decoder_sq/total_sq:.6g}; cross {cross_sq/total_sq:.6g}.',
          'Shares add to1 including cross; negative cross means cancellation. They are fixed-anchor descriptive diagnostics.',
          '', '| Quantity | median | p90 | p99 |','|---|---:|---:|---:|']
        for k in ('delta_phi_exact','delta_phi_latent','delta_phi_material','exact_response_norm','maxwell_remainder_norm','decoder_secant_error_norm','saturation_response_norm','property_saturation_interaction_norm','saturation_regime_changed'):
            s=result['quantiles'][k];lines.append(f"| {k} | {s['median']:.7g} | {s['p90']:.7g} | {s['p99']:.7g} |")
        lines+=['','Decoder contribution separates saturation finite secants from the property/fraction interaction. Maxwell remainder holds the actual material secant fixed and measures anchor material-response failure.','']+['- '+v for v in result['limitations']]
        (dest/'REPORT.md').write_text('\n'.join(lines)+'\n')
        print(json.dumps({'path':str(dest),'n_endpoint':len(rows),'pooled':result['pooled_squared_error'],
                          'selected_quantiles':{k:result['quantiles'][k] for k in ('delta_phi_exact','delta_phi_latent','delta_phi_material','saturation_regime_changed')}}))
    except BaseException as exc:
        status='ERROR';error=repr(exc);raise
    finally: book.close(status,error)

if __name__=='__main__': main()
