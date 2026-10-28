"""Controlled discrete-fibre coupling experiment, original Full66 decoder.

The 4-point matched curves are deterministic controls, not prior-conditional
or posterior samples. Their equal-weight product is an explicit intervention.
No Maxwell adapter, Book, observation, truth data, or neural model is used.
"""
from __future__ import annotations
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import sys,json,time,resource
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
import decoder_geometry as dg

P=Path(__file__).resolve().parents[1]
OUT=P/'audit/copula_full66'
PARENTS=['R7_E1_092','R7_E1_201','R7_E1_233','R7_E1_251']
T=np.array([.001,.01,.1,.5])
OFF_DIAGONAL=np.array([[0,3],[3,0]],dtype=int)

def statistics(material,axes):
    mean=np.mean(material,axis=axes)
    reshaped=material.reshape(-1,576)
    residual=reshaped-mean
    variance=np.mean(abs(residual)**2,axis=0)
    rr=np.mean(residual.real**2,axis=0);ii=np.mean(residual.imag**2,axis=0)
    ri=np.mean(residual.real*residual.imag,axis=0)
    return mean,variance,np.stack([rr,ri,ii]),float(np.mean(variance))

def main():
    start=time.perf_counter();cpu=time.process_time();OUT.mkdir(parents=True,exist_ok=True)
    protocol=dict(experiment='CONTROLLED_DISCRETE_FIBRE',status='FROZEN_BEFORE_DECODER_CALLS',
      parents=PARENTS,curve_t=T,correct_distribution='matched diagonal points (z_i,xi_i), each weight 1/4',
      intervention='independent product of the two discrete marginals: (z_i,xi_j), each weight 1/16',
      grid='all 4x4 combinations; rows are spatial field coordinates, columns are property coordinates',
      predeclared_off_diagonal_pairs=OFF_DIAGONAL,predeclared_diagonal_control=[0,0],
      source_curves='audit/E1_FOUR_RADIAL_CURVES.npz',source_alpha='raw/decoder/NONLINEAR_EQUAL_MATERIAL_CORRECTIONS.json',
      additional_draws=0,new_Maxwell_RHS=0,new_forward=0,new_LU=0,NN_training=0,
      qualifications=['Deterministic equal-weight discrete fibre experiment; not original-prior conditional samples.',
        'Not posterior samples, not a new posterior estimate, and not evidence that this caused R6 errors.',
        'All crossed points are checked for legality under the original Full66 decoder.',
        'Finite same-material equality is numerical; no global exact gauge theorem is inferred.',
        'No measurement, external ground truth, or physics selection is used.'])
    protocol=dg.plain(protocol)
    if (OUT/'PROTOCOL.json').exists():
        if json.loads((OUT/'PROTOCOL.json').read_text())!=protocol:raise RuntimeError('incompatible frozen protocol')
    else:dg.write(OUT/'PROTOCOL.json',protocol)
    prior,_,_=dg.load_prior()
    with np.load(P/'audit/E1_FOUR_RADIAL_CURVES.npz',allow_pickle=False) as f:curves={k:f[k] for k in f.files}
    with np.load(P/'raw/decoder/FULL66_DECODER_ARRAYS.npz',allow_pickle=False) as f:
        center_x=f['x'][:256].copy();center_material=f['material'][:256].copy()
    alpha_records=json.loads((P/'raw/decoder/NONLINEAR_EQUAL_MATERIAL_CORRECTIONS.json').read_text())['records']
    records=[];allx=[];allq=[];allchi=[];allfield=[];allfraction=[];allalpha=[];allcenter=[];correct_means=[];product_means=[];correct_vars=[];product_vars=[];covars=[];endpoints=[]
    for parent in PARENTS:
        ids=np.flatnonzero(curves['parent_id']==parent)
        ids=ids[np.argsort(curves['t'][ids])]
        assert np.array_equal(curves['t'][ids],T)
        x=curves['x'][ids].copy();oldchi=curves['material'][ids]
        pi=int(parent.rsplit('_',1)[1]);base=center_material[pi];q0=dg.from_probit(center_x[pi]);_,aux0=prior.decode(q0)
        delta0=aux0['inclusion_increment'];increment=base-prior.host
        alpha=np.array([next(r['alpha'] for r in alpha_records if r['parent_id']==parent and r['t']==t) for t in T])
        mixed_x=np.empty((4,4,66));mixed_q=np.empty_like(mixed_x);material=np.empty((4,4,576),complex)
        field=np.empty((4,576));fraction=np.empty((4,576));legality=[];alpha_error=[]
        for i in range(4):
            _,auxi=prior.decode(dg.from_probit(x[i]));field[i]=auxi['field'];fraction[i]=auxi['fraction']
            for j in range(4):
                xx=np.r_[x[i,:64],x[j,64:]];q=dg.from_probit(xx);chi,aux=prior.decode(q)
                mixed_x[i,j]=xx;mixed_q[i,j]=q;material[i,j]=chi
                legality.append(prior.legal(q) and bool(np.all((q[64:]>0)&(q[64:]<1))))
                delta=aux['inclusion_increment'];alpha_error.append(max(abs(delta.real/delta0.real-alpha[j]),abs(delta.imag/delta0.imag-alpha[j])))
        assert all(legality)
        diagonal=material[np.arange(4),np.arange(4)]
        cm,cv,cc,cvar=statistics(diagonal,0);pm,pv,pc,pvar=statistics(material,(0,1))
        # Long-double finite sums and the nonnegative pairwise expression
        # avoid cancellation when the prescribed alpha range is very narrow.
        al=alpha.astype(np.longdouble);ea=np.mean(al);ei=np.mean(1/al)
        raw_beta=ea*ei
        stable_bias=np.sum((al[:,None]-al[None,:])**2/(al[:,None]*al[None,:]))/(2*len(al)**2)
        ratio=al[None,:]/al[:,None];ratio_mean=np.mean(ratio);ratio_var=np.mean((ratio-ratio_mean)**2)
        predicted_bias=float(stable_bias)*increment
        predicted_mean=prior.host+float(raw_beta)*increment
        amplitude_RMS=float(np.linalg.norm(increment)/24)
        perpoint_prediction=prior.host+np.asarray(ratio,dtype=float)[:,:,None]*increment[None,None,:]
        obs_bias=pm-cm;pred_var=amplitude_RMS**2*float(ratio_var)
        point_records=[]
        for i,j in OFF_DIAGONAL:
            row=dict(parent_id=parent,field_index=int(i),property_index=int(j),field_t=float(T[i]),property_t=float(T[j]),
                alpha_ratio=float(alpha[j]/alpha[i]),material_RMS_from_correct_mean=float(np.linalg.norm(material[i,j]-cm)/24),
                material_RMS_predicted=abs(float(alpha[j]/alpha[i])-1)*amplitude_RMS,
                x=mixed_x[i,j],candidate_only=True,Maxwell_status='NOT_RUN')
            endpoints.append(row);point_records.append(row)
        endpoints.append(dict(parent_id=parent,field_index=0,property_index=0,field_t=float(T[0]),property_t=float(T[0]),
            alpha_ratio=1.,material_RMS_from_correct_mean=float(np.linalg.norm(material[0,0]-cm)/24),x=mixed_x[0,0],
            candidate_only=True,Maxwell_status='NOT_RUN'))
        rec=dict(parent_id=parent,experiment='CONTROLLED_DISCRETE_FIBRE',alpha=alpha,
            E_alpha=float(ea),E_inverse_alpha=float(ei),E_alpha_E_inverse_alpha=float(raw_beta),
            stable_E_alpha_E_inverse_alpha_minus_one=float(stable_bias),
            cancellation_diagnostic_beta_minus_one=float(raw_beta-1),
            correct_material_variance_RMS_squared=cvar,product_material_variance_RMS_squared=pvar,
            predicted_product_material_variance_RMS_squared=pred_var,
            product_material_variance_relative_error=abs(pvar-pred_var)/pred_var,
            observed_mean_bias_material_RMS=float(np.linalg.norm(obs_bias)/24),
            predicted_mean_bias_material_RMS=float(stable_bias)*amplitude_RMS,
            mean_bias_vector_prediction_error_RMS=float(np.linalg.norm(obs_bias-predicted_bias)/24),
            mean_material_prediction_error_RMS=float(np.linalg.norm(pm-predicted_mean)/24),
            pointwise_alpha_ratio_prediction_max_material_RMS=float(np.max(np.linalg.norm(material-perpoint_prediction,axis=2)/24)),
            matched_curve_max_material_RMS_from_center=float(np.max(np.linalg.norm(diagonal-base,axis=1)/24)),
            matched_curve_archive_redecode_max_absolute_error=float(np.max(abs(diagonal-oldchi))),
            maximum_property_alpha_recovery_error=max(alpha_error),all_16_crossed_points_legal=all(legality),
            fixed_off_diagonal=point_records,
            qualification='Finite artificial matched pairs versus their equal-weight product. This is not a prior-conditional or posterior distribution.')
        records.append(rec);allx.append(mixed_x);allq.append(mixed_q);allchi.append(material);allfield.append(field);allfraction.append(fraction);allalpha.append(alpha);allcenter.append(base)
        correct_means.append(cm);product_means.append(pm);correct_vars.append(cv);product_vars.append(pv);covars.append(np.stack([cc,pc]))
    X=np.array(allx);Q=np.array(allq);CHI=np.array(allchi)
    dg.save(OUT/'CONTROLLED_DISCRETE_FIBRE_ARRAYS.npz',parent_id=np.array(PARENTS),t=T,alpha=np.array(allalpha),
        x_crossed=X,q_crossed=Q,material_crossed=CHI,
        x_correct_curve=X[:,np.arange(4),np.arange(4)],material_correct_curve=CHI[:,np.arange(4),np.arange(4)],
        x_fixed_off_diagonal=X[:,OFF_DIAGONAL[:,0],OFF_DIAGONAL[:,1]],material_fixed_off_diagonal=CHI[:,OFF_DIAGONAL[:,0],OFF_DIAGONAL[:,1]],
        curve_field=np.array(allfield),curve_fraction=np.array(allfraction),
        material_reference=np.array(allcenter),correct_material_mean=np.array(correct_means),product_material_mean=np.array(product_means),
        correct_material_variance=np.array(correct_vars),product_material_variance=np.array(product_vars),
        material_real_covariances_rr_ri_ii=np.array(covars),correct_weights=np.full(4,.25),product_weights=np.full((4,4),1/16),
        fixed_off_diagonal_pairs=OFF_DIAGONAL,spatial_indices=np.arange(64),property_indices=np.array([64,65]))
    dg.write(OUT/'FIXED_ENDPOINT_CANDIDATES.json',dict(experiment='CONTROLLED_DISCRETE_FIBRE',selection='frozen (0,3),(3,0) plus (0,0) for every parent; no outcome selection',
        Maxwell_status='NOT_RUN',new_Maxwell_RHS=0,endpoints=endpoints))
    result=dict(status='PASS_CONTROLLED_DISCRETE_FIBRE_DECODER_EXPERIMENT_COMPLETED',experiment='CONTROLLED_DISCRETE_FIBRE',
      n_parents=4,matched_points_per_parent=4,product_points_per_parent=16,total_decoder_grid_points=64,
      new_Maxwell_RHS=0,new_forward=0,new_LU=0,NN_training=0,GT_access=False,observation_access=False,
      records=records,qualifications=protocol['qualifications'],wall_seconds=time.perf_counter()-start,CPU_seconds=time.process_time()-cpu,
      peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    dg.write(OUT/'RESULTS.json',result)
    with (OUT/'RUN_RECEIPTS.jsonl').open('a') as stream:
        stream.write(json.dumps(dg.plain({k:result[k] for k in ('experiment','new_Maxwell_RHS','new_forward','new_LU','NN_training','wall_seconds','CPU_seconds','peak_RSS_bytes')}))+'\n')
    write_report(result);figure(result,np.array(allchi),np.array(allcenter),prior)
    print(json.dumps(dg.plain([{k:r[k] for k in ('parent_id','stable_E_alpha_E_inverse_alpha_minus_one','observed_mean_bias_material_RMS','product_material_variance_RMS_squared','mean_bias_vector_prediction_error_RMS')} for r in records]),indent=2))

def write_report(result):
    rows=['# Controlled discrete nonuniform Full66 fibre coupling experiment','',
      '**CONTROLLED_DISCRETE_FIBRE:** four deterministic matched points per parent, each weighted 1/4, versus the 4×4 product of their field/property marginals, each weighted 1/16. These are not original-prior conditional draws or posterior samples. All 64 crossed coordinates are legal under the unchanged original Full66 decoder. New Maxwell RHS / states / LU / NN training = 0.',
      '', 'Rows use spatial coordinates from curve point i; columns use Gaussianized property coordinates from curve point j. The diagonal retains the numerically same-material coupling. Off-diagonal crossings are fixed before decoding. The finite decoder prediction is chi_ij = host + (alpha_j/alpha_i)(chi_reference-host), with the finite equality residual measured separately. Thus the elementary finite-sum mean factor is E[alpha] E[1/alpha]. No new general theorem is claimed.',
      '', '| Parent | Ealpha E1/alpha − 1 | Mean bias, material RMS | Product variance, RMS² | Bias prediction error, RMS |',
      '|---|---:|---:|---:|---:|']
    for r in result['records']:
        rows.append(f"| {r['parent_id']} | {r['stable_E_alpha_E_inverse_alpha_minus_one']:.9g} | {r['observed_mean_bias_material_RMS']:.9g} | {r['product_material_variance_RMS_squared']:.9g} | {r['mean_bias_vector_prediction_error_RMS']:.3g} |")
    rows+=['','The matched-curve material variance remains at numerical roundoff. The independent product produces nonzero material variance and the predicted mean bias in these finite constructions. This demonstrates what the specified destruction of coupling does in this controlled ensemble; it does not establish original-prior conditional weights, a posterior effect, a causal diagnosis of R6 results, or any scattering consequence.',
      '', 'Full x, native q, 576-cell complex materials, field/fraction arrays, means, per-cell complex variance, Re/Im covariance, weights, and prescribed pairs are in `CONTROLLED_DISCRETE_FIBRE_ARRAYS.npz`. `FIXED_ENDPOINT_CANDIDATES.json` contains the predeclared first/last crossings and one matched control per parent for possible later separately billed physics. Their Maxwell status is NOT_RUN.',
      '', 'The beta−1 calculation uses a nonnegative pairwise finite-sum expression to avoid cancellation for narrow alpha ranges. All raw alpha values, direct Ealpha E1/alpha, pointwise ratio residuals and covariance checks remain in `RESULTS.json`. The source curves and original E0/E1 evidence are read-only.']
    (OUT/'REPORT.md').write_text('\n'.join(rows)+'\n')

def figure(result,materials,reference,prior):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    f,axes=plt.subplots(1,4,figsize=(15,3.8))
    for k,(ax,r) in enumerate(zip(axes,result['records'])):
        distances=np.linalg.norm(materials[k]-reference[k],axis=2)/24
        im=ax.imshow(distances,cmap='magma');ax.set(xlabel='Property point j (t)',ylabel='Field point i (t)',title=r['parent_id'])
        ax.set_xticks(range(4),['.001','.01','.1','.5']);ax.set_yticks(range(4),['.001','.01','.1','.5'])
        f.colorbar(im,ax=ax,fraction=.046,pad=.04)
    f.suptitle('CONTROLLED_DISCRETE_FIBRE: material RMS change under crossed field/property combinations',fontsize=12)
    f.tight_layout();f.savefig(OUT/'CONTROLLED_FIBRE_COPULA.png',dpi=180);f.savefig(OUT/'CONTROLLED_FIBRE_COPULA.pdf');plt.close(f)

if __name__=='__main__':main()
