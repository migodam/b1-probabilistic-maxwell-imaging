"""Independent saved-weight KKT/FW, morphology and actual L-estimator audit."""
from pathlib import Path
import csv,importlib.util,json
import numpy as np
from .data import release_root,safe,case_metadata,offline_targets
from .frozen_a26_metrics import image_metrics
from .topology_supplement import topology
from .statistics import paired,paired_difference,COMPLEX


def csv_write(path,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)


def topology_metrics(image,truth,bg):
    results={}
    for v in (.15,.25,.35):
        threshold=v*float(np.max(np.abs(truth-bg)))
        pc,ph=topology((np.abs(image-bg)>threshold).reshape(24,24))
        tc,th=topology((np.abs(truth-bg)>threshold).reshape(24,24))
        key={.15:'15',.25:'25',.35:'35'}[v]
        results.update({f'components_{key}':pc,f'holes_{key}':ph,f'component_abs_error_{key}':abs(pc-tc),
                        f'hole_abs_error_{key}':abs(ph-th),f'truth_components_{key}':tc,f'truth_holes_{key}':th})
    magnitude=np.abs(image-bg).reshape(24,24);scale=float(np.max(np.abs(truth-bg)))
    results['intermediate_pixel_fraction']=float(np.mean((magnitude>.15*scale)&(magnitude<.85*scale)))
    results['magnitude_TV']=float(np.abs(np.diff(magnitude,axis=0)).sum()+np.abs(np.diff(magnitude,axis=1)).sum())
    return results


def audit(root):
    root=Path(root).resolve();r=release_root(root);output=root/'results/a30/hull'
    spec=importlib.util.spec_from_file_location('frozen_hull',root/'theory/package/a29_1_convex_diagnostic.py')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    metadata=case_metadata(root);ids=[x['case_id'] for x in metadata]
    truth,bg,n3=offline_targets(root,ids)
    weights=json.loads((output/'A29_1_WEIGHTS.json').read_text())
    saved={x['case_id']:x for x in weights}
    if set(saved)!=set(ids):raise ValueError('MISSING_OR_DUPLICATE_SAVED_WEIGHTS')
    records=[];cert=[];arrays={'case_id':np.asarray(ids),'truth':truth,'background':bg,'N3':n3}
    names=['N3','full_hard','small_hard','full_hull','small_hull','full_L_MAP','small_L_MAP','full_L_mean','small_L_mean']
    collected={k:[] for k in names}
    for i,meta in enumerate(metadata):
        case=meta['case_id']
        with np.load(safe(r,'results/a29/pools/'+case+'.npz'),allow_pickle=False) as z:p={k:z[k].copy() for k in z.files}
        chi=np.asarray(p['chi'],np.complex128);small=np.asarray(p['small_indices'],int)
        den=float(np.linalg.norm(truth[i]-bg[i]));norm=lambda x:float(np.linalg.norm(x-truth[i])/den)
        lf=m.actual_linear_decision(p,np.arange(len(chi)),truth[i],den);ls=m.actual_linear_decision(p,small,truth[i],den)
        images=dict(N3=chi[0],full_hard=chi[np.argmin(np.linalg.norm(chi-truth[i],axis=1))],
                    small_hard=chi[small[np.argmin(np.linalg.norm(chi[small]-truth[i],axis=1))]],
                    full_L_MAP=chi[lf['hard_index']],small_L_MAP=chi[ls['hard_index']],full_L_mean=lf['mean'],small_L_mean=ls['mean'])
        for kind,indices in [('full',np.arange(len(chi))),('small',small)]:
            w=np.asarray(saved[case][kind+'_oracle'],np.float64);a=chi[indices]
            if w.shape!=(len(a),):raise ValueError('WEIGHTS_INDEX_LAYOUT_MISMATCH')
            image=w@a;residual=(image-truth[i])/den
            # Direct vector form independent of the script's Gram arithmetic.
            gradient=np.real(((a-a[0])/den).conj()@residual)
            rawgap=float(w@gradient-np.min(gradient));gap=max(0.,rawgap)
            dual=gradient-gradient.min();complementarity=float(np.max(np.abs(w*dual)))
            e=norm(image);lo=float(np.sqrt(max(0.,e*e-2*gap)))
            row=dict(**meta,pool=kind,simplex_sum_residual=abs(float(w.sum())-1),simplex_min_weight=float(w.min()),
                direct_FW_gap=gap,raw_FW_gap=rawgap,dual_feasibility_min=float(dual.min()),
                complementarity_max=complementarity,active_gradient_spread=float(np.ptp(gradient[w>1e-6])) if np.any(w>1e-6) else None,
                achieved_nrmse=e,optimal_nrmse_lower=lo,optimal_nrmse_upper=e,interval_width=e-lo,
                material_real_min=float(image.real.min()),material_real_max=float(image.real.max()),
                material_imag_min=float(image.imag.min()),material_imag_max=float(image.imag.max()),
                interpretation='OFFLINE_CONVEX_HULL_NUMERICAL_CERTIFICATE_NOT_INTERVAL_ARITHMETIC')
            row['certified']=bool(gap<=1e-8 and rawgap>=-1e-12 and row['simplex_sum_residual']<=1e-10 and w.min()>=-1e-12)
            if not row['certified']:raise ValueError('INDEPENDENT_HULL_CERTIFICATE_FAILED:'+case+kind)
            cert.append(row);images[kind+'_hull']=image
        for method in names:
            image=images[method];collected[method].append(image)
            metrics=image_metrics(image,truth[i],bg[i]);top=topology_metrics(image,truth[i],bg[i])
            records.append(dict(**meta,method=method,N3_nrmse=norm(n3[i]),N3_edge=image_metrics(n3[i],truth[i],bg[i])['edge_f1'],
                N3_iou=image_metrics(n3[i],truth[i],bg[i])['support_iou'],**metrics,**top,
                mass_squared_error=float(np.vdot(image-truth[i],image-truth[i]).real/576),
                ESS=lf['ess'] if method.startswith('full_L') else ls['ess'] if method.startswith('small_L') else None,
                actual_L_weighted_estimator='NOT_CALIBRATED_POSTERIOR' if '_L_' in method else 'OFFLINE_OR_BASELINE'))
    csv_write(output/'HULL_INDEPENDENT_CERTIFICATES.csv',cert);csv_write(output/'HULL_MATERIAL_AND_STRUCTURE.csv',records)
    for key,values in collected.items():arrays[key]=np.asarray(values)
    np.savez_compressed(output/'HULL_RECONSTRUCTIONS.npz',**arrays)
    summary={}
    n3rows={x['case_id']:x for x in records if x['method']=='N3'}
    for method in names:
        rr=[dict(x,N3_component=n3rows[x['case_id']]['component_abs_error_25'],N3_hole=n3rows[x['case_id']]['hole_abs_error_25'])
            for x in records if x['method']==method]
        summary[method]=dict(complex=paired(rr,'delta_nrmse'),family={f:paired(rr,'delta_nrmse',families={f}) for f in sorted(COMPLEX|{'smooth'})},
            edge=paired_difference(rr,'edge_f1','N3_edge'),iou=paired_difference(rr,'support_iou','N3_iou'),
            component_MAE_increment=paired_difference(rr,'component_abs_error_25','N3_component'),
            hole_MAE_increment=paired_difference(rr,'hole_abs_error_25','N3_hole'))
    summary['certificates']=dict(status='PASS',count=len(cert),max_direct_FW_gap=max(x['direct_FW_gap'] for x in cert),
        max_complementarity=max(x['complementarity_max'] for x in cert),new_Maxwell_states=0)
    (output/'HULL_AUDIT_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k=='certificates'}));return summary


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path('.'));a=p.parse_args();audit(a.root)
