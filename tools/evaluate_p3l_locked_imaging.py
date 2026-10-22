"""Separate frozen locked-imaging evaluation; authoring grants no execution.

No Target, solver, trainer or Book import. Main orchestrator must verify the
online processes exited and wrap execution in its analysis-only receipt.
Only this dedicated offline phase reads locked truth, after all128 estimates
and frozen density/context/draws have been loaded and recorded first.
"""
from pathlib import Path
import argparse
import json
import sys
import zipfile

import numpy as np
from scipy.special import ndtri

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from p3l.learned_images import load_selected_density

DEFAULT_OUT=ROOT/'research/p3l_local_posterior_imaging'
PARENTS=[f'SBI_locked_{i:04d}' for i in range(128)]
METHODS=['cheap16','independent_native96','probit64_coordinate_mode',
         'model_times_strong_Gaussian_mean','conditional_mode_plugin_mean','original_prior_material_mean']
FALLBACK=['independent_native96','cheap16','original_prior_material_mean']
CAP=128*1024*1024


def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,ensure_ascii=False,allow_nan=False,indent=2)+'\n')
    temp.replace(path)


def inside(path,out):
    path=Path(path).resolve()
    if not path.is_relative_to(out):raise ValueError('path escapes campaign root')
    return path


def load(path,names,*,out):
    path=inside(path,out)
    with zipfile.ZipFile(path) as archive:
        if any(archive.getinfo(name+'.npy').file_size>CAP for name in names):
            raise MemoryError('individual array exceeds frozen128MiB cap')
    with np.load(path,allow_pickle=False) as archive:
        result={name:archive[name].copy() for name in names}
    if any(not np.isfinite(value).all() for value in result.values()):raise ValueError('nonfinite saved array')
    return result


def material(path,key,out):
    value=load(path,[key],out=out)[key]
    if value.shape!=(576,):raise ValueError('material vector must have576pixels')
    return value


def freeze_config(out,spec):
    optional=spec.get('full_corrected_IS_evaluation')
    if optional is not None:
        if (optional.get('predeclared') is not True or
                not all(isinstance(optional.get(k),str) for k in ('array_template','mean_key')) or
                optional['array_template'].count('{parent}')!=1):
            raise RuntimeError('optional IS must be explicitly preregistered with one parent placeholder')
        for key in ('array_template','metadata_template'):
            if key in optional:
                template=optional[key]
                if not isinstance(template,str) or template.count('{parent}')!=1 or 'offline_only' in template:
                    raise RuntimeError('invalid saved online IS evidence template')
                for parent in PARENTS:inside(out/template.format(parent=parent),out)
    config=dict(revision='P3L_LOCKED_IMAGING_OFFLINE_EVALUATION_V1',parents=PARENTS,
        methods=METHODS+(['full_corrected_IS_unqualified_mean'] if optional is not None else []),
        baselines=['cheap16','independent_native96','probit64_coordinate_mode','original_prior_material_mean'],
        fallback_order=FALLBACK,fallback_scope='each missing/failed method; retain raw failure and effective fallback separately',
        prior_mean_source=dict(path='B/gaussian_images/DEV_000.npz',key='original_prior_material_mean',
            interpretation='same original prior material expectation for every parent; no GT selection'),
        losses=['complex_MSE','complex_MAE'],masks=['full_grid','ROI'],
        bootstrap=dict(replicates=2000,seed=2026102299,quantiles=[.025,.975],unit='paired128parents',
            relative_reduction='1-mean_method_risk/mean_baseline_risk; never average per-parent ratios',
            missing_truth='incomplete128cohort gives null aggregate; never remove parents'),
        calibration=dict(coordinates='[Vw.T@(zGT-qcheap[:64]),ndtri(uGT)]',draws=2048,
            empirical_PIT='(count(draw<truth)+0.5count(draw==truth)+0.5)/(2048+1)',
            central_interval=[.05,.95],components=10,parent_units=128,
            qualification='WEAK_COORDINATE_GENERATIVE_CALIBRATION_NOT_FULL_OR_MATERIAL_CREDIBLE_COVERAGE',
            correlated_components_not_independent_parent_units=True,no_property_clipping=True),
        gallery_parents=PARENTS[:4],gallery_scales='sharedReepsilon acrossallparents/methods; sharedImepsilon separately',
        optional_full_corrected_IS=optional,
        selected_density=spec['selected_model'],
        probit_mode_interpretation='coordinate density mode in x=(z,eta); NOT_NATIVE_MAP_NOT_BAYES_ACTION',
        image_MCSE='finite IID model-reference draw error only; NOT posterior MCSE, approximation bias or calibrated uncertainty',
        qualification='NO_AUTOMATIC_POSTERIOR_OR_TEACHER_QUALIFICATION',new_physics_RHS=0)
    path=out/'configs/LOCKED_IMAGING_EVALUATION.json'
    if path.exists():
        if json.loads(path.read_text())!=config:raise RuntimeError('locked evaluation differs from prior frozen config')
    else:
        with path.open('x') as stream:stream.write(json.dumps(config,ensure_ascii=False,allow_nan=False,indent=2)+'\n')
    return config


def gather(out,config,online_rows):
    """All online arrays are materialized before caller starts truth phase."""
    prior=material(out/config['prior_mean_source']['path'],config['prior_mean_source']['key'],out)
    estimates={};public={};rows=[]
    for parent in PARENTS:
        folder=out/f'D/locked_imaging/{parent}'
        values={'original_prior_material_mean':prior.copy()}
        row=dict(parent=parent,online_status=online_rows[parent].get('status'),
            online_phases=online_rows[parent].get('phases',{}),raw_failures={},fallbacks={},
            risks={},raw_risks={},calibration=dict(status='NOT_EVALUATED'))
        paths={'cheap16':(folder/'cheap16_ESTIMATE.npz','chi'),
            'independent_native96':(folder/'independent_native96_ESTIMATE.npz','chi'),
            'probit64_coordinate_mode':(folder/'PROBIT64_ESTIMATE.npz','chi'),
            'model_times_strong_Gaussian_mean':(folder/'LEARNED_REFERENCE_IMAGES.npz','full_learned_gaussian_material_mean'),
            'conditional_mode_plugin_mean':(folder/'LEARNED_REFERENCE_IMAGES.npz','conditional_mode_plugin_material_mean')}
        optional=config['optional_full_corrected_IS']
        if optional is not None:
            paths['full_corrected_IS_unqualified_mean']=(out/optional['array_template'].format(parent=parent),optional['mean_key'])
            row['IS_interpretation']='UNQUALIFIED_FINITE_IMPORTANCE_IMAGE; no posterior qualification'
            if optional.get('metadata_template'):
                try:
                    source=inside(out/optional['metadata_template'].format(parent=parent),out)
                    meta=json.loads(source.read_text());key=optional.get('ESS_key','importance_ESS')
                    row['IS_saved_ESS']=meta.get(key);row['IS_metadata']=meta
                except Exception as exc:row['IS_metadata_error']=repr(exc)
        for method,(path,key) in paths.items():
            try:
                phase_names={'cheap16':('cheap16','SAVED'),
                    'independent_native96':('independent_native96','SAVED'),
                    'probit64_coordinate_mode':('probit_refinement64','SAVED_PROBIT_COORDINATE_MODE_NOT_NATIVE_MAP')}
                if method in phase_names:
                    phase_name,expected_status=phase_names[method]
                    phase=online_rows[parent].get('phases',{}).get(phase_name,{})
                    if phase.get('status')!=expected_status or phase.get('report',{}).get('failure') is not None:
                        raise RuntimeError('online phase failed or missing; partial array cannot bypass fallback')
                if method in ('model_times_strong_Gaussian_mean','conditional_mode_plugin_mean'):
                    if online_rows[parent].get('learned_images',{}).get('status')!='SAVED_UNQUALIFIED_MODEL_REFERENCE':
                        raise RuntimeError('online learned image phase failed or missing')
                values[method]=material(path,key,out)
            except Exception as exc:row['raw_failures'][method]=dict(status='MISSING_OR_FAILED_RETAINED',error=repr(exc))
        available=dict(values)
        for method in config['methods']:
            if method not in values:
                replacement=next(name for name in FALLBACK if name in available)
                values[method]=available[replacement]
                row['fallbacks'][method]=replacement
        info={}
        try:
            context=load(folder/'CHEAP16_CONTEXT.npz',['scalar','image','q'],out=out)
            if context['scalar'].shape!=(1154,) or context['image'].shape!=(4,24,24) or context['q'].shape!=(66,):
                raise ValueError('cheap-context layout mismatch')
            if not np.array_equal(context['scalar'][1024:1090],context['q']):raise ValueError('context q disagreement')
            info['context']=context
        except Exception as exc:row['calibration_context_error']=repr(exc)
        try:
            draws=load(folder/'LEARNED_REFERENCE_IMAGES.npz',['outer_a8_eta2'],out=out)['outer_a8_eta2']
            if draws.shape!=(2048,10):raise ValueError('outer2048x10draw layout mismatch')
            info['draws']=draws
        except Exception as exc:row['calibration_draws_error']=repr(exc)
        row['model_MC_error_diagnostic']={}
        for prefix,name in [('full_learned_gaussian','model_times_strong_Gaussian_mean'),('conditional_mode_plugin','conditional_mode_plugin_mean')]:
            try:
                mcse=material(folder/'LEARNED_REFERENCE_IMAGES.npz',prefix+'_iid_model_mcse_modulus',out)
                if np.iscomplexobj(mcse) or np.any(mcse<0):raise ValueError('invalid saved modelMCSE')
                row['model_MC_error_diagnostic'][name]=dict(mean_pixel_MCSE=float(mcse.mean()),
                    RMS_pixel_MCSE=float(np.sqrt(np.mean(mcse**2))),maximum_pixel_MCSE=float(mcse.max()),
                    status='FINITE_IID_MODEL_DRAW_ERROR_NOT_POSTERIOR_MCSE')
            except Exception as exc:row['model_MC_error_diagnostic'][name]=dict(status='UNAVAILABLE',error=repr(exc))
        public[parent]=info
        estimates[parent]=dict(effective=values,raw=available)
        rows.append(row)
    return estimates,public,rows


def risk(value,truth,masks):
    return {name:dict(complex_MSE=float(np.mean(np.abs(value[mask]-truth[mask])**2)),
        complex_MAE=float(np.mean(np.abs(value[mask]-truth[mask])))) for name,mask in masks.items()}


def density_logp(loaded,label,context):
    model=loaded['model']
    if hasattr(model,'parameters'):
        import torch
        parameter=next(model.parameters(),None)
        dtype=torch.float32 if parameter is None else parameter.dtype
        with torch.no_grad():
            value=model.log_prob(torch.as_tensor(label[None,None],dtype=dtype),
                torch.as_tensor(context['scalar'][None],dtype=dtype),torch.as_tensor(context['image'][None],dtype=dtype))
        result=float(value.detach().cpu().numpy().reshape(-1)[0])
    else:
        result=float(np.asarray(model.log_prob(label[None,None],context['scalar'][None],context['image'][None])).reshape(-1)[0])
    if not np.isfinite(result):raise ValueError('nonfinite proper weak-coordinate log density')
    return result


def calibrate(qtruth,info,Vw,loaded):
    if 'context' not in info or 'draws' not in info:raise RuntimeError('saved calibration context/draws unavailable')
    u=qtruth[64:]
    if np.any((u<=0)|(u>=1)):raise ValueError('strict0<uGT<1required; no clipping')
    label=np.r_[Vw.T@(qtruth[:64]-info['context']['q'][:64]),ndtri(u)]
    draws=info['draws']
    ranks=np.sum(draws<label,axis=0)+.5*np.sum(draws==label,axis=0)
    pit=(ranks+.5)/(len(draws)+1)
    lower,upper=np.quantile(draws,[.05,.95],axis=0)
    return dict(status='WEAK_COORDINATE_GENERATIVE_CALIBRATION_ONLY',true_a8_eta2=label.tolist(),
        proper_raw_weak_eta_NLL=-density_logp(loaded,label,info['context']),
        empirical_component_rank=ranks.tolist(),PIT=pit.tolist(),
        central90_lower=lower.tolist(),central90_upper=upper.tolist(),
        component90_covered=((label>=lower)&(label<=upper)).tolist(),
        component90_width=(upper-lower).tolist(),draws=2048,
        full_posterior_or_material_coverage=False)


def aggregate(rows,config):
    rng=np.random.default_rng(config['bootstrap']['seed'])
    indices=rng.integers(0,128,size=(2000,128))
    result=dict(parent_count=128,missing_truth=[r['parent'] for r in rows if r.get('evaluation_status')!='OFFLINE_LOCKED_EVALUATED'],
        mean_risks={},paired_comparisons={},raw_mean_risks={},failure_counts={},fallback_counts={})
    for method in config['methods']:
        result['failure_counts'][method]=sum(method in r['raw_failures'] for r in rows)
        result['fallback_counts'][method]={name:sum(r['fallbacks'].get(method)==name for r in rows) for name in FALLBACK}
        for mask in config['masks']:
            for metric in config['losses']:
                key=f'{method}/{mask}/{metric}'
                for source,destination in [('risks','mean_risks'),('raw_risks','raw_mean_risks')]:
                    vals=[r.get(source,{}).get(method,{}).get(mask,{}).get(metric) for r in rows]
                    missing=[r['parent'] for r,v in zip(rows,vals) if v is None]
                    result[destination][key]=dict(mean=None if missing else float(np.mean(vals)),missing_parents=missing,
                        status='INCOMPLETE128_NO_PARENT_REMOVAL' if missing else 'COMPLETE128')
                for baseline in config['baselines']:
                    if method==baseline:continue
                    vals=[r.get('risks',{}).get(method,{}).get(mask,{}).get(metric) for r in rows]
                    base=[r.get('risks',{}).get(baseline,{}).get(mask,{}).get(metric) for r in rows]
                    missing=[r['parent'] for r,a,b in zip(rows,vals,base) if a is None or b is None]
                    item=dict(mean_difference=None,paired_bootstrap95_difference=None,relative_reduction=None,
                        paired_bootstrap95_relative_reduction=None,missing_parents=missing)
                    if not missing:
                        a,b=np.asarray(vals),np.asarray(base);diff=a-b
                        bboot=b[indices].mean(1);aboot=a[indices].mean(1)
                        item.update(mean_difference=float(diff.mean()),paired_bootstrap95_difference=np.quantile(diff[indices].mean(1),[.025,.975]).tolist(),
                            relative_reduction=None if b.mean()<=0 else float(1-a.mean()/b.mean()))
                        if np.all(bboot>0):item['paired_bootstrap95_relative_reduction']=np.quantile(1-aboot/bboot,[.025,.975]).tolist()
                    result['paired_comparisons'][f'{method}_minus_{baseline}/{mask}/{metric}']=item
    calibrated=[r for r in rows if r['calibration'].get('status')=='WEAK_COORDINATE_GENERATIVE_CALIBRATION_ONLY']
    calibrated_names={r['parent'] for r in calibrated}
    result['calibration']=dict(parent_units=128,assessed_parents=len(calibrated),missing_parents=[r['parent'] for r in rows if r['parent'] not in calibrated_names],
        correlated_components_not_independent_parent_units=True,full_material_credible_coverage=False)
    if len(calibrated)==128:
        nll=np.asarray([r['calibration']['proper_raw_weak_eta_NLL'] for r in calibrated])
        pit=np.asarray([r['calibration']['PIT'] for r in calibrated]);covered=np.asarray([r['calibration']['component90_covered'] for r in calibrated])
        result['calibration'].update(mean_proper_NLL=float(nll.mean()),parent_NLL_standard_error=float(nll.std(ddof=1)/np.sqrt(128)),
            component90_coverage=covered.mean(0).tolist(),
            component_PIT_histogram10=[np.histogram(pit[:,j],bins=np.linspace(0,1,11))[0].tolist() for j in range(10)],
            component_PIT_mean=pit.mean(0).tolist())
    result['MC_error_unavailable_counts']={method:sum(r['model_MC_error_diagnostic'][method]['status']=='UNAVAILABLE' for r in rows)
        for method in ['model_times_strong_Gaussian_mean','conditional_mode_plugin_mean']}
    return result


def gallery(out,estimates,truth,config):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder=out/'D/locked_evaluation/figures';folder.mkdir(parents=True,exist_ok=True)
    columns=['truth']+config['methods'];result={}
    for component in ('real','imag'):
        values=[getattr(value,component)+(1 if component=='real' else 0)
            for parent in PARENTS[:4] if parent in truth
            for value in [truth[parent]]+list(estimates[parent]['effective'].values())]
        if not values:result[component]=dict(status='NO_TRUTH_FIGURE_PARENT_RETAINED');continue
        lo=min(float(v.min()) for v in values);hi=max(float(v.max()) for v in values)
        fig,axes=plt.subplots(4,len(columns),figsize=(3*len(columns),12),squeeze=False,constrained_layout=True)
        image=None
        for i,parent in enumerate(PARENTS[:4]):
            for j,method in enumerate(columns):
                ax=axes[i,j]
                value=truth.get(parent) if method=='truth' else estimates[parent]['effective'][method]
                if value is None:ax.text(.5,.5,'MISSING TRUTH',ha='center',transform=ax.transAxes)
                else:image=ax.imshow((getattr(value,component)+(1 if component=='real' else 0)).reshape(24,24),origin='lower',vmin=lo,vmax=hi,cmap='viridis')
                if i==0:ax.set_title(method.replace('_','\n'),fontsize=8)
                ax.set_xticks([]);ax.set_yticks([])
                if j==0:ax.set_ylabel(parent)
        if image is not None:fig.colorbar(image,ax=axes.ravel().tolist(),shrink=.7,label='epsilon '+component+'; shared scale')
        fig.suptitle('First four frozen locked parents; effective fallback actions; reference means NOT qualified posterior')
        path=folder/f'fixed_first4_epsilon_{component}.png';fig.savefig(path,dpi=150,bbox_inches='tight');plt.close(fig)
        result[component]=dict(path=path.relative_to(out).as_posix(),shared_min=lo,shared_max=hi)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,default=DEFAULT_OUT)
    args=parser.parse_args();out=args.out.resolve()
    result_path=out/'D/LOCKED_IMAGE_RISK.json'
    if result_path.exists():raise FileExistsError('locked evaluation already exists; no implicit replay')
    closure=json.loads((out/'D/jobs/locked_imaging_frozen_panel/CLOSED.json').read_text())
    if closure.get('status')!='COMPLETE':raise RuntimeError('final COMPLETE online closure required before locked truth')
    frozen=json.loads((out/'METHODS_FROZEN_FOR_LOCKED.json').read_text());spec=frozen['locked_imaging']
    optional=spec.get('full_corrected_IS_evaluation')
    if optional is not None:
        is_closure=json.loads(inside(out/optional['online_job_closure'],out).read_text())
        if is_closure.get('status')!='COMPLETE' or is_closure.get('error') is not None:
            raise RuntimeError('registered Full-IS job must close before any offline locked truth')
    if spec.get('parents')!=PARENTS or spec.get('heldout_risk_requirement_registered') is not True:
        raise RuntimeError('exact frozen128locked identities required')
    manifest=json.loads((out/'configs/PARENTS.json').read_text())
    if [r['parent'] for r in manifest if r['group']=='locked']!=PARENTS:raise RuntimeError('manifest locked identities/order differ')
    online=json.loads((out/'D/LOCKED_ONLINE_IMAGING.json').read_text())
    if [r['parent'] for r in online['rows']]!=PARENTS or online.get('GT_read') is not False:
        raise RuntimeError('exact128online output rows/noGT contract required')
    if any(r.get('status') not in ('ONLINE_ESTIMATES_RECORDED_FAILURES_RETAINED','PARENT_FAILURE_RETAINED') for r in online['rows']):
        raise RuntimeError('online panel contains unrun/interrupted parents')
    if online.get('model_selection',{}).get('selected')!=spec['selected_model']:
        raise RuntimeError('online selected model differs from method freeze')
    config=freeze_config(out,spec) # immutable before any locked GT read
    frame=load(out/'D/FIXED_FRAME.npz',['Vw','Vs'],out=out)
    frame_metadata=json.loads((out/'D/FIXED_FRAME.json').read_text())
    if frame_metadata.get('training_only') is not True or frame_metadata.get('GT_used') is not False:
        raise RuntimeError('training-only frame provenance required')
    Vw,Vs=frame['Vw'],frame['Vs'];Q=np.column_stack((Vw,Vs))
    if Vw.shape!=(64,8) or Vs.shape!=(64,56) or not np.allclose(Q.T@Q,np.eye(64),atol=1e-10,rtol=0):
        raise RuntimeError('frozen complete orthogonal frame invalid')
    model=load_selected_density(spec['selected_model'])
    estimates,public,rows=gather(out,config,{r['parent']:r for r in online['rows']})
    write(out/'D/LOCKED_EVALUATION_INPUTS.json',dict(status='ALL128_ONLINE_ESTIMATES_LOADED_AND_FALLBACKS_FROZEN_BEFORE_TRUTH',
        parents=PARENTS,rows=rows,GT_read=False,online_closure=closure,config='configs/LOCKED_IMAGING_EVALUATION.json'))
    target=json.loads((out/'configs/TARGET.json').read_text())
    points=np.asarray(json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points'])
    masks=dict(full_grid=np.ones(576,bool),ROI=np.linalg.norm(points,axis=1)<=target['roi_radius'])
    # Dedicated offline GT phase starts only here. No online estimator update.
    truth_for_gallery={}
    for row in rows:
        parent=row['parent']
        try:
            truth=load(out/f'inputs/offline_only/{parent}.npz',['chi','q'],out=out)
            if truth['chi'].shape!=(576,) or truth['q'].shape!=(66,):raise ValueError('truth layout mismatch')
            row['risks']={name:risk(value,truth['chi'],masks) for name,value in estimates[parent]['effective'].items()}
            row['raw_risks']={name:risk(value,truth['chi'],masks) for name,value in estimates[parent]['raw'].items()}
            row['evaluation_status']='OFFLINE_LOCKED_EVALUATED'
            try:row['calibration']=calibrate(truth['q'],public[parent],Vw,model)
            except Exception as exc:row['calibration']=dict(status='FAILED_PARENT_RETAINED',error=repr(exc))
            if parent in PARENTS[:4]:truth_for_gallery[parent]=truth['chi'].copy()
        except Exception as exc:row.update(evaluation_status='TRUTH_READ_FAILED_PARENT_RETAINED',truth_error=repr(exc))
    report=dict(status='OFFLINE_LOCKED_IMAGING_EVALUATED_NO_POSTERIOR_QUALIFICATION',
        config='configs/LOCKED_IMAGING_EVALUATION.json',all128_retained=True,parents=PARENTS,
        rows=rows,aggregation=aggregate(rows,config),new_physics_RHS=0,online_GT_use=False,
        selected_density=spec['selected_model'],probit_mode_interpretation=config['probit_mode_interpretation'],
        calibration_scope=config['calibration']['qualification'],model_MC_error_scope=config['image_MCSE'],
        IS_qualification='UNQUALIFIED_ONLY' if config['optional_full_corrected_IS'] is not None else 'NOT_REGISTERED',
        final_scientific_verdict='PENDING_PARENT_REVIEW')
    try:report['figures']=gallery(out,estimates,truth_for_gallery,config)
    except Exception as exc:report['figures']=dict(status='FIGURE_FAILURE_RISKS_RETAINED',error=repr(exc))
    write(result_path,report)
    print(json.dumps(dict(status=report['status'],parents=128,report=str(result_path),new_physics_RHS=0)))


if __name__=='__main__':main()
