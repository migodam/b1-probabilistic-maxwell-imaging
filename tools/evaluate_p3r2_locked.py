#!/usr/bin/env python3
"""Separate offline locked evaluator. Truth opens only after all completion gates."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import json
import sys
import time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import torch
from scipy.special import ndtri
from p3r2.density import load_density,labels_from_q,LOG2PI
from infer_p3r2_locked import DEFAULT_OUT,POINT_METHODS,frozen_contract,read_online,write,save


def closed(path):
    value=json.loads(Path(path).read_text())
    if value.get('status')!='COMPLETE':raise ValueError('required job not COMPLETE: '+str(path))


def bootstrap(candidate,baseline,indices,relative_risk=True):
    delta=candidate-baseline
    means=delta[indices].mean(1)
    denominator=baseline[indices].mean(1)
    relative=-means/np.maximum(denominator,1e-300)
    result=dict(mean_difference=float(delta.mean()),interpretation='negative difference favors candidate')
    if relative_risk:result['relative_risk_reduction']=float(1-candidate.mean()/baseline.mean())
    for confidence in (.95,.975):
        tail=(1-confidence)/2;label='95' if confidence==.95 else '97_5'
        result['difference_CI_'+label]=np.quantile(means,[tail,1-tail]).tolist()
        if relative_risk:result['relative_reduction_CI_'+label]=np.quantile(relative,[tail,1-tail]).tolist()
    return result


def gallery(out,parent,truth,predictions,cfg):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    names=('truth','anchor','fused','direct','residual')
    images={'truth':truth,**predictions}
    fig,axes=plt.subplots(2,5,figsize=(14,5.8),constrained_layout=True)
    limits=((cfg['host_chi'][0],cfg['host_chi'][0]+cfg['inclusion_increment_real_interval'][1]),
            (cfg['host_chi'][1],cfg['host_chi'][1]+cfg['inclusion_increment_imag_interval'][1]))
    for row,component in enumerate(('real','imag')):
        for col,name in enumerate(names):
            values=getattr(images[name],component).reshape(24,24)
            image=axes[row,col].imshow(values,origin='lower',cmap='viridis',vmin=limits[row][0],vmax=limits[row][1])
            axes[row,col].set_title(name+' '+component);axes[row,col].set_xticks([]);axes[row,col].set_yticks([])
        fig.colorbar(image,ax=axes[row,:],shrink=.85)
    fig.suptitle(parent+' — shared material color scales')
    path=out/'gallery'/f'{parent}.png';path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path,dpi=150);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=DEFAULT_OUT);parser.add_argument('--methods',type=Path)
    args=parser.parse_args();out=args.root.resolve();methods=args.methods or out/'METHODS_FROZEN_FOR_LOCKED.json'
    torch.set_num_threads(1)
    if torch.get_num_interop_threads()!=1:torch.set_num_interop_threads(1)
    manifest,parents,paths,fusion=frozen_contract(out,methods)
    estimates=out/'L/estimates';complete=json.loads((estimates/'COMPLETE.json').read_text())
    if complete.get('status')!='COMPLETE' or complete.get('parents')!=parents or complete.get('count')!=256:
        raise ValueError('all256 completed inference parents required before truth access')
    if complete.get('source_revision')!=manifest['source_revision'] or float(complete.get('frozen_unix',-1))!=float(manifest['frozen_unix']):
        raise ValueError('inference belongs to another freeze')
    generator=Path(manifest.get('generator_closed_json',out/'L/jobs/data_locked_0_256/CLOSED.json'))
    if not generator.is_absolute():generator=out/generator
    closed(generator)
    generated=json.loads((generator.parent/'DATA_MANIFEST.json').read_text())
    if generated.get('group')!='locked' or set(generated.get('parents',[]))!=set(parents):
        raise ValueError('completed locked generator must cover the exact frozen256 parents')
    context_jobs=manifest.get('context_closed_jsons',[str(out/'L/jobs/anchors_locked_0_256_96_64/CLOSED.json')])
    if not context_jobs:raise ValueError('context job closures required')
    for item in context_jobs:
        path=Path(item);closed(path if path.is_absolute() else out/path)
    for item in manifest.get('additional_online_closures',[]):
        path=Path(item);closed(path if path.is_absolute() else out/path)
    for parent in parents:
        if not (estimates/f'{parent}.npz').exists() or not (out/'L/contexts'/f'{parent}.npz').exists():
            raise ValueError('all context and estimate files required before truth access')
        sidecar=json.loads((out/'L/contexts'/f'{parent}.json').read_text())
        if sidecar.get('native_cap')!=96 or sidecar.get('probit_cap')!=64 or sidecar.get('GT_read') is not False:
            raise ValueError('required truth-free native96 plus probit64 context: '+parent)
        if manifest.get('include_matched_point176',False):
            with np.load(estimates/f'{parent}.npz',allow_pickle=False) as estimate:
                if 'matched_point176' not in estimate.files:raise ValueError('all256 saved matched estimates required before truth access')
    dest=out/'L/evaluation'
    if (dest/'COMPLETE.json').exists():raise FileExistsError('locked evaluation already complete; no implicit repeat')
    dest.mkdir(parents=True,exist_ok=True)
    start_wall,start_cpu=time.monotonic(),time.process_time()
    cfg=json.loads(paths['target_config'].read_text())
    with np.load(estimates/'DENSITY_FRAMES.npz',allow_pickle=False) as z:
        frames={name:z[name+'_frame'] for name in ('density','weak_density')}
    density_models={};density_load_failures={}
    for name,key in (('density','density_checkpoint'),('weak_density','weak_density_checkpoint')):
        try:density_models[name]=load_density(paths[key],device='cpu')
        except Exception as error:density_models[name]=None;density_load_failures[name]=repr(error)
    ablation_names=['ablation__'+entry['name'] for entry in manifest.get('ablation_checkpoints',[])]
    names=list(POINT_METHODS)+ablation_names
    if 'classical_action_json' in paths:names.append('classical_postprocess')
    matched_present=manifest.get('include_matched_point176',False)
    if matched_present:names.append('matched_point176')
    metrics={key:[] for key in ('fullgrid_MSE','fullgrid_MAE','ROI_MSE','ROI_MAE')}
    density_rows={name:dict(nll=[],unit_prior_nll=[],coverage=[],fallback=[]) for name in ('density','weak_density')}
    statuses=[];gallery_data=[];truth_opened=[]
    for index,parent in enumerate(parents):
        # This is the first and only labelled archive access in either script.
        truth_path=out/'inputs/offline_only'/f'{parent}.npz'
        if truth_path.stat().st_mtime<float(manifest['frozen_unix']):
            raise ValueError('locked truth predates frozen methods: '+parent)
        with np.load(truth_path,allow_pickle=False) as truth_file:
            truth=truth_file['chi'].copy();q=truth_file['q'].copy()
        truth_opened.append(parent)
        if truth.shape!=(576,) or q.shape!=(66,) or not np.isfinite(truth).all() or not np.isfinite(q).all():
            raise ValueError('invalid locked truth; no parent exclusion permitted')
        with np.load(estimates/f'{parent}.npz',allow_pickle=False) as archive:
            if str(archive['source_revision'])!=manifest['source_revision'] or float(archive['frozen_unix'])!=float(manifest['frozen_unix']):
                raise ValueError('per-parent estimate belongs to another freeze')
            roi=archive['roi_indices'];status=json.loads(str(archive['status_json']))
            predictions={name:archive[name].copy() for name in names}
            intervals={name:archive[name+'_interval90'].copy() for name in ('density','weak_density')}
            fallbacks={name:bool(archive[name+'_used_prior_fallback']) for name in ('density','weak_density')}
        matrix=np.stack([predictions[name] for name in names])
        if matrix.shape!=(len(names),576) or not np.isfinite(matrix).all():raise ValueError('all finite material estimates required')
        error=np.abs(matrix-truth)
        metrics['fullgrid_MSE'].append(np.mean(error**2,axis=1));metrics['fullgrid_MAE'].append(np.mean(error,axis=1))
        metrics['ROI_MSE'].append(np.mean(error[:,roi]**2,axis=1));metrics['ROI_MAE'].append(np.mean(error[:,roi],axis=1))
        try:observation=read_online(out/'L/contexts'/f'{parent}.npz')['observation']
        except Exception:observation=None
        for name in ('density','weak_density'):
            label=labels_from_q(q[None,:],frames[name])[0]
            prior_nll=float(.5*np.sum(label.astype('float64')**2+LOG2PI))
            fallback=fallbacks[name]
            if not fallback and density_models[name] is not None and observation is not None:
                try:nll=float(-density_models[name].log_prob_labels(observation,label[None,:])[0])
                except Exception as exception:
                    # A score failure is retained explicitly and is never silently
                    # replaced by a different evaluated learned method.
                    raise RuntimeError('frozen density scoring failed for '+parent+'/'+name) from exception
            elif fallback:nll=prior_nll
            else:raise RuntimeError('successful frozen density unavailable for scoring: '+name)
            if not np.isfinite(nll):raise ValueError('nonfinite density score; no exclusions')
            coverage=(label>=intervals[name][0])&(label<=intervals[name][1])
            density_rows[name]['nll'].append(nll);density_rows[name]['unit_prior_nll'].append(prior_nll)
            density_rows[name]['coverage'].append(coverage);density_rows[name]['fallback'].append(fallback)
        statuses.append(status)
        if index<4:gallery_data.append((parent,truth,predictions))
    metric_arrays={key:np.asarray(value,dtype='float64') for key,value in metrics.items()}
    indices=np.random.default_rng(2026102292).integers(0,256,size=(2000,256))
    baseline_index=names.index('probit');fused_index=names.index('fused')
    summaries={};comparisons={}
    for metric,values in metric_arrays.items():
        summaries[metric]={name:dict(mean=float(values[:,j].mean()),parent_sem=float(values[:,j].std(ddof=1)/16)) for j,name in enumerate(names)}
        comparisons[metric]={name:bootstrap(values[:,j],values[:,baseline_index],indices) for j,name in enumerate(names) if name!='probit'}
    additional={}
    for baseline in ('native','matched_point176'):
        if baseline in names:additional['fused_vs_'+baseline]=bootstrap(metric_arrays['fullgrid_MSE'][:,fused_index],metric_arrays['fullgrid_MSE'][:,names.index(baseline)],indices)
    raw=dict(parent=np.asarray(parents),method=np.asarray(names),**metric_arrays)
    fallback_flags=np.asarray([[name in status['failures'] for name in names] for status in statuses],dtype=bool)
    raw['method_used_fallback']=fallback_flags
    fallback_counts={name:int(fallback_flags[:,j].sum()) for j,name in enumerate(names)}
    fallback_counts['context']=sum('context' in status['failures'] for status in statuses)
    baseline_contains_fallbacks=any(fallback_counts[key] for key in ('native','probit','context'))
    density_summary={}
    for name,rows in density_rows.items():
        nll=np.asarray(rows['nll']);prior=np.asarray(rows['unit_prior_nll']);coverage=np.asarray(rows['coverage'],dtype=bool)
        raw[name+'_nll']=nll;raw[name+'_unit_prior_nll']=prior;raw[name+'_coverage90']=coverage
        raw[name+'_used_prior_fallback']=np.asarray(rows['fallback'])
        density_summary[name]=dict(dimension=10,mean_nll=float(nll.mean()),unit_prior_mean_nll=float(prior.mean()),
            paired_nll_comparison=bootstrap(nll,prior,indices,relative_risk=False),component_coverage90=coverage.mean(0).tolist(),
            mean_component_coverage90=float(coverage.mean()),fallback_count=int(np.sum(rows['fallback'])),
            interpretation='proper NLL and marginal fixed-coordinate interval coverage only; no full-posterior calibration or material Bayes-mean claim',
            interval_samples=1024,coordinate_definition='frozen orthonormal Gaussian z8 + probit u2')
    save(dest/'RAW_RISKS.npz',**raw)
    result=dict(status='COMPLETE',parents=parents,count=256,no_exclusions=True,primary_metric='parent mean full-grid complex material MSE',
        primary_comparison=comparisons['fullgrid_MSE']['fused'],primary_baseline='strong probit point: native96 plus incremental probit64',
        per_method_fallback_counts=fallback_counts,primary_baseline_contains_fallbacks=baseline_contains_fallbacks,primary_baseline_nominal_contexts_complete=not baseline_contains_fallbacks,
        primary_baseline_qualification='COMPARISON_WITH_FALLBACKS: native/probit/context failures present; comparisons retain all256 and include explicitly recorded fallback policies' if baseline_contains_fallbacks else 'all256 nominal native96+probit64 baseline contexts available',
        include_matched_point176=matched_present,
        classical_postprocess='frozen validation-selected ordinary smoothing/shrink action; secondary attribution baseline' if 'classical_action_json' in paths else 'NOT_INCLUDED_IN_FROZEN_MANIFEST',
        co_primary_adjustment='two planned comparisons: fused and separately integrated curvature-arc mean; two-sided97.5% intervals reported',
        co_primary_second_comparison='NOT_COMPUTED_HERE; curvature-arc mean comparison is integrated separately by the main analysis',
        bootstrap=dict(draws=2000,seed=2026102292,unit='paired independent parent',confidence_levels=[.95,.975]),
        summaries=summaries,secondary_comparisons_vs_probit=comparisons,additional_baseline_comparisons=additional,
        frozen_ablation_checkpoints=manifest.get('ablation_checkpoints',[]),
        ablation_interpretation='prespecified secondary comparisons; no locked-test reselection of primary fused estimator',
        density=density_summary,density_load_failures=density_load_failures,parent_status=statuses,
        source_revision=manifest['source_revision'],frozen_unix=manifest['frozen_unix'],truth_opened_only_after_all_gates=True,
        physics_calls=0,scientific_status='locked risk and fixed-coordinate density evidence; no posterior-reference equivalence asserted')
    write(dest/'RESULTS.json',result)
    for parent,truth,predictions in gallery_data:gallery(dest,parent,truth,predictions,cfg)
    write(dest/'COMPLETE.json',dict(status='COMPLETE',count=256,parents=parents,truth_archives_opened=truth_opened,
        source_revision=manifest['source_revision'],frozen_unix=manifest['frozen_unix'],physics_calls=0,
        cpu_seconds=time.process_time()-start_cpu,wall_seconds=time.monotonic()-start_wall))
    print(json.dumps(dict(primary=result['primary_comparison'],density=density_summary)),flush=True)


if __name__=='__main__':main()
