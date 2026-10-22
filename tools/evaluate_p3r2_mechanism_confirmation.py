"""Offline confirmation scorer: GT opens only after all new-cohort gates pass.

The three primary contrasts are frozen. Each parent/family risk is the arithmetic
mean of three seed-specific image risks, never risk of an averaged image.
No old256-parent score is opened and no family is selected on the new test.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,csv,json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r2.core import OUT,write,save
from p3r2_mechanism.contract import VERSION,MODEL_COUNT,PRIMARY_CONFIDENCE,PARENTS,FAMILIES,PRIMARY,validate_freeze,validate_online_inputs,completed_jobs,checkpoint_metadata,check_record,file_record,resolve


def bootstrap_comparison(candidate,baseline,indices,confidence):
    delta=candidate-baseline;means=delta[indices].mean(axis=1);tail=(1-confidence)/2
    ci=np.quantile(means,[tail,1-tail]);denominator=baseline[indices].mean(axis=1)
    relative=None;relative_ci=None
    if baseline.mean()>0 and np.all(denominator>0):
        relative=float(1-candidate.mean()/baseline.mean());relative_ci=np.quantile(-means/denominator,[tail,1-tail]).tolist()
    return dict(candidate_mean_risk=float(candidate.mean()),baseline_mean_risk=float(baseline.mean()),
        mean_risk_difference=float(delta.mean()),difference_CI=ci.tolist(),CI_confidence=confidence,
        relative_risk_reduction=relative,relative_reduction_CI=relative_ci,
        direction='Negative difference favors candidate',
        descriptive_verdict='CANDIDATE_LOWER_RISK' if ci[1]<0 else 'BASELINE_LOWER_RISK' if ci[0]>0 else 'UNRESOLVED',
        bootstrap='Paired nonparametric parent bootstrap, all256 including every fallback; finite-sample approximate interval')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=OUT);p.add_argument('--manifest',type=Path);a=p.parse_args()
    out=a.root.resolve();manifest=a.manifest or out/'MECHANISM_CONFIRMATION_FROZEN.json'
    m,paths,frozen_record=validate_freeze(out,manifest)
    closures=completed_jobs(out,m);validate_online_inputs(out,m)
    estimates=out/'L/mechanism_confirmation/estimates';receipt_path=estimates/'COMPLETE.json';receipt=json.loads(receipt_path.read_text())
    if receipt.get('status')!='COMPLETE' or receipt.get('version')!=VERSION or receipt.get('parents')!=PARENTS or receipt.get('count')!=256 or receipt.get('model_count')!=MODEL_COUNT or receipt.get('GT_read') is not False:
        raise ValueError('256 completed truth-free online inference receipts required before GT')
    if receipt.get('manifest_record')!=frozen_record or receipt.get('source_revision')!=m['source_revision'] or float(receipt.get('frozen_unix',-1))!=float(m['frozen_unix']):raise ValueError('inference belongs to another explicit freeze')
    names=[x['name'] for x in m['models']]
    if receipt.get('estimator_action')!='raw_unprojected_no_blend' or receipt.get('primary_metric')!='complex_fullgrid_MSE':raise ValueError('raw full-grid confirmation receipt required')
    if receipt.get('model_names')!=names:raise ValueError('all21 frozen models required')
    metadata={entry['name']:checkpoint_metadata(paths['models'][entry['name']],entry) for entry in m['models']}
    if metadata!=receipt.get('checkpoint_metadata'):raise ValueError('checkpoint metadata/version changed since inference')
    expected_estimates={str((estimates/f'{parent}.npz').resolve()) for parent in PARENTS}
    if {str(resolve(r['path'],out)) for r in receipt.get('estimate_files',[])}!=expected_estimates or len(receipt.get('estimate_files',[]))!=256:raise ValueError('all256 estimate file metadata required')
    expected_online={str(path.resolve()) for parent in PARENTS for path in (out/'L/contexts'/f'{parent}.npz',out/'L/contexts'/f'{parent}.json',out/'inputs/online'/f'{parent}.npz')}
    if {str(resolve(r['path'],out)) for r in receipt.get('online_files',[])}!=expected_online or len(receipt.get('online_files',[]))!=768:raise ValueError('all source-isolated online metadata required')
    for record in receipt['estimate_files']+receipt['online_files']:check_record(record,out)
    cfg=json.loads(paths['target'].read_text());points=np.array(json.loads(paths['geometry'].read_text())['points']);roi=np.flatnonzero(np.linalg.norm(points,axis=1)<=cfg['roi_radius'])
    if len(roi)!=256:raise ValueError('frozen P1 ROI required')
    # Fully validate every saved online image before the first labelled read.
    predictions=np.empty((256,MODEL_COUNT,576),complex);anchors=np.empty((256,576),complex);fallbacks=np.empty((256,MODEL_COUNT),bool)
    statuses=[]
    for i,parent in enumerate(PARENTS):
        with np.load(estimates/f'{parent}.npz',allow_pickle=False) as z:
            if int(z['version'])!=VERSION or str(z['parent'])!=parent or str(z['source_revision'])!=m['source_revision'] or float(z['frozen_unix'])!=float(m['frozen_unix']) or json.loads(str(z['manifest_record_json']))!=frozen_record or not np.array_equal(z['roi_indices'],roi):raise ValueError('per-parent frozen receipt mismatch')
            if str(z['estimator_action'])!='raw_unprojected_no_blend':raise ValueError('per-parent raw action receipt mismatch')
            anchors[i]=z['anchor'];statuses.append(json.loads(str(z['status_json'])))
            if statuses[-1].get('GT_read') is not False:raise ValueError('online receipt must be truth-free')
            for j,name in enumerate(names):predictions[i,j]=z[name];fallbacks[i,j]=bool(z[name+'__used_fallback'])
        if not np.isfinite(predictions[i]).all() or not np.isfinite(anchors[i]).all():raise ValueError('all parents must have finite estimates/fallbacks')
    # No old result is consulted; paths are generated solely from strict00256..00511 IDs.
    truth_paths=[out/'inputs/offline_only'/f'{parent}.npz' for parent in PARENTS]
    for path in truth_paths:
        if not path.is_file() or path.stat().st_mtime<float(m['frozen_unix']):raise ValueError('new post-freeze GT archive required')
    _,_,again=validate_freeze(out,manifest)
    if again!=frozen_record:raise ValueError('freeze changed during gating')
    dest=out/'L/mechanism_confirmation/evaluation'
    if (dest/'COMPLETE.json').exists():raise FileExistsError('new confirmation already evaluated; no implicit repeat')
    dest.mkdir(parents=True,exist_ok=True)
    write(dest/'GT_OPEN_GATES.json',dict(status='ALL_GATES_PASSED_BEFORE_FIRST_GT_READ',version=VERSION,manifest_record=frozen_record,
        parents=PARENTS,count=256,model_count=MODEL_COUNT,completed_job_closures=[file_record(v) for v in closures],inference_receipt=file_record(receipt_path),
        all_frozen_source_and_checkpoint_metadata_unchanged=True,all256_online_estimates_and_inputs_validated=True,old_cohort_scores_read=False,GT_opened=False))
    start=time.monotonic();ROI_MSE=np.empty((256,MODEL_COUNT));fullgrid_MSE=np.empty((256,MODEL_COUNT));anchor_risk=np.empty(256);anchor_ROI_risk=np.empty(256)
    with (dest/'GT_ACCESS_LOG.jsonl').open('x') as access:
        for i,(parent,path) in enumerate(zip(PARENTS,truth_paths)):
            # Recheck declared immutable source/checkpoint metadata immediately
            # before each GT open. No GT-dependent model or action is selected.
            _,_,current=validate_freeze(out,manifest)
            if current!=frozen_record:raise ValueError('freeze changed before GT open')
            with np.load(path,allow_pickle=False) as z:truth=z['chi'].copy()
            if truth.shape!=(576,) or not np.isfinite(truth).all():raise ValueError('invalid GT; no exclusion or replacement permitted')
            error=abs(predictions[i]-truth)**2;ROI_MSE[i]=error[:,roi].mean(axis=1);fullgrid_MSE[i]=error.mean(axis=1)
            anchor_risk[i]=float(np.mean(abs(anchors[i]-truth)**2));anchor_ROI_risk[i]=float(np.mean(abs(anchors[i,roi]-truth[roi])**2))
            access.write(json.dumps(dict(parent=parent,path=str(path),opened_after_all_gates=True,utc_unix=time.time()))+'\n');access.flush()
    family_risk=np.empty((256,len(FAMILIES)));family_ROI_risk=np.empty((256,len(FAMILIES)))
    for k,family in enumerate(FAMILIES):
        ix=[j for j,entry in enumerate(m['models']) if entry['family']==family]
        family_risk[:,k]=fullgrid_MSE[:,ix].mean(axis=1);family_ROI_risk[:,k]=ROI_MSE[:,ix].mean(axis=1)
    rng=np.random.default_rng(2026102393);indices=rng.integers(0,256,size=(2000,256))
    save(dest/'PARENT_RISKS.npz',parent=np.array(PARENTS),model=np.array(names),family=np.array(FAMILIES),ROI_MSE_model=ROI_MSE,fullgrid_MSE_model=fullgrid_MSE,
        fullgrid_MSE_three_seed_average=family_risk,ROI_MSE_three_seed_average=family_ROI_risk,anchor_fullgrid_MSE=anchor_risk,anchor_ROI_MSE=anchor_ROI_risk,fallback=fallbacks)
    save(dest/'BOOTSTRAP_PARENT_INDICES.npz',indices=indices,seed=2026102393)
    primary=[]
    for candidate,baseline in PRIMARY:
        row=bootstrap_comparison(family_risk[:,FAMILIES.index(candidate)],family_risk[:,FAMILIES.index(baseline)],indices,PRIMARY_CONFIDENCE)
        primary.append(dict(candidate=candidate,baseline=baseline,preregistered_primary=True,**row))
    secondary=[]
    for baseline in ('local_anchor','capacity','permuted'):
        secondary.append(dict(candidate='full_geometry',baseline=baseline,preregistered_primary=False,
            **bootstrap_comparison(family_risk[:,0],family_risk[:,FAMILIES.index(baseline)],indices,.95)))
    for k,family in enumerate(FAMILIES):
        secondary.append(dict(candidate=family,baseline='anchor',preregistered_primary=False,
            **bootstrap_comparison(family_risk[:,k],anchor_risk,indices,.95)))
    with (dest/'PER_PARENT_RISKS.csv').open('w') as f:
        fields=['parent','anchor_fullgrid_MSE']+[family+'_three_seed_mean_fullgrid_MSE' for family in FAMILIES]+['anchor_ROI_MSE']+[family+'_three_seed_mean_ROI_MSE' for family in FAMILIES]+['fallback_models']
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for i,parent in enumerate(PARENTS):writer.writerow(dict(parent=parent,anchor_fullgrid_MSE=anchor_risk[i],anchor_ROI_MSE=anchor_ROI_risk[i],fallback_models=int(fallbacks[i].sum()),**{family+'_three_seed_mean_fullgrid_MSE':family_risk[i,k] for k,family in enumerate(FAMILIES)},**{family+'_three_seed_mean_ROI_MSE':family_ROI_risk[i,k] for k,family in enumerate(FAMILIES)}))
    _,_,final_record=validate_freeze(out,manifest)
    if final_record!=frozen_record:raise ValueError('freeze changed during scoring')
    result=dict(status='COMPLETE',version=VERSION,protocol_kind=m['protocol_kind'],source_revision=m['source_revision'],manifest_record=frozen_record,
        parents=PARENTS,count=256,model_count=MODEL_COUNT,objective='SUPERVISED_RISK',primary_metric='complex_fullgrid_MSE',primary_metric_definition='Per-parent complex MSE over all576 cells, averaged over three seed-specific risks',estimator_action='raw_unprojected_no_blend',
        seed_average_is_not_ensemble_image=True,primary_comparisons=primary,secondary_comparisons=secondary,
        bootstrap_replicates=2000,bootstrap_seed=2026102393,co_primary_CI=PRIMARY_CONFIDENCE,
        all256_including_fallbacks_retained=True,no_test_based_winner_selection=True,old256_scores_read=False,
        family_mean_risk={family:float(family_risk[:,k].mean()) for k,family in enumerate(FAMILIES)},
        secondary_family_mean_ROI_MSE={family:float(family_ROI_risk[:,k].mean()) for k,family in enumerate(FAMILIES)},
        model_mean_risk={name:float(fullgrid_MSE[:,j].mean()) for j,name in enumerate(names)},
        fallback_counts={name:int(fallbacks[:,j].sum()) for j,name in enumerate(names)},checkpoint_metadata=metadata,
        GT_opened_only_after_all_gates=True,GT_parents_opened=256,physics_calls=0,training=False,posterior_claim=False,wall_seconds=time.monotonic()-start)
    write(dest/'RESULT.json',result);write(dest/'COMPLETE.json',dict(status='COMPLETE',version=VERSION,count=256,parents=PARENTS,manifest_record=frozen_record,GT_parents_opened=256))
    print(json.dumps(dict(status='COMPLETE',parents=256,primary_comparisons=primary),indent=2),flush=True)


if __name__=='__main__':main()
