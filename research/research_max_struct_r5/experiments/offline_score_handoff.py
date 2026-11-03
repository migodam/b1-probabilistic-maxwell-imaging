"""Seal and score two natural-stage handoffs offline, without physics calls."""
from pathlib import Path
from datetime import datetime,timezone
import importlib.util,json,time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'RESULTS'
AUDIT=ROOT/'audit'
REPO=ROOT.parents[1]
R2=REPO/'research/research_max_struct_r2'
CASE='obs_b20d5a55677d4eca'
SEAL=AUDIT/'HANDOFF_SEAL.json'
RECEIPT=AUDIT/'OFFLINE_HANDOFF_SCORING_RECEIPT.json'
spec=importlib.util.spec_from_file_location('material_scoring_only',ROOT/'experiments/offline_score_stage_i.py')
score=importlib.util.module_from_spec(spec);spec.loader.exec_module(score)
metrics=score.metrics
def utc():return datetime.now(timezone.utc).isoformat()
def same(entries):return all(metrics.unchanged(e) for e in entries)

def from_saved(m,relative,objective,public,truth,support):
    result=score.material(m,truth,support)
    whitened=float(relative*np.linalg.norm(public['Y'])/float(public['sigma']))
    data=whitened<=float(public['noise_radius'])
    result.update(relative_residual=float(relative),physical_objective=float(objective),whitened_residual=whitened,
        noise_radius99=float(public['noise_radius']),noise_gate=bool(data),
        material_and_data_gate=bool(result['material_gate'] and data and result['material_legal']),
        joint_gate=bool(result['image_gate'] and data and result['material_legal']))
    return result

def run():
    wall,cpu=time.perf_counter(),time.process_time()
    summary=json.loads((OUT/'HANDOFF_SUMMARY.json').read_text())
    if summary['status']!='COMPLETED_TWO_FROZEN_HANDOFFS' or len(summary['runs'])!=2:raise ValueError('Incomplete handoff summary')
    original=[metrics.stamp(AUDIT/n) for n in ['STAGE_I_SEAL.json','FACTORIAL_EXTENSION_SEAL.json']]
    files=[metrics.stamp(OUT/'HANDOFF_SUMMARY.json')]
    for run in summary['runs']:
        folder=OUT/'handoff'/run['handoff_id']
        rec=json.loads((folder/'RECEIPT.json').read_text())
        if rec['status']!='COMPLETED_DEVELOPMENT_STAGE_HANDOFF':raise ValueError('Incomplete handoff')
        files.extend(metrics.stamp(p) for p in sorted(folder.iterdir()) if p.is_file())
    seal=dict(status='TWO_NATURAL_STAGE_HANDOFFS_SEALED',sealed_utc=utc(),case_id=CASE,
        classification='PUBLIC_DEVELOPMENT_DIAGNOSTIC_NOT_BLIND',runs=2,unique_cases=1,
        truth_arrays_opened_by_this_scoring_run=False,files=files,hashes='not computed')
    score.dump(SEAL,seal)
    receipt=dict(status='SCORING',model_assignment='Sol medium',seal_path=str(SEAL.resolve()),sealed_utc=seal['sealed_utc'],
        first_GT_read_utc=utc(),classification=seal['classification'],physical_LU=0,physical_complex_RHS=0,
        physics_imported=False,parameters_changed=False,selection_changed=False,metrics_contract=metrics.frozen_contract(),
        truth_use='OFFLINE_DIAGNOSTIC only; never online selection or algorithm parameter tuning')
    score.dump(RECEIPT,receipt)
    truth_path=R2/'sealed_truth'/(CASE+'.npz')
    with np.load(truth_path,allow_pickle=False) as f:
        truth=f['m'].copy().ravel();support=f['mask'].copy().ravel().astype(bool)
    receipt['truth_reads']=[dict(path=str(truth_path.resolve()),arrays_read=['m','mask'],opened_utc=receipt['first_GT_read_utc'])]
    RECEIPT.write_text(json.dumps(receipt,indent=2)+'\n')
    public=metrics.public_metrics(R2/'data_public'/(CASE+'.npz'))
    rows=[]
    for run in summary['runs']:
        handoff=run['handoff_id'];folder=OUT/'handoff'/handoff
        rec=json.loads((folder/'RECEIPT.json').read_text())
        common=dict(case=CASE,handoff_id=handoff,source='new_R5_natural_stage_handoff',natural_stage_split='after original mf3 stage0; continue original mf6/mf10 stages',
            corrector_stop=rec['corrector']['stop'],corrector_iterations=rec['corrector']['iterations'],
            convergence_qualified=False,convergence_note='1000-iteration budget reached; no convergence certificate')
        with np.load(folder/'AUXILIARY_CANDIDATES.npz',allow_pickle=False) as f:
            cm=f['m'].copy();phi=f['original_Phi'].copy();relative=f['physical_relative_residual'].copy()
        for i,candidate in enumerate(rec['auxiliary_candidate_scores']):
            stage=int(candidate['stage'])
            label='reused_host' if stage==-2 else 'handoff_initial_original_stage0' if stage==-1 else 'remaining_original_stage'+str(stage+1)
            row=dict(**common,state_kind='auxiliary_candidate',candidate_index=i,stage_local=stage,state_label=label,
                prediction_path=str((folder/'AUXILIARY_CANDIDATES.npz').resolve()),
                physical_metrics_source='Already saved exact physical candidate validation; not recomputed',
                **from_saved(cm[i],relative[i],phi[i],public,truth,support))
            row['historical_gate_status']='PASS' if row['joint_gate'] else 'FAIL';rows.append(row)
        selected=rec['auxiliary_selected_stage']
        idx=next(i for i,x in enumerate(rec['auxiliary_candidate_scores']) if x['stage']==selected)
        with np.load(folder/'SELECTED_BEFORE_CORRECTOR.npz',allow_pickle=False) as f:seed=f['m'].copy()
        if not np.array_equal(seed,cm[idx]):raise ValueError('Saved selection mismatch')
        row=dict(**common,state_kind='selected_before_corrector',candidate_index=idx,stage_local=selected,
            state_label='Frozen original-Phi selection; not GT-selected',prediction_path=str((folder/'SELECTED_BEFORE_CORRECTOR.npz').resolve()),
            physical_metrics_source='Saved selected candidate physical validation',
            **from_saved(seed,relative[idx],phi[idx],public,truth,support))
        row['historical_gate_status']='PASS' if row['joint_gate'] else 'FAIL';rows.append(row)
        m,yp=metrics.prediction(folder/'prediction.npz',(24,24))
        row=dict(**common,state_kind='final_prediction',candidate_index=None,stage_local=None,state_label='Sealed final original-Phi incumbent',
            prediction_path=str((folder/'prediction.npz').resolve()),physical_metrics_source='Saved Y_prediction, no new forward',
            **metrics.metrics(m,yp,public,truth,support,(24,24)))
        row['historical_gate_status']='PASS' if row['joint_gate'] else 'FAIL';rows.append(row)
    result=dict(status='COMPLETE',classification=seal['classification'],seal_path=str(SEAL.resolve()),sealed_utc=seal['sealed_utc'],
        first_GT_read_utc=receipt['first_GT_read_utc'],metrics_contract=metrics.frozen_contract(),new_physical_LU=0,new_complex_RHS=0,
        denominators=dict(new_handoff_runs=2,public_development_cases=1,new_blind_cases=0,final_predictions=2,auxiliary_candidates=8,
            selected_before_corrector_references=2,total_table_rows=12),
        interpretation='Natural stage split with original PR reset and frequency zero-padding. Both correctors reached 1000 budget; no convergence claim. Candidate/seed rows are not additional independent runs.',rows=rows)
    score.dump(OUT/'OFFLINE_HANDOFF_SCORES.json',result)
    metrics.write_csv_new(OUT/'OFFLINE_HANDOFF_SCORES.csv',rows)
    lines=['# R5 natural-stage handoff：离线评分','',
        '两次handoff均来自同一个已公开旧弯条开发对象；独立新盲测对象0。自然分界为原mf3/stage0之后，续跑原mf6/mf10阶段；阶段边界PR重置和Fourier零填充保持原规则。',
        '',f'独立封存：{seal["sealed_utc"]}；本次首次GT读取：{receipt["first_GT_read_utc"]}。本次评分新增物理LU/RHS为0，原封存未修改。',
        '', '两个corrector均达到1000迭代预算上限，未认证收敛。候选和选中seed仅作轨迹参考，不计为额外独立运行；不按GT改选。',
        '', '| handoff | 状态 | RMSE/12 | IoU | 相对残差 | Phi | 历史联合门 |','|---|---|---:|---:|---:|---:|---|']
    for r in rows:lines.append(f'|{r["handoff_id"]}|{r["state_kind"]}: {r["state_label"]}|{r["amplitude_normalized_nrmse"]:.9f}|{r["otsu_iou"]:.9f}|{r["relative_residual"]:.9f}|{r["physical_objective"]:.9g}|{r["historical_gate_status"]}|')
    lines+=['','沿用历史RMSE/12≤.15、预测自身Otsu256 IoU≥.7、公共99%噪声球与材料盒，未增加门。GT仅用于封存后描述，禁止在线择根或参数调整。','']
    with (OUT/'OFFLINE_HANDOFF_SCORES.md').open('x') as f:f.write('\n'.join(lines))
    receipt.update(status='COMPLETE',finished_utc=utc(),wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,
        denominators=result['denominators'],sealed_files_unchanged=same(files),previous_seals_unchanged=same(original),
        output_paths=[str((OUT/n).resolve()) for n in ['OFFLINE_HANDOFF_SCORES.json','OFFLINE_HANDOFF_SCORES.csv','OFFLINE_HANDOFF_SCORES.md']])
    if not receipt['sealed_files_unchanged'] or not receipt['previous_seals_unchanged']:raise ValueError('Sealed inputs changed')
    RECEIPT.write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status='COMPLETE',sealed_files=len(files),final_predictions=[{k:r[k] for k in ['handoff_id','amplitude_normalized_nrmse','otsu_iou','relative_residual','physical_objective','historical_gate_status']} for r in rows if r['state_kind']=='final_prediction'])))
if __name__=='__main__':run()
