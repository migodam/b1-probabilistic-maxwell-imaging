"""Seal then score three disclosed development extensions; no physics imports."""
from pathlib import Path
from datetime import datetime,timezone
import importlib.util,json,time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
R4=ROOT.parent/'research_max_struct_r4'
OUT=ROOT/'RESULTS'
AUDIT=ROOT/'audit'
CASES=['obs_42cb233e4f224bb3','obs_63abafa750fc483c','obs_100bd65498644322']
SEAL=AUDIT/'FACTORIAL_EXTENSION_SEAL.json'
RECEIPT=AUDIT/'OFFLINE_FACTORIAL_EXTENSION_SCORING_RECEIPT.json'
spec=importlib.util.spec_from_file_location('metrics_only',R4/'experiments/offline_score_common.py')
metrics=importlib.util.module_from_spec(spec);spec.loader.exec_module(metrics)
def utc():return datetime.now(timezone.utc).isoformat()
def new(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')
def same(stamps):return all(metrics.unchanged(x) for x in stamps)

def run():
    wall,cpu=time.perf_counter(),time.process_time()
    original=metrics.stamp(AUDIT/'STAGE_I_SEAL.json')
    stamps=[]
    for case in CASES:
        folder=OUT/'factorial_extension'/(case+'__material_only')
        receipt=json.loads((folder/'RECEIPT.json').read_text())
        if receipt.get('status')!='COMPLETED_DEVELOPMENT_FACTORIAL':raise ValueError('Incomplete extension '+case)
        for name in ('prediction.npz','RECEIPT.json','FROZEN_CONFIG.json'):stamps.append(metrics.stamp(folder/name))
    seal=dict(status='THREE_FACTORIAL_EXTENSION_PREDICTIONS_SEALED',sealed_utc=utc(),
        classification='PUBLIC_DISCLOSED_DEVELOPMENT_NOT_NEW_BLIND',case_ids=CASES,new_runs=3,
        truth_arrays_opened_by_this_scoring_run=False,files=stamps,hashes='not computed',
        original_stage_i_seal='Unchanged separate seal; not extended or rewritten')
    new(SEAL,seal)
    receipt=dict(status='SCORING',model_assignment='Sol medium',classification=seal['classification'],
        seal_path=str(SEAL.resolve()),sealed_utc=seal['sealed_utc'],first_GT_read_utc=utc(),
        truth_reads=[],physical_LU=0,physical_complex_RHS=0,physics_imported=False,
        parameters_changed=False,online_selection_changed=False,metrics_contract=metrics.frozen_contract(),
        truth_use='OFFLINE_DIAGNOSTIC only; no subsequent online GT selection or parameter tuning',
        original_stage_i_seal_stamp=original)
    new(RECEIPT,receipt)
    rows=[]
    for case in CASES:
        truth_path=R4/'sealed_truth/confirmation'/(case+'.npz')
        opened=utc()
        with np.load(truth_path,allow_pickle=False) as f:
            truth=f['m'].copy().ravel();support=f['mask'].copy().ravel().astype(bool)
        receipt['truth_reads'].append(dict(path=str(truth_path.resolve()),opened_utc=opened,arrays_read=['m','mask']))
        RECEIPT.write_text(json.dumps(receipt,indent=2)+'\n')
        public=metrics.public_metrics(R4/'data_public/confirmation'/(case+'.npz'))
        definitions=[('new_R5_extension','material_only',OUT/'factorial_extension'/(case+'__material_only')/'prediction.npz')]
        definitions += [('historical_R4_reference',arm,R4/'results/confirmation_remote/predictions'/case/(arm+'.npz')) for arm in ['global','onehop','native_extended']]
        for source,arm,path in definitions:
            m,yp=metrics.prediction(path,(24,24))
            score=metrics.metrics(m,yp,public,truth,support,(24,24))
            row=dict(case=case,arm=arm,source=source,prediction_path=str(path.resolve()),
                platform_comparability='Historical different-platform endpoint; not same execution' if source.startswith('historical') else 'New local R5 material-only extension',
                **score,historical_gate_status='PASS' if score['joint_gate'] else 'FAIL')
            rows.append(row)
    result=dict(status='COMPLETE',classification=seal['classification'],seal_path=str(SEAL.resolve()),sealed_utc=seal['sealed_utc'],
        first_GT_read_utc=receipt['first_GT_read_utc'],metrics_contract=metrics.frozen_contract(),
        denominators=dict(new_material_only_runs=3,new_development_cases=3,new_blind_cases=0,
            historical_references=9,historical_methods=3,historical_case_count=3),
        interpretation='Frozen disclosed-development extension; history references from another platform are not same-run controls. No equal-compute or fresh-blind claim.',
        new_physical_LU=0,new_complex_RHS=0,rows=rows)
    new(OUT/'OFFLINE_FACTORIAL_EXTENSION_SCORES.json',result)
    metrics.write_csv_new(OUT/'OFFLINE_FACTORIAL_EXTENSION_SCORES.csv',rows)
    lines=['# R5 material-only 扩展：离线评分','',
        '新增3个已公开开发对象各1次material-only；历史R4引用为同3例的global/onehop/native三个原端点，共9条。新盲测对象0。历史参考来自另一平台，不称同运行或等算力对照。',
        '',f'独立封存时间：{seal["sealed_utc"]}；本次首次GT读取：{receipt["first_GT_read_utc"]}。原STAGE_I_SEAL未修改；评分新增物理LU/RHS为0。',
        '', '历史合同：RMSE/12≤.15，预测自身Otsu256 IoU≥.7，公共99%噪声球与材料盒；未增加门或改变参数/选择。',
        '', '| 对象 | 来源 | 臂 | RMSE/12 | IoU | 相对残差 | Phi | 历史联合门 |', '|---|---|---|---:|---:|---:|---:|---|']
    for r in rows:lines.append(f'|{r["case"]}|{r["source"]}|{r["arm"]}|{r["amplitude_normalized_nrmse"]:.9f}|{r["otsu_iou"]:.9f}|{r["relative_residual"]:.9f}|{r["physical_objective"]:.9g}|{r["historical_gate_status"]}|')
    lines+=['','GT只用于封存后的离线评分，不能用于在线择根或后续参数调整。三个开发对象的描述性结果不能解释为总体成功率或新盲测优势。','']
    with (OUT/'OFFLINE_FACTORIAL_EXTENSION_SCORES.md').open('x') as f:f.write('\n'.join(lines))
    receipt.update(status='COMPLETE',finished_utc=utc(),wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,
        denominators=result['denominators'],sealed_files_unchanged=same(stamps),original_stage_i_seal_unchanged=metrics.unchanged(original),
        output_paths=[str((OUT/name).resolve()) for name in ['OFFLINE_FACTORIAL_EXTENSION_SCORES.json','OFFLINE_FACTORIAL_EXTENSION_SCORES.csv','OFFLINE_FACTORIAL_EXTENSION_SCORES.md']])
    if not receipt['sealed_files_unchanged'] or not receipt['original_stage_i_seal_unchanged']:raise ValueError('Input seal changed')
    RECEIPT.write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status='COMPLETE',new_rows=[{k:r[k] for k in ['case','amplitude_normalized_nrmse','otsu_iou','relative_residual','physical_objective','historical_gate_status']} for r in rows if r['source']=='new_R5_extension'],sealed_files=len(stamps))))
if __name__=='__main__':run()
