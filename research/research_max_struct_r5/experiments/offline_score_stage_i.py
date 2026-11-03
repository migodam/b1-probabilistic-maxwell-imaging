"""OFFLINE_DIAGNOSTIC only: score sealed saved materials; no physics imports."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import importlib.util
import json
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
R2=REPO/'research/research_max_struct_r2'
R4=REPO/'research/research_max_struct_r4'
OUT=ROOT/'RESULTS'
SEAL=ROOT/'audit/STAGE_I_SEAL.json'
RECEIPT=ROOT/'audit/OFFLINE_SCORING_RECEIPT.json'
DEV='obs_b20d5a55677d4eca'
CASES=[DEV,'obs_42cb233e4f224bb3','obs_63abafa750fc483c','obs_100bd65498644322']
spec=importlib.util.spec_from_file_location('metrics_only',R4/'experiments/offline_score_common.py')
metrics=importlib.util.module_from_spec(spec)
spec.loader.exec_module(metrics)

def utc():return datetime.now(timezone.utc).isoformat()
def dump(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2,allow_nan=False,ensure_ascii=False);f.write('\n')
def flat(d,prefix=''):
    out={}
    for k,v in d.items():
        key=prefix+k
        if isinstance(v,dict):out.update(flat(v,key+'_'))
        elif isinstance(v,list):out[key]=json.dumps(v,allow_nan=False)
        else:out[key]=v
    return out
def csv_write(path,rows):metrics.write_csv_new(path,[flat(r) for r in rows])
def stamps_match(seal):
    failures=[]
    for e in seal['files']:
        p=ROOT/e['path'];s=p.stat()
        if s.st_size!=e['bytes'] or s.st_mtime_ns!=e['mtime_ns']:failures.append(e['path'])
    return dict(entries=len(seal['files']),matched=len(seal['files'])-len(failures),failures=failures,method='size+mtime_ns; no hashes')

def material(m,t,mask):
    m=np.asarray(m,float).ravel();t=t.ravel()
    if m.shape!=t.shape or not np.isfinite(m).all():raise ValueError('Invalid material vector')
    predicted,threshold=metrics.otsu(m)
    e=m-t;l2=float(np.linalg.norm(e));rmse=l2/np.sqrt(t.size);nrmse=rmse/12
    overlap=metrics.iou(predicted,mask);violation=float(max(0.,.02-float(m.min()),float(m.max())-20.))
    return dict(amplitude_normalized_nrmse=nrmse,RMSE=rmse,error_l2=l2,error_squared_l2=float(e@e),
        otsu_iou=overlap,otsu_threshold=threshold,material_gate=bool(nrmse<=.15),
        image_gate=bool(nrmse<=.15 and overlap>=.7),material_legal=bool(violation<=1e-12),
        material_constraint_violation=violation,material_min=float(m.min()),material_max=float(m.max()))

def changes(before,after,t,mask):
    a=material(before,t,mask);b=material(after,t,mask);step=after-before;direction=t-before
    denom=float(np.linalg.norm(step)*np.linalg.norm(direction))
    return dict(before=a,after=b,delta_nrmse=b['amplitude_normalized_nrmse']-a['amplitude_normalized_nrmse'],
        delta_RMSE=b['RMSE']-a['RMSE'],delta_error_l2=b['error_l2']-a['error_l2'],
        delta_error_squared_l2=b['error_squared_l2']-a['error_squared_l2'],delta_otsu_iou=b['otsu_iou']-a['otsu_iou'],
        actual_step_norm=float(np.linalg.norm(step)),step_dot_GT_minus_before=float(step@direction),
        cosine_step_GT_minus_before=float(np.clip((step@direction)/denom,-1.,1.)) if denom>0 else None)

def run():
    wall,cpu=time.perf_counter(),time.process_time()
    receipt=json.loads(RECEIPT.read_text())
    seal=json.loads(SEAL.read_text());receipt['seal_check_before']=stamps_match(seal)
    if receipt['seal_check_before']['failures']:raise ValueError('Sealed inputs changed')
    receipt.update(status='SCORING',metrics_contract=metrics.frozen_contract(),
        scoring_source=str(Path(__file__).resolve()),truth_read_scope='Only m and mask arrays; no clean Y/noise/seed; already public development cases',
        downstream_use='Truth directions are offline descriptive only, forbidden for online root selection, algorithm parameters or modifying this sealed batch')
    RECEIPT.write_text(json.dumps(receipt,indent=2)+'\n')
    truths={}
    for case in CASES:
        p=(R2/'sealed_truth'/(case+'.npz')) if case==DEV else R4/'sealed_truth/confirmation'/(case+'.npz')
        opened_utc=utc()
        with np.load(p,allow_pickle=False) as f:
            if receipt.get('first_GT_read_utc') is None:
                receipt['first_GT_read_utc']=opened_utc
                RECEIPT.write_text(json.dumps(receipt,indent=2)+'\n')
            keys=list(f.files);t=f['m'].copy().ravel();mask=f['mask'].copy().ravel().astype(bool)
        truths[case]=(t,mask)
        receipt['truth_reads'].append(dict(path=str(p.resolve()),opened_utc=opened_utc,available_keys=keys,arrays_read=['m','mask'],m_shape=list(t.shape),mask_shape=list(mask.shape)))
        RECEIPT.write_text(json.dumps(receipt,indent=2)+'\n')
    pub=metrics.public_metrics(R2/'data_public'/(DEV+'.npz'))
    factorial=[]
    defs=[('new_R5_run','current_only',OUT/'factorial'/(DEV+'__current_only')/'prediction.npz'),
          ('new_R5_run','material_only',OUT/'factorial'/(DEV+'__material_only')/'prediction.npz'),
          ('historical_R4_reference','global',R4/'results/global_common_corrector/prediction.npz'),
          ('historical_R4_reference','onehop',R4/'results/onehop_absorption/prediction.npz'),
          ('historical_R4_reference','self',R4/'results/self_absorption/prediction.npz'),
          ('historical_R4_reference','phase_pi',R4/'results/onehop_phase_pi/prediction.npz')]
    for source,arm,p in defs:
        m,yp=metrics.prediction(p,(24,24));t,mask=truths[DEV]
        row=dict(case=DEV,arm=arm,source=source,prediction_path=str(p.resolve()),**metrics.metrics(m,yp,pub,t,mask,(24,24)))
        row['historical_gate_status']='PASS' if row['joint_gate'] else 'FAIL'
        factorial.append(row)
    context=dict(classification='OFFLINE_DIAGNOSTIC_PUBLIC_DEVELOPMENT_NOT_BLIND',seal_path=str(SEAL.resolve()),seal_utc=seal['sealed_utc'],
        first_GT_read_utc=receipt['first_GT_read_utc'],metrics_contract=metrics.frozen_contract(),new_physical_LU=0,new_complex_RHS=0,
        interpretation='Descriptive fixed states only; no GT selection or algorithm tuning; historical success contract retained without new gates')
    fact=dict(**context,denominators=dict(new_runs=2,new_case_count=1,historical_references=4,historical_case_count=1,
        independent_new_blind_cases=0),rows=factorial)
    dump(OUT/'OFFLINE_FACTORIAL_SCORES.json',fact);csv_write(OUT/'OFFLINE_FACTORIAL_SCORES.csv',factorial)
    lines=['# R5 factorial 离线评分','',
        '仅旧开发弯条 obs_b20d5a55677d4eca：新增2臂、历史引用4臂；独立新盲测对象0。材料/数据阈值完全沿用历史合同：RMSE/12≤.15、预测自身Otsu256 IoU≥.7、公共99%噪声球与材料盒。',
        '',f'全部输入先封存：{seal["sealed_utc"]}。首次GT读取：{receipt["first_GT_read_utc"]}。本评分新增物理LU/RHS为0。',
        '', '| 来源 | 臂 | RMSE/12 | IoU | 相对残差 | Phi | 历史联合门 |','|---|---|---:|---:|---:|---:|---|']
    for r in factorial:lines.append(f'|{r["source"]}|{r["arm"]}|{r["amplitude_normalized_nrmse"]:.9f}|{r["otsu_iou"]:.9f}|{r["relative_residual"]:.9f}|{r["physical_objective"]:.9g}|{r["historical_gate_status"]}|')
    lines+=['','新增2臂中material_only满足历史合同、current_only不满足。这个单一公开开发对象上的位置消融不等于新分布成功率或普遍优势。GT方向仅供离线解释，不能用于后续在线择根或算法参数。历史引用没有重跑，也不计作新增运行。','', '详细输入绝对路径、逐项噪声/材料/联合门见 JSON/CSV。','']
    with (OUT/'OFFLINE_FACTORIAL_SCORES.md').open('x') as f:f.write('\n'.join(lines))
    rows=[];cf_count=0
    for p in sorted((OUT/'trajectory_diagnostics').glob('obs*/checkpoint_s*_i*.npz')):
        runid=p.parent.name;case,arm=runid.split('__');t,mask=truths[case]
        summary=json.loads(p.with_suffix('.json').read_text())
        with np.load(p,allow_pickle=False) as f:
            before=f['m_before'].copy();after=f['m_after'].copy()
            row=dict(case=case,arm=arm,source='sealed_replay_checkpoint',run_id=runid,checkpoint=p.stem,
                stage=summary['stage'],iteration=summary['iteration'],diagnostic_path=str(p.resolve()),
                actual_step_available=True,**changes(before,after,t,mask))
            counter={}
            for kind in ('ordinary','self','plus','minus'):
                key='counterfactual_m_'+kind
                if key in f.files:
                    score=changes(before,f[key].copy(),t,mask)
                    counter[kind]=dict(**score,scope='Frozen-J/Q material readout from same m_before; not actual trajectory',
                        saved_physical=summary['counterfactual_readouts'][kind]['Phi'])
                    cf_count+=1
                else:counter[kind]=None
            row['counterfactuals']=counter
            row['counterfactual_missing_reason']=None if any(v is not None for v in counter.values()) else 'Not generated at nonendpoint under frozen design'
        rows.append(row)
    actual_rows=len(rows)
    for p in sorted((OUT/'trajectory_diagnostics').glob('legacy*/stage_*.npz')):
        t,mask=truths[DEV];summary=json.loads(p.with_suffix('.json').read_text())
        with np.load(p,allow_pickle=False) as f:
            state=f['m_state'].copy();counter={}
            for kind in ('ordinary','self','plus','minus'):
                key='counterfactual_m_'+kind
                if key in f.files:
                    counter[kind]=dict(**changes(state,f[key].copy(),t,mask),scope='Counterfactual from legacy fixed m_state, not missing actual step');cf_count+=1
                else:counter[kind]=None
            row=dict(case=DEV,arm=p.parent.name.removeprefix('legacy_'),source='historical_fixed_state',run_id=p.parent.name,
                checkpoint=p.stem,stage=summary['stage'],iteration=500,diagnostic_path=str(p.resolve()),
                actual_step_available=False,before=None,after=None,state=material(state,t,mask),
                delta_nrmse=None,delta_RMSE=None,delta_error_l2=None,delta_error_squared_l2=None,delta_otsu_iou=None,
                actual_step_norm=None,step_dot_GT_minus_before=None,cosine_step_GT_minus_before=None,
                missing_reason='No saved adjacent single-step before state; cross-stage subtraction is forbidden',counterfactuals=counter)
            rows.append(row)
    if actual_rows!=168 or len(rows)-actual_rows!=6 or cf_count!=120:raise ValueError('Frozen denominator mismatch')
    output=dict(**context,denominators=dict(actual_trajectory_checkpoints=168,replay_runs=8,case_count=4,legacy_fixed_states=6,
        actual_single_steps_available=168,legacy_actual_steps_available=0,counterfactual_materials_scored=cf_count),
        sign_convention='delta=after-before; negative material-error delta means improvement. Cosine uses actual single-step delta_m against GT-m_before. Undefined zero-norm cosine is null.',rows=rows)
    dump(OUT/'OFFLINE_TRAJECTORY_SCORES.json',output);csv_write(OUT/'OFFLINE_TRAJECTORY_SCORES.csv',rows)
    receipt.update(status='COMPLETE',finished_utc=utc(),wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,
        denominators=output['denominators'],factorial_denominators=fact['denominators'],seal_check_after=stamps_match(seal),
        physical_LU=0,physical_complex_RHS=0,physics_imported=False,tests_run=0,
        output_paths=[str((OUT/n).resolve()) for n in ['OFFLINE_TRAJECTORY_SCORES.json','OFFLINE_TRAJECTORY_SCORES.csv','OFFLINE_FACTORIAL_SCORES.json','OFFLINE_FACTORIAL_SCORES.csv','OFFLINE_FACTORIAL_SCORES.md']])
    if receipt['seal_check_after']['failures']:raise ValueError('Sealed inputs changed during scoring')
    RECEIPT.write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status='COMPLETE',denominators=output['denominators'],factorial=fact['denominators'],seal_unchanged=True)))

if __name__=='__main__':run()
