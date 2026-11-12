"""Aggregate existing post-seal score tables; no material NPZ, truth or solver calls.

Independent campaigns require a parent authorization JSON with
  authorized_aggregation=true, campaigns=[...]
and each source SCORING_RECEIPT must point to an authorized valid prediction seal.
Optional family labels additionally require family_labels_authorized=true.
"""
import argparse,csv,itertools,json,re
from pathlib import Path

METRICS=['RMSE','NRMSE','boxrange_NRMSE','relativeL2','IoU','mean_amplitude_bias','mean_amplitude_bias_normalized','D','R','Phi','residual_over_delta']
VIEWS=['raw','corrected','primary','secondary','feasible_raw_retain','min_change_paid_pool']


def load(path):return json.loads(Path(path).read_text())


def decode(value):
    if value in ['',None]:return None
    try:return json.loads(value)
    except (ValueError,TypeError):return value


def source_table(folder,name):
    path=folder/(name+'.json')
    if path.is_file():return load(path)
    path=folder/(name+'.csv')
    if path.is_file():
        with path.open(newline='') as f:return [{k:decode(v) for k,v in r.items()} for r in csv.DictReader(f)]
    return []


def table(out,name,rows):
    (out/(name+'.json')).write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (out/(name+'.csv')).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in rows)


def development(name):return name.startswith(('development_','representation_'))


def paired(block,records,method_a,method_b,expected_cases=None,views=VIEWS):
    lookup={(r['case_id'],r['method'],r['view']):r for r in records}
    cases=expected_cases or sorted({r['case_id'] for r in records if r['method'] in [method_a,method_b]});rows=[]
    for view,cid in itertools.product(views,cases):
        aa=lookup.get((cid,method_a,view),{});bb=lookup.get((cid,method_b,view),{});available=bool(aa.get('available') and bb.get('available'))
        row=dict(block=block,view=view,case_id=cid,object_id=aa.get('object_id',bb.get('object_id',cid)),method=method_a,reference=method_b,paired_available=available,method_available=bool(aa.get('available')),reference_available=bool(bb.get('available')),host=aa.get('host',bb.get('host')),domain=aa.get('domain',bb.get('domain')))
        for key in METRICS:
            a,b=aa.get(key),bb.get(key);row['method_'+key]=a;row['reference_'+key]=b;row['delta_'+key]=a-b if available and a is not None and b is not None else None
        rows.append(row)
    return rows


def paired_summary(rows):
    import numpy as np
    from scipy.stats import wilcoxon
    groups={}
    for r in rows:groups.setdefault((r['block'],r['view'],r['method'],r['reference'],r.get('host'),r.get('domain')),[]).append(r)
    output=[]
    for key,items in groups.items():
        objects={}
        for r in items:
            if r.get('delta_NRMSE') is not None:objects.setdefault(r['object_id'],[]).append(r['delta_NRMSE'])
        values=np.array([np.mean(v) for v in objects.values()]);rng=np.random.default_rng(131313)
        boot=np.mean(values[rng.integers(0,len(values),size=(4000,len(values)))],axis=1) if len(values) else None
        try:p=float(wilcoxon(values).pvalue) if len(values) and np.any(values!=0) else 1. if len(values) else None
        except ValueError:p=None
        output.append(dict(block=key[0],view=key[1],method=key[2],reference=key[3],host=key[4],domain=key[5],expected_case_pairs=len(items),available_case_pairs=sum(r.get('delta_NRMSE') is not None for r in items),material_objects=len(values),NRMSE_paired_mean=float(values.mean()) if len(values) else None,bootstrap95_low=float(np.quantile(boot,.025)) if len(values) else None,bootstrap95_high=float(np.quantile(boot,.975)) if len(values) else None,error_decreases=int(np.count_nonzero(values<0)),error_increases=int(np.count_nonzero(values>0)),unchanged=int(np.count_nonzero(values==0)),wilcoxon_descriptive_uncorrected_p=p,bootstrap_seed=131313,bootstrap_replicates=4000,statistical_unit='material object; repeated conditions averaged within object',qualification='descriptive unadjusted; no new success threshold'))
    return output


def summaries(records):
    import numpy as np
    groups={}
    for r in records:
        strata=[('all','all')]
        for key in ['grid_contract','acquisition_noise_stratum','family']:
            if r.get(key) is not None:strata.append((key,str(r[key])))
        for s,v in strata:groups.setdefault((r['campaign'],r['host'],r['domain'],r['method'],r['view'],s,v),[]).append(r)
    result=[]
    for key,items in groups.items():
        row=dict(campaign=key[0],host=key[1],domain=key[2],method=key[3],view=key[4],stratum_type=key[5],stratum=key[6],expected_rows=len(items),available_rows=sum(bool(r.get('available')) for r in items),unique_material_objects=len({r.get('object_id',r['case_id']) for r in items}),missing_or_unavailable_case_ids=[r['case_id'] for r in items if not r.get('available')])
        for metric in METRICS:
            values=[r[metric] for r in items if r.get('available') and r.get(metric) is not None];row[metric+'_mean_available']=float(np.mean(values)) if values else None;row[metric+'_count']=len(values)
        result.append(row)
    return result


def conditional(root,campaigns,records,costs):
    if 'confirmation_core48' not in campaigns or 'confirmation_rescue48' not in campaigns:return [],[]
    config=load(root/'configs/confirmation_rescue48.json');policy=load(root/'configs/conditional_rescue_policy.json');parent={(r['case_id'],r['view']):r for r in records if r['campaign']=='confirmation_core48' and r['method']==policy['parent_method']};rescues={(r['case_id'],r['method'],r['view']):r for r in records if r['campaign']=='confirmation_rescue48'}
    clook={(r['campaign'],r['case_id'],r['method']):r for r in costs};rows=[];costrows=[]
    declared=config['all48_deployment_rows']
    if len(declared)!=48:raise ValueError('Conditional deployment must account for all48 declared parents')
    for strategy in policy['methods']:
        for case in declared:
            cid=case['case_id'];trigger=bool(case['trigger']);base=parent.get((cid,'primary'),{});pcost=clook.get(('confirmation_core48',cid,policy['parent_method']),{});rcost=clook.get(('confirmation_rescue48',cid,strategy),{}) if trigger else{}
            pcpu=pcost.get('external_child_CPU_seconds');rcpu=rcost.get('external_child_CPU_seconds') if trigger else 0.
            cpu=pcpu+rcpu if pcpu is not None and rcpu is not None else None
            costrow=dict(campaign='conditional_deployment48',case_id=cid,method=strategy,host=base.get('host','i9_remote'),domain=base.get('domain','scalar_tm_2d/fresh48'),triggered=trigger,parent_external_CPU_seconds=pcpu,rescue_external_CPU_seconds=rcpu,deployment_external_CPU_seconds=cpu,bypass=not trigger,cost_rule='charge each48 parent once per deployed strategy; add only actual triggered rescue. Do not sum competing strategy totals as experimental consumption.')
            for metric in ['physical_LU','physical_complex_RHS']:
                aa=pcost.get(metric);bb=rcost.get(metric) if trigger else 0;costrow['deployment_'+metric]=aa+bb if aa is not None and bb is not None else None
            costrows.append(costrow)
            for view in VIEWS:
                chosen=rescues.get((cid,strategy,view),{}) if trigger else base
                row={k:v for k,v in chosen.items() if k!='metadata'}
                row.update(campaign='conditional_deployment48',physics_kind=base.get('physics_kind','scalar_tm_2d'),cohort='fresh48_conditional_deployment',case_id=cid,object_id=base.get('object_id',cid),method=strategy,view=view,host=costrow['host'],domain=costrow['domain'],available=bool(chosen.get('available')),triggered=trigger,bypass=not trigger,source_campaign='confirmation_rescue48' if trigger else'confirmation_core48',source_method=strategy if trigger else policy['parent_method'],source_view=view if trigger else'primary',external_child_CPU_seconds=cpu,parent_primary_NRMSE=base.get('NRMSE'),availability_reason=chosen.get('availability_reason','TRIGGERED_RESCUE_MISSING' if trigger else'PARENT_MISSING'))
                rows.append(row)
    return rows,costrows


def figures(out,records,pairs,stats):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fd=out/'FIGURES';fd.mkdir()
    def save(fig,name):fig.savefig(fd/(name+'.png'),dpi=160,bbox_inches='tight');fig.savefig(fd/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
    for block in sorted({r['block'] for r in stats}):
        items=[r for r in stats if r['block']==block and r['view']=='primary' and r['material_objects']]
        for page,start in enumerate(range(0,len(items),18),1):
            batch=items[start:start+18];fig,axes=plt.subplots(1,2,figsize=(14,max(4,.4*len(batch)+1)));y=np.arange(len(batch));labels=[r['method']+' vs '+r['reference'] for r in batch]
            means=np.array([r['NRMSE_paired_mean'] for r in batch]);lo=np.array([r['bootstrap95_low'] for r in batch]);hi=np.array([r['bootstrap95_high'] for r in batch]);axes[0].errorbar(means,y,xerr=np.vstack([means-lo,hi-means]),fmt='o',ms=3);axes[0].axvline(0,color='gray',ls=':');axes[0].set_yticks(y,labels,fontsize=8);axes[0].set_xlabel('Paired amplitude-NRMSE difference; object bootstrap95%')
            axes[1].barh(y,[-r['error_decreases'] for r in batch],label='error decreased');axes[1].barh(y,[r['error_increases'] for r in batch],label='error increased');axes[1].set_yticks(y,labels,fontsize=8);axes[1].set_xlabel('Objects with signed changes; no success gate');axes[1].legend(fontsize=8);fig.suptitle(block+' · primary; descriptive unadjusted');fig.tight_layout();save(fig,'paired_changes_'+re.sub('[^A-Za-z0-9_-]+','_',block)+f'_{page:02d}')
    groups={}
    for r in records:
        if r['view']=='primary' and r.get('available') and r.get('NRMSE') is not None:groups.setdefault((r['campaign'],r['host'],r['domain']),[]).append(r)
    for (campaign,host,domain),items in groups.items():
        fig,ax=plt.subplots(figsize=(8,5))
        for method in sorted({r['method'] for r in items}):
            vals=[r for r in items if r['method']==method and r.get('external_child_CPU_seconds') is not None]
            if vals:ax.scatter([r['external_child_CPU_seconds'] for r in vals],[r['NRMSE'] for r in vals],s=20,alpha=.7,label=method)
        ax.set_xlabel('Actual whole-child CPU (s); same listed host only');ax.set_ylabel('Primary amplitude NRMSE');ax.set_title(campaign+' / '+host+'\n'+domain);handles,_=ax.get_legend_handles_labels()
        if handles:ax.legend(fontsize=7,ncol=2)
        ax.grid(alpha=.2);fig.tight_layout();save(fig,'cost_material_'+re.sub('[^A-Za-z0-9_-]+','_',campaign+'_'+host+'_'+domain))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--campaign',action='append',required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--authorization',type=Path);p.add_argument('--family-metadata',type=Path);p.add_argument('--no-figures',action='store_true');a=p.parse_args();root=a.root.resolve();campaigns=a.campaign
    
    if len(set(campaigns))!=len(campaigns):raise ValueError('Duplicate source campaign names would duplicate runs/material rows')
    independent=[name for name in campaigns if not development(name)];authorization=load(a.authorization) if a.authorization else{}
    # This guard precedes every source score/config/private-label read.
    if independent and not(authorization.get('authorized_aggregation') is True and set(independent)<=set(authorization.get('campaigns',[]))):raise PermissionError('Parent aggregation authorization listing every independent campaign required before any independent score read.')
    if a.family_metadata and (authorization.get('family_labels_authorized') is not True or 'confirmation_core48' not in independent):raise PermissionError('Explicit parent family-label authorization and included scored core48 campaign are both required.')
    from seal_results import validate_seal
    for name in campaigns:
        folder=root/'analysis'/name;receipt=load(folder/'SCORING_RECEIPT.json')
        if receipt['campaign']!=name:raise ValueError('Source scoring campaign mismatch')
        validate_seal(load(receipt['authorized_seal']),require_authorization=True)
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);records=[];costs=[];sources=[];cost_context=[];attempt_history=[]
    for name in campaigns:
        folder=root/'analysis'/name;source_rows=source_table(folder,'SCORES');source_costs=source_table(folder,'COST_LEDGER');host='Apple_local' if 'local' in name else'i9_remote';clook={r['task_id']:r for r in source_costs}
        if(folder/'COST_TOTALS.json').is_file():cost_context.append(dict(campaign=name,host=host,source_cost_totals=load(folder/'COST_TOTALS.json'),scope='Original totals preserved; not added to per-run counters or other hosts'))
        attempt_history.extend(dict(r,campaign=name,host=host) for r in source_table(folder,'QUEUE_ATTEMPT_HISTORY'))
        sources.append(dict(campaign=name,source_analysis=str(folder),scoring_receipt=str(folder/'SCORING_RECEIPT.json'),score_rows=len(source_rows),host=host))
        for original in source_rows:
            row={k:v for k,v in original.items() if k!='metadata'};row.update(campaign=name,source_analysis=str(folder),host=host,domain=row.get('domain',row.get('physics_kind','unknown')+'/'+row.get('cohort','unknown')),external_child_CPU_seconds=clook.get(row['task_id'],{}).get('external_child_CPU_seconds'));records.append(row)
        costs.extend(dict(r,campaign=name,host=host) for r in source_costs)
    fresh=[r for r in records if r['case_id'].startswith('r13_fresh_')]
    if fresh:
        acquisition={r['case_id']:r for r in load(root/'data_public/fresh48/MANIFEST.json')['case_files']}
        for row in fresh:
            info=acquisition[row['case_id']];row['grid_contract']='matched' if info['forward_side']==info['inverse_side'] else'crossgrid_model_mismatch';row['acquisition_noise_stratum']=f"sources{info['sources']}_receivers{info['receivers']}_noise{info['relative_noise_RMS']}"
    if a.family_metadata:
        families={r['case_id']:r['family'] for r in load(a.family_metadata)['cases']}
        for row in fresh:row['family']=families.get(row['case_id'])
    deployment,deployment_costs=conditional(root,campaigns,records,costs);records.extend(deployment)
    allpairs=[]
    for campaign,host,domain in sorted({(r['campaign'],r['host'],r['domain']) for r in records}):
        items=[r for r in records if(r['campaign'],r['host'],r['domain'])==(campaign,host,domain)];methods=sorted({r['method'] for r in items})
        for aa,bb in itertools.combinations(methods,2):allpairs.extend(paired(campaign+'/'+domain,items,aa,bb))
    fixed=[f'r13_fresh_{i:04d}' for i in range(0,48,3)]
    if 'confirmation_core48' in campaigns and 'confirmation_mechanism16' in campaigns:
        core=[r for r in records if r['campaign']=='confirmation_core48' and r['case_id'] in fixed];mech=[r for r in records if r['campaign']=='confirmation_mechanism16' and r['case_id'] in fixed]
        for method in ['dr','g19','direct']:allpairs.extend(paired('NN_common_predeclared16',core+mech,'nn_'+method,method,expected_cases=fixed))
        online_sources={'confirmation_core48','confirmation_mechanism16','confirmation_online_controls16'}
        online_methods=['dr','online_event','online_random','online_double','online_noevent','online_reset']
        online=[dict(r,source_campaign=r['campaign'],campaign='online_controls_common16',domain='scalar_tm_2d/online_controls_common16',cohort='online_controls_common16',comparison_view_only_not_new_run=True) for r in records if r['campaign'] in online_sources and r['case_id'] in fixed and r['method'] in online_methods]
        for method in online_methods[1:]:
            allpairs.extend(paired('online_controls_common16_vs_DR',online,method,'dr',expected_cases=fixed))
            if method!='online_event':allpairs.extend(paired('online_controls_common16_vs_event',online,method,'online_event',expected_cases=fixed))
        present=sorted({r['method'] for r in online});table(out,'ONLINE_COMMON16_ROWS',online)
        table(out,'ONLINE_CONTROL_COVERAGE',[dict(requested_controls=online_methods[1:],source_campaigns_selected=sorted(online_sources & set(campaigns)),present_methods=present,missing_controls=[m for m in online_methods[1:] if m not in present],expected_objects=16,policy='Join actual core+mechanism+C3 score rows. Every declared comparison keeps16cases and all stage views; missing rows stay unavailable, never fabricated.')])
    local_common=[]
    local_sources={'confirmation_local_reference16','confirmation_local_trf16'}
    if local_sources<=set(campaigns):
        local_methods=['cie','dual_cie','g19','dr','direct','physical_globaldual','physical_trf']
        local_common=[dict(r,source_campaign=r['campaign'],campaign='local_reference_common16',domain='scalar_tm_2d/local_reference_common16',cohort='local_reference_common16',host='Apple_local',comparison_view_only_not_new_run=True) for r in records if r['campaign'] in local_sources and r['case_id'] in fixed and r['method'] in local_methods]
        for aa,bb in itertools.combinations(local_methods,2):allpairs.extend(paired('local_reference_common16_360CPU',local_common,aa,bb,expected_cases=fixed))
        table(out,'LOCAL_REFERENCE_COMMON16_ROWS',local_common);table(out,'LOCAL_REFERENCE_COMMON16_SUMMARY',summaries(local_common))
        table(out,'LOCAL_REFERENCE_COVERAGE',[dict(source_campaigns=sorted(local_sources),expected_objects=16,expected_methods=local_methods,expected_case_ids=fixed,excluded_exposed_case='B__r12_confirm24_0010',budget_CPU_seconds=360,host='Apple_local',policy='Same original public observations/uniform initialization; comparison views reuse source rows, not additional objects or CPU consumption.')])
    if deployment:
        originals=[]
        for r in records:
            if r['campaign']=='confirmation_core48' and r['method']=='g19' and r['view']=='primary':
                for view in VIEWS:originals.append(dict(r,method='original_g19_primary',view=view))
        for method in sorted({r['method'] for r in deployment}):allpairs.extend(paired('conditional_ALL48_vs_original_g19',deployment+originals,method,'original_g19_primary',expected_cases=[f'r13_fresh_{i:04d}' for i in range(48)]))
    stats=paired_summary(allpairs);summary=summaries(records)
    for name,rows in [('ALL_ESTIMATOR_ROWS',records),('STRATIFIED_SUMMARY',summary),('PAIRED_OBJECT_CHANGES',allpairs),('PAIRED_OBJECT_BOOTSTRAP',stats),('CONDITIONAL_DEPLOYMENT_ALL48',deployment),('CONDITIONAL_DEPLOYMENT_COSTS',deployment_costs),('SOURCE_RUN_COSTS',costs),('SOURCE_COST_CONTEXT',cost_context),('SOURCE_QUEUE_ATTEMPT_HISTORY',attempt_history),('SOURCE_ANALYSES',sources)]:table(out,name,rows)
    if 'confirmation_mechanism16' in campaigns and(root/'delegated/evidence/nnprep/PREPROCESS_RECEIPT.json').is_file():
        nn=load(root/'delegated/evidence/nnprep/PREPROCESS_RECEIPT.json');(out/'NN_PREPROCESS_COST_CONTEXT.json').write_text(json.dumps(dict(Apple_preprocessing_CPU_subtotal=nn.get('process_CPU_seconds'),Apple_preprocessing_wall_subtotal=nn.get('process_wall_seconds'),batch_objects=nn['objects'],old_training_CPU_seconds=None,old_training_GPU_hours=None,scope='Charge preprocessing once to research batch; each independently deployed16objectNN strategy must disclose applicable setup. Do not add Apple CPU seconds to i9 CPU as a same-host speed figure.'),indent=2)+'\n')
    if not a.no_figures:figures(out,records+local_common,allpairs,stats)
    identity={(r.get('physics_kind',r['domain']),r.get('object_id',r['case_id'])) for r in records};receipt=dict(status='TABLE_AGGREGATION_COMPLETE_PARENT_INTERPRETATION_PENDING',campaigns=campaigns,estimator_rows=len(records),unique_material_identities=len(identity),paired_rows=len(allpairs),new_physical_forward=0,new_adjoint=0,new_optimization=0,material_or_truth_NPZ_opened=False,family_labels_read=bool(a.family_metadata),source_campaign_counts_are_not_independent_material_counts=True,host_CPU_not_pooled=True)
    (out/'AGGREGATION_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');(out/'REPORT.md').write_text('# R13 表格聚合附件\n\n来源：'+', '.join(campaigns)+'。所有数值来自已授权封存后评分表；本脚本未打开材料/真值NPZ或调用求解器。\n\n保留全部估计器与缺失分母，按物理域、host、网格及采集条件分层。样本单位为材料对象，重复运行/初值/策略不能当作新增独立对象。配对变化与bootstrap是描述性统计，不新增成功门。\n\n条件策略使用全48原G19分母：未触发保留原primary，触发使用实际rescue估计器；缺失触发结果保持不可用。每种部署策略收费包括全48parent加实际触发rescue；竞争策略总费用不能相加冒充实验消费。Apple本地参考与i9分开。\n\n原始表：ALL_ESTIMATOR_ROWS、STRATIFIED_SUMMARY、PAIRED_OBJECT_CHANGES、PAIRED_OBJECT_BOOTSTRAP及SOURCE_RUN_COSTS。图只绘实际评分/费用记录，没有重构或模拟材料图。最终科学解释、负结果讨论与发布由父线程完成。\n');print(json.dumps(receipt))
if __name__=='__main__':main()
