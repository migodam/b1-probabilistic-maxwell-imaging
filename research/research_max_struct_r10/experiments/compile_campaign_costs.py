"""R10 cost-only union of canonical receipts. No arrays, GT, solve or scorer.

Whole-process CPU rows are additive; operation counters/nested timers are
separate views. Refresh branch producers first, then rerun this script.
"""
from pathlib import Path
from collections import defaultdict
from datetime import datetime,timezone
import json,time
ROOT=Path(__file__).resolve().parents[1]


def read(p):return json.loads(Path(p).read_text())
def utc():return datetime.now(timezone.utc).isoformat()
def finite_cpu(x):return float(x) if isinstance(x,(int,float)) and x>=0 else None
def rel(p):
    p=Path(p)
    try:return str(p.relative_to(ROOT))
    except ValueError:return str(p)
def platform_name(x):
    if isinstance(x,dict):x=x.get('system') or x.get('platform') or x.get('host')
    x=str(x or '')
    if x.lower() in ('darwin','mac','macos'):return 'Mac'
    if x.lower() in ('windows','xinan','win32','xinan_windows_i9'):return 'XINAN'
    return x if x else 'UNKNOWN'


def main():
    started,cpu=time.perf_counter(),time.process_time()
    rows=[];views=[];unknown=[];producers=[];covered=set();uids=set()
    def add(branch,identity,value,source,ledger,scope,status=None,platform=None,activity='compute',**extra):
        uid=branch+':'+str(identity)
        if uid in uids:
            unknown.append(dict(source=source,reason='Duplicate producer identity excluded; root must resolve',identity=uid));return
        uids.add(uid);value=finite_cpu(value)
        rows.append(dict(identity=uid,branch=branch,activity=activity,CPU_seconds=value,
            platform=platform_name(platform),source=source,source_ledger=ledger,
            status=status or 'RECORDED',scope=scope,additive=True,**extra))
        if value is None:unknown.append(dict(source=source,identity=uid,reason='CPU not recorded; not zero'))
        covered.add(str(source).replace('\\','/'))
    def producer(path,branch):
        p=ROOT/path
        if not p.exists():
            producers.append(dict(branch=branch,source=path,status='MISSING_OR_NOT_YET_COMPLETE'))
            unknown.append(dict(source=path,reason='Producer absent; rerun after branch completion'));return None
        d=read(p);producers.append(dict(branch=branch,source=path,status=d.get('status','RECORDED'),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns));return d
    # This audit already unions historical energy, DEV attempts/failures and
    # old scoring. Do NOT add its energy ledger or raw attempt receipts again.
    source='delegated/cost_audit/COMPLETED_COSTS.json';d=producer(source,'historical_DEV')
    if d:
        for x in d['rows']:
            subgroup=x['group'];activity='energy_or_scoring' if subgroup in ('historical_energy','historical_postprocessing') else 'DEV'
            add('historical_DEV',x['id'],x.get('recorded_CPU_seconds'),x['source'],source,x.get('scope'),x.get('status'),
                (x.get('extra') or {}).get('execution_platform'),activity,producer_role=x.get('cost_aggregation_role'))
        for key in ('physical_counts_by_matrix_dimension','shifted_auxiliary_counts_by_matrix_dimension','denominator_small_counts_by_dimension'):
            views.append(dict(branch='historical_DEV',source=source,kind=key,counts=d.get(key),additive=False))
        unknown.extend(dict(branch='historical_DEV',**u) for u in d.get('unknowns',[]))
        calls=ROOT/'delegated/cost_audit/BUILD_COSTS.jsonl'
        if calls.exists():
            a=[json.loads(line) for line in calls.read_text().splitlines() if line]
            add('historical_DEV','cost_audit_compiler_calls',sum(x['CPU_seconds'] for x in a),rel(calls),source,'Disjoint compiler invocations, not already in DEV compute sum',platform='Mac',activity='metadata')
    # Canonical 2D ledgers include generation once (baseline only), complete
    # child GetProcessTimes and engineering components. Internal R8/aux timers
    # remain views, never additional CPU charges.
    for branch,path in [('twod_baseline','delegated/twod/COST_LEDGER.json'),('twod_matrix_candidates','delegated/twod/DENOMINATOR_COST_LEDGER.json'),
                        ('right_exposed_diagnostic','delegated/twod/RIGHT_COST_LEDGER.json'),('conditional_confirmation','delegated/twod/CONDITIONAL_COST_LEDGER.json')]:
        d=producer(path,branch)
        if not d:continue
        for x in d.get('rows',[]):
            if x.get('additive') is not True:
                views.append(dict(branch=branch,source=path,row=x,additive=False));continue
            ident=x.get('id',x.get('identity'))
            activity=('generation' if 'generation' in str(ident) else 'packaging_transport' if any(w in str(ident).lower() for w in ('package','transport','ssh','scp','archive','fetch')) else 'inversion' if x.get('LU') is not None else 'control_or_metadata')
            add(branch,ident,x.get('CPU_seconds'),x.get('artifact',path),path,x.get('scope'),x.get('status'),x.get('platform'),activity)
            if x.get('LU') is not None:
                views.append(dict(branch=branch,identity=ident,kind='original_physical_Maxwell',matrix_dimension=576,LU=x.get('LU'),complex_RHS=x.get('complex_RHS'),source=path,additive=False))
        views.append(dict(branch=branch,source=path,kind='nested_internal_method_timers',rows=d.get('nested_views',d.get('deployment_views',[])),additive=False))
        unknown.extend(dict(branch=branch,**u) for u in d.get('unknown_costs',[]))
    # Vector producer owns process/transport/probe/failure accounting; never
    # re-add raw per-arm costs. FINAL supersedes CURRENT, not an extra bill.
    vector='RESULTS/vector/ACCOUNTING_FINAL.json' if (ROOT/'RESULTS/vector/ACCOUNTING_FINAL.json').exists() else 'RESULTS/vector/ACCOUNTING_CURRENT.json'
    d=producer(vector,'vector')
    if d:
        for x in d['rows']:
            if x.get('additive') is not True:
                views.append(dict(branch='vector',source=vector,row=x,additive=False));continue
            activity='packaging_transport' if any(w in x['identity'].lower() for w in ('transport','ssh','scp')) else 'vector_process'
            add('vector',x['identity'],x.get('CPU_seconds'),x['source'],vector,x.get('scope'),d['status'],x.get('platform'),activity)
            if x.get('counts') is not None:views.append(dict(branch='vector',source=x['source'],kind='producer_dimensioned_counts',counts=x['counts'],additive=False))
        views.append(dict(branch='vector',kind='producer_nested_views',rows=d.get('nonadditive_views',[]),counter_scope=d.get('counter_scope'),additive=False))
        unknown.extend(dict(branch='vector',**x) for x in d.get('unknowns',[]))
    continuation='RESULTS/vector_continuation_v04/ACCOUNTING_FINAL.json'
    if not (ROOT/continuation).exists() and (ROOT/'RESULTS/vector_continuation_v04/ACCOUNTING_CURRENT.json').exists():
        continuation='RESULTS/vector_continuation_v04/ACCOUNTING_CURRENT.json'
    d=producer(continuation,'vector_exposed_continuation')
    if d:
        for x in d.get('rows',[]):
            if x.get('additive') is True:
                add('vector_exposed_continuation',x.get('identity',x.get('id')),x.get('CPU_seconds'),x.get('source',continuation),continuation,x.get('scope'),d.get('status'),x.get('platform'),'vector_process')
                if x.get('counts') is not None:views.append(dict(branch='vector_exposed_continuation',source=x.get('source',continuation),kind='producer_dimensioned_counts',matrix_dimension=5184 if x['counts'].get('physical_LU') else None,dimension_source='fixed grid12 vector field3*12^3, producer protocol',counts=x['counts'],additive=False))
            else:views.append(dict(branch='vector_exposed_continuation',source=continuation,row=x,additive=False))
        views.append(dict(branch='vector_exposed_continuation',kind='producer_nested_views',rows=d.get('nonadditive_views',[]),counter_scope=d.get('counter_scope'),additive=False))
        unknown.extend(dict(branch='vector_exposed_continuation',**u) for u in d.get('unknowns',[]))
    # Canonical postseal indexes supersede their raw receipt windows.
    for folder in (ROOT/'delegated').iterdir():
        if not folder.is_dir() or not any(w in folder.name for w in ('_score','_energy')):continue
        for p in sorted(folder.rglob('*POSTSEAL*COST*.json')):
            name=rel(p);index=producer(name,'postseal_analysis')
            if not index:continue
            for x in index.get('rows',[]):
                src=x.get('source',name)
                if x.get('additive') is not True:
                    views.append(dict(branch='postseal_analysis',source=name,row=x,additive=False));continue
                if src in covered:continue
                add('postseal_analysis',src,x.get('recorded_CPU_seconds'),src,name,x.get('scope'),'RECORDED',x.get('platform'),'offline_postprocessing')
                counter_source=ROOT/src
                if counter_source.is_file() and counter_source.name.endswith('RECEIPT.json') and not any(w in counter_source.parts for w in ('sealed_truth','data_private')):
                    metadata=read(counter_source)
                    counters={k:v for k,v in metadata.items() if any(w in k for w in ('LU','RHS','SVD','dimension','new_forward','new_adjoint')) and isinstance(v,(int,float,str,type(None)))}
                    views.append(dict(branch='postseal_analysis',source=src,kind='offline_operation_counters',counts=counters,additive=False))
    # Score/energy process RECEIPTs only (never SCORES/GT containers). Already
    # covered historical audit sources are skipped. Each receipt's whole timer
    # includes its nested geometry, figures and offline GT-forward windows.
    extras={ROOT/'delegated/input_output_projection/RECEIPT.json',ROOT/'delegated/vector_score/scored/NOISE_MODEL_COMPATIBILITY.json'}
    for folder in (ROOT/'delegated').iterdir():
        if folder.is_dir() and any(w in folder.name for w in ('_score','_energy')):
            extras.update(folder.rglob('RECEIPT.json'))
    for p in sorted(extras):
        name=rel(p)
        if name in covered:continue
        d=producer(name,'scoring_energy')
        if not d:continue
        value=next((d[k] for k in ('CPU_seconds','cpu_seconds','scoring_CPU_seconds') if k in d),None)
        add('scoring_energy',name,value,name,name,d.get('cost_scope',d.get('scope','Whole recorded scoring/energy process; no nested timers re-added')),d.get('status'),d.get('execution_platform'),'offline_postprocessing')
        counters={k:v for k,v in d.items() if any(w in k for w in ('LU','RHS','SVD','dimension','new_forward','new_adjoint')) and isinstance(v,(int,float,str,type(None)))}
        if name=='delegated/input_output_projection/RECEIPT.json':counters['auxiliary_small_matrix_dimension']=19
        views.append(dict(branch='scoring_energy',source=name,kind='offline_operation_counters',counts=counters,additive=False))
    # Explicit missing known producers remain pending, never free or failed.
    for name in ['delegated/twod_score/scored/RECEIPT.json','delegated/twod_energy/scored/RECEIPT.json','delegated/vector_score/scored/RECEIPT.json']:
        if not (ROOT/name).exists():producers.append(dict(branch='scoring_energy',source=name,status='PENDING'))
    for name in ['delegated/twod_score/PREPARATION.json','delegated/conditional_confirmation_score/PREPARATION.json']:
        if (ROOT/name).exists():unknown.append(dict(source=name,reason='Preparation CPU not instrumented; not zero or another scientific inverse'))
    delivery='audit/DELIVERY_CHECK.json';d=producer(delivery,'delivery_metadata')
    if d:
        add('delivery_metadata','final_delivery_check',d.get('CPU_seconds'),delivery,delivery,'Root recorded syntax/local-link delivery check window; no rerun by compiler',d.get('status','RECORDED'),d.get('execution_platform'),'metadata')
        unknown.append(dict(source=delivery,reason='Earlier uninstrumented editing CPU remains null; no invented charge or zero'))
    path='audit/PACKAGING_RECEIPT.json';d=producer(path,'publication_packaging')
    if d:add('publication_packaging','publication_archive',d.get('CPU_seconds'),path,path,'Recorded evidence archive process, separate from scientific inversions',d.get('status','RECORDED'),d.get('execution_platform'),'packaging_transport')
    literature='audit/literature/cie2019_error_formula_check.json'
    if (ROOT/literature).exists():
        add('external_retrieval','cie2019_error_formula_check',None,literature,literature,'Root primary-formula tool check; CPU/tool billing unavailable, not zero',platform='UNKNOWN',activity='literature')
    # Detailed scalar small-LU/RHS counts are dimensions, not extra process CPU
    # or original Maxwell factorizations. Read only cost/setup JSON artifacts.
    for root in [ROOT/'RESULTS/twod',ROOT/'RESULTS/right_projection_diagnostic',ROOT/'RESULTS/conditional_confirmation']:
        if not root.exists():continue
        for p in root.rglob('AUXILIARY_COST.json'):
            c=read(p)
            if 'denominator_small_LU' not in c:continue
            setup=p.parent/'DENOMINATOR_SETUP.json';s=read(setup) if setup.exists() else {}
            dim=s.get('actual_denominator_rank',s.get('twofold_rank') if s.get('basis')=='twofold' else s.get('green_strong_rank'))
            views.append(dict(source=rel(p),kind='auxiliary_denominator_small_LU_NOT_physical_Maxwell',matrix_dimension=dim,
                LU=c['denominator_small_LU'],LU_dimension_sum=c.get('denominator_small_LU_dimension_sum'),complex_RHS=c.get('denominator_small_triangular_complex_RHS'),
                projection=s.get('projection'),additive=False))
    generation=ROOT/'delegated/twod/GENERATION_RECEIPT.json'
    if generation.exists():views.append(dict(source=rel(generation),kind='primary24_generation_only_once',matrix_dimension=576,counts=read(generation).get('physical_cost'),additive=False))
    # Invocation log is bookkeeping CPU only, not another physics test.
    audit=ROOT/'audit';audit.mkdir(exist_ok=True);log=audit/'CAMPAIGN_COST_COMPILER.jsonl'
    old=[json.loads(x) for x in log.read_text().splitlines() if x] if log.exists() else []
    current=dict(utc=utc(),CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-started,new_physical_calls=0,GT_arrays_read=0)
    with log.open('a') as fh:fh.write(json.dumps(current)+'\n')
    add('campaign_bookkeeping','compiler_invocations',sum(x['CPU_seconds'] for x in old+[current]),rel(log),rel(log),'Disjoint compiler calls; timer excludes final output serialization',platform='Mac',activity='metadata')
    grouped=defaultdict(lambda:dict(CPU_seconds_known=0.,known_rows=0,unknown_rows=0))
    for x in rows:
        g=grouped[(x['branch'],x['platform'],x['activity'])]
        if x['CPU_seconds'] is None:g['unknown_rows']+=1
        else:g['CPU_seconds_known']+=x['CPU_seconds'];g['known_rows']+=1
    summary=[dict(branch=k[0],platform=k[1],activity=k[2],**v) for k,v in sorted(grouped.items())]
    result=dict(status='KNOWN_RECORDED_COST_SNAPSHOT',created_utc=utc(),rows=rows,group_summaries=summary,
        CPU_seconds_known_lower_bound=sum(x['CPU_seconds'] for x in rows if x['CPU_seconds'] is not None),
        producer_statuses=producers,operation_counter_views=views,unknowns=unknown,
        no_double_count='Canonical producer additive rows once; old DEV includes historical energy/scoring; vector FINAL supersedes CURRENT; production CIE prefix shared once; deployment/nested timers/counters are nonadditive',
        dimension_policy='Never pool scalar576, vector5184/17496, shifted auxiliary matrices or small19/30 denominators as one LU/RHS speed metric',
        wall_policy='No sum of overlapping process walls is a campaign makespan',new_physical_calls=0,GT_arrays_read=0,scorer_invocations=0,
        refresh='Update branch cost compilers after completion; rerun this script. Missing/current producers remain explicit lower bounds, not algorithm failures or zero costs.')
    (audit/'CAMPAIGN_COSTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    lines=['# R10 已记录费用账本','',f"当前已知CPU下界：**{result['CPU_seconds_known_lower_bound']:.6f} 秒**。这是成本记录汇总，不是科学结论，也不是包含全部未知开销的精确总额。",'',
        '|分支|平台|类别|已知CPU秒|未知行|','|---|---|---|---:|---:|']
    for x in summary:lines.append(f"|{x['branch']}|{x['platform']}|{x['activity']}|{x['CPU_seconds_known']:.6f} ({x['known_rows']}行)|{x['unknown_rows']}|")
    lines+=['','whole-child/wait4/GetProcessTimes与controller各计一次；R8/R10内部、auxiliary、corrector、部署视图均不再次收费。原24例数据生成仅在baseline计一次；conditional4为独立新数据，另计一次。两种条件部署共享CIE前缀的费用视图不能相加为生产总账。','',
        '打包/传输/归档在表中单列类别。DEV成本审计已包含其历史energy/评分及失败尝试，不能再加其raw receipt或energy总和。vector只读取ACCOUNTING_FINAL（若有）或CURRENT的additive行，compiler/transport已由该producer收费，不再加每臂receipt。','',
        '原Maxwell物理LU/RHS、移位辅助矩阵、低秩分母小LU/RHS保留在JSON独立计数视图，维度不合并；fresh Twofold实际rank30，不能沿用历史rank28标签。GPU事件时间不加CPU，且不是active-kernel时间。并行wall不相加为日历工期。','',
        f"显式unknown记录：{len(unknown)}条。缺失/进行中producer见CAMPAIGN_COSTS.json；无记录不填0。",'',
        '最终刷新：先更新vector与各分支cost ledger，再运行 `python research/research_max_struct_r10/experiments/compile_campaign_costs.py`。本脚本只读取已知cost/setup JSON，不加载GT/NPZ、不运行物理或评分。']
    (ROOT/'COST_LEDGER.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(status=result['status'],CPU_seconds_known_lower_bound=result['CPU_seconds_known_lower_bound'],additive_rows=len(rows),unknowns=len(unknown),GT_arrays_read=0,new_physical_calls=0)))


if __name__=='__main__':main()
