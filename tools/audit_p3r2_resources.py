"""Stream R2 resource evidence; no experiment, array, truth or model loading.

Main schedules under its D analysis wrapper:
  python tools/audit_p3r2_resources.py --accounted-job FINAL_RESOURCE_AUDIT
The wrapper must use that exact job name. Without a wrapper, active_wall must
be empty. Writes new AUDIT/resources/{RESULTS.json,RESOURCE_SUMMARY.md} only.
"""
from pathlib import Path
from collections import Counter, defaultdict
from datetime import datetime, timezone
import argparse
import fcntl
import json
import math
import statistics
import sys
import time
from scipy.stats import chi2

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'research/p3l_r2_risk_posterior_exploration'
LIMIT=8*1024*1024
FIELDS=('RHS','states','LU','controlled_wall_seconds','physics_wall_seconds','GPU_seconds')
COMPONENTS=('native','probit','final_context_J','matched_extra')
TRAINING=('training','training32768','residual_training','residual_controls',
          'residual_local_controls','density_training','weak_density_training',
          'mechanism_controls','random_sensitivity_controls')
BACKEND=('full_forward_calls','a30_fresh_full_state_calls','full_forward_source_rhs',
         'full_forward_RHS','linear_solve_rhs','LU_factorizations','full_LU_factorizations',
         'receiver_transfer_builds','receiver_transfer_adjoint_RHS','likelihood_gradient_RHS',
         'saved_background_file_reads','saved_background_imports','paid_cache_loads',
         'paid_background_array_setup_seconds','a30_background_file_read_seconds',
         'a30_saved_background_import_seconds','a30_cached_background_equation_audit_seconds')


def stats(values):
    values=[float(x) for x in values if isinstance(x,(float,int)) and math.isfinite(x)]
    if not values:return dict(status='MISSING',count=0)
    values.sort();position=(len(values)-1)*.95;lo=math.floor(position);hi=math.ceil(position)
    return dict(count=len(values),total=sum(values),median=statistics.median(values),
                P95=values[lo]+(position-lo)*(values[hi]-values[lo]),minimum=values[0],maximum=values[-1])


def close(a,b):
    if isinstance(a,(int,float)) and isinstance(b,(int,float)):
        return math.isfinite(a) and math.isfinite(b) and math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-7)
    return a==b


class Audit:
    def __init__(self,out):
        self.out=out;self.errors=[];self.error_count=0;self.missing=Counter();self.sources=set()

    def error(self,message):
        self.error_count+=1
        if len(self.errors)<100:self.errors.append(message)

    def read(self,relative):
        path=self.out/relative
        if path.is_symlink() or self.out not in path.resolve().parents:
            raise ValueError('External/symlink evidence refused: '+str(relative))
        if not path.is_file():self.missing[str(relative)]+=1;return None
        if path.stat().st_size>LIMIT:raise ValueError('Metadata exceeds bounded limit: '+str(relative))
        self.sources.add(str(relative));return json.loads(path.read_text())

    def ledger(self,snapshot):
        previous=None;sequence=0;jobs={};roles=defaultdict(Counter);events=defaultdict(Counter)
        backend=defaultdict(Counter);backend_coverage=defaultdict(Counter)
        max_active_nonanalysis=0;last_utc=None;gpu_absent=0
        path=self.out/'RESOURCE_LEDGER.jsonl'
        if path.is_symlink():raise ValueError('Symlink ledger refused')
        self.sources.add('RESOURCE_LEDGER.jsonl')
        with path.open() as stream:
            while True:
                line=stream.readline(LIMIT+1)
                if not line:break
                if len(line)>LIMIT:raise ValueError('Oversize ledger row; bounded reader stopped')
                if not line.endswith('\n'):self.error('Truncated last ledger row');break
                row=json.loads(line);n=row.get('sequence');after=row['after'];stage=row['stage'];job=row['job'];event=row['event']
                if n!=sequence+1:self.error(f'Sequence gap: expected {sequence+1}, observed {n}')
                sequence=n
                if after.get('sequence')!=n:self.error(f'Row/after sequence differs at {n}')
                if previous is None:
                    previous=dict(sequence=0,stages={s:{k:0 for k in FIELDS} for s in after['stages']},active_wall={})
                if set(previous['stages'])!=set(after['stages']):self.error(f'Stage set changed at {n}')
                if 'GPU_seconds' not in row:gpu_absent+=1
                for s,current in after['stages'].items():
                    for key in FIELDS:
                        if key not in current or key not in previous['stages'].get(s,{}):
                            self.error(f'Missing cumulative field {s}/{key} at {n}');continue
                        delta=current[key]-previous['stages'][s][key]
                        if s!=stage:
                            expected=0
                        elif key in row:
                            expected=row[key]
                        elif key=='GPU_seconds':
                            # Book currently omits gpu from row-level metadata;
                            # only a zero cumulative change is verifiable here.
                            expected=0
                            if delta:self.error(f'Unverifiable nonzero GPU delta at {n}')
                        else:
                            self.missing['ledger_row_field/'+key]+=1
                            self.error(f'Missing row increment {key} at {n}');continue
                        if not close(delta,expected):self.error(f'Increment mismatch seq{n} {s}/{key}: delta={delta}, row={expected}')
                before_active=previous.get('active_wall',{});active=after.get('active_wall',{})
                oldkeys=set(before_active);newkeys=set(active)
                if event=='JOB_START':
                    if job in jobs:self.error(f'Reused job identifier {job} at {n}')
                    if newkeys!=oldkeys|{job} or job in oldkeys:self.error(f'JOB_START active_wall mismatch at {n}')
                    jobs[job]=dict(stage=stage,start_sequence=n,analysis_only=row.get('analysis_only'),
                                   reserved_wall_seconds=active.get(job),reserved={k:0 for k in FIELDS})
                elif event=='JOB_END':
                    if job not in oldkeys or newkeys!=oldkeys-{job}:self.error(f'JOB_END active_wall mismatch at {n}')
                elif active!=before_active:self.error(f'Unexpected active_wall transition at {n}')
                for key in oldkeys&newkeys:
                    if active[key]!=before_active[key]:self.error(f'Existing active reservation changed at {n}: {key}')
                info=jobs.setdefault(job,dict(stage=stage,reserved={k:0 for k in FIELDS}))
                for key in FIELDS:
                    if key in row:info['reserved'][key]+=row[key]
                events[stage][event]+=1
                if row.get('RHS',0):
                    role='full_forward_attempt_reservation' if event=='FULL_FORWARD_ATTEMPT' else row.get('role','UNSPECIFIED_'+event)
                    roles[stage][role]+=row['RHS']
                if event=='JOB_END':
                    info.update(end_sequence=n,status=row.get('status'),error=row.get('error'),
                                unused_forward_credit=row.get('unused_forward_credit'),
                                counters=row.get('counts'),stage_cumulative_at_end=after['stages'][stage])
                    counts=row.get('counts')
                    if isinstance(counts,dict):
                        # JOB_END and CLOSED repeat the same counters; sum once.
                        for key,value in counts.items():
                            if isinstance(value,(float,int)):
                                backend[stage][key]+=value;backend_coverage[stage][key]+=1
                    else:self.missing['JOB_END/counts']+=1
                nonanalysis=sum(jobs.get(name,{}).get('analysis_only') is False for name in active)
                max_active_nonanalysis=max(max_active_nonanalysis,nonanalysis)
                previous=after;last_utc=row.get('utc_unix')
        if previous is None:raise ValueError('Empty ledger')
        if previous!=snapshot:self.error('Final streamed after does not exactly match current RESOURCE_LEDGER.json snapshot')
        stage_summary={}
        for stage,cumulative in previous['stages'].items():
            actual={k:backend[stage][k] for k in backend[stage]}
            stage_summary[stage]=dict(ledger_cumulative=cumulative,reserved_RHS_by_role=dict(roles[stage]),
                ledger_event_counts=dict(events[stage]),actual_backend_closed_job_counters=actual,
                actual_counter_job_coverage=dict(backend_coverage[stage]),
                missing_requested_backend_keys=[k for k in BACKEND if k not in actual],
                reservation_minus_actual=dict(
                    state_attempts_minus_reported_forward_calls=cumulative['states']-actual['full_forward_calls'] if 'full_forward_calls' in actual else None,
                    reserved_LU_minus_reported_factorizations=cumulative['LU']-actual['LU_factorizations'] if 'LU_factorizations' in actual else None,
                    reserved_RHS_minus_reported_solved_RHS=cumulative['RHS']-actual['linear_solve_rhs'] if 'linear_solve_rhs' in actual else None))
        return dict(last_sequence=sequence,last_ledger_utc_unix=last_utc,snapshot_matches=previous==snapshot,
                    stages=stage_summary,jobs=jobs,active_wall=previous['active_wall'],
                    row_GPU_increment_field_absent_count=gpu_absent,
                    GPU_row_audit_limit='Book rows omit GPU increment field; zero after-deltas checked; nonzero unrecorded deltas fail',
                    max_declared_nonanalysis_active_reservations=max_active_nonanalysis)

    def closures(self,ledger):
        failed=[];credits=[]
        for job,info in ledger['jobs'].items():
            if 'end_sequence' not in info:
                info['closure_scope']='ACTIVE_OR_NO_JOB_END; counters and wall incomplete';continue
            rel=f"{info['stage']}/jobs/{job}/CLOSED.json";closed=self.read(rel)
            if closed is None:continue
            for key,expected in (('status',info.get('status')),('error',info.get('error')),('counters',info.get('counters')),('stage_cumulative',info.get('stage_cumulative_at_end'))):
                if key not in closed:self.missing['CLOSED/'+key]+=1
                elif closed[key]!=expected:self.error(f'{rel}: {key} differs from ledger JOB_END')
            wall=closed.get('controlled_wall_seconds');expected=info['reserved']['controlled_wall_seconds']
            if wall is None:self.missing['CLOSED/controlled_wall_seconds']+=1
            elif not close(wall,expected):self.error(f'{rel}: controlled wall differs from JOB_END accumulation')
            if info.get('status')!='COMPLETE':failed.append(dict(job=job,stage=info['stage'],status=info.get('status'),error=info.get('error'),reserved=info['reserved']))
            if info.get('unused_forward_credit') is not None:
                credits.append(dict(job=job,stage=info['stage'],unused_forward_credit=info['unused_forward_credit']))
            else:self.missing['JOB_END/unused_forward_credit']+=1
        return dict(failed_jobs=failed,unused_forward_credits=credits,
                    interpretation='Reserved attempts and credits remain billed. Backend counters are reported completed operations, not reservations; gaps are retained, not silently reconciled.')

    def context_cohort(self,start,count,label,matched_extra):
        parents=[];component_values={c:defaultdict(list) for c in COMPONENTS}
        native=[];probit=[];phis=[];optimizer_status=Counter();threshold=float(chi2.ppf(.95,512)/2)
        for index in range(start,start+count):
            parent=f'R2_locked_{index:05d}';cost=self.read(f'L/contexts/costs/{parent}.json');station=self.read(f'L/contexts/{parent}.json')
            entry=dict(parent=parent,cohort=label,components={})
            if not matched_extra:
                entry['components']['matched_extra']=dict(status='NOT_RUN',reason='Second-cohort native96+probit64 design has no matched16 control')
            if cost:
                if cost.get('parent')!=parent:self.error('Cost parent mismatch: '+parent)
                if cost.get('GT_used') is not False:self.error('Cost receipt missing GT_used=false: '+parent)
                for component in COMPONENTS:
                    if component=='matched_extra' and not matched_extra:continue
                    row=cost.get('components',{}).get(component)
                    if row is None:self.missing['context_cost/'+label+'/'+component]+=1;continue
                    selected={k:row[k] for k in ('new_RHS','new_states','new_LU','wall_seconds') if k in row}
                    for key in ('new_RHS','new_states','new_LU','wall_seconds'):
                        if key in selected:component_values[component][key].append(selected[key])
                        else:self.missing['context_cost/'+label+'/'+component+'/'+key]+=1
                    counters=row.get('counters')
                    selected['actual_backend_counters']={k:counters[k] for k in BACKEND if k in counters} if isinstance(counters,dict) else None
                    entry['components'][component]=selected
            if station:
                if station.get('parent')!=parent:self.error('Stationarity parent mismatch: '+parent)
                if station.get('GT_read') is not False:self.error('Stationarity receipt missing GT_read=false: '+parent)
                ns=station.get('native_stationarity',{});ps=station.get('probit_stationarity',{})
                entry['native_stationarity']={k:ns[k] for k in ('status','scaled_KKT','raw_projected_gradient_inf') if k in ns}
                entry['probit_stationarity']={k:ps[k] for k in ('gradient_inf','scaled_gradient_inf','coordinate_mode','qualification','optimizer') if k in ps}
                if 'scaled_KKT' in ns:native.append(ns['scaled_KKT'])
                else:self.missing[label+'/native_stationarity/scaled_KKT']+=1
                if 'scaled_gradient_inf' in ps:probit.append(ps['scaled_gradient_inf'])
                else:self.missing[label+'/probit_stationarity/scaled_gradient_inf']+=1
                opt=ps.get('optimizer',{})
                optimizer_status[str(opt.get('status','MISSING'))+'/success='+str(opt.get('success','MISSING'))]+=1
                phi=station.get('Phi');entry['final_context_Phi']=phi
                if isinstance(phi,(int,float)) and math.isfinite(phi):phis.append(phi);entry['descriptive_noise_ball95']=phi<=threshold
                else:self.missing[label+'/stationarity/Phi']+=1
            parents.append(entry)
        return dict(cohort=label,expected_parents=count,parents=parents,
            component_distributions={c:(dict(status='NOT_RUN',reason='Not in second-cohort design') if c=='matched_extra' and not matched_extra
                else {k:stats(v) for k,v in fields.items()} if fields else dict(status='MISSING')) for c,fields in component_values.items()},
            native_scaled_KKT=stats(native),probit_scaled_gradient_inf=stats(probit),probit_optimizer_status_counts=dict(optimizer_status),
            final_context_Phi=stats(phis),noise_ball95=dict(threshold=threshold,definition='Phi <= chi2.ppf(.95,512)/2; 256 complex observations',
                observed_contexts=len(phis),within_count=sum(phi<=threshold for phi in phis),missing_contexts=count-len(phis)),
            stationarity_scope='Native projected KKT and probit-coordinate gradients are distinct numerical diagnostics. No claim that all points are MAP converged. Noise-ball membership is descriptive and not a fit or posterior certificate.',
            measurement_scope='Saved per-parent incremental measurements from one serial, paid-background/cached worker execution; no repeated cold/warm benchmarks. First-parent background setup can be included in native time. Final-J reuse may legitimately have zero additional RHS.')

    def contexts(self):
        cohorts={'first256':self.context_cohort(0,256,'first256',True)}
        # Presence only: no model/config access, and no extra status gate needed
        # for a resource snapshot. Active jobs are still rejected by main.
        if (self.out/'MECHANISM_CHECKPOINT.json').exists() or (self.out/'L/contexts/costs/R2_locked_00256.json').exists():
            cohorts['second256_confirmation']=self.context_cohort(256,256,'second256_confirmation',False)
        return dict(expected_total_parents=sum(c['expected_parents'] for c in cohorts.values()),cohorts=cohorts,
            cohorts_not_pooled=True,second_cohort_matched16='NOT_RUN_BY_DESIGN',
            scope='All registered first/confirmation parent costs and stationarity are reported separately; no cross-cohort MAP convergence assertion')

    def training(self):
        rows=[]
        for folder in TRAINING:
            rel=f'A/{folder}/RESULTS.json';result=self.read(rel)
            if result is None:continue
            for run in result.get('results',[]):
                rows.append(dict(family=folder,kind=run.get('kind'),seed=run.get('seed'),train_size=run.get('train_size'),
                    selected_update=run.get('selected_update'),cpu_seconds=run.get('cpu_seconds'),wall_seconds=run.get('wall_seconds'),
                    device=run.get('training_device',run.get('device')),source=rel))
                for key in ('cpu_seconds','wall_seconds'):
                    if key not in run:self.missing['training/'+folder+'/'+key]+=1
        return dict(runs=rows,CPU_seconds=stats([r['cpu_seconds'] for r in rows]),wall_seconds=stats([r['wall_seconds'] for r in rows]),
                    interpretation='Child per-run CPU/wall is a subset of accounted subprocess wall; do not add it again to ledger controlled wall. Validation selection, packing and orchestration overhead may be outside per-run timings but inside wrapper time.')


def markdown(result):
    ledger=result['ledger'];lines=['# R2 资源审计','',f"审计状态：**{result['status']}**；发现 {result['error_count']} 个账本/回执不一致。",
        '',f"快照范围：sequence {ledger['last_sequence']}，最后账本时间戳 {ledger['last_ledger_utc_unix']}；读取时活动预留 `{ledger['active_wall']}`。",
        '本次自身分析作业的结束时间及 JOB_END 尚未进入此快照；本报告不是最终封账。','',
        '| stage | reserved RHS | state attempts | reserved LU | controlled wall s | physics wall s | GPU s |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for stage,value in ledger['stages'].items():
        c=value['ledger_cumulative'];lines.append('| '+stage+' | '+' | '.join(str(c[k]) for k in FIELDS)+' |')
    lines.extend(['','RHS 角色、后端实际 forward/LU/receiver/gradient 计数及其覆盖作业数见 RESULTS.json。预算预留和后端完成操作分别列出，不把两者当成同一个计数；JOB_END 与 CLOSED 重复的 counters 只累计一次。',
        '',f"失败作业数：{len(result['closures']['failed_jobs'])}；未消耗 forward credit 回执逐作业保留。",'',
        '## 每父样本上下文实测成本',''])
    for cohort,context in result['contexts']['cohorts'].items():
        lines.extend(['### '+cohort,'',
            '| component | RHS total | RHS median | RHS P95 | wall total s | wall median s | wall P95 s |',
            '|---|---:|---:|---:|---:|---:|---:|'])
        for component,value in context['component_distributions'].items():
            r=value.get('new_RHS',{});w=value.get('wall_seconds',{})
            marker='NOT_RUN' if value.get('status')=='NOT_RUN' else 'MISSING'
            lines.append('| '+cohort+'/'+component+' | '+' | '.join(str(d.get(k,marker)) for d,k in ((r,'total'),(r,'median'),(r,'P95'),(w,'total'),(w,'median'),(w,'P95')))+' |')
        n=context['noise_ball95']
        lines.extend(['',f"{cohort} 最终 probit 上下文的描述性 noise-ball：{n['within_count']} / {n['observed_contexts']} 个有 Phi 回执的父样本；缺失 {n['missing_contexts']}。阈值 {n['threshold']}。",''])
    lines.extend(['第二批为 native96 + probit64；matched16 为 NOT_RUN（设计未执行），不计为缺失失败。',
        '两批成本与 stationarity 分开报告；native KKT 与 probit 梯度分别统计。优化器未成功、预算到限或缺失回执均保留，不能宣称所有点已达到 MAP 收敛。',
        '', '## 测量边界','',
        '已付背景数组导入、缓存读取、setup 与 CPU 训练/分析的受控墙钟由所属账本作业包含。已有付费缓存不是从零免费建立；每轮训练的 CPU/wall 只作分解，不重复加到总账。',
        f"GPU 累计：{result['GPU_total_seconds']} 秒；这是快照值，不是未结束作业的最终资源证明。GPU 本行增量缺字段的范围见 JSON。",
        '单物理 worker 与 receiver transfer cache≤4 是源码策略；本工具只审计活动预留，不能据此证明操作系统层面没有并发。',
        '成本是原执行中保存的增量测量，不包含新冷启动/热启动重复实验；首父样本可能承担背景 setup，后续复用缓存。未记录的 cold/warm 分拆标为不可用，不能由总量反推。',
        '只读 R2 ledger/CLOSED、上下文 costs/stationarity JSON 与训练 RESULTS；未打开 GT、观测数组、LU 或模型，未重复实验。',
        '',f"缺字段/缺文件类别：{len(result['missing_fields_or_files'])}；详见 RESULTS.json。"])
    if result['errors_first100']:lines.extend(['','## 不一致（最多前100条）','']+['- '+x for x in result['errors_first100']])
    return '\n'.join(lines)+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=OUT)
    parser.add_argument('--accounted-job',help='Exact own stage D analysis wrapper job name')
    args=parser.parse_args();out=args.root.resolve();destination=out/'AUDIT/resources'
    if destination.exists():raise FileExistsError('Preserve prior resource audit; output already exists')
    audit=Audit(out);start=time.monotonic()
    lockpath=out/'budget.lock'
    if lockpath.is_symlink():raise ValueError('Symlink budget lock refused')
    with lockpath.open('r+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        snapshot=audit.read('RESOURCE_LEDGER.json')
        if snapshot is None:raise ValueError('RESOURCE_LEDGER.json required')
        expected={args.accounted_job} if args.accounted_job else set()
        if set(snapshot.get('active_wall',{}))!=expected:raise ValueError('Only the explicitly named own audit job may be active')
        ledger=audit.ledger(snapshot)
        if args.accounted_job:
            own=ledger['jobs'].get(args.accounted_job,{})
            if own.get('stage')!='D' or own.get('analysis_only') is not True or 'end_sequence' in own:
                raise ValueError('Own accounted job must be active stage D analysis_only')
        closures=audit.closures(ledger);contexts=audit.contexts();training=audit.training()
        result=dict(status=('FAIL_INTERNAL_ACCOUNTING' if audit.error_count else
                    'CHECKED_WITH_MISSING_FIELDS' if audit.missing else 'PASS_INTERNAL_ACCOUNTING'),
            scientific_verdict=None,captured_utc=datetime.now(timezone.utc).isoformat(),ledger=ledger,closures=closures,
            contexts=contexts,training=training,GPU_total_seconds=sum(s['ledger_cumulative']['GPU_seconds'] for s in ledger['stages'].values()),
            missing_fields_or_files=dict(audit.missing),error_count=audit.error_count,errors_first100=audit.errors,sources=sorted(audit.sources),
            scope='Sequence-timestamped in-progress resource snapshot, not final closure; own wrapper reservation allowed, own JOB_END wall not yet included',
            cache_policy=dict(paid_workers=1,receiver_transfer_cache_entries=4,policy_source='src/p3r2/core.py; not re-executed'),
            no_GT_arrays_LU_checkpoints_physics_training_or_benchmarks=True,wall_seconds=time.monotonic()-start)
        destination.mkdir(parents=True,exist_ok=False)
        with (destination/'RESULTS.json').open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False)
        with (destination/'RESOURCE_SUMMARY.md').open('x') as stream:stream.write(markdown(result))
    print(json.dumps(dict(status=result['status'],errors=audit.error_count,output=str(destination))))
    if audit.error_count:sys.exit(2)


if __name__=='__main__':main()
