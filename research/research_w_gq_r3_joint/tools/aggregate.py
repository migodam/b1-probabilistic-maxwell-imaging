"""Aggregate sealed/scored artifacts and exact accounting without opening private truth."""
from pathlib import Path
import json,csv,argparse,datetime,fcntl
ROOT=Path(__file__).resolve().parents[1]
FIELDS=['parent','noise','method','cohort','status','position_error','shape_IoU','epsilon_real_relative_error','loss_absolute_error','minor_semiaxis_abs_error','major_semiaxis_abs_error','minor_semiaxis_relative_error','major_semiaxis_relative_error','epsilon_complex_relative_error_f1','IoU_fine1601x1201','IoU_raster_change','center_distance_over_true_equivalent_radius','boundary_symmetric_Hausdorff_sampled','boundary_symmetric_mean_arclength','shape_mask_relative_L2','target_material_contrast_relative_L2','exact_projected_gradient','stationarity_scope','stop','converged','engineering_gate','full_objective','data_chi2','branches','full_fe','rhs','CPU_seconds','CPU_missing','wall_seconds','peak_RSS_bytes','polled_peak_RSS_bytes','pid','runtime_missing_status','cost_authority','output']
def read(path,default=None):return json.loads(path.read_text()) if path.exists() else default
def jsonl(path):return [json.loads(x) for x in path.read_text().splitlines() if x] if path.exists() else []
def write_json(path,data):path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
def write_csv(path,rows,fields):
    with path.open('w') as f:w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
def norm(path):
    p=Path(path);return str(p if p.is_absolute() else ROOT/p)
def fmt(x,digits=4):
    if x is None:return 'NOT_RECORDED'
    if isinstance(x,bool):return '是' if x else '否'
    if isinstance(x,int):return str(x)
    if isinstance(x,float):return f'{x:.{digits}g}'
    return str(x)
def table(rows,columns):
    return '| '+' | '.join(t for k,t in columns)+' |\n| '+' | '.join('---' for _ in columns)+' |\n'+'\n'.join('| '+' | '.join(('NOT_RUN_TIME' if r.get(k) is None and r.get('status')=='NOT_RUN_TIME' else fmt(r.get(k))) for k,t in columns)+' |' for r in rows)+'\n'
def continuous_metrics(score):
    z=score.get('estimated_parameters');t=score.get('truth_parameters')
    if z is None or t is None:return {}
    axes=sorted(z[14:16]);truth=sorted(t[14:16]);absolute=[abs(a-b) for a,b in zip(axes,truth)];relative=[a/b for a,b in zip(absolute,truth)]
    complex_eps_error=abs(complex(z[16]-t[16],z[17]-t[17]))/abs(complex(1+t[16],t[17]))
    return dict(minor_semiaxis_abs_error=absolute[0],major_semiaxis_abs_error=absolute[1],minor_semiaxis_relative_error=relative[0],major_semiaxis_relative_error=relative[1],epsilon_complex_relative_error_f1=complex_eps_error)
def main():
    a=argparse.ArgumentParser();a.add_argument('--out',type=Path,default=ROOT/'results/aggregation');args=a.parse_args();out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    with (ROOT/'.ledger.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_SH);ledger=jsonl(ROOT/'PHYSICS_LEDGER.jsonl');fcntl.flock(lock,fcntl.LOCK_UN)
    physics={}
    for r in ledger:
        if r.get('kind')=='reservation':
            p=physics.setdefault(r.get('pid'),dict(full_fe=0,rhs=0));p['full_fe']+=r.get('factorizations',0);p['rhs']+=r.get('rhs',0)
    wrappers=jsonl(ROOT/'JOB_LEDGER.jsonl');wrapper_by_output={}
    for r in wrappers:
        cmd=r.get('command',[])
        if not any(Path(str(t)).name in ['controls.py','deploy_r3.py','run_r3.py','ordinary_controls.py','oracle_control.py'] for t in cmd):continue
        for i,t in enumerate(cmd[:-1]):
            if t in ('--out','--output'):wrapper_by_output[norm(cmd[i+1])]=r
    rows=[];branches=[];nonprimary=[];unscored=[]
    seals=sorted((ROOT/'results').glob('**/ONLINE_SEALED.json'))
    for sealpath in seals:
        d=sealpath.parent;s=read(d/'SCORES.json');seal=read(sealpath,{});summary=seal.get('summary',read(d/'SUMMARY.json',{}))
        if s is None:unscored.append(dict(output=str(d.relative_to(ROOT)),method=summary.get('method'),status='SEALED_NOT_SCORED'));continue
        q=s['selected'];runtime=read(d/'runtime_cost.json',read(d/'RECEIPT.json',{}));wrap=wrapper_by_output.get(str(d));pid=runtime.get('pid')
        access=d/'FILE_ACCESS.jsonl'
        if access.exists():
            with access.open() as f:first=json.loads(next(f))
            pid=first.get('pid',pid)
        if pid is None and wrap:pid=wrap.get('pid')
        actual=physics.get(pid,runtime.get('job_physical_cost',{}));cpu=runtime.get('process_cpu_seconds');wall=runtime.get('wall_seconds');rss=runtime.get('peak_rss_bytes');authority='runtime_cost/RECEIPT plus exact per-PID reservations'
        if wrap:
            cpu=wrap.get('process_cpu_seconds');wall=wrap.get('wall_seconds');authority='run_suite inclusive child CPU once; exact per-PID reservations'
        if not actual and runtime.get('method') in ['PRIOR_ONLY','WALL_FIRST_FIXED','KNOWN_WALL_ORACLE']:actual=runtime.get('physical_cost',{})
        exact_pg=q.get('final_exact_projected_gradient');stop=q.get('stop');selected_rec=next((r for r in summary.get('starts',[]) if r.get('stage')==summary.get('selected_start')),{})
        exact_pg=selected_rec.get('final_exact_projected_gradient') if exact_pg is None else exact_pg;stop=selected_rec.get('stop') if stop is None else stop
        gate='PASS' if q['position_error']<=.1 and q['shape_IoU']>=.7 and q['epsilon_real_relative_error']<=.2 else 'FAIL'
        row=dict(parent=s['case'],noise=s['repeat'],method=s.get('method',summary.get('method',d.name)),cohort=s.get('cohort','DEVELOPMENT_EXPOSED' if s['case']=='dev' else 'INDEPENDENT_FROZEN_BLIND'),status=s['status'],position_error=q['position_error'],shape_IoU=q['shape_IoU'],epsilon_real_relative_error=q['epsilon_real_relative_error'],loss_absolute_error=q['loss_absolute_error'],shape_mask_relative_L2=q.get('shape_mask_relative_L2'),target_material_contrast_relative_L2=q['target_material_contrast_relative_L2'],exact_projected_gradient=exact_pg,stop=stop,converged=q.get('converged',summary.get('converged')),engineering_gate=gate,full_objective=summary.get('objective',q.get('objective')),data_chi2=summary.get('data_chi2'),branches=len(s['branches']),full_fe=actual.get('full_fe'),rhs=actual.get('rhs'),CPU_seconds=cpu,CPU_missing=cpu is None,wall_seconds=wall,peak_RSS_bytes=rss,polled_peak_RSS_bytes=(wrap.get('polled_peak_RSS_bytes') if wrap else runtime.get('polled_peak_RSS_bytes')),pid=pid,runtime_missing_status=runtime.get('missing_metrics_status','NOT_RECORDED' if cpu is None else 'RECORDED'),cost_authority=authority,output=str(d.relative_to(ROOT)))
        row['stationarity_scope']='seven target coordinates conditional on fixed wall' if row['method'] in ['WALL_FIRST_FIXED','KNOWN_WALL_ORACLE'] else 'not optimized' if row['method']=='PRIOR_ONLY' else 'all 19 joint coordinates'
        row.update(continuous_metrics(q))
        bad=any(any(k in part for k in ['INVALID','REPAIR','RANK8','UNDERRESOLVED','SUPERSEDED']) for part in d.parts)
        (nonprimary if bad else rows).append(row)
        for i,b in enumerate(s['branches']):
            br=dict(parent=s['case'],noise=s['repeat'],method=row['method'],nonprimary=bad,branch=i,endpoint=b['endpoint'])
            for k in ['objective','stop','evaluations','accepted_steps','position_error','shape_IoU','epsilon_real_relative_error','loss_absolute_error','shape_mask_relative_L2','target_material_contrast_relative_L2','final_exact_projected_gradient']:br[k]=b.get(k)
            br.update(continuous_metrics(b));branches.append(br)
    expected=[(c,n,m) for c in ['blind_normal','blind_cavity','blind_strong'] for n in [0,1] for m in ['FULL_GN','FIXED_CURRENT_CHART','ADAPTIVE_CURRENT_CHART']]+[(c,0,m) for c in ['blind_normal','blind_cavity','blind_strong'] for m in ['PRIOR_ONLY','WALL_FIRST_FIXED','KNOWN_WALL_ORACLE']]
    present={(r['parent'],r['noise'],r['method']) for r in rows}
    for parent,noise,method in expected:
        if (parent,noise,method) in present:continue
        d=ROOT/'results/online'/f'{parent}_r{noise}'/method;runtime=read(d/'runtime_cost.json',read(d/'RECEIPT.json',{}));wrap=wrapper_by_output.get(str(d));access=d/'FILE_ACCESS.jsonl';pid=runtime.get('pid')
        if access.exists():
            with access.open() as f:first=json.loads(next(f))
            pid=first.get('pid',pid)
        attempted=bool(pid or (d/'TRACE.jsonl').exists() or (d/'endpoint.npz').exists());sealed=(d/'ONLINE_SEALED.json').exists();status='SEALED_NOT_SCORED' if sealed else 'ATTEMPTED_NOT_SEALED' if attempted else 'NOT_RUN_TIME';cost=physics.get(pid,runtime.get('job_physical_cost',{}));cpu=wrap.get('process_cpu_seconds') if wrap else runtime.get('process_cpu_seconds')
        row={k:None for k in FIELDS};row.update(parent=parent,noise=noise,method=method,cohort='ORACLE_DIAGNOSTIC' if method=='KNOWN_WALL_ORACLE' else 'FROZEN_BLIND',status=status,engineering_gate=status,branches=0,full_fe=cost.get('full_fe',0 if not attempted else None),rhs=cost.get('rhs',0 if not attempted else None),CPU_seconds=cpu,CPU_missing=attempted and cpu is None,wall_seconds=wrap.get('wall_seconds') if wrap else runtime.get('wall_seconds'),peak_RSS_bytes=runtime.get('peak_rss_bytes'),polled_peak_RSS_bytes=wrap.get('polled_peak_RSS_bytes') if wrap else runtime.get('polled_peak_RSS_bytes'),pid=pid,runtime_missing_status='NOT_APPLICABLE_NOT_RUN' if not attempted else runtime.get('missing_metrics_status','NOT_RECORDED' if cpu is None else 'RECORDED'),cost_authority='Not executed' if not attempted else 'Attempt receipts and per-PID reservations; no scored endpoint',output=str(d.relative_to(ROOT)));rows.append(row)
    geometry=read(ROOT/'results/geometry_resolution/GEOMETRY_RESOLUTION.json',{})
    geom_by_parent={r['parent']:r for r in geometry.get('rows',[])}
    for r in rows:
        if r['method']=='FULL_GN' and r['noise']==0 and r['parent'] in geom_by_parent:
            g=geom_by_parent[r['parent']];r.update(IoU_fine1601x1201=g['IoU_1601x1201'],IoU_raster_change=g['IoU_raster_change'],center_distance_over_true_equivalent_radius=g['center_distance_over_true_equivalent_radius'],boundary_symmetric_Hausdorff_sampled=g['boundary_symmetric_Hausdorff_sampled'],boundary_symmetric_mean_arclength=g['boundary_symmetric_mean_arclength'])
    rows.sort(key=lambda r:(r['parent'],r['noise'],r['method']));blind=[r for r in rows if r['parent'].startswith('blind')];scored_blind=[r for r in blind if r['status']=='POST_SEAL_SCORED'];parents=sorted(set(r['parent'] for r in scored_blind));repeats={c:sorted(set(r['noise'] for r in scored_blind if r['parent']==c)) for c in parents}
    total=dict(full_fe=sum(x['full_fe'] for x in physics.values()),rhs=sum(x['rhs'] for x in physics.values()),authority='All R3 reservations, including invalid/development/rejected/interrupted/failed attempts; not just scored endpoints')
    cpu_known=sum(r['CPU_seconds'] for r in blind if r['CPU_seconds'] is not None);missing=sum(r['CPU_missing'] for r in blind)
    primary=['FULL_GN','FIXED_CURRENT_CHART','ADAPTIVE_CURRENT_CHART']
    paired_parents=[c for c in parents if all(any(r['parent']==c and r['noise']==0 and r['method']==m for r in scored_blind) for m in primary)]
    paired=[]
    for method in primary:
        rr=[next(r for r in scored_blind if r['parent']==c and r['noise']==0 and r['method']==method) for c in paired_parents]
        if not rr:continue
        n=len(rr);rec=dict(method=method,n_independent_parents=n,noise=0,position_RMSE=(sum(r['position_error']**2 for r in rr)/n)**.5,mean_IoU=sum(r['shape_IoU'] for r in rr)/n,mean_epsilon_relative_error=sum(r['epsilon_real_relative_error'] for r in rr)/n,mean_loss_absolute_error=sum(r['loss_absolute_error'] for r in rr)/n,mean_epsilon_complex_relative_error_f1=sum(r['epsilon_complex_relative_error_f1'] for r in rr)/n,mean_minor_semiaxis_relative_error=sum(r['minor_semiaxis_relative_error'] for r in rr)/n,mean_major_semiaxis_relative_error=sum(r['major_semiaxis_relative_error'] for r in rr)/n,mean_material_map_relative_L2=sum(r['target_material_contrast_relative_L2'] for r in rr)/n,engineering_passes=sum(r['engineering_gate']=='PASS' for r in rr),full_fe=sum(r['full_fe'] or 0 for r in rr),rhs=sum(r['rhs'] or 0 for r in rr),known_CPU_subtotal=sum(r['CPU_seconds'] for r in rr if r['CPU_seconds'] is not None),missing_CPU_rows=sum(r['CPU_missing'] for r in rr),total_CPU_seconds=None if any(r['CPU_missing'] for r in rr) else sum(r['CPU_seconds'] for r in rr),wall_seconds=sum(r['wall_seconds'] or 0 for r in rr))
        paired.append(rec)
    result=dict(updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),registered_primary_slots=18,registered_control_slots=9,n_completed_primary=sum(r['status']=='POST_SEAL_SCORED' and r['method'] in primary for r in blind),n_completed_controls=sum(r['status']=='POST_SEAL_SCORED' and r['method'] not in primary for r in blind),n_independent_blind_parents=len(parents),noise_repeats_by_parent=repeats,statistical_scope='Three frozen target/wall parent strata; repeated noises within parent are dependent descriptive repeats, not six independent parents. No population posterior/calibration claim.',engineering_gate='position<=.1 AND shapeIoU>=.7 AND real_epsilon_relative_error<=.2; loss and mapL2 reported separately',boundary_air_defect='Historical blind_cavity is one boundary-connected air macrocell, not an enclosed internal cavity.',geometry_resolution_review=geometry,geometry_review_scope='Only r0 selected FullGN endpoints; no new threshold or branch selection; original gates retained',paired_primary_r0_parentset=paired_parents,paired_primary_r0_summary=paired,rows=rows,nonprimary_development_rows=nonprimary,all_branches=branches,sealed_not_scored=unscored,global_physics_reservations=total,known_primary_or_control_CPU_seconds=cpu_known,n_missing_CPU_rows=missing,CPU_total_scope='Known CPU is a subtotal only when any row missing; never substitute wall time, ps or kernel timing for missing CPU.')
    write_json(out/'AGGREGATE.json',result)
    write_csv(out/'METHOD_SUMMARIES_PAIRED_R0.csv',paired,list(paired[0]) if paired else ['method','n_independent_parents'])
    write_csv(out/'IMAGING_RESULTS.csv',rows,FIELDS);write_csv(out/'NONPRIMARY_DEVELOPMENT.csv',nonprimary,FIELDS);write_csv(out/'BRANCH_RESULTS.csv',branches,list(branches[0]) if branches else ['parent','noise','method','branch']);write_csv(out/'SEALED_NOT_SCORED.csv',unscored,['output','method','status'])
    cols=[('parent','Parent / 历史ID'),('noise','噪声'),('method','方法'),('status','运行状态'),('position_error','位置误差'),('shape_IoU','IoU'),('epsilon_real_relative_error','ε实部相对误差'),('loss_absolute_error','损耗绝对误差'),('epsilon_complex_relative_error_f1','复ε相对误差f1'),('minor_semiaxis_relative_error','短半轴相对误差'),('major_semiaxis_relative_error','长半轴相对误差'),('target_material_contrast_relative_L2','材料图相对L2'),('engineering_gate','三门合取')]
    costcols=[('parent','Parent'),('noise','噪声'),('method','方法'),('exact_projected_gradient','完整J投影梯度'),('stationarity_scope','活动坐标范围'),('stop','停止原因'),('converged','收敛'),('full_objective','完整目标'),('data_chi2','数据χ²'),('full_fe','FE'),('rhs','RHS'),('CPU_seconds','CPU秒'),('CPU_missing','CPU缺测')]
    report='# 独立评估数据表（自动生成，科学解释由主线程负责）\n\n'+f"更新：{result['updated_utc']}。当前有 {len(parents)} 个独立盲测 parent；各 parent 噪声编号为 {repeats}。不能把同一 parent 的两次噪声当成两个独立世界。\n\n"+'主工程门严格为：位置≤0.1、IoU≥0.7、目标实部介电常数相对误差≤0.2，三者同时满足。损耗、f=1复介电常数相对误差、排序后两个连续半轴的绝对/相对误差与图像L2单独报告；排序只移除椭圆参数的轴交换对称，不使用像素分类标签。所有端点先按完整测量目标选择、封存，再打开私有真值评分。\n\n历史ID `blind_cavity` 是连边界空气缺口，不是封闭内腔。非线性前向 F 在各臂均为完整 Maxwell；forcing-current tangent chart 只近似材料切向 J。\n\n'+table(blind,cols)+'\n完整目标、收敛与费用：\n\n'+table(blind,costcols)+'\n只有三主方法都完成的同一 r0 parent 才进入以下配对描述汇总；当前 parent 集为 '+str(paired_parents)+'。不作总体置信推断。\n\n'+table(paired,[('method','方法'),('n_independent_parents','独立parent数'),('position_RMSE','位置RMSE'),('mean_IoU','平均IoU'),('mean_epsilon_relative_error','平均ε误差'),('engineering_passes','合取通过数'),('full_fe','FE'),('rhs','RHS'),('total_CPU_seconds','完整CPU秒'),('missing_CPU_rows','CPU缺测行')])+'\n连续轮廓/分辨率复核（0 PDE，不改冻结门）：\n\n'+table([r for r in blind if r['method']=='FULL_GN' and r['noise']==0 and r['status']=='POST_SEAL_SCORED'],[('parent','Parent'),('shape_IoU','原评分IoU'),('IoU_fine1601x1201','细栅格IoU'),('center_distance_over_true_equivalent_radius','中心距/真值等效半径'),('boundary_symmetric_mean_arclength','平均边界距λ₀'),('boundary_symmetric_Hausdorff_sampled','最大边界距λ₀')])+'\n轮廓放大图使用同一0.5λ₀物理视窗和等比例。细化栅格只复核显示/采样，不提高测量信息；FEM独立前向核验另列，不能从高分辨率图推出统计分辨率或连续PDE证书。\n\n开发案例（不计为盲测）：\n\n'+table([r for r in rows if r['parent']=='dev'],cols)+'\n保留的非主开发尝试：\n\n'+table(nonprimary,cols)+f"\n全部 R3 已预约费用（含失败和拒绝，不只是表中端点）：{total['full_fe']} FE / {total['rhs']} RHS。主表及控制已记录 CPU 小计 {cpu_known:.6f} 秒；{missing} 行 CPU 缺测，缺失值保持 NOT_RECORDED，不能将小计称为完整 CPU。\n\n四起点的逐分支指标见 BRANCH_RESULTS.csv；所有已封存但未评分项目见 SEALED_NOT_SCORED.csv。未运行或未完成部分不能以现有成功端点替代。没有全局最优、任意未知墙、连续 PDE 或全局 posterior 主张。\n"
    (out/'EVALUATION_TABLES_ZH.md').write_text(report)
    print(json.dumps(dict(rows=len(rows),blind_rows=len(blind),nonprimary=len(nonprimary),branch_rows=len(branches),sealed_not_scored=len(unscored),CPU_missing=missing,global_physics=total)))
if __name__=='__main__':main()
