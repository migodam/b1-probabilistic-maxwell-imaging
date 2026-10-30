"""Post-score cost and paired reports. No truth access, solver import or physical solves."""
import argparse, collections, csv, json, pathlib, statistics

COMPARATORS=('native_trf','reciprocal_trf','data_lbfgs')
ROOT=pathlib.Path(__file__).resolve().parent

def read(path):return json.loads(pathlib.Path(path).read_text(encoding='utf-8-sig'))
def mean(values):
    values=[x for x in values if isinstance(x,(int,float))]
    return statistics.fmean(values) if values else None
def write_csv(path,rows):
    columns=list(dict.fromkeys(k for row in rows for k in row))
    with pathlib.Path(path).open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader()
        writer.writerows({k:json.dumps(v,sort_keys=True) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in rows)
def outcome(left,right):
    tolerance=max(.02,.10*min(left,right));delta=left-right
    return ('tie' if abs(delta)<=tolerance else 'win' if delta<0 else 'loss'),delta,tolerance

def ledger_row(receipt,score_row,source):
    result=receipt.get('result') or {};top=result.get('top1') or {};counts=result.get('costcounts')
    accounting=result.get('accounting') or {};complete=counts is not None
    counts=counts or {}
    row={'category':'inference','task_id':score_row['task_id'],'backend':score_row['backend'],
         'method':score_row['method'],'observation_id':score_row['observation_id'],
         'status':receipt.get('status',score_row.get('transport_status','NOT_RUN')),
         'solver_status':result.get('status',score_row.get('solver_status')),
         'costs_status':receipt.get('costs_status','COMPLETE' if complete else 'UNKNOWN'),
         'source':str(source) if source else None,
         'wall_seconds':result.get('wall_seconds',score_row.get('wall_seconds')),
         'cpu_seconds':result.get('cpu_seconds',score_row.get('cpu_seconds')),
         'runner_solve_wall_seconds':receipt.get('solve_wall_s'),
         'worker_total_seconds':receipt.get('total_worker_s'),
         'top1_stage':top.get('stage'),'top1_origin':top.get('origin'),
         'all_unknowns_released':result.get('all_unknowns_released'),
         'stop_reason':result.get('stop_reason',score_row.get('stop_reason')),
         'LU':counts.get('LU',0) if complete else score_row.get('LU'),
         'actual_complex_RHS':counts.get('actual_complex_RHS',0) if complete else score_row.get('actual_complex_RHS'),
         'source_RHS':counts.get('source_RHS',0) if complete else None,
         'adjoint_RHS':counts.get('adjoint_RHS',0) if complete else None,
         'receiver_transpose_RHS':counts.get('receiver_transpose_RHS',0) if complete else None,
         'Gram_RHS':counts.get('Gram_RHS',0) if complete else None,
         'FullJ':counts.get('FullJ',0) if complete else None,
         'Gram':counts.get('Gram',0) if complete else None,
         'Gram_factorizations':counts.get('Gram_factorizations',0) if complete else None,
         'failed_operations':counts.get('failed_operations',0) if complete else None,
         'budget_overshoot_seconds':result.get('budget_overshoot_seconds'),
         'full_costcounts':counts if complete else None,
         'wall_by_operation':accounting.get('wall_by_operation'),
         'cpu_by_operation':accounting.get('cpu_by_operation')}
    for key in ('RMSE','NRMSE','relative_data_residual','data_gate','material_gate','joint_gate'):
        row[key]=score_row.get(key)
    if row['status']=='NOT_RUN':row['costs_status']='NOT_RUN'
    # Gram_RHS are small 12x12 Gram-system solves; never add them to full 576 RHS.
    return row

def extra_costs(audit,cost_notes):
    rows=[];generation=audit/'GENERATION_RECEIPT.json'
    if generation.exists():
        data=read(generation);counts=data.get('generation_costcounts',{})
        rows.append({'category':'generation','status':data.get('status'),'source':str(generation),
            'wall_seconds':data.get('generation_wall_seconds'),'LU':counts.get('LU'),
            'actual_complex_RHS':counts.get('actual_complex_RHS'),'full_costcounts':counts,
            'costs_status':'COMPLETE_PHYSICAL_COUNTERS','note':'36 full clean states; noise/Born diagnostic adds no LU'})
    failure=audit/'DEV_GPU_IMPORT_FAILURE.json'
    if failure.exists():
        data=read(failure);rows.append({'category':'import','backend':'cuda','status':data.get('status'),
            'source':str(failure),'LU':data.get('physical_LU'),'actual_complex_RHS':data.get('physical_RHS'),
            'costs_status':'PHYSICAL_KNOWN_TIME_UNKNOWN','note':data.get('reason')})
    correction=audit/'CPU_LAUNCH_CORRECTION.json'
    if correction.exists():
        data=read(correction);rows.append({'category':'launch','backend':'cpu','status':'FAILED_INITIAL_DETACHED_LAUNCH',
            'source':str(correction),'costs_status':'UNKNOWN','note':data.get('unrecovered_physical_cost'),
            'completed_predictions':data.get('initial_receipts')})
    for path in (audit/'DEV_CPU.json',audit.parent/'hardware/dev_gpu_actual/dev_cuda_globaldual.json'):
        if not path.exists():continue
        data=read(path);receipt=data if 'result' in data else {'result':data,'status':'DEV_RESULT_RETURNED'}
        result=receipt.get('result') or {}
        scalar={'task_id':receipt.get('task_id',path.stem),'backend':result.get('backend'),
                'method':result.get('method'),'observation_id':'OLD_PUBLIC_Y_DEV'}
        row=ledger_row(receipt,scalar,path);row['category']='DEV'
        row['note']='Old exposed public-Y development; cost outside production comparisons, not another independent object'
        rows.append(row)
    for path in cost_notes:
        data=read(path);rows.extend(data if isinstance(data,list) else data['rows'])
    # Explicit unknown rows preserve unmeasured overhead. Unknown never means free.
    for category in ('setup','DEV','import','launch'):
        if not any(r.get('category')==category for r in rows):
            rows.append({'category':category,'status':'UNMEASURED','costs_status':'UNKNOWN',
                         'note':'Supply actual scalar overhead/cost receipt using --cost-notes; no estimate invented'})
    return rows

def paired_rows(rows):
    lookup={(r['backend'],r['observation_id'],r['method']):r for r in rows}
    observations={(r['backend'],r['observation_id']):r for r in rows};out=[]
    for (backend,ident),case in observations.items():
        left=lookup.get((backend,ident,'globaldual'))
        for comparator in COMPARATORS:
            right=lookup.get((backend,ident,comparator));finite=bool(left and right and left.get('finite_endpoint') and right.get('finite_endpoint'))
            pair={k:case.get(k) for k in ('shape','kind','amplitude','noise')}
            pair.update(backend=backend,observation_id=ident,comparator=comparator,finite_pair=finite,
                        outcome='MISSING',globaldual_status=left.get('transport_status') if left else 'NOT_RUN',
                        comparator_status=right.get('transport_status') if right else 'NOT_RUN')
            for name,row in (('globaldual',left),('comparator',right)):
                for metric in ('RMSE','NRMSE','joint_gate','data_gate','material_gate'):
                    pair[name+'_'+metric]=row.get(metric) if row else None
            if finite:
                status,delta,tolerance=outcome(left['NRMSE'],right['NRMSE'])
                pair.update(outcome=status,delta_NRMSE=delta,delta_RMSE=left['RMSE']-right['RMSE'],practical_tolerance=tolerance)
            out.append(pair)
    return out

def grouped_pairs(pairs,fields):
    groups=collections.defaultdict(list)
    for pair in pairs:groups[tuple(pair.get(f) for f in fields)].append(pair)
    out=[]
    for key,group in groups.items():
        row=dict(zip(fields,key));row['planned_pairs']=len(group)
        for status in ('win','tie','loss','MISSING'):row[status]=sum(p['outcome']==status for p in group)
        row['distinct_shapes']=len({p['shape'] for p in group})
        for key2 in ('globaldual_NRMSE','comparator_NRMSE','globaldual_RMSE','comparator_RMSE','delta_NRMSE','delta_RMSE'):
            # Endpoint means use only matched finite pairs.
            row[key2+'_mean']=mean([p.get(key2) for p in group if p['finite_pair']])
        for method in ('globaldual','comparator'):
            for state in ('PASS','FAIL','NOT_RUN'):
                row[method+'_joint_'+state]=sum(p.get(method+'_joint_gate')==state for p in group)
        out.append(row)
    return out

def shape_clusters(pairs):
    groups=collections.defaultdict(list)
    for pair in pairs:groups[(pair['backend'],pair['comparator'],pair['shape'])].append(pair)
    out=[]
    for (backend,comparator,shape),group in groups.items():
        finite=[p for p in group if p['finite_pair']]
        row={'backend':backend,'comparator':comparator,'shape':shape,'planned_variants':len(group),
             'finite_paired_variants':len(finite),'missing_variants':len(group)-len(finite),
             'cluster_complete':len(group)==6 and len(finite)==6,'outcome':'MISSING'}
        for metric in ('NRMSE','RMSE'):
            for method in ('globaldual','comparator'):
                row[method+'_'+metric]=mean([p[method+'_'+metric] for p in finite])
            row['delta_'+metric]=mean([p['delta_'+metric] for p in finite])
        for method in ('globaldual','comparator'):
            row[method+'_joint_PASS_fraction_matched']=mean([int(p[method+'_joint_gate']=='PASS') for p in finite])
        if finite:
            status,_,tol=outcome(row['globaldual_NRMSE'],row['comparator_NRMSE'])
            row['outcome']=status;row['practical_tolerance_on_shape_means']=tol
        out.append(row)
    return out

def analyze(summary_path,receipts_dir,out,audit,cost_notes):
    summary=read(summary_path);rows=summary['rows'];out=pathlib.Path(out);out.mkdir(parents=True,exist_ok=True)
    paths={} if receipts_dir else {r['task_id']:r.get('path') for r in read(summary['seal'])['receipts']}
    ledger=[]
    for row in rows:
        path=pathlib.Path(receipts_dir)/(row['task_id']+'.json') if receipts_dir else paths.get(row['task_id'])
        exists=path and pathlib.Path(path).exists()
        receipt=read(path) if exists else {'status':row.get('transport_status','NOT_RUN')}
        ledger.append(ledger_row(receipt,row,path if exists else None))
    ledger.extend(extra_costs(pathlib.Path(audit),cost_notes));write_csv(out/'COST_LEDGER.csv',ledger)
    pairs=paired_rows(rows);write_csv(out/'PAIRED_VARIANTS.csv',pairs)
    write_csv(out/'PAIRED_BY_BACKEND_AMPLITUDE_NOISE.csv',grouped_pairs(pairs,('backend','amplitude','noise','comparator')))
    clusters=shape_clusters(pairs);write_csv(out/'SHAPE_CLUSTERS.csv',clusters)
    overview=[]
    for backend in ('cpu','cuda'):
        for comparator in COMPARATORS:
            group=[r for r in clusters if r['backend']==backend and r['comparator']==comparator]
            complete=[r for r in group if r['cluster_complete']]
            row={'backend':backend,'comparator':comparator,'planned_shapes':len(group),'complete_shapes':len(complete),
                 'partial_shapes':sum(not r['cluster_complete'] and r['finite_paired_variants']>0 for r in group)}
            for metric in ('NRMSE','RMSE'):row['equal_shape_mean_delta_'+metric]=mean([r['delta_'+metric] for r in complete])
            for status in ('win','tie','loss','MISSING'):row[status+'_complete_shapes']=sum(r['outcome']==status for r in complete)
            overview.append(row)
    write_csv(out/'SHAPE_CLUSTER_SUMMARY.csv',overview)
    cost_totals=[]
    for category in sorted({r['category'] for r in ledger}):
        group=[r for r in ledger if r['category']==category]
        for backend in sorted({r.get('backend') or 'unspecified' for r in group}):
            sub=[r for r in group if (r.get('backend') or 'unspecified')==backend]
            item={'category':category,'backend':backend,'rows':len(sub)}
            for key in ('LU','actual_complex_RHS','Gram_RHS','wall_seconds','cpu_seconds'):
                values=[r[key] for r in sub if isinstance(r.get(key),(float,int))]
                item[key+'_known_sum']=sum(values) if values else None;item[key+'_unknown_rows']=len(sub)-len(values)
            cost_totals.append(item)
    report={'scope':'Descriptive paired reporting only; no significance tests or superiority guarantees',
        'tie_rule':'abs(deltaNRMSE)<=max(.02,.10*min(leftNRMSE,rightNRMSE))',
        'independent_objects':'12 CPU shapes; CUDA is a paired 3-shape subset, never extra independent objects',
        'variants':'6 paired contrast/noise variants per shape; shape means use matched variants; overview requires complete shapes',
        'counter_scope':'Gram_RHS counts small Gram solves, separate from full-system actual_complex_RHS',
        'cost_time_scope':'Solver, worker-total, operation timing overlap; do not sum these time columns together',
        'no_new_physics':True,'cost_totals':cost_totals,'shape_cluster_summary':overview}
    (out/'ANALYSIS.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--summary',required=True)
    ap.add_argument('--receipts',help='Optional local directory overriding seal receipt paths')
    ap.add_argument('--out',default='analysis');ap.add_argument('--audit',default=ROOT/'audit')
    ap.add_argument('--cost-notes',action='append',default=[],help='JSON array or {rows:[...]} of scalar setup/DEV/import/launch costs')
    a=ap.parse_args();analyze(a.summary,a.receipts,a.out,a.audit,a.cost_notes)
