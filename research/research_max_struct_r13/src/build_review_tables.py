"""Compact authorized aggregate tables only; no pixels, solvers or new statistics gates."""
import argparse,csv,io,json,math,re,statistics,time
from collections import defaultdict
from pathlib import Path
FIELDS='campaign host domain physics_kind cohort family case_id object_id task_id method view available availability_reason run_status partial_failed_run eligible legal material_legal RMSE NRMSE boxrange_NRMSE relativeL2 IoU mean_amplitude_bias Phi D R residual_over_delta noise_feasible noise_certificate_eligible noise_qualification external_child_CPU_seconds CPU_seconds triggered bypass source_campaign source_method source_view parent_primary_NRMSE'.split()

def read(path):return json.loads(path.read_text())
def finite(v):return isinstance(v,(float,int)) and not isinstance(v,bool) and math.isfinite(v)
def write(out,name,rows,manifest):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row));buf=io.StringIO();writer=csv.DictWriter(buf,fieldnames=keys);writer.writeheader();chunks=[];n=0
    for row in rows:
        encoded={k:json.dumps(v,separators=(',',':')) if isinstance(v,(dict,list)) else v for k,v in row.items()}
        trial=io.StringIO();csv.DictWriter(trial,fieldnames=keys).writerow(encoded);line=trial.getvalue()
        if len((buf.getvalue()+line).encode())>1_800_000 and n:
            chunks.append((buf.getvalue(),n));buf=io.StringIO();csv.DictWriter(buf,fieldnames=keys).writeheader();n=0
        buf.write(line);n+=1
    chunks.append((buf.getvalue(),n))
    for i,(text,count) in enumerate(chunks):
        filename=name+(f'_part{i+1:02}' if len(chunks)>1 else '')+'.csv';f=out/filename;f.write_text(text);assert f.stat().st_size<2_000_000;manifest.append(dict(path=filename,rows=count,bytes=f.stat().st_size))

def main():
    cpu,wall=time.process_time(),time.perf_counter();p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--authorization',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();root=a.root.resolve();auth=read(a.authorization)
    if not auth.get('authorized_aggregation'):raise PermissionError('Parent authorized aggregation required')
    source=root/'analysis/combined';receipt=read(source/'AGGREGATION_RECEIPT.json')
    if not set(receipt['campaigns'])<=set(auth['campaigns']):raise PermissionError('Campaign not authorized')
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);manifest=[];rows=read(source/'ALL_ESTIMATOR_ROWS.json');deployment=read(source/'CONDITIONAL_DEPLOYMENT_ALL48.json');costs=read(source/'CONDITIONAL_DEPLOYMENT_COSTS.json');pairs=read(source/'PAIRED_OBJECT_CHANGES.json')
    compact=lambda rr:[{k:r.get(k) for k in FIELDS} for r in rr]
    for campaign in receipt['campaigns']:write(out,'CAMPAIGN_'+campaign,compact([r for r in rows if r['campaign']==campaign]),manifest)
    originals=[dict(r,method='original_g19_primary') for r in rows if r['campaign']=='confirmation_core48' and r['method']=='g19' and r['view']=='primary']
    write(out,'CONDITIONAL_ALL48_WITH_ORIGINAL_G19',compact(originals+deployment),manifest);write(out,'CONDITIONAL_PARENT_PLUS_RESCUE_COSTS',costs,manifest)
    fixed={f'r13_fresh_{i:04d}' for i in range(0,48,3)};nn=[r for r in rows if r['case_id'] in fixed and ((r['campaign']=='confirmation_core48' and r['method'] in ['dr','g19','direct']) or (r['campaign']=='confirmation_mechanism16' and r['method'] in ['nn_dr','nn_g19','nn_direct']))]
    write(out,'NN_COMMON16_ROWS',compact(nn),manifest)
    for name in ['ONLINE_COMMON16_ROWS','LOCAL_REFERENCE_COMMON16_ROWS']:write(out,name,compact(read(source/(name+'.json'))),manifest)
    for block in sorted({r['block'] for r in pairs}):write(out,'PAIRED_'+re.sub('[^A-Za-z0-9_-]','_',block),[r for r in pairs if r['block']==block],manifest)
    write(out,'PAIRED_OBJECT_BOOTSTRAP',read(source/'PAIRED_OBJECT_BOOTSTRAP.json'),manifest)
    grouped=defaultdict(list)
    for r in rows+deployment+originals:grouped[tuple(r.get(k) for k in ['campaign','host','domain','method','view'])].append(r)
    summary=[]
    for key,rr in grouped.items():
        item=dict(zip(['campaign','host','domain','method','view'],key));item.update(expected_rows=len(rr),available_rows=sum(bool(r.get('available')) for r in rr),missing_rows=sum(not r.get('available') for r in rr),ball_known_rows=sum(isinstance(r.get('noise_feasible'),bool) for r in rr),ball_feasible_rows=sum(r.get('noise_feasible') is True for r in rr))
        for metric in ['RMSE','NRMSE','boxrange_NRMSE','relativeL2','IoU','external_child_CPU_seconds']:
            vals=[r[metric] for r in rr if finite(r.get(metric))];item[metric+'_count']=len(vals);item[metric+'_mean']=statistics.mean(vals) if vals else None;item[metric+'_median']=statistics.median(vals) if vals else None
        summary.append(item)
    write(out,'REVIEW_MEANS_MEDIANS_BY_HOST_DOMAIN',summary,manifest)
    common_summary=[]
    common_sets={'NN_COMMON16_ROWS':nn,'ONLINE_COMMON16_ROWS':read(source/'ONLINE_COMMON16_ROWS.json'),'LOCAL_REFERENCE_COMMON16_ROWS':read(source/'LOCAL_REFERENCE_COMMON16_ROWS.json')}
    for name,source_rows in common_sets.items():
        groups=defaultdict(list)
        for row in source_rows:groups[tuple(row.get(k) for k in ['host','domain','method','view'])].append(row)
        for key,rr in groups.items():
            item=dict(comparison=name,host=key[0],domain=key[1],method=key[2],view=key[3],expected_rows=len(rr),available_rows=sum(bool(x.get('available')) for x in rr),ball_feasible_rows=sum(x.get('noise_feasible') is True for x in rr),ball_known_rows=sum(isinstance(x.get('noise_feasible'),bool) for x in rr))
            for metric in ['NRMSE','RMSE','IoU','external_child_CPU_seconds']:
                vals=[x[metric] for x in rr if finite(x.get(metric))];item[metric+'_count']=len(vals);item[metric+'_mean']=statistics.mean(vals) if vals else None;item[metric+'_median']=statistics.median(vals) if vals else None
            common_summary.append(item)
    write(out,'COMMON16_REVIEW_MEANS_MEDIANS',common_summary,manifest)
    (out/'REVIEW_TABLES_RECEIPT.json').write_text(json.dumps(dict(status='AUTHORIZED_TABLE_SERIALIZATION_COMPLETE',source_aggregation=str(source),host='Apple_local',source_estimator_rows=len(rows),conditional_rows=len(deployment),original_G19_primary_rows=len(originals),paired_rows=len(pairs),files=manifest,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,new_Maxwell_LU=0,new_RHS=0,truth_pixels_read=False,new_statistical_gate=False,scope='All source rows retained including missing; file chunks below2MB. Means/medians use recorded finite metrics and explicit counts; hosts/domains separate. Candidate CPU differs from external whole-run CPU; conditional latter includes parent+rescue. Repeated view/strategy CPU must not be added as research consumption.'),indent=2)+'\n')
    print(json.dumps(dict(files=len(manifest),source_rows=len(rows),paired_rows=len(pairs),max_file_bytes=max(x['bytes'] for x in manifest))))
if __name__=='__main__':main()
