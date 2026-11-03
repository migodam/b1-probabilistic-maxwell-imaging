"""Read-only inherited-text/array/history coverage. No physical solver import.

R1 per-field inventory is reused only after direct byte equality to its already
audited copy, not a hash. R2 closed numerical arrays are read anew. Long JSON
histories are parsed and summarized without dumping them into chat context.
"""
from pathlib import Path
import ast,json,time
from collections import Counter
import numpy as np
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
PREVIOUS=REPO.parent/'r2_worktree'

def read():
    start,cpu=time.perf_counter(),time.process_time()
    previous_inventory=json.loads((REPO/'research/research_max_struct_r2/audit/R1_ARRAY_READ_AUDIT.json').read_text())
    oldarrays={x['path']:x for x in previous_inventory['archives']}
    texts=[];archives=[];histories=[];errors=[];json_count=0
    for study in ('research_max_struct_r1','research_max_struct_r2'):
        path=REPO/'research'/study
        for p in sorted(path.rglob('*')):
            if not p.is_file()or'__pycache__'in p.parts:continue
            rel=str(p.relative_to(REPO))
            if p.suffix=='.npz':
                old=PREVIOUS/rel
                if study=='research_max_struct_r1'and rel in oldarrays and old.exists()and p.read_bytes()==old.read_bytes():
                    rec=dict(oldarrays[rel]);rec['inspection']='Prior full R1 inventory reused after byte equality to previously read file; no hashes'
                    archives.append(rec);continue
                with np.load(p,allow_pickle=False)as f:
                    fields={}
                    for k in f.files:
                        a=f[k];z=dict(shape=list(a.shape),dtype=str(a.dtype),size=int(a.size))
                        if np.issubdtype(a.dtype,np.number):
                            z.update(all_finite=bool(np.isfinite(a).all()),l2_norm=float(np.linalg.norm(a.ravel())),
                                real_min=float(a.real.min())if a.size else None,real_max=float(a.real.max())if a.size else None)
                            if np.iscomplexobj(a):z.update(imag_min=float(a.imag.min())if a.size else None,imag_max=float(a.imag.max())if a.size else None)
                        fields[k]=z
                archives.append(dict(path=rel,inspection='All fields opened read-only; no forward solve',fields=fields))
            elif p.suffix in ('.md','.txt','.tex','.py','.csv','.json','.jsonl'):
                # Author-provided old hashes are not interpreted or verified.
                if 'HASHES'in p.name or'FILE_HASHES'in p.name:continue
                raw=p.read_bytes();s=raw.decode('utf-8',errors='replace');old=PREVIOUS/rel
                prior_equal=old.exists()and raw==old.read_bytes()
                rec=dict(path=rel,bytes=len(raw),lines=len(s.splitlines()),
                    byte_identical_prior_worktree=prior_equal,
                    handling='Full content opened; JSON histories programmatically parsed; prose/source semantic findings recorded in inheritance report')
                if p.suffix=='.py':
                    try:rec['functions']=[n.name for n in ast.walk(ast.parse(s))if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))]
                    except SyntaxError as e:errors.append(dict(path=rel,error=str(e)))
                if p.suffix=='.json':
                    try:
                        data=json.loads(s);json_count+=1
                        rec['top_keys']=list(data)if isinstance(data,dict)else None
                        def walk(x,at=''):
                            if isinstance(x,dict):
                                for k,v in x.items():
                                    here=at+'.'+k if at else k
                                    if k=='history'and isinstance(v,list):
                                        status=Counter(str(a.get('status',a.get('event','iteration')))for a in v if isinstance(a,dict))
                                        stages=Counter(str(a.get('stage','NA'))for a in v if isinstance(a,dict))
                                        scalar_keys=sorted(set(kk for a in v if isinstance(a,dict)for kk,vv in a.items()if isinstance(vv,(int,float,str,bool))or vv is None))
                                        histories.append(dict(source=rel,json_path=here,entries=len(v),status_counts=dict(status),stage_counts=dict(stages),scalar_keys=scalar_keys,
                                            note='Nested copies can duplicate the same trajectory. Counts describe parsed storage, not independent runs or chargeable work.'))
                                    walk(v,here)
                            elif isinstance(x,list):
                                for i,v in enumerate(x):
                                    if isinstance(v,(dict,list)):walk(v,f'{at}[{i}]')
                        walk(data)
                    except (ValueError,TypeError)as e:errors.append(dict(path=rel,error=str(e)))
                texts.append(rec)
    result=dict(classification='INHERITED_DISCLOSED_R1_R2_READ_ONLY',fresh_R3_truth_opened=False,
      text_records=texts,arrays=archives,parsed_history_records=histories,parse_errors=errors,
      text_count=len(texts),json_parsed=json_count,array_archive_count=len(archives),
      array_field_count=sum(len(x['fields'])for x in archives),R1_archives_reused=sum('reused'in x['inspection']for x in archives),
      all_numeric_arrays_finite=all(v.get('all_finite',True)for x in archives for v in x['fields'].values()),
      new_physical_LU=0,new_forward_RHS=0,no_physics_reason='Only files, AST, JSON and saved-array arithmetic read; no solver imported',
      wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,
      timer_scope='Read audit after imports; own JSON output excluded',hash_checks=False)
    (ROOT/'audit/INHERITED_READ_COVERAGE.json').write_text(json.dumps(result,indent=2,allow_nan=False))
    print(json.dumps({k:result[k]for k in ['text_count','json_parsed','array_archive_count','array_field_count','R1_archives_reused','all_numeric_arrays_finite','parse_errors','wall_seconds']}))

if __name__=='__main__':
    with threadpool_limits(limits=1):read()
