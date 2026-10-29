"""Evidence-based data-flow audit. This is not an OS inaccessible sandbox."""
from pathlib import Path
import argparse,ast,json,time
import numpy as np
from generate import ROOT,write_json
FORBIDDEN=('inputs/private','inputs\\private','truth.npz','truth.json','SCORES.json','score_endpoint','from generate','import generate','from score','import score','observation.json','shape_stratum')
def main():
    p=argparse.ArgumentParser();p.add_argument('--online-script',type=Path,action='append',default=[]);p.add_argument('--access-log',type=Path,action='append',default=[]);p.add_argument('--output',type=Path,default=ROOT/'results/ISOLATION_AUDIT.json');a=p.parse_args()
    rows=[]
    for path in sorted((ROOT/'inputs/online').glob('*')):
        npz=path/'observation.npz'
        if not npz.exists():continue
        bundle=np.load(npz);keys=list(bundle.keys());cfg=json.loads((path/'config.json').read_text())
        forbidden_config=[k for k in cfg if any(s in k.lower() for s in ('truth','clean','sampled','cavity_cell','target_pose','wall_eps','wall_loss'))]
        metadata=json.loads((path/'observation.json').read_text());shape_keys=[k for k in metadata if k in ('stratum','shape','shape_stratum','target')]
        rows.append(dict(kind='online_input',path=str(path.relative_to(ROOT)),keys=keys,data_shape=list(bundle['y'].shape) if 'y' in keys else None,pass_y_only=keys==['y'],forbidden_config_keys=forbidden_config,forbidden_shape_metadata_keys=shape_keys,passed=keys==['y'] and not forbidden_config and not shape_keys))
    configs=[json.loads(p.read_text()) for p in sorted((ROOT/'inputs/online').glob('*/config.json'))]
    if configs:rows.append(dict(kind='public_config_uniformity',n_cases=len(configs),passed=all(c==configs[0] for c in configs),scope='No case-specific truth coefficient or shape in public configuration'))
    for path in a.online_script:
        text=path.read_text();ast.parse(text)
        hits=[token for token in FORBIDDEN if token in text]
        rows.append(dict(kind='online_source_static',path=str(path),forbidden_tokens=hits,passed=not hits,scope='literal path/import scan only; dynamic reads not proven absent'))
    for path in a.access_log:
        events=[json.loads(l) for l in path.read_text().splitlines() if l]
        opens=[e for e in events if e.get('event')=='file_open'];denied=[e for e in opens if e.get('denied')]
        forbidden_reads=[e for e in opens if not e.get('write') and ('/inputs/private/' in e['path'] or e['path'].endswith('/SCORES.json'))]
        metadata_reads=[e['path'] for e in opens if not e.get('write') and e['path'].endswith('/observation.json')]
        rows.append(dict(kind='dynamic_python_access_log',path=str(path),open_events=len(opens),denied_reads=denied,forbidden_reads=forbidden_reads,observation_metadata_reads=metadata_reads,finished=any(e.get('event')=='guard_finished' for e in events),passed=not forbidden_reads and not metadata_reads,scope='Python open audit; native library and child process reads outside scope'))
    for seal in sorted((ROOT/'results').glob('**/ONLINE_SEALED.json')):
        out=seal.parent;scores=out/'SCORES.json';endpoint=out/'endpoint.npz'
        fullseal=json.loads(seal.read_text());summary=fullseal.get('summary',{});compactseal={k:fullseal.get(k) for k in ['time_unix','selected_by','truth_access','diagnostic_only']};compactseal.update({k:summary.get(k) for k in ['method','selected_start','objective','converged']})
        row=dict(kind='endpoint_seal',path=str(out.relative_to(ROOT)),endpoint_exists=endpoint.exists(),score_exists=scores.exists(),seal_summary=compactseal)
        if endpoint.exists():row['endpoint_before_seal']=endpoint.stat().st_mtime<=seal.stat().st_mtime+.01
        if scores.exists():row['score_after_seal']=scores.stat().st_mtime>=seal.stat().st_mtime
        row['passed']=row['endpoint_exists'] and row.get('endpoint_before_seal',False) and row.get('score_after_seal',True);rows.append(row)
    for scope in sorted((ROOT/'results/online').glob('*/KNOWN_WALL_ORACLE/ORACLE_SCOPE.json')):
        d=scope.parent;case=d.parent.name.rsplit('_r',1)[0];wall=np.load(ROOT/'inputs/private'/case/'oracle_wall.npz');seal=json.loads((d/'ONLINE_SEALED.json').read_text());keys=list(wall.keys());ok=set(keys)=={'eps','loss'} and wall['eps'].shape==(6,) and wall['loss'].shape==(6,) and bool(seal.get('diagnostic_only')) and 'wall' in str(seal.get('truth_access',''))
        rows.append(dict(kind='isolated_known_wall_diagnostic',path=str(d.relative_to(ROOT)),private_input_keys=keys,private_input_shapes={k:list(wall[k].shape) for k in keys},true_target_in_input=False,diagnostic_only=seal.get('diagnostic_only'),passed=ok,scope='Expected labeled wall truth; separate process, target initialized from unchanged public seeds. Source/data schema evidence, not OS file-open proof.'))
    result=dict(status='PASS_DATA_FLOW' if rows and all(r['passed'] for r in rows) else 'FAIL_OR_NOT_READY',checked_unix=time.time(),rows=rows,truth_exposure='Development and each blind case truth only consumed by post-seal scorer/plotter; known-wall oracle must be explicitly labeled and separate',limitations=['Filesystem remains shared and readable; no OS access restriction','Literal static scan cannot certify arbitrary dynamic file access','Audit establishes artifact order and input contents, not human memory absence'])
    write_json(a.output,result);print(json.dumps(dict(status=result['status'],rows=len(rows),failed=sum(not r['passed'] for r in rows))))
if __name__=='__main__':main()
