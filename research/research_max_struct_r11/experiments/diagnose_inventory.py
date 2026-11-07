"""Read-only inventory of root-approved exposed R10 endpoints; no numerical replay."""
from pathlib import Path
import json, zipfile, struct, ast
ROOT=Path(__file__).absolute().parents[1]
R10=ROOT.parents[2]/'r10_worktree/research/research_max_struct_r10'
OUT=ROOT/'delegated/diagnostics'

def schema(path):
    result={}
    with zipfile.ZipFile(path) as z:
        for name in z.namelist():
            with z.open(name) as b:
                magic=b.read(8); width=2 if magic[6]==1 else 4
                n=struct.unpack('<H' if width==2 else '<I',b.read(width))[0]
                h=ast.literal_eval(b.read(n).decode())
                result[name.removesuffix('.npy')]=dict(shape=h['shape'],dtype=h['descr'])
    return result

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    rows=[]
    for scope,ids in [('twod',(1,4,22)),('conditional_confirmation',(0,3))]:
        for index in ids:
            cid=f'r10_twod_{index:04d}' if scope=='twod' else f'r10_conditional_{index:04d}'
            folder=R10/'RESULTS'/scope/(f'{cid}__matrix_green19__noise_ball_99' if scope=='twod' else f'fallback/{cid}__conditional_fallback_G19_PGP')
            files={name:str(folder/name) for name in ('PUBLIC_TASK.json','FROZEN_CONFIG.json','NUMERICAL_SCHEDULE.json','SOURCE_METRICS.npz','DENOMINATOR_BASIS.npz','DENOMINATOR_SETUP.json')}
            stages=[]
            for stage in (0,3):
                path=folder/f'stage_{stage:02d}_joint.npz'; meta=folder/f'stage_{stage:02d}_joint.json'
                info=json.loads(meta.read_text()) if meta.exists() else {}
                stages.append(dict(stage=stage,path=str(path),available=path.is_file(),schema=schema(path) if path.is_file() else None,
                    metadata_path=str(meta),metadata=info,next_true_material_array_available=False,
                    next_true_material_note='Stage endpoint is after final material readout; no per-inner next material saved. Later stage material includes intervening current updates and is not same-state next step.'))
            pubdir=R10/'data_public'/scope
            rows.append(dict(case_id=cid,scope=scope,folder=str(folder),files=files,file_available={k:Path(v).is_file() for k,v in files.items()},
                portable_public_files=[str(pubdir/'operators.npz'),str(pubdir/f'{cid}.npz')],
                portable_truth_path_OFFLINE_ONLY=str(R10/'sealed_truth'/scope/f'{cid}.npz'),stages=stages,
                metric_schema=schema(folder/'SOURCE_METRICS.npz'),basis_schema=schema(folder/'DENOMINATOR_BASIS.npz')))
    payload=dict(status='COMPLETE',truth_access='NONE; paths inventoried only',historical_data_status='EXPOSED_PUBLIC_R10',
      scope='Root approved five objects, G19 PGP trajectories, stage0/mf3 and stage3/mf12 endpoints only',
      interpretation='Full denotes full derivative of saved low-rank metric, not full-rank Green operator',objects=rows)
    (OUT/'STATE_INVENTORY.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(dict(objects=len(rows),states=sum(s['available'] for r in rows for s in r['stages']),output=str(OUT/'STATE_INVENTORY.json'))))
if __name__=='__main__':main()
