"""Explicit truth-free context_x export from proven online reconstruction x.

Never accept a generic joint NPZ's x: in that format x is the supervised label.
The source here is only the canonical per-parent online reconstruction folder,
with exact schema and Gaussianized-anchor/material consistency checks.
"""
from pathlib import Path
import argparse,sys,json,numpy as np
from scipy.special import ndtr
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,Book,Target,save,write
from p1_material.prior import MaterialPrior
from p3r4.evaluation import checkpoint_identity

def main():
    p=argparse.ArgumentParser();p.add_argument('group',choices=['dev','train','validation','locked','confirmation'])
    p.add_argument('--count',type=int,required=True);a=p.parse_args()
    prefix=dict(dev='R4DEV',train='R2_train',validation='R2_validation',locked='R4L1',confirmation='R4L2')[a.group]
    stage='IV' if a.group in ['locked','confirmation'] else 'I'
    if stage=='IV' and not (OUT/f'configs/FROZEN_{a.group.upper()}.json').exists():raise ValueError('test protocol freeze required')
    dest=OUT/'inputs/anchor_conditions'/a.group
    if dest.exists():raise FileExistsError('do not overwrite exported conditions')
    book=Book(stage,'export_anchor_conditions_'+a.group,analysis_only=True,max_wall=180)
    status='COMPLETE';error=None
    try:
        cfg=json.loads((OUT/'configs/TARGET.json').read_text());points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
        prior=MaterialPrior(cfg,points);dest.mkdir(parents=True);rows=[]
        expected={'observation','q0','chi0','J','prediction','residual','phi','x','H','parent'}
        for k in range(a.count):
            parent=f'{prefix}_{k:05d}';source=OUT/'results/contexts'/a.group/f'{parent}.npz'
            before=checkpoint_identity(source)
            with np.load(source,allow_pickle=False) as z:
                if set(z.files)!=expected or str(z['parent'].item())!=parent:raise ValueError('canonical online-only context schema/identity required')
                values={key:z[key] for key in expected}
            q,x=values['q0'],values['x']
            if q.shape!=(66,) or x.shape!=(66,) or not np.isfinite(x).all() or not np.isfinite(q).all():raise ValueError('finite proven online x/q0 required')
            if np.max(abs(x[:64]-q[:64]))>1e-12 or np.max(abs(ndtr(x[64:])-q[64:]))>1e-12:raise ValueError('x is not Gaussianized online anchor q0')
            chi=prior.decode(q)[0]
            if not np.array_equal(chi,values['chi0']):raise ValueError('online anchor material inconsistent with original decoder')
            if values['observation'].shape!=(32,8) or values['J'].shape!=(512,66):raise ValueError('frozen observation/J shape')
            if np.max(abs(values['observation']-values['prediction']-values['residual']))>1e-12:raise ValueError('online prediction/residual identity')
            if checkpoint_identity(source)!=before:raise ValueError('online context changed during export')
            # Aliases are explicit only after the above provenance checks.
            output={key:values[key] for key in ['parent','observation','q0','chi0','J','residual','prediction','phi']}
            output.update(context_x=x,context_J=values['J'],base_mean=x,target_id=np.array(cfg['target_id']),sigma=np.array(cfg['noise']['sigma']),
                          source_context_identity=np.array(json.dumps(before,sort_keys=True)),context_role=np.array('proven_online_reconstruction_not_GT'))
            save(dest/f'{parent}.npz',**output);rows.append(dict(parent=parent,source=before,export=checkpoint_identity(dest/f'{parent}.npz')))
        write(dest/'EXPORT_COMPLETE.json',dict(status='COMPLETE',group=a.group,parents=a.count,rows=rows,
            truth_read=False,new_physics_RHS=0,source_scope='canonical isolated online reconstruction only',
            ambiguous_joint_x_never_accepted=True,source_context_cost='paid96+64 reconstruction; reused not free cold setup'))
    except BaseException as e:status='FAILED';error=repr(e);raise
    finally:book.close(status,error)
if __name__=='__main__':main()
