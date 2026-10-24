#!/usr/bin/env python3
"""Explicitly released arbitrary material-action Full diagnostics, no GT/prior law.

Uses the unchanged backend directly, never Target.evaluate(q); no material
projection or decoder-realizability inference is performed.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import argparse,json,sys,time,traceback,resource
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.representative_eval import action_properties
from p3r4.evaluation import checkpoint_identity,write_new_json
STAGE=ROOT/'research/p3l_r4_conditional_transport';EXPERIMENT='raw_material_forward_diagnostic_v1'


def freeze(a):
    prepared=json.loads(a.prepared_request.read_text())
    if prepared.get('status')!='FROZEN_READY' or prepared.get('GT_read') is not False:raise ValueError('observation-only frozen representative request required')
    if any(p.startswith(('R4L1_','R4L2_')) for p in prepared['request']['parents']) and a.stage!='IV':raise ValueError('fresh raw diagnostic freeze requires explicit StageIV')
    rows=[]
    for item in prepared['qualified']:
        for key in item['outputs']:rows.append(dict(parent=item['parent'],action=key,frozen_mu=item['frozen_mu'],mu_identity=checkpoint_identity(item['frozen_mu']),observation=item['online_observation'],observation_identity=checkpoint_identity(item['online_observation'])))
    spec=dict(status='FROZEN_ARBITRARY_MATERIAL_FORWARD',stage=a.stage,experiment=EXPERIMENT,rows=rows,planned_RHS=8*len(rows),GT_read=False,
              decoder_realizability='NOT_PROVEN',P1_legal_world='NOT_CLAIMED',arbitrary_materials_preserved_no_clipping=True,declared_parents=prepared['request']['parents'],predeclared_unqualified_skips=prepared['skips'],scope='qualified frozen_mu snapshots only; other declared parents explicit NOT_EVALUATED')
    write_new_json(a.out,spec);print(json.dumps(dict(status=spec['status'],actions=len(rows),planned_RHS=spec['planned_RHS'])))


def run(a):
    spec=json.loads(a.freeze.read_text());release=json.loads(a.release.read_text())
    request_stage=spec.get('stage','II')
    if a.stage!=request_stage or release.get('stage','II')!=request_stage:raise ValueError('explicit raw request/release/CLI stage mismatch')
    if any(p.startswith(('R4L1_','R4L2_')) for p in spec.get('declared_parents',[])) and request_stage!='IV':raise ValueError('fresh locked/confirmation raw diagnostics require StageIV')
    if spec['status']!='FROZEN_ARBITRARY_MATERIAL_FORWARD' or release.get('released') is not True or release.get('experiment')!=EXPERIMENT:raise ValueError('explicit matching release/freeze required')
    if not a.job or a.max_rhs is None or a.max_wall is None or a.max_rhs<spec['planned_RHS'] or a.max_wall<=15:raise ValueError('explicit bounded job/RHS/wall required')
    if a.out.exists():raise FileExistsError('no repeated paidrun')
    for row in spec['rows']:
        if 'offline_only' in Path(row['frozen_mu']).parts or 'offline_only' in Path(row['observation']).parts:raise ValueError('no offline data in physics')
        if checkpoint_identity(row['frozen_mu'])!=row['mu_identity'] or checkpoint_identity(row['observation'])!=row['observation_identity']:raise ValueError('frozen online action/observation changed')
    # Only here does the explicitly released physical driver import Book/backend.
    from p3r4.core import Book,BudgetStop,write,save,BoundedTransfers
    from p1_material.prior import MaterialPrior
    cfg=json.loads((STAGE/'configs/TARGET.json').read_text());points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points'];roi=np.flatnonzero(MaterialPrior(cfg,points).mask)
    sys.path.insert(0,str(ROOT.parents[1]/'A30/a30-maxwell-material-manifold/src'));from a30.physics import FrozenMaxwellAdapter
    book=Book(request_stage,a.job,max_wall=a.max_wall);original=book._update;costs=dict(RHS=0,states=0,LU=0);status='COMPLETE';error=None;rows=[];start=time.monotonic()
    def limited(event,**kw):
        if costs['RHS']+kw.get('rhs',0)>a.max_rhs:raise BudgetStop('raw material diagnostic hardRHS cap')
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)>3*1024**3:raise BudgetStop('3GiB cap')
        result=original(event,**kw)
        for k,x in [('RHS','rhs'),('states','states'),('LU','lu')]:costs[k]+=kw.get(x,0)
        return result
    book._update=limited
    try:
        a.out.mkdir(parents=True);write(a.out/'FROZEN_REQUEST.json',spec)
        model=FrozenMaxwellAdapter.from_background_npz(STAGE/'inputs/background/weak.npz',book,purpose_prefix='R4_ARBITRARY_ACTION_DIAGNOSTIC');model._transfers=BoundedTransfers(model._transfers,limit=2)
        for index,item in enumerate(spec['rows']):
            record=dict(parent=item['parent'],action=item['action'],status='FAILED',decoder_realizability='NOT_PROVEN',P1_legal_world='NOT_CLAIMED');before=costs['RHS']
            try:
                with np.load(item['frozen_mu'],allow_pickle=False) as z:chi=z[item['action']]
                with np.load(item['observation'],allow_pickle=False) as z:y=z['observation']
                properties=action_properties(chi,cfg,roi);book.context=dict(parent=item['parent'],action=item['action'],role='arbitrary_material_action_not_P1_world')
                save(a.out/'attempts'/f'a{index:05d}.npz',chi=chi,observation=y)
                with book.physical('arbitrary_material_Full8sources'):state=model.full_state(chi,purpose=f'raw_action_{index:05d}')
                if costs['RHS']-before!=8:raise AssertionError('raw Fullmustcharge8sourceRHS')
                if state.residual>1e-10:raise RuntimeError('raw-action backward residual contract')
                phi=float(np.linalg.norm(state.data-y)**2/cfg['noise']['sigma']**2)
                save(a.out/'endpoints'/f'a{index:05d}.npz',chi=chi,current=state.current,field=state.field,prediction=state.data,Phi=np.array(phi),backward_residual=np.array(state.residual))
                record.update(status='COMPLETE',Phi=phi,backward_residual=float(state.residual),properties=properties,new_RHS=costs['RHS']-before)
                del state;model._transfers.clear()
            except BudgetStop:raise
            except Exception as exc:record.update(error=repr(exc),new_RHS=costs['RHS']-before)
            rows.append(record);write(a.out/'PROGRESS.json',dict(rows=rows,costs=costs));print(json.dumps(record),flush=True)
        write(a.out/'RESULTS.json',dict(status='COMPLETE',stage=request_stage,rows=rows,costs=costs,GT_read=False,no_clipping=True,decoder_realizability='NOT_PROVEN',posterior_sample=False,
             scope='nonlinear Full response/passivity/ROI/rectangle diagnostics of arbitrary image only; not prior legal world',predeclared_unqualified_skips=spec.get('predeclared_unqualified_skips',[]),declared_parents=spec.get('declared_parents',[]),wall_seconds=time.monotonic()-start))
    except BaseException as exc:status='ERROR';error=repr(exc);write(a.out/'ERROR.json',dict(error=error,traceback=traceback.format_exc(),costs=costs));raise
    finally:book.close(status,error)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('freeze','run'))
    for key in ('prepared-request','freeze','release','out'):p.add_argument('--'+key,type=Path,required=key=='out')
    p.add_argument('--stage',choices=('II','IV'),default='II');p.add_argument('--job');p.add_argument('--max-rhs',type=int);p.add_argument('--max-wall',type=float);a=p.parse_args()
    if a.command=='freeze':freeze(a)
    else:run(a)

if __name__=='__main__':main()
