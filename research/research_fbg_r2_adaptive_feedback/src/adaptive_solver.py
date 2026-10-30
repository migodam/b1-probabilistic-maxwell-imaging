"""Finite, counted FBG research executor over the existing resident CUDA physics.

No truth/labels are read here. Every optimization is ranked by original data
loss plus the same R4 spatial regularizer. Feedback is frozen within a phase.
"""
from __future__ import annotations
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import importlib.util
import json
from pathlib import Path
import sys
import time
import traceback
import numpy as np
import torch

HERE=Path(__file__).resolve().parents[1]
PIPELINE=Path(os.environ.get('CUDA_PIPELINE_ROOT',str(HERE.parent/'research_z_alg_cuda_pipeline')))
def _module(name,file):
    spec=importlib.util.spec_from_file_location(name,file)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
    return module
base=_module('fbg_cuda_device',PIPELINE/'src/device.py')
engine=_module('fbg_cuda_executor',PIPELINE/'src/executor.py')
from .feedback_fit import FeedbackOperator, fit_feedback, feedback_value_gradient


class FeedbackLane(base.ObjectiveLane):
    def __init__(self,*args,feedback=None,task_context=None,**kwargs):
        super().__init__(*args,**kwargs)
        self.feedback=feedback
        self.task_context=task_context
        self.budget_enabled=True

    def _budget(self):
        if not self.budget_enabled or self.task_context is None:return
        counts=Counter()
        for lane in self.task_context['lanes']:counts.update(lane.counts)
        cfg=self.task_context['config']
        rhs=sum(counts[k] for k in ('source_RHS','adjoint_RHS','receiver_RHS','tangent_RHS'))
        if counts['LU']>=cfg.get('max_lu',float('inf')):raise engine.BudgetStop('LU_CAP')
        if rhs>=cfg.get('max_rhs',float('inf')):raise engine.BudgetStop('RHS_CAP')

    def _forward(self,x):
        self._budget()
        return super()._forward(x)

    def close(self,x):
        self.budget_enabled=False
        result=super().close(x)
        self.budget_enabled=True
        return result

    def evaluate(self,x):
        if self.feedback is None and self.objective!='fixedgram':
            return super().evaluate(x)
        s=self.shared
        with torch.cuda.stream(self.stream),torch.no_grad():
            z=self._forward(x)
            if self.objective=='fixedgram' and self.feedback is None:
                qr=s.feedback_reference_q
                chol,gram_info=torch.linalg.cholesky_ex(qr.mH@qr,check_errors=False)
                weighted=torch.cholesky_solve(z['r'].mH,chol).mH
                factor=.5*s.reference_field_mean/self.scale2
                value=factor*(z['r'].conj()*weighted).sum().real
                adj=torch.linalg.lu_solve(z['lu'],z['pivots'],s.SH@weighted,adjoint=True)
                g=2*factor*(adj.conj()*z['q']).sum(1).real
            else:
                factor=.5*s.reference_field_mean*self.feedback.reference_scale/self.scale2
                value,g,gram_info=feedback_value_gradient(s,z,self.feedback,factor)
            d=z['d']
            gd=(2*self.config['smooth_weight']*d+self.config['tv_weight']*d/z['tvq'])/(s.edges*19.5)
            greg=torch.zeros_like(z['m']);greg.index_add_(0,s.right,gd);greg.index_add_(0,s.left,-gd)
            g=g+greg
            if self.stage!=s.side:g=self.B.T@g
            packed=torch.cat(((value+z['reg']).reshape(1),g,z['data'].reshape(1),z['reg'].reshape(1),z['info'].double().reshape(1),gram_info.double().reshape(1)))
            result=self._pull(packed,self.output_host)
            self.counts.update({'adjoint_RHS':s.ns,'gradient':1,'Gram':1,'feedback_evaluations':1})
        if result[-1]!=0 or result[-2]!=0:raise ArithmeticError('LU/feedback Gram failed')
        self.last_physical_objective=float(result[-4]+result[-3])
        return float(result[0]),result[1:1+self.n_coeffs]


def _spatial_directions(side,m):
    """Public cosine/region modes, including amplitude; no GT geometry."""
    yy,xx=np.meshgrid(np.linspace(0,1,side),np.linspace(0,1,side),indexing='ij')
    hs=[np.ones(side*side)]
    for fy,fx in ((1,0),(0,1),(1,1),(2,0),(0,2)):
        hs.append((np.cos(np.pi*fy*yy)*np.cos(np.pi*fx*xx)).ravel())
    hs.append((np.asarray(m)-np.mean(m)))
    return np.asarray([h/(np.max(np.abs(h))+1e-12) for h in hs])


def make_feedback(lane,x,kind,cfg):
    """Anchor-only exact current fitting; all source RHS and fitting charged."""
    s=lane.shared;t0=time.perf_counter()
    with torch.cuda.stream(lane.stream),torch.no_grad():
        z=lane._forward(x)
        m=lane._pull(z['m'])
        needs_tangent=kind in ('diag','block','shuffle') or cfg.get('record_tangent_metrics',False)
        if needs_tangent and (getattr(lane,'feedback_anchor_x',None) is None or not np.array_equal(x,lane.feedback_anchor_x)):
            hs=torch.as_tensor(_spatial_directions(s.side,m),dtype=torch.float64,device=s.device)
            rhs=(hs[:,:,None]*z['q'][None,:,:]).permute(1,0,2).reshape(s.n,-1)
            tangents=torch.linalg.lu_solve(z['lu'],z['pivots'],rhs).reshape(s.n,len(hs),s.ns).permute(1,0,2)
            lane.feedback_anchor_x=np.asarray(x).copy();lane.feedback_anchor_tangents=tangents
            lane.counts.update({'tangent_RHS':int(len(hs)*s.ns)})
        tangents=lane.feedback_anchor_tangents if needs_tangent else None
        lane.counts['feedback_fits']+=1
        calibration=None
        if kind=='qbeta':
            if 'calibrated_beta' not in cfg:
                qr=s.feedback_reference_q;rr=s.S@s.feedback_reference_j-lane.target
                cf=torch.linalg.cholesky(qr.mH@qr)
                wr=torch.cholesky_solve(rr.mH,cf).mH;D=wr.mH@wr
                adj=torch.linalg.lu_solve(s.feedback_reference_lu,s.feedback_reference_pivots,
                                        s.SH@wr-s.GH@(qr@D),adjoint=True)
                g0=(adj.conj()*qr).sum(1).real;d=((qr@D).conj()*qr).sum(1).real
                gd=lane._pull(torch.stack((g0,d)));g,d=gd[0],gd[1]
                valid=(g>0)&(d>np.median(d));t=g[valid]/d[valid];t=t[t<10.]
                thresholds=t/(1-.1*t)
                beta=float(np.clip(np.median(thresholds),.05,.4)) if len(thresholds) else 0.
                cfg['calibrated_beta']=beta
                calibration={'beta':beta,'eligible_pixels':len(thresholds),
                             'rule':'median exact uniform-reference positive-gradient flip threshold in upper half d; clip .05..4; empty0',
                             'truth_access':'NONE'}
                lane.counts.update({'calibration_adjoint_RHS':s.ns,'adjoint_RHS':s.ns,'Gram':1,'beta_calibrations':1})
            beta=cfg['calibrated_beta']
            feedback=FeedbackOperator.diagonal(torch.full((s.n,),beta,dtype=torch.complex128,device=s.device),name=kind)
        elif kind.startswith('beta'):
            beta=float(kind.split(':')[1])
            feedback=FeedbackOperator.diagonal(torch.full((s.n,),beta,dtype=torch.complex128,device=s.device),name=kind)
        elif kind=='self':
            feedback=FeedbackOperator.diagonal(-torch.diagonal(s.G),name=kind)
        elif kind=='zero':
            feedback=FeedbackOperator.diagonal(torch.zeros(s.n,dtype=torch.complex128,device=s.device),name=kind)
        else:
            variant='block' if kind=='block' else 'complex'
            feedback=fit_feedback(s.G,z['q'],tangents,kind=variant,side=s.side,ridge=float(cfg.get('feedback_ridge',.05)))
            feedback.name=kind
            if kind=='shuffle':
                feedback=FeedbackOperator.diagonal(torch.roll(feedback.diag,37),name=kind)
        feedback,box_receipt=feedback.box_certify(lower_m=.02,upper_m=20.,perturbation_cap=.8)
        # One common public reference; scale cannot be chosen from the target.
        ref=torch.full((s.n,),.1,dtype=torch.float64,device=s.device)
        ref_q=s.feedback_reference_q
        ref_v=ref_q+feedback.apply(ref[:,None]*ref_q)
        feedback.reference_scale=(ref_v.abs().square().sum()/ref_q.abs().square().sum()).real
        safety=feedback.safety(z['m'],lower=.02,upper=20.)
        lane.counts.update({'feedback_reference':1})
        v=z['q']+feedback.apply(z['m'][:,None]*z['q'])
        energy=None
        if tangents is not None:
            norm0=(s.G@tangents).abs().square().sum()
            norm1=(s.G@tangents+feedback.apply(tangents)).abs().square().sum()
            energy=float(lane._pull((norm1/norm0).reshape(1))[0])
        details={'kind':kind,'tangent_energy_ratio':energy,'calibration':calibration,
                 'safety':safety,'whole_box':box_receipt,'fit_wall_s':time.perf_counter()-t0,'directions':0 if tangents is None else len(tangents),
                 'normalization':'fixed public material .1; not current auxiliary value'}
        feedback.details=details
    return feedback,details


def select_feedback(lane,x,cfg):
    """Finite physical pilot: auxiliary scale is never the ranking criterion."""
    records=[];best=None;best_score=float('inf')
    candidates=cfg.get('feedback_candidates',['zero','beta:0.2','self','diag','block','shuffle'])
    for kind in candidates:
        feedback,detail=make_feedback(lane,x,kind,cfg)
        lane.feedback=feedback
        value,g=lane.evaluate(x)
        h=-np.asarray(g);h*=.5/(np.max(np.abs(h))+1e-15)
        trial=np.clip(np.asarray(x)+h,.02,20.)
        with torch.cuda.stream(lane.stream),torch.no_grad():
            z=lane._forward(trial)
            score=float(lane._pull((z['data']+z['reg']).reshape(1))[0])
        records.append({**detail,'pilot_physical_score':score,'auxiliary_value_at_anchor':value})
        if score<best_score:best,best_score=feedback,score
    lane.feedback=best
    return best,{'selected':best.name,'pilot_candidates':records,'selection':'minimum original physical score after .5-contrast bounded gradient pilot'}


def solve(shared,public,task,context=None):
    cfg=dict(task.get('config',{}));method=task['method']
    wall=float(cfg.get('wall_seconds',45));started=time.perf_counter()
    fractions=cfg.get('stage_fractions',[.2,.25,.55]);stages=cfg.get('stages',[8,12,24])
    old=None;coeff=None;phases=[];fit_records=[]
    task_calibrated_beta=None
    for stage,fraction in zip(stages,fractions):
        remaining=wall-(time.perf_counter()-started)
        if remaining<=0:break
        c={**cfg,'stage':stage,'wall_seconds':min(remaining,wall*fraction)}
        if task_calibrated_beta is not None:c['calibrated_beta']=task_calibrated_beta
        if coeff is not None:c['initial_coefficients']=engine._resample(coeff,old,stage)
        elif 'initial_pixels' in task:
            c['initial_coefficients']=engine._resample(task['initial_pixels'],24,stage)
        else:c['initial_material']=task.get('initial_material',.1)
        x=np.asarray(c.get('initial_coefficients',np.full(stage*stage,c.get('initial_material',.1))))
        is_trf=method in ('native_trf','reciprocal_trf')
        feedback_kind=None
        if method.startswith('feedback_'):feedback_kind=method[len('feedback_'):]
        use_dual=method in ('gd_all','globaldual','fixedgram') or (method=='gd_late' and stage==24) or feedback_kind is not None
        objective='fixedgram' if method=='fixedgram' else ('dual' if use_dual else 'data')
        lane=FeedbackLane(shared,public,stage=stage,objective=objective,config=c,task_context=context)
        if context is not None:context['lanes'].append(lane)
        try:
            lane._budget()
            if feedback_kind:
                if feedback_kind=='selected':fb,record=select_feedback(lane,x,c)
                else:fb,record=make_feedback(lane,x,feedback_kind,c)
                lane.feedback=fb;fit_records.append({'stage':stage,**record})
                if 'calibrated_beta' in c:task_calibrated_beta=c['calibrated_beta']
        except engine.BudgetStop as exc:
            result=lane.close(x)
            result.update({'stage':stage,'phase':'feedback_budget_stop','wrapper_status':'BUDGET_STOP',
                           'wrapper_stop_reason':str(exc),'evaluation_trace':[]})
            phases.append(result)
            break
        # Same proposal/correction split for data-LBFGS and all weighted arms.
        if is_trf:
            plan=[(method,1.)]
        else:plan=[('dual' if use_dual else 'data',.55),('data',.45)]
        phase_left=c['wall_seconds']
        for phase_i,(optimization,pfrac) in enumerate(plan):
            rem=min(wall-(time.perf_counter()-started),phase_left)
            if rem<=0:break
            pc={**c,'wall_seconds':rem if phase_i==len(plan)-1 else rem*pfrac,
                'initial_coefficients':x,'maxiter':cfg.get('maxiter',2000),'maxfun':cfg.get('maxfun',4000)}
            if phase_i:
                previous_counts=Counter(lane.counts)
                lane.feedback=None;lane.objective='data'
            else:previous_counts=Counter()
            pstart=time.perf_counter()
            result=engine.solve_task(lane,{'method':optimization,'stage':stage,'config':pc})
            x=result.pop('_coefficients');phase_left-=time.perf_counter()-pstart
            result['counts']=dict(Counter(result['counts'])-previous_counts)
            result['phase']='corrector' if phase_i else 'proposal'
            phases.append(result)
        coeff=x;old=stage
    if not phases:raise RuntimeError('No completed physical phase')
    # Prolongation preserves every saved coarse material: coarse candidates are
    # eligible and are identified; do not pretend every selected image is fine.
    chosen=min(phases,key=lambda p:p['ranking_objective'])
    counts=Counter()
    for row in phases:counts.update(row['counts'])
    return {**chosen,'method':method,'counts':dict(counts),'phases':phases,'feedback_records':fit_records,
            'task_wall_s':time.perf_counter()-started,'stages_reached':sorted(set(p['stage'] for p in phases)),
            'online_selection':'minimum original full-data loss + common regularizer over paid candidates',
            'truth_access':'NONE','budget':cfg}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('manifest');a=ap.parse_args()
    manifest_path=Path(a.manifest).resolve();root=manifest_path.parent
    cfg=json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    if cfg.get('vendor_path'):sys.path.insert(0,cfg['vendor_path'])
    from threadpoolctl import threadpool_limits
    torch.set_num_threads(1);torch.cuda.set_device(0)
    started=time.perf_counter();out=Path(cfg['output_dir']);out.mkdir(parents=True,exist_ok=True)
    shared_cache={}
    setup_counts=Counter()
    for operator_path in dict.fromkeys(task.get('operators_path',cfg['operators_path']) for task in cfg['tasks']):
        with np.load(operator_path,allow_pickle=False) as z:operators=dict(z)
        shared=base.SharedOperators(operators)
        ref_lu=torch.linalg.lu_factor_ex(shared.eye-.1*shared.G,check_errors=False)
        shared.feedback_reference_j=torch.linalg.lu_solve(ref_lu.LU,ref_lu.pivots,.1*shared.E0)
        shared.feedback_reference_q=shared.E0+shared.G@shared.feedback_reference_j
        shared.feedback_reference_lu,shared.feedback_reference_pivots=ref_lu.LU,ref_lu.pivots
        shared.setup_counts.update({'LU':2,'source_RHS':2*shared.ns})
        shared.ready.record()
        shared_cache[operator_path]=shared;setup_counts.update(shared.setup_counts)
    cache={}
    for task in cfg['tasks']:
        if task['public_path'] not in cache:
            with np.load(task['public_path'],allow_pickle=False) as z:cache[task['public_path']]=dict(z)
    setup=time.perf_counter()-started
    hybrid=_module('fbg_hybrid_svd',PIPELINE/'src/hybrid_svd.py').HybridSVD();hybrid.install()
    receipts=[]
    def worker(task):
        began=time.perf_counter()
        context={'lanes':[],'config':task.get('config',{})}
        try:
            shared=shared_cache[task.get('operators_path',cfg['operators_path'])]
            result=solve(shared,cache[task['public_path']],task,context)
            row={'task':task,'status':'OK','result':result}
        except Exception as exc:
            partial=Counter()
            for lane in context['lanes']:partial.update(lane.counts)
            row={'task':task,'status':'FAIL','error':str(exc),'traceback':traceback.format_exc(),
                 'partial_costs':dict(partial),'costs_scope':'completed operations; interrupted operation may be uncounted'}
        row['wall_s']=time.perf_counter()-began
        row['finished_s']=time.perf_counter()-started
        return row
    try:
        with threadpool_limits(limits=1),ThreadPoolExecutor(max_workers=cfg.get('lanes',4)) as pool:
            futures=[pool.submit(worker,t) for t in cfg['tasks']]
            for future in as_completed(futures):
                row=future.result();receipts.append(row)
                engine.write_json(out/(row['task']['task_id']+'.json'),row)
                print(row['task']['task_id'],row['status'],round(row['wall_s'],2),flush=True)
    finally:hybrid.restore()
    engine.write_json(out/'SUMMARY.json',{'utc':datetime.now(timezone.utc).isoformat(),'tasks':len(receipts),
          'ok':sum(r['status']=='OK' for r in receipts),'wall_s':time.perf_counter()-started,
          'setup_wall_s':setup,'setup_counts':dict(setup_counts),'manifest':cfg,'svd':hybrid.summary()})

if __name__=='__main__':main()
