"""A0 measurement-only three-baseline runner. Importing this file runs nothing.

Reads inputs/online, the frozen contracts, E0 gate and the preparation noise
rule only. Holdout and physical labels belong exclusively to a separate scorer.
"""
from pathlib import Path
import argparse
import json
import time
from unittest.mock import patch

import numpy as np
from scipy.optimize import least_squares as scipy_least_squares
from scipy.special import hankel1

from wslam.accounting import BudgetBook, BudgetExceeded, json_default
from wslam.forward import Evaluator
from wslam.geometry import ROI, validate_geometry
from wslam import optimize

METHODS = ('FULL_J_TR', 'MULTISTART_JOINT_MAP',
           'HYBRID_MULTISCALE_BIRTH_DEATH_EXCHANGE')


def write_json(path, data):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data,indent=2,default=json_default,allow_nan=False)+'\n')
    temporary.replace(path)


def save_npz(path, data):
    temporary = path.with_suffix('.tmp.npz')
    np.savez_compressed(temporary,**data)
    temporary.replace(path)


def hierarchical_score(fit):
    # Same physical-coordinate prior constants in both levels; material remains
    # latent for empty wall. Each level has probability .5, conditional masks
    # carry normalized Bernoulli(.65) mass. No density from proposal-bank weights.
    return float(fit['objective']-np.log(optimize.support_mass(fit['support']))+np.log(2.))


def born_prediction(support, theta, degree=4):
    """One common analytic free-space Born score for both public grid levels."""
    s = np.asarray(support,bool)
    if s.shape not in ((2,2),(3,3)) or not validate_geometry(s)['legal']:
        raise ValueError('Illegal public occupancy in initializer pool')
    tx = np.array([[-1.5,0.],[1.5,0.]])
    rx = np.column_stack((np.linspace(-2,2,16),np.zeros(16)))
    k = 2*np.pi
    def green(points,sources):
        r = np.linalg.norm(points[:,None,:]-sources[None,:,:],axis=-1)
        if np.any(r==0):
            raise ValueError('Singular analytic quadrature')
        return .25j*hankel1(0,k*r)
    def integral(points, weights):
        return k*k*np.einsum('qt,qr,q->tr',green(points,tx),green(points,rx),weights)
    nodes,weights = np.polynomial.legendre.leggauss(degree)
    xmin,xmax,ymin,ymax = ROI
    wall = np.zeros((2,16),complex)
    ny,nx = s.shape
    for row,col in np.argwhere(s):
        dx,dy = (xmax-xmin)/nx,(ymax-ymin)/ny
        x = xmin+(col+.5)*dx+nodes*dx/2
        yy = ymin+(row+.5)*dy+nodes*dy/2
        xx,zz = np.meshgrid(x,yy)
        wall += integral(np.column_stack((xx.ravel(),zz.ravel())),np.outer(weights,weights).ravel()*dx*dy/4)
    radius = (nodes+1)*.12/2
    angles = np.arange(4*degree)*2*np.pi/(4*degree)
    rr,aa = np.meshgrid(radius,angles,indexing='ij')
    points = theta[2:4]+np.column_stack((rr.ravel()*np.cos(aa.ravel()),rr.ravel()*np.sin(aa.ravel())))
    disk_weights = np.repeat(weights*.12/2*radius,4*degree)*2*np.pi/(4*degree)
    return (green(rx,tx).T+(theta[0]-1+1j*theta[1])*wall+(3.+.02j)*integral(points,disk_weights)).ravel()


def initializer(y, support, theta, origin):
    theta = np.array(theta,float,copy=True)
    profile = optimize.profile_gain_phase(y,born_prediction(support,theta))
    theta[4:] = profile['gain_offset'],profile['phase']
    return dict(support=np.asarray(support,bool).copy(),theta=theta,
                support_id=optimize.support_id(support),origin=origin,
                cheap_score=float(profile['objective']-np.log(optimize.support_mass(support))+np.log(2.)),
                cheap_model='COMMON_ANALYTIC_FREE_SPACE_BORN_PROPOSAL_ONLY',
                cheap_profile=profile,physics_rhs=0)


def common_pool(y, bankdata, grid):
    ranking = optimize.bank_rank(y,bankdata,grid)
    coarse = sorted(ranking['best_by_support'],key=lambda x:(x['ranking_objective'],x['support_id'],x['row_id']))[:2]
    pool = [initializer(y,row['support'],row['theta'],'measurement_bank_top_two_distinct_supports') for row in coarse]
    mapped = optimize.coarse_support_to3x3(coarse[0]['support'])
    if not validate_geometry(mapped['support'])['legal']:
        raise ValueError('Center-sampled coarse-to-fine mapping is illegal')
    neighbors = optimize.legal_support_neighbors(mapped['support'],include_noop=True)
    screened = optimize.born_region_screen(y,coarse[0]['theta'],[r['support'] for r in neighbors])
    # Profile all analytic predictions with the same proper nuisance prior.
    fine = [initializer(y,row['support'],coarse[0]['theta'],'public_3x3_neighbor_Born_proposal') for row in screened]
    fine.sort(key=lambda r:(r['cheap_score'],r['support_id']))
    pool.append(fine[0])
    pool.sort(key=lambda r:(r['cheap_score'],r['support_id']))
    return pool,coarse[0],dict(bank_all_supports_retained=True,
        bank_support_count=len(ranking['best_by_support']),representation_change=mapped['representation_change'],
        fine_neighbor_ids=[r['support_id'] for r in fine],all_methods_receive_same_pool=True)


class RecordedEvaluator:
    """Record paid data/J even if a subsequent optimization trial fails."""
    def __init__(self, solver, path):
        self.solver,self.path = solver,path
        self.book,self.stage = solver.book,solver.stage
        self.records = []
    def evaluate(self,world,**kwargs):
        ev = self.solver.evaluate(world,**kwargs)
        self.records.append(dict(theta=np.array([world.epsinf,world.ell,*world.pose]),
            support=world.support.copy(),data=np.asarray(ev.data).copy(),
            J=None if ev.jacobian is None else np.asarray(ev.jacobian).copy(),metadata=dict(ev.metadata)))
        self.save()
        return ev
    def save(self, trials=None):
        arrays = {}
        for i,row in enumerate(self.records):
            prefix=f'paid_{i:03d}_'
            arrays[prefix+'physical_theta']=row['theta']
            arrays[prefix+'support']=row['support']
            arrays[prefix+'data']=row['data']
            arrays[prefix+'metadata_json']=np.array(json.dumps(row['metadata'],default=json_default))
            if row['J'] is not None:
                arrays[prefix+'J']=row['J']
        for i,row in enumerate(trials or []):
            prefix=f'trial_{i:03d}_'
            for key,value in row.items():
                arrays[prefix+key]=value
        save_npz(self.path,arrays)


def run_fit(y,start,limit,book,directory,name):
    path = directory/(name+'.npz')
    solver=Evaluator(book,'A0')
    recorder=RecordedEvaluator(solver,path)
    before=book.totals('A0');began=time.perf_counter();cpu=time.process_time()
    trials=[]
    def traced_least_squares(fun,x0,**kwargs):
        actual_jac=kwargs['jac']
        def residual(theta):
            row=dict(kind=np.array('residual_trial'),theta=theta.copy())
            trials.append(row)
            recorder.save(trials)
            value=fun(theta)
            row['residual']=value.copy()
            row['rawtrainwhite']=value[:-2].copy()
            key=theta[:4]
            matching=[r for r in recorder.records if np.array_equal(r['theta'],key)]
            if matching:
                row['pred']=(1+theta[4])*np.exp(1j*theta[5])*matching[-1]['data']
            recorder.save(trials)
            return value
        def jac(theta):
            row=dict(kind=np.array('jacobian_trial'),theta=theta.copy())
            trials.append(row);recorder.save(trials)
            value=actual_jac(theta)
            row['residual_J']=value.copy();recorder.save(trials)
            return value
        kwargs['jac']=jac
        return scipy_least_squares(residual,x0,**kwargs)
    try:
        with patch.object(optimize,'least_squares',traced_least_squares):
            fit=optimize.fit_joint(y,start['support'],start['theta'],evaluator=recorder,max_nfev=limit)
        fit['hierarchical_map_objective']=hierarchical_score(fit)
        fit['fit_status']='CONVERGED' if fit['success'] else 'MAX_NFEV_OR_NONCONVERGENCE'
        fit['trace_path']=str(path)
        fit['initial_origin']=start['origin']
        fit['max_nfev_allocated']=limit
        arrays={}
        for key in ('theta','support','pred','rawtrainwhite'):
            arrays['final_'+key]=fit[key]
        for i,row in enumerate(fit['saved_evaluations']):
            arrays[f'cache_{i:03d}_physical_theta']=row['physical_theta']
            arrays[f'cache_{i:03d}_data']=row['F']
            if row['J'] is not None: arrays[f'cache_{i:03d}_J']=row['J']
        recorder.save(trials)
        # Final/cache arrays have their own companion; paid/trial receipt remains
        # durable throughout failures and uninterrupted optimization alike.
        save_npz(directory/(name+'_final.npz'),arrays)
        return {k:v for k,v in fit.items() if k not in ('traces','saved_evaluations')}
    except BudgetExceeded:
        recorder.save(trials)
        raise
    except Exception as exc:
        recorder.save(trials)
        return dict(fit_status='FAILED',success=False,error_type=type(exc).__name__,error=str(exc),
                    trace_path=str(path),max_nfev_allocated=limit,initial_origin=start['origin'],
                    nfev=sum('residual' in t for t in trials),
                    njev=sum('residual_J' in t for t in trials),
                    residual_trials_attempted=sum(t['kind'].item()=='residual_trial' for t in trials),
                    stage_cost_diffs={k:book.totals('A0')[k]-before[k] for k in before})
    finally:
        solver.clear()
        write_json(directory/(name+'_receipt.json'),dict(trace_path=str(path),
            elapsed_wall_seconds=time.perf_counter()-began,elapsed_process_cpu_seconds=time.process_time()-cpu,
            stage_cost_diffs={k:book.totals('A0')[k]-before[k] for k in before},
            allocated_max_nfev=limit,paid_records=len(recorder.records),trials=len(trials)))


def run_method(method,y,pool,coarse,book,directory):
    before=book.totals('A0')
    fits=[];selection=None;acceptance=None
    if method=='FULL_J_TR':
        fits=[run_fit(y,pool[0],6,book,directory,'single')]
    elif method=='MULTISTART_JOINT_MAP':
        starts=[]
        for row in pool:
            if row['support_id'] not in {r['support_id'] for r in starts}:starts.append(row)
            if len(starts)==2:break
        for i,start in enumerate(starts):fits.append(run_fit(y,start,3,book,directory,f'start_{i}'))
    elif method=='HYBRID_MULTISCALE_BIRTH_DEATH_EXCHANGE':
        start=dict(support=coarse['support'],theta=coarse['theta'],origin='measurement_bank_best_coarse')
        first=run_fit(y,start,3,book,directory,'coarse')
        fits.append(first)
        if first.get('fit_status')!='FAILED':
            mapped=optimize.coarse_support_to3x3(first['support'])
            neighbors=optimize.legal_support_neighbors(mapped['support'],include_noop=True)
            screen=optimize.born_region_screen(y,first['theta'],[r['support'] for r in neighbors])
            proposed=[initializer(y,r['support'],first['theta'],'reversible_public_3x3_neighbor') for r in screen]
            proposed.sort(key=lambda r:(r['cheap_score'],r['support_id']))
            second=run_fit(y,proposed[0],3,book,directory,'fine_proposal')
            fits.append(second)
            acceptance=dict(rule='FULL_PHYSICAL_HIERARCHICAL_MAP_DECREASE_OR_RETAIN_COARSE_NOOP',
                proposal_support_id=proposed[0]['support_id'],neighbor_ids=[r['support_id'] for r in proposed],
                accepted=second.get('fit_status')!='FAILED' and second['hierarchical_map_objective']<first['hierarchical_map_objective'])
            selection=second if acceptance['accepted'] else first
    else:raise ValueError(method)
    usable=[f for f in fits if f.get('fit_status')!='FAILED']
    if selection is None and usable:
        selection=min(usable,key=lambda f:(f['hierarchical_map_objective'],f['support_id']))
    return dict(method=method,status='CANDIDATE' if selection is not None else 'FAILED',fits=fits,
        selection=selection,acceptance=acceptance,max_nfev_total=6,
        actual_nfev=sum(f.get('nfev',0) for f in fits),actual_njev=sum(f.get('njev',0) for f in fits),
        stage_cost_diffs={k:book.totals('A0')[k]-before[k] for k in before},
        all_methods_support_access='COMMON_POOL_2x2_AND_3x3',performance_claim=False,
        posterior_claim=False,holdout_status='NOT_READ_SEPARATE_SCORER_REQUIRED')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    args=parser.parse_args();root=args.root.resolve()
    gate=json.loads((root/'results/e0/gate.json').read_text())
    if gate.get('decision')!='PASS':raise RuntimeError('E0 PASS required')
    # Read-only frozen preregistration guards; do not manufacture missing fields.
    json.loads((root/'EXPERIMENT_CONTRACT.json').read_text())
    supplement=json.loads((root/'A0_ALGORITHM_SUPPLEMENT.json').read_text())
    levels=supplement['hierarchical_support_prior']['level_probability']
    if levels!={'2':.5,'3':.5}:raise RuntimeError('Required preregistered .5/.5 grid-level prior absent')
    rule=json.loads((root/'results/a0/preparation/DISCREPANCY_RULE.json').read_text())
    upper=float(rule['upper_noise_energy'])
    if not np.isfinite(upper) or upper<=0:raise ValueError('Invalid frozen discrepancy envelope')
    online=root/'inputs/online'
    with np.load(online/'COARSE_BANK.npz',allow_pickle=False) as raw:
        if int(raw['completed'])!=int(raw['total']):raise RuntimeError('Proposal bank incomplete')
        bankdata=raw['data'][:,:2].copy()
    with np.load(online/'COARSE_GRID.npz',allow_pickle=False) as raw:grid={k:raw[k].copy() for k in raw.files}
    with np.load(online/'A0_OBSERVATIONS.npz',allow_pickle=False) as raw:
        observations=raw['base'].copy();parent_ids=raw['parent_id'].copy();repeats=raw['noise_repeat'].copy()
    if observations.shape!=(24,2,16):raise ValueError('Frozen A0 panel requires 24 base observations')
    out=root/'results/a0/optimization';out.mkdir(parents=True,exist_ok=False)
    book=BudgetBook(root/'results/budget_ledger.json')
    records=[dict(observation_id=i,parent_id=int(parent_ids[i]),noise_repeat=int(repeats[i]),
                  method=m,status='NOT_RUN') for i in range(24) for m in METHODS]
    report=dict(scope='MEASUREMENT_ONLY_A0_CANDIDATES_NOT_ACCEPTANCE',status='RUNNING',
        records=records,noise_discrepancy_rule=rule,holdout_read=False,
        performance_claim=False,posterior_claim=False,level_prior=levels,
        budget_start=book.totals('A0'),common_pool_model='ANALYTIC_FREE_SPACE_BORN_PROPOSAL_ONLY')
    write_json(out/'RESULTS.json',report)
    current=None
    try:
        for i,y in enumerate(observations):
            pool,coarse,pool_metadata=common_pool(y,bankdata,grid)
            obsdir=out/f'observation_{i:02d}';obsdir.mkdir()
            write_json(obsdir/'INITIALIZER_POOL.json',dict(pool=pool,metadata=pool_metadata))
            for method in METHODS:
                current=next(r for r in records if r['observation_id']==i and r['method']==method)
                directory=obsdir/method;directory.mkdir()
                result=run_method(method,y,pool,coarse,book,directory)
                if result['selection'] is not None:
                    energy=float(result['selection']['rawtrainwhite']@result['selection']['rawtrainwhite'])
                    result['training_noise_envelope_compatible']=energy<=upper
                    result['training_noise_energy']=energy
                    result['training_qualification']='ONE_SIDED_NOISE_ENVELOPE_ONLY_NOT_A0_PASS'
                current.update(result)
                write_json(out/'RESULTS.json',report)
                print(json.dumps(dict(observation=i,method=method,status=current['status'],cost=book.totals('A0'))),flush=True)
        report['status']='CANDIDATES_COMPLETE_AWAIT_SEPARATE_HOLDOUT_AND_PARENT_JUDGMENT'
    except BudgetExceeded as exc:
        if current is not None:current.update(status='STOPPED_BUDGET',error=str(exc))
        report.update(status='STOPPED_A0_BUDGET',error=str(exc))
        for row in records:
            if row['status']=='NOT_RUN':row['reason']='A0_SHARED_CAMPAIGN_BUDGET_STOP'
    except Exception as exc:
        if current is not None:current.update(status='FAILED',error_type=type(exc).__name__,error=str(exc))
        report.update(status='STOPPED_IMPLEMENTATION_OR_INPUT_FAILURE',error_type=type(exc).__name__,error=str(exc))
        for row in records:
            if row['status']=='NOT_RUN':row['reason']='STAGE_STOP_AFTER_FAILURE'
    finally:
        report['budget_end']=book.totals('A0')
        write_json(out/'RESULTS.json',report)
        write_json(out/'gate.json',dict(decision='INCONCLUSIVE',execution_status=report['status'],
            reason='Only training candidates generated; independent holdout and scientific review remain separate',holdout_read=False))


if __name__=='__main__':
    main()
