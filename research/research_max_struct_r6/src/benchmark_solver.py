"""Four frozen R6 arms and observation-only original R4 physical corrector.

Numerical provenance:
* cie_single / cie_dual: FBG-R3 mature run_cie, locally observed copy;
* global_dual: research_max_struct_r2.src.cie_global.run_cie;
* a_material: research_max_struct_r5.src.cie_factorial.run_factorial_cie;
* corrector: research_max_struct_r4.experiments.cooperative_paths.refine.
No truth, scoring, parallel queue, or parameter search is implemented here.
"""
from collections import Counter
from copy import deepcopy
import time
import numpy as np
from scipy import linalg as la
from scipy.optimize import minimize
from research.research_max_struct_r6.src.mature_cie_observed import run_cie as mature_cie
from research.research_max_struct_r2.src.cie_global import run_cie as global_cie
from research.research_max_struct_r5.src.cie_factorial import run_factorial_cie

ARMS = ('cie_single', 'cie_dual', 'global_dual', 'a_material')
SOURCE_PATHS = {
    'cie_single': 'research/research_fbg_r3_experiments/src/cie_tsom_baseline.py::run_cie',
    'cie_dual': 'research/research_fbg_r3_experiments/src/cie_tsom_baseline.py::run_cie',
    'global_dual': 'research/research_max_struct_r2/src/cie_global.py::run_cie',
    'a_material': 'research/research_max_struct_r5/src/cie_factorial.py::run_factorial_cie',
    'corrector': 'research/research_max_struct_r4/experiments/cooperative_paths.py::refine',
}


def difference(after, before):
    return {k: after.get(k, 0)-before.get(k, 0) for k in set(after)|set(before)}


def arm_config(arm, common, Y, E, continuation_stages=None):
    """Freeze the shared source metrics once in original public source coordinates.

    Continuation must be explicit in the config. It is appended to the same
    run_cie invocation so Fourier coefficients, material and J are transported.
    This function adds no wall/CPU/function-evaluation caps.
    """
    if arm not in ARMS:
        raise ValueError('Unknown fixed arm: '+str(arm))
    cfg = deepcopy(common)
    cfg.pop('host')
    for key, value in {'coefficient_tol': 0., 'miniter': 50,
                       'material_inner_maxiter': 60,
                       'material_inner_ftol': 1e-12,
                       'material_inner_gtol': 1e-8}.items():
        if key in cfg and cfg[key] != value:
            raise ValueError('R6 frozen setting mismatch: '+key)
        cfg[key] = value
    if tuple(cfg.get('material_box', ())) != (.02, 20.):
        raise ValueError('R6 original material box must be [.02,20]')
    if float(common['host']) != .1:
        raise ValueError('R6 shared original host must be .1')
    if ('dpcs_rank' in cfg) == ('dpcs_rtol' in cfg):
        raise ValueError('Specify exactly one public-frozen DPCS rule')
    # These are the original paper normalizations, identical across all arms.
    for key in ('data_metric', 'state_metric', 'data_weights', 'state_weights',
                'stages', 'sourcewise_cg', 'current_transform', 'material_transform',
                'current_solver', 'method'):
        if key in cfg:
            raise ValueError('Arm-specific setting belongs to fixed R6 arm: '+key)
    cfg.update(data_weights=1/np.maximum(np.sum(abs(Y)**2, axis=0), 1e-300),
               state_weights=1/np.maximum(np.sum(abs(E)**2, axis=0), 1e-300),
               method='cie', current_solver='pr', sourcewise_cg=arm.startswith('cie_'))
    stages = []
    for mf, bm in ((3, 1.), (6, .2), (10, .2)):
        stages.append(dict(mf=mf, beta_material=bm,
                           beta_current=bm if arm == 'cie_single' else 1., maxiter=500))
    for stage in continuation_stages or []:
        stage = dict(stage)
        if set(stage) != {'mf', 'maxiter'} or stage['mf'] != 12 or int(stage['maxiter']) <= 0:
            raise ValueError('Explicit natural continuation accepts mf12/maxiter only')
        stages.append(dict(mf=12, maxiter=int(stage['maxiter']), beta_material=.2,
                           beta_current=.2 if arm == 'cie_single' else 1.))
    cfg['stages'] = stages
    if arm == 'a_material':
        cfg.update(current_transform='scalar', material_transform='onehop',
                   cache_material_value_gradient=True)
    return cfg


def run_auxiliary(model, Y, m0, arm, config):
    fn = mature_cie if arm.startswith('cie_') else global_cie if arm == 'global_dual' else run_factorial_cie
    wall, cpu, before = time.perf_counter(), time.process_time(), dict(model.cost)
    result = fn(model, Y, m0, config)
    result.setdefault('execution_receipt', dict(status='COMPLETED',
        wall_seconds=time.perf_counter()-wall, cpu_seconds=time.process_time()-cpu,
        operation_counts={k: v for k, v in result['cost'].items() if k not in ('model', 'cpu_seconds', 'wall_seconds')},
        physical_model_delta=difference(dict(model.cost), before)))
    result['execution_receipt']['source_entrypoint'] = SOURCE_PATHS[arm]
    delta=result['execution_receipt']['physical_model_delta']
    result['execution_receipt'].setdefault('physical_LU',delta.get('LU',0))
    result['execution_receipt'].setdefault('physical_complex_RHS',delta.get('actual_complex_RHS',0))
    return result


class PaidStateObserver:
    """Record predictions already evaluated by the solver, never call forward.

    MaxwellModel._current is the state retained by its real value-gradient
    function. Its arrays are copied without retaining the dense LU/L matrix.
    """
    def __init__(self, model, Y, save_snapshot, gradient_checkpoints=(100, 500, 1000)):
        self.model, self.Y, self.save_snapshot = model, Y, save_snapshot
        self.gradient_checkpoints = set(gradient_checkpoints)
        self.validation_records = []
        self.gradient_records = []
        self.best = None
        self.last_evaluation = None
        self.phase = 'auxiliary'
        self.wall, self.cpu = time.perf_counter(), time.process_time()

    def __getattr__(self, name):
        return getattr(self.model, name)

    def record(self, m, state, value=None, gradient=None):
        residual = state['Y']-self.Y
        value = (self.model.physical_score(m, self.Y) if value is None else float(value))
        row = dict(m=np.asarray(m).copy(), j=state['j'].copy(), Q=state['Q'].copy(),
                   Y_prediction=state['Y'].copy(), objective=value,
                   physical_relative_residual=float(la.norm(residual)/la.norm(self.Y)),
                   state_relative_residual=float(state['state_relative_residual']))
        if gradient is not None:
            row['gradient'] = np.asarray(gradient).copy()
        if self.best is None or value < self.best['objective']:
            self.best = row
        return row

    def state(self, m):
        state = self.model.state(m)
        if self.phase == 'auxiliary':
            row = self.record(m, state)
            self.validation_records.append(row)
            self.save_snapshot('initial' if len(self.validation_records) == 1 else
                               f'stage_{len(self.validation_records)-2:02d}_physical', row)
        return state

    def physical_score(self, m, Y):
        # Same original physical objective; model.state's own cache applies.
        state = self.model.state(m)
        value = self.model.physical_score(m, Y)
        self.last_score = self.record(m, state, value=value)
        return value

    def data_value_gradient(self, m, Y):
        value, gradient = self.model.data_value_gradient(m, Y)
        state = self.model._current
        if state is None or not np.array_equal(state['m'], m):
            raise RuntimeError('Original value-gradient no longer retains the paid state')
        row = self.record(m, state, value=value, gradient=gradient)
        self.last_evaluation = row
        index = len(self.gradient_records)+1
        self.gradient_records.append(dict(evaluation=index, objective=float(value),
            best_objective=float(self.best['objective']), gradient_inf=float(la.norm(gradient, np.inf)),
            wall_seconds=time.perf_counter()-self.wall, cpu_seconds=time.process_time()-self.cpu,
            physical_model_cost=dict(self.model.cost)))
        if index in self.gradient_checkpoints:
            self.save_snapshot(f'corrector_evaluation_{index:04d}', row)
        return value, gradient


def refine_observed(model, Y, m, maximum=1000, accepted_checkpoints=(100, 500, 1000)):
    """Original R4 refine with a callback observing already paid evaluations.

    Only additions to the original numerical body are this callback and timing/
    state records. maxfun=4*maximum, scale=1e6, ftol=1e-12, gtol=.01,
    maxls=30 and the real box are the original R4 settings.
    """
    before=Counter(model.cost);history=[];scale=1e6;start=time.perf_counter()
    cpu=time.process_time();accepted=[];wanted=set(accepted_checkpoints)
    model.phase = 'corrector'
    best=[m.copy(),model.physical_score(m,Y)]
    def objective(x):
        value,g=model.data_value_gradient(x,Y)
        history.append([float(value),float(la.norm(g,np.inf))])
        if value<best[1]:best[:]=[x.copy(),float(value)]
        return scale*value,scale*g
    def callback(x):
        row=model.last_evaluation
        if row is None or not np.array_equal(x,row['m']):
            raise RuntimeError('Accepted corrector point is absent from last paid evaluation')
        iteration=len(accepted)+1
        accepted.append(dict(iteration=iteration, evaluation=len(history), objective=row['objective'],
            best_objective=float(best[1]), gradient_inf=float(la.norm(row['gradient'],np.inf)),
            wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,
            physical_model_cost=dict(model.cost)))
        if iteration in wanted:model.save_snapshot(f'corrector_accepted_{iteration:04d}',row)
    result=minimize(objective,m,method='L-BFGS-B',jac=True,bounds=[(.02,20.)]*model.n,
        callback=callback,
        options=dict(maxiter=maximum,maxfun=4*maximum,ftol=1e-12,gtol=.01,maxls=30))
    return best[0],dict(stop=str(result.message),iterations=int(result.nit),evaluations=int(result.nfev),
        success=bool(result.success),history=history,accepted_history=accepted,
        evaluation_history=model.gradient_records,cost=difference(dict(model.cost),dict(before)),
        wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,
        checkpoint_semantics='Snapshots copy already paid state/value/gradient; no extra forward or adjoint',
        best_objective=float(best[1]),source_entrypoint=SOURCE_PATHS['corrector'])
