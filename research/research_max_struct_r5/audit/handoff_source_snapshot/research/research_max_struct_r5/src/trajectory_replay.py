"""R5 observation-only callbacks around the unchanged historical optimizers."""
from pathlib import Path
import json
import time
import numpy as np

CHECKPOINTS = (1, 2, 10, 50, 150, 300, 500)


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False))


class TrajectoryRecorder:
    def __init__(self, directory, model, m0):
        self.directory = Path(directory)
        self.model = model
        self.previous = dict(m=np.asarray(m0).copy())
        self.checkpoints = []
        self.callback_wall_seconds = 0.
        self.callback_cpu_seconds = 0.
        self.save_wall_seconds = 0.
        self.extra_green_wall_seconds = 0.
        self.extra_green_cpu_seconds = 0.
        self.extra_green_columns = 0
        self.callback_calls = 0

    def iteration(self, state):
        wall, cpu = time.perf_counter(), time.process_time()
        self.callback_calls += 1
        if state['iteration'] in CHECKPOINTS:
            start, ct = time.perf_counter(), time.process_time()
            q = self.model.E + self.model.G @ state['J']
            self.extra_green_wall_seconds += time.perf_counter()-start
            self.extra_green_cpu_seconds += time.process_time()-ct
            self.extra_green_columns += state['J'].shape[1]
            name = f"checkpoint_s{state['stage']}_i{state['iteration']:03d}.npz"
            arrays = dict(stage=state['stage'], iteration=state['iteration'],
                m_before=self.previous['m'], m_after=state['m'],
                delta_m=state['m']-self.previous['m'],
                J_after_current=state['J'], A=state['coefficients'], Q=q)
            if 'J' in self.previous:
                arrays['J_previous_callback'] = self.previous['J']
                arrays['A_previous_callback'] = self.previous['coefficients']
                arrays['previous_stage'] = self.previous['stage']
                arrays['previous_iteration'] = self.previous['iteration']
            start = time.perf_counter()
            np.savez_compressed(self.directory/name, **arrays)
            self.save_wall_seconds += time.perf_counter()-start
            self.checkpoints.append(name)
        # Solver callback supplies private copies. Retaining these references
        # performs no forward/Green operation on non-checkpoint iterations.
        self.previous = state
        self.callback_wall_seconds += time.perf_counter()-wall
        self.callback_cpu_seconds += time.process_time()-cpu

    def stage(self, state):
        wall, cpu = time.perf_counter(), time.process_time()
        np.savez_compressed(self.directory/f"stage_{state['stage']}.npz",
            **{key:state[key] for key in ('m','J','Q','coefficients')})
        self.save_wall_seconds += time.perf_counter()-wall
        self.callback_wall_seconds += time.perf_counter()-wall
        self.callback_cpu_seconds += time.process_time()-cpu
        print(json.dumps(dict(event='STAGE_SAVED', stage=state['stage'])), flush=True)

    def receipt(self):
        return dict(callback_calls=self.callback_calls, checkpoints=self.checkpoints,
            callback_wall_seconds=self.callback_wall_seconds,
            callback_cpu_seconds=self.callback_cpu_seconds,
            callback_save_wall_seconds=self.save_wall_seconds,
            extra_checkpoint_Green_columns=self.extra_green_columns,
            extra_checkpoint_Green_calls=len(self.checkpoints),
            extra_Green_wall_seconds=self.extra_green_wall_seconds,
            extra_Green_cpu_seconds=self.extra_green_cpu_seconds,
            physical_LU=0, physical_complex_RHS=0,
            scope='Callbacks only; solver-created callback array copies stay in auxiliary time. Q=E+GJ exactly once per frozen checkpoint; no other extra physics.')


def difference(a, b):
    a, b = np.asarray(a), np.asarray(b)
    delta = a-b
    return dict(l2=float(np.linalg.norm(delta)), max_abs=float(np.max(abs(delta))),
        relative_l2=float(np.linalg.norm(delta)/max(np.linalg.norm(b), 1e-300)),
        exact_array_equal=bool(np.array_equal(a,b)))


def compare_reference(repo, config, auxiliary):
    reference = config['reference']
    stages = []
    if reference['kind'] == 'stage_npz':
        for state in auxiliary['stage_states']:
            path = repo/reference['directory']/f"stage_{state['stage']}.npz"
            with np.load(path,allow_pickle=False) as old:
                stages.append(dict(stage=state['stage'], reference=str(path.relative_to(repo)),
                    differences={key:difference(state[key],old[key]) for key in ('m','J','Q','coefficients')}))
    else:
        path = repo/reference['prediction']
        with np.load(path,allow_pickle=False) as old:
            for state in auxiliary['stage_states']:
                stages.append(dict(stage=state['stage'], reference=str(path.relative_to(repo)),
                    differences=dict(m=difference(state['m'],old['candidate_m'][state['stage']+1])),
                    J_comparison='NOT_AVAILABLE: historical auxiliary stage J was not saved',
                    Q_comparison='NOT_AVAILABLE: historical auxiliary stage Q was not saved',
                    coefficients_comparison='NOT_AVAILABLE: historical auxiliary stage A was not saved'))
    return dict(status='COMPARED_WITH_REAL_DIFFERENCES', stages=stages,
        interpretation='CPU/LAPACK/platform trajectory drift is measured; no claim of bitwise identity. Historical final physical j/Q belong to corrected material and are not auxiliary stage J/Q.')
